/*
 * Implementation of Goldberg-Tarjan's parallel push-relabel algorithm.
 * Description can be found in Goldberg, Andrew and Tarjan, Robert, A New
 * Approach to the Maximum-Flow Problem, J. ACM, 1988.
 *
 * This implementation is also based on detailed pseudocode presented in
 * Baumstark, Niklas, Speeding up Maximum Flow Computations on Shared-Memory
 * Platforms, KIT, Karlsruhe, 2014.
 */

#ifndef MAXFLOW_PARALLEL_PUSH_RELABEL_H
#define MAXFLOW_PARALLEL_PUSH_RELABEL_H

#include "util.h"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <omp.h>

#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif

namespace parallel_push_relabel {
class max_flow_instance {
    // this actually has a big perf impact
    struct alignas(CACHE_LINE_SIZE) vertex {
        uint32_t excess{0};
        std::atomic<uint32_t> new_excess{0};
        uint32_t label;
        uint32_t new_label;
        std::atomic_flag discovered = ATOMIC_FLAG_INIT;
    };

    std::vector<std::vector<edge>> residual;
    std::unique_ptr<vertex[]> vertices;
    std::unique_ptr<uint32_t[]> active{};
    buffer_pool<uint32_t> pool;
    uint32_t s, t, relabel_threshold, nactive;
    std::size_t relabel_progress;
    const uint32_t nthreads;

  public:
    max_flow_instance(
        std::vector<std::vector<edge>> graph, uint32_t source, uint32_t sink,
        std::size_t thread_count = static_cast<size_t>(omp_get_max_threads()))
        : residual(std::move(graph)),
          vertices(std::make_unique<vertex[]>(residual.size())),
          active(std::make_unique<uint32_t[]>(residual.size())),
          pool(buffer_pool<uint32_t>(thread_count,
                                                      residual.size())),
          s(source), t(sink), nactive(0), relabel_progress(0),
          nthreads(thread_count) {
        omp_set_num_threads(static_cast<int>(nthreads));
        init();
    }

    uint64_t nphase = 0;
    uint64_t npush = 0;
    uint64_t nupdate = 0;

    uint32_t find_max_flow() noexcept {
        find_max_flow_inner();

#ifdef DEBUG
        std::cout << "global updates:\t" << nupdate << std::endl;
        std::cout << "phase cnt: " << nphase << std::endl;
        std::cout << "pushes: " << npush << std::endl;
#endif
        return vertices[t].new_excess + vertices[t].excess;
    }

    void preflow_to_flow() {
        std::swap(s, t);
        find_max_flow_inner();
        std::swap(s, t);
#ifdef DEBUG
        for (std::size_t i = 0; i < residual.size(); ++i)
            if (i != s && i != t)
                if (vertices[i].excess > 0)
                    std::cerr << "Excess violation: vertex " << i << ", excess "
                              << vertices[i].excess << '\n';
#endif
    }

    auto steal_network() { return std::move(residual); }

  private:
    static constexpr uint32_t ALPHA = 6, BETA = 12;
    static constexpr double GLOBAL_RELABEL_FREQ = 0.5;

    void init() noexcept {
#pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < residual[s].size(); ++i) {
            auto &edge = residual[s][i];
            auto &rev_edge = residual[edge.dst][edge.rev_index];

            vertices[edge.dst].excess = edge.cap;
            edge.rev_cap += edge.cap;
            rev_edge.cap += edge.cap;
            rev_edge.rev_cap -= edge.cap;
            edge.cap = 0;
        }

        uint32_t m = 0;
        for (std::size_t i = 0; i < residual.size(); ++i)
            m += residual[i].size();
        relabel_threshold = residual.size() * ALPHA + m / 2;
    }

    void find_max_flow_inner() {
        global_relabel();

        for (;;) {
            if (nactive == 0)
                return;

            ++nphase;
            uint64_t push_per_phase = 0;

#pragma omp parallel
            {
#pragma omp for schedule(static) reduction(+ : push_per_phase)
                for (uint32_t i = 0; i < nactive; ++i) {
                    auto thr_id = omp_get_thread_num();
                    auto vid = active[i];
                    auto &vertex = vertices[vid];
                    if (vertex.label == residual.size())
                        continue;
                    push(vid, vertex.label, thr_id, push_per_phase);
                }
// stage 2
#pragma omp for schedule(static) reduction(+ : relabel_progress)
                for (uint32_t i = 0; i < nactive; ++i) {
                    auto thr_id = omp_get_thread_num();
                    auto vertex = active[i];
                    relabel(vertex, thr_id, relabel_progress);
                }
// stage 3
#pragma omp for schedule(static)
                for (uint32_t i = 0; i < nactive; ++i) {
                    auto &vertex = vertices[active[i]];
                    vertex.label = vertex.new_label;
                    vertex.discovered.clear(std::memory_order_relaxed);
                }
// stage 4
#pragma omp single
                nactive = pool.swap_data(active);

#pragma omp for schedule(static)
                for (uint32_t i = 0; i < nactive; ++i) {
                    auto &vertex = vertices[active[i]];
                    vertex.excess +=
                        vertex.new_excess.load(std::memory_order_relaxed);
                    vertex.new_excess.store(0, std::memory_order_relaxed);
                    vertex.discovered.clear(std::memory_order_relaxed);
                }
            }

            if (relabel_progress * GLOBAL_RELABEL_FREQ >= relabel_threshold ||
                push_per_phase == 0) {
                relabel_progress = 0;
                global_relabel();
            }

            npush += push_per_phase;
        }
    }

    inline void push(const uint32_t vertex, const uint32_t label, int thr_id,
                     uint64_t &push_cnt) noexcept {
        const auto target_label = label - 1;
        for (auto &edge : residual[vertex]) {
            if (edge.cap > 0 && vertices[edge.dst].label == target_label) {
                auto flow = std::min(vertices[vertex].excess, edge.cap);
                if (edge.dst != s && edge.dst != t)
                    if (!vertices[edge.dst].discovered.test_and_set(
                            std::memory_order_relaxed))
                        pool.push_back(edge.dst, static_cast<size_t>(thr_id));
                ++push_cnt;
                vertices[vertex].excess -= flow;
                vertices[edge.dst].new_excess.fetch_add(
                    flow, std::memory_order_relaxed);
                edge.cap -= flow;
                edge.rev_cap += flow;
                residual[edge.dst][edge.rev_index].rev_cap -= flow;
                residual[edge.dst][edge.rev_index].cap += flow;
                if (vertices[vertex].excess == 0)
                    return;
            }
        }
    }

    inline void relabel(const uint32_t vertex, const int thr_id,
                        std::size_t &relabel_progress) noexcept {
        if (vertices[vertex].excess > 0 ||
            vertices[vertex].label == residual.size()) {
            relabel_progress += BETA;
            vertices[vertex].new_label = calculate_new_label(vertex);
            relabel_progress += residual[vertex].size();
            if (vertices[vertex].new_label == residual.size()) {
                vertices[vertex].excess += vertices[vertex].new_excess;
                vertices[vertex].new_excess = 0;
                return;
            }

            if (!vertices[vertex].discovered.test_and_set(
                    std::memory_order_relaxed))
                pool.push_back(vertex, static_cast<size_t>(thr_id));
        } else
            vertices[vertex].new_label = vertices[vertex].label;
    }

    inline uint32_t calculate_new_label(const uint32_t vertex) noexcept {
        uint32_t increase_to = residual.size() - 1;
        for (auto &edge : residual[vertex]) {
            if (edge.cap == 0)
                continue;

            increase_to = std::min(increase_to, vertices[edge.dst].label);
        }
        return increase_to + 1;
    }

    void global_relabel() noexcept {
        ++nupdate;
        const auto not_reached = residual.size();

#pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < residual.size(); ++i)
            vertices[i].label = not_reached;

        vertices[t].label = 0;
        vertices[t].discovered.test_and_set();
        assert(pool.empty());
        active[0] = t;
        std::size_t current_queue_size = 1;
        uint32_t current_distance = 0;

        while (current_queue_size > 0) {
#pragma omp parallel for schedule(static)
            for (std::size_t i = 0; i < current_queue_size; ++i) {
                auto thr_id = omp_get_thread_num();
                auto current_vertex = active[i];

                for (auto edge : residual[current_vertex]) {
                    if (edge.rev_cap > 0) {
                        if (!vertices[edge.dst].discovered.test_and_set(
                                std::memory_order_relaxed)) {
                            vertices[edge.dst].label = current_distance + 1;
                            pool.push_back(edge.dst,
                                           static_cast<std::size_t>(thr_id));
                        }
                    }
                }
            }
            current_queue_size = pool.swap_data(active);
            ++current_distance;
        }

#pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < residual.size(); ++i) {
            auto thr_id = omp_get_thread_num();
            if (vertices[i].label != not_reached && vertices[i].excess > 0 &&
                i != t)
                pool.push_back(i, static_cast<size_t>(thr_id));
            vertices[i].discovered.clear(std::memory_order_relaxed);
        }

        nactive = pool.swap_data(active);
    }
};
} // namespace parallel_push_relabel

#endif // MAXFLOW_PARALLEL_PUSH_RELABEL_H
