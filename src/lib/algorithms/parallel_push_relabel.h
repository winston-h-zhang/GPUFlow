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

#include "../common_types.h"
#include "../data_structures/thread_local_buffer_pool.h"
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
    struct alignas(CACHE_LINE_SIZE) vertex {
        uint32_t excess{0};
        std::atomic<uint32_t> new_excess{0};
        uint32_t label;
        uint32_t new_label;
        std::atomic_flag discovered = ATOMIC_FLAG_INIT;
    };

    std::vector<std::vector<edge>> residual;
    std::unique_ptr<vertex[]> _vertices;
    std::unique_ptr<uint32_t[]> _active{};
    data_structures::thread_local_buffer_pool<uint32_t> _pool;
    uint32_t _source, _sink, _relabel_threshold, nactive;
    std::size_t _relabel_progress;
    const uint32_t nthreads;

  public:
    max_flow_instance(
        std::vector<std::vector<edge>> graph, uint32_t source, uint32_t sink,
        std::size_t thread_count = static_cast<size_t>(omp_get_max_threads()))
        : residual(std::move(graph)),
          _vertices(std::make_unique<vertex[]>(residual.size())),
          _active(std::make_unique<uint32_t[]>(residual.size())),
          _pool(data_structures::thread_local_buffer_pool<uint32_t>(
              thread_count, residual.size())),
          _source(source), _sink(sink), nactive(0), _relabel_progress(0),
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
        return _vertices[_sink].new_excess + _vertices[_sink].excess;
    }

    void preflow_to_flow() {
        std::swap(_source, _sink);
        find_max_flow_inner();
        std::swap(_source, _sink);
#ifdef DEBUG
        for (std::size_t i = 0; i < residual.size(); ++i)
            if (i != _source && i != _sink)
                if (_vertices[i].excess > 0)
                    std::cerr << "Excess violation: vertex " << i << ", excess "
                              << _vertices[i].excess << '\n';
#endif
    }

    auto steal_network() { return std::move(residual); }

  private:
    static constexpr uint32_t ALPHA = 6, BETA = 12;
    static constexpr double GLOBAL_RELABEL_FREQ = 0.5;

    void init() noexcept {
#pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < residual[_source].size(); ++i) {
            auto &edge = residual[_source][i];
            auto &rev_edge = residual[edge.dst][edge.rev_index];

            _vertices[edge.dst].excess = edge.cap;
            edge.rev_cap += edge.cap;
            rev_edge.cap += edge.cap;
            rev_edge.rev_cap -= edge.cap;
            edge.cap = 0;
        }

        uint32_t m = 0;
        for (std::size_t i = 0; i < residual.size(); ++i)
            m += residual[i].size();
        _relabel_threshold = residual.size() * ALPHA + m / 2;
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
                    auto vid = _active[i];
                    auto &vertex = _vertices[vid];
                    if (vertex.label == residual.size())
                        continue;
                    push(vid, vertex.label, thr_id, push_per_phase);
                }
// stage 2
#pragma omp for schedule(static) reduction(+ : _relabel_progress)
                for (uint32_t i = 0; i < nactive; ++i) {
                    auto thr_id = omp_get_thread_num();
                    auto vertex = _active[i];
                    relabel(vertex, thr_id, _relabel_progress);
                }
// stage 3
#pragma omp for schedule(static)
                for (uint32_t i = 0; i < nactive; ++i) {
                    auto &vertex = _vertices[_active[i]];
                    vertex.label = vertex.new_label;
                    vertex.discovered.clear(std::memory_order_relaxed);
                }
// stage 4
#pragma omp single
                nactive = _pool.swap_data(_active);

#pragma omp for schedule(static)
                for (uint32_t i = 0; i < nactive; ++i) {
                    auto &vertex = _vertices[_active[i]];
                    vertex.excess +=
                        vertex.new_excess.load(std::memory_order_relaxed);
                    vertex.new_excess.store(0, std::memory_order_relaxed);
                    vertex.discovered.clear(std::memory_order_relaxed);
                }
            }

            if (_relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold ||
                push_per_phase == 0) {
                _relabel_progress = 0;
                global_relabel();
            }

            npush += push_per_phase;
        }
    }

    inline void push(const uint32_t vertex, const uint32_t label, int thr_id,
                     uint64_t &push_cnt) noexcept {
        const auto target_label = label - 1;
        for (auto &edge : residual[vertex]) {
            if (edge.cap > 0 && _vertices[edge.dst].label == target_label) {
                auto flow = std::min(_vertices[vertex].excess, edge.cap);
                if (edge.dst != _source && edge.dst != _sink)
                    if (!_vertices[edge.dst].discovered.test_and_set(
                            std::memory_order_relaxed))
                        _pool.push_back(edge.dst, static_cast<size_t>(thr_id));
                ++push_cnt;
                _vertices[vertex].excess -= flow;
                _vertices[edge.dst].new_excess.fetch_add(
                    flow, std::memory_order_relaxed);
                edge.cap -= flow;
                edge.rev_cap += flow;
                residual[edge.dst][edge.rev_index].rev_cap -= flow;
                residual[edge.dst][edge.rev_index].cap += flow;
                if (_vertices[vertex].excess == 0)
                    return;
            }
        }
    }

    inline void relabel(const uint32_t vertex, const int thr_id,
                        std::size_t &relabel_progress) noexcept {
        if (_vertices[vertex].excess > 0 ||
            _vertices[vertex].label == residual.size()) {
            relabel_progress += BETA;
            _vertices[vertex].new_label = calculate_new_label(vertex);
            relabel_progress += residual[vertex].size();
            if (_vertices[vertex].new_label == residual.size()) {
                _vertices[vertex].excess += _vertices[vertex].new_excess;
                _vertices[vertex].new_excess = 0;
                return;
            }

            if (!_vertices[vertex].discovered.test_and_set(
                    std::memory_order_relaxed))
                _pool.push_back(vertex, static_cast<size_t>(thr_id));
        } else
            _vertices[vertex].new_label = _vertices[vertex].label;
    }

    inline uint32_t calculate_new_label(const uint32_t vertex) noexcept {
        uint32_t increase_to = residual.size() - 1;
        for (auto &edge : residual[vertex]) {
            if (edge.cap == 0)
                continue;

            increase_to = std::min(increase_to, _vertices[edge.dst].label);
        }
        return increase_to + 1;
    }

    void global_relabel() noexcept {
        ++nupdate;
        const auto not_reached = residual.size();

#pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < residual.size(); ++i)
            _vertices[i].label = not_reached;

        _vertices[_sink].label = 0;
        _vertices[_sink].discovered.test_and_set();
        assert(_pool.empty());
        _active[0] = _sink;
        std::size_t current_queue_size = 1;
        uint32_t current_distance = 0;

        while (current_queue_size > 0) {
#pragma omp parallel for schedule(static)
            for (std::size_t i = 0; i < current_queue_size; ++i) {
                auto thr_id = omp_get_thread_num();
                auto current_vertex = _active[i];

                for (auto edge : residual[current_vertex]) {
                    if (edge.rev_cap > 0) {
                        if (!_vertices[edge.dst].discovered.test_and_set(
                                std::memory_order_relaxed)) {
                            _vertices[edge.dst].label = current_distance + 1;
                            _pool.push_back(edge.dst,
                                            static_cast<std::size_t>(thr_id));
                        }
                    }
                }
            }
            current_queue_size = _pool.swap_data(_active);
            ++current_distance;
        }

#pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < residual.size(); ++i) {
            auto thr_id = omp_get_thread_num();
            if (_vertices[i].label != not_reached && _vertices[i].excess > 0 &&
                i != _sink)
                _pool.push_back(i, static_cast<size_t>(thr_id));
            _vertices[i].discovered.clear(std::memory_order_relaxed);
        }

        nactive = _pool.swap_data(_active);
    }
};
} // namespace parallel_push_relabel

#endif // MAXFLOW_PARALLEL_PUSH_RELABEL_H
