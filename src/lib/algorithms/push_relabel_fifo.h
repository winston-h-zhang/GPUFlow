/*
 * Push-relabel, FIFO active vertex selection.
 */

#ifndef MAXFLOW_PUSH_RELABEL_FIFO_H
#define MAXFLOW_PUSH_RELABEL_FIFO_H

#include "../common_types.h"
#include <chrono>
#include <iostream>
#include <memory>
#include <queue>

namespace push_relabel_fifo {

class max_flow_instance {
    struct vertex {
        uint32_t excess{0};
        uint32_t label;
    };
    using pair = std::pair<uint32_t, uint32_t>;
    std::vector<std::vector<edge>> residual;
    std::unique_ptr<vertex[]> _vertices;
    std::queue<uint32_t> _q;
    std::queue<pair> _distance_q;
    uint32_t _source, _sink, _relabel_progress{0}, _relabel_threshold;

    // statistics
    uint64_t npush{0}, nrelabel{0}, nglobal{0};

  public:
    max_flow_instance(std::vector<std::vector<edge>> graph, uint32_t source,
                      uint32_t sink)
        : residual(std::move(graph)),
          _vertices(std::make_unique<vertex[]>(residual.size())),
          _q(std::queue<uint32_t>{}), _distance_q(std::queue<pair>{}),
          _source(source), _sink(sink) {
        init();
    }

    uint32_t find_max_flow() {
        find_max_flow_inner();
#ifdef DEBUG
        std::cout << "pushes:\t\t" << npush << std::endl;
        std::cout << "relabels:\t" << nrelabel << std::endl;
        std::cout << "global updates:\t" << nglobal << std::endl;
#endif
        return _vertices[_sink].excess;
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

    void init() {
        for (auto &edge : residual[_source]) {
            if (edge.cap > 0) {
                _vertices[edge.dst].excess = edge.cap;
                edge.rev_cap += edge.cap;
                residual[edge.dst][edge.rev_index].cap += edge.cap;
                residual[edge.dst][edge.rev_index].rev_cap -= edge.cap;
                edge.cap = 0;
                ++npush;
            }
        }

        uint32_t m = 0;
        for (auto &vec : residual)
            m += vec.size();
        _relabel_threshold = residual.size() * ALPHA + m / 2;
    }

    void find_max_flow_inner() {
        global_relabel();
        for (;;) {
            if (_q.empty())
                return;

            auto vertex = _q.front();
            _q.pop();
            auto label = _vertices[vertex].label;
            discharge(vertex, label);

            if (_relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold) {
                _relabel_progress = 0;
                global_relabel();
            }
        }
    }

    void discharge(const uint32_t vertex, uint32_t label) {
        while (label < residual.size()) {
            if (push(vertex, label))
                return;
            label = relabel(vertex);
        }
    }

    bool push(const uint32_t vertex, const uint32_t label) {
        const auto target_label = label - 1;
        for (auto &edge : residual[vertex]) {
            if (edge.cap > 0 && _vertices[edge.dst].label == target_label) {
                ++npush;
                auto flow = std::min(_vertices[vertex].excess, edge.cap);
                if (_vertices[edge.dst].excess == 0 && edge.dst != _sink)
                    _q.push(edge.dst);
                _vertices[vertex].excess -= flow;
                _vertices[edge.dst].excess += flow;
                edge.cap -= flow;
                edge.rev_cap += flow;
                residual[edge.dst][edge.rev_index].rev_cap -= flow;
                residual[edge.dst][edge.rev_index].cap += flow;
                if (_vertices[vertex].excess == 0)
                    return true;
            }
        }
        return false;
    }

    uint32_t relabel(const uint32_t vertex) {
        ++nrelabel;
        _relabel_progress += BETA;
        _vertices[vertex].label = calculate_new_label(vertex);
        return _vertices[vertex].label;
    }

    uint32_t calculate_new_label(const uint32_t vertex) {
        uint32_t increase_to = residual.size() - 1;
        for (auto &edge : residual[vertex]) {
            if (edge.cap == 0)
                continue;
            increase_to = std::min(increase_to, _vertices[edge.dst].label);
        }
        _relabel_progress += residual[vertex].size();
        return increase_to + 1;
    }

    void global_relabel() {
        ++nglobal;
        auto not_reached = residual.size();

        for (std::size_t i = 0; i < residual.size(); ++i)
            _vertices[i].label = not_reached;

        _q = std::queue<uint32_t>();
        _distance_q = std::queue<pair>();
        _distance_q.push(std::make_pair(_sink, 0));
        _vertices[_sink].label = 0;

        while (!_distance_q.empty()) {
            auto current_elem = _distance_q.front();
            _distance_q.pop();

            auto current_vertex = current_elem.first;
            auto current_distance = current_elem.second;
            for (auto &edge : residual[current_vertex]) {
                if (edge.rev_cap > 0 &&
                    _vertices[edge.dst].label == not_reached) {
                    _vertices[edge.dst].label = current_distance + 1;
                    _distance_q.push(
                        std::make_pair(edge.dst, current_distance + 1));
                    if (_vertices[edge.dst].excess > 0)
                        _q.push(edge.dst);
                }
            }
        }
    }
};
} // namespace push_relabel_fifo

#endif // MAXFLOW_PUSH_RELABEL_FIFO_H
