//
// Created by Jan Groschaft on 28.10.18.
//

#ifndef MAXFLOW_MEASURE_H
#define MAXFLOW_MEASURE_H

#include "common_types.h"
#include <cassert>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <optional>
#include <vector>

struct measurement_result {
    uint32_t max_flow;
    std::chrono::milliseconds time_read;
    std::chrono::milliseconds time_init;
    std::chrono::milliseconds time_solve;
};

template <typename S>
auto measure_single(S &solver, measurement_result &result) {
    auto start = std::chrono::high_resolution_clock::now();
    auto res = solver.find_max_flow();
    auto end = std::chrono::high_resolution_clock::now();
    result.max_flow = res;
    result.time_solve =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
}

template <typename S>
auto measure_parallel(std::vector<std::vector<edge>> graph, uint32_t source,
                      uint32_t sink, std::size_t thread_count) {
    measurement_result measurement;
    auto start = std::chrono::high_resolution_clock::now();
    S solver(std::move(graph), source, sink, thread_count);
    auto end = std::chrono::high_resolution_clock::now();
    measurement.time_init =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    measure_single(solver, measurement);
    return measurement;
}

template <typename S>
auto measure_sequential(std::vector<std::vector<edge>> graph, uint32_t source,
                        uint32_t sink) {
    measurement_result measurement;
    auto start = std::chrono::high_resolution_clock::now();
    S solver(std::move(graph), source, sink);
    auto end = std::chrono::high_resolution_clock::now();
    measurement.time_init =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    measure_single(solver, measurement);
    return measurement;
}

#endif // MAXFLOW_MEASURE_H
