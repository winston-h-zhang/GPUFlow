#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <cstdint>
#include <limits>
#include <utility>

struct edge {
    edge() = default;

    edge(uint32_t dst_vertex, uint32_t capacity,
                uint32_t reverse_edge_index)
        : dst_vertex(dst_vertex), reverse_edge_index(reverse_edge_index),
          r_capacity(capacity) {}

    uint32_t dst_vertex;
    uint32_t reverse_edge_index;
    uint32_t r_capacity;
    uint32_t reverse_r_capacity;
};

#endif // COMMON_TYPES_H
