#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <cstdint>
#include <limits>
#include <utility>

struct edge {
    edge() = default;

    edge(uint32_t dst, uint32_t capacity, uint32_t rev_index)
        : dst(dst), rev_index(rev_index), cap(capacity) {}

    uint32_t dst;
    uint32_t rev_index;
    uint32_t cap;
    uint32_t rev_cap;
};

#endif // COMMON_TYPES_H
