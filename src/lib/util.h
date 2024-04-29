#ifndef UTIL_H
#define UTIL_H

#include <cstdint>
#include <limits>
#include <utility>
#include <mutex>
#include <memory>
#include <cstring>
#include <vector>

#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif

struct edge {
    edge() = default;

    edge(uint32_t dst, uint32_t capacity, uint32_t rev_index)
        : dst(dst), rev_index(rev_index), cap(capacity) {}

    uint32_t dst;
    uint32_t rev_index;
    uint32_t cap;
    uint32_t rev_cap;
};

/*
 * Thread safe buffer, accesses are protected by a mutex.
 */
template <typename T> class thread_safe_buffer {
    std::unique_ptr<T[]> _buffer;
    std::size_t _pos;
    std::mutex mutex;

  public:
    explicit thread_safe_buffer(std::size_t size)
        : _buffer(std::make_unique<T[]>(size)), _pos(0) {}

    void clear() noexcept { 
        std::lock_guard guard(mutex);
        _pos = 0; 
    }

    bool empty() noexcept { return _pos == 0; }

    void append(const T *const data, std::size_t size) noexcept {
        std::lock_guard guard(mutex);
        std::memcpy(_buffer.get() + _pos, data, size * sizeof(T));
        _pos += size;
    }

    void push_back(const T &val) noexcept {
        std::lock_guard guard(mutex);
        _buffer[_pos] = val;
        _pos++;
    }

    std::size_t size() noexcept { return _pos; }

    void swap_data(std::unique_ptr<T[]> &other) noexcept {
        std::lock_guard guard(mutex);
        _buffer.swap(other);
        _pos = 0;
    }
};

/*
 * Buffer which provides concurrent access to its underlying container, without
 * the need of synchronization for every insert. Each thread stores data to its
 * own, smaller buffer, which is later copied to the shared thread_safe_buffer.
 */
template <typename T> class buffer_pool {
    // this actually has a big perf impact
    struct alignas(CACHE_LINE_SIZE) aligned_vector {
        std::vector<T> data;
    };

    thread_safe_buffer<T> _buffer;
    std::unique_ptr<aligned_vector[]> pool;
    const std::size_t nthreads;

  public:
    buffer_pool(std::size_t thread_count,
                             std::size_t total_max_size)
        : _buffer(total_max_size),
          pool(std::make_unique<aligned_vector[]>(thread_count)),
          nthreads(thread_count) {
        for (std::size_t i = 0; i < nthreads; ++i)
            pool[i].data.reserve(
                std::max(1ul << 4, total_max_size / nthreads));
    }

    void push_back(const T &value, std::size_t thread_id) noexcept {
        auto &target_buffer = pool[thread_id].data;
        if (target_buffer.size() == target_buffer.capacity()) {
            _buffer.append(target_buffer.data(), target_buffer.size());
            target_buffer.clear();
        }
        pool[thread_id].data.push_back(value);
    }

    bool empty() const noexcept {
        if (!_buffer.empty())
            return false;
        for (std::size_t i = 0; i < nthreads; ++i)
            if (!pool[i].data.empty())
                return false;
        return true;
    }

    auto swap_data(std::unique_ptr<T[]> &other) noexcept {
#pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < nthreads; ++i) {
            auto &target_buffer = pool[i].data;
            if (target_buffer.size() > 0)
                _buffer.append(target_buffer.data(), target_buffer.size());
            target_buffer.clear();
        }
        auto prev_size = _buffer.size();
        _buffer.swap_data(other);
        return prev_size;
    }
};

#endif // UTIL_H
