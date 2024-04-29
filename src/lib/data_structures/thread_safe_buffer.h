/*
 * Thread safe array, data is retrieved using the swap method.
 */

#ifndef MAXFLOW_THREAD_SAFE_BUFFER_H
#define MAXFLOW_THREAD_SAFE_BUFFER_H

#include <cstring>
#include <memory>
#include <mutex>

namespace data_structures {
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
} // namespace data_structures

#endif // MAXFLOW_THREAD_SAFE_BUFFER_H
