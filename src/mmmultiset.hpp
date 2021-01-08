#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <functional>
#include <unordered_set>
#include <utility>
#include <thread>

#include <fcntl.h>
#include <cstdio>
#include <cstdlib>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <mio/mmap.hpp>

#include "sdsl/bit_vectors.hpp"
#include "ips4o.hpp"
#include "atomic_queue.h"

namespace mmmulti {

/*
'mmmulti::set' is a disk-backed multiset values are stored
in a binary file. The key space is assumed to be numeric, but values
may be of arbitrary size.  To build the multiset we first append
values. We then sort the values. We can now iterate over the unique values
of the multiset, such as to obtain their counts.
*/

template <typename Value> class set {

    struct ValueLess {
        bool operator()(const Value &a, const Value &b) const { return a < b; }
	};


    mio::mmap_source reader;
    std::string filename;
    std::string index_filename;
    bool sorted = false;
    // key information
    uint64_t n_records = 0;
    bool indexed = false;
    std::thread writer_thread;
    atomic_queue::AtomicQueue2<Value, 2 << 16> value_queue;
    std::atomic<bool> work_todo;
    
public:

    // forward declaration for iterator types
    class iterator;
    class const_iterator;

    // make sure we and our friend the compiler don't try anything silly
    set(void) = delete;
    set(const set& m) = delete;
    set(set&& m) = delete;
    set& operator=(set&& m) = delete;

    set(const std::string& f) : filename(f) {
        work_todo.store(false);
    }

    ~set(void) {
        close_writer();
        close_reader();
    }

    void set_base_filename(const std::string& f) {
        filename = f;
    }

    void writer_func(void) {
        std::ofstream writer(filename.c_str(), std::ios::binary | std::ios::trunc);
        if (writer.fail()) {
            throw std::ios_base::failure(std::strerror(errno));
        }
        Value value;
        while (work_todo.load() || !value_queue.was_empty()) {
            if (value_queue.try_pop(value)) {
                do {
                    writer.write((char*)&value, sizeof(Value));
                } while (value_queue.try_pop(value));
            } else {
                std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            }
        }
        writer.close();
    }

    // open backing file and start writer_thread
    void open_writer(void) {
        if (!work_todo.load()) {
            work_todo.store(true);
            writer_thread = std::thread(&set::writer_func, this);
        }
    }

    void close_writer(void) {
        if (work_todo.load()) {
            work_todo.store(false);
            while (!value_queue.was_empty()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            if (writer_thread.joinable()) {
                writer_thread.join();
            }
        }
    }

    // è mio
    void open_reader(void) {
        if (reader.is_mapped()) return;
        std::error_code error;
        reader = mio::make_mmap_source(
            filename.c_str(), 0, mio::map_entire_file, error);
        if (error) { assert(false); }
        // set up for sequential access
        madvise((void*)reader.begin(),
                reader.mapped_length(),
                POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);
    }

    void close_reader(void) {
        if (reader.is_mapped()) reader.unmap();
    }

    /// write the pair to end of backing file
    /// open_writer() must be called first to set up our buffer and writer
    void append(const Value& v) {
        value_queue.push(v);
    }

    /// return the number of records, which will only work after indexing
    size_t size(void) const {
        return n_records;
    }
    
    /// iterator to first value
    iterator begin(void) {
        return iterator((Value*)reader.begin());
    }
    
    /// iterator to one past end
    iterator end(void) {
        return iterator((Value*)reader.end());
    }

    /// const iterator to first value
    const_iterator begin(void) const {
        return const_iterator((Value*)reader.begin());
    }

    /// const iterator to one past end
    const_iterator end(void) const {
        return const_iterator((Value*)reader.end());
    }

    /// return the size of each combined record
    size_t get_record_size(void) const {
        return sizeof(Value);
    }

    /// get the record count
    size_t record_count(void) {
        return reader.size() / sizeof(Value);
    }

    std::ifstream::pos_type filesize(const char* filename) {
        std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
        return in.tellg();
    }

    void sort(int num_threads) {
        if (sorted) return;
        //std::cerr << "sorting!" << std::endl;
        std::error_code error;
        mio::mmap_sink buffer = mio::make_mmap_sink(
            filename.c_str(), 0, mio::map_entire_file, error);
        if (error) { assert(false); }
        // sort in parallel (uses OpenMP if available, std::thread otherwise)
        ips4o::parallel::sort((Value*)buffer.begin(),
                              (Value*)buffer.end(),
                              ValueLess(),
                              num_threads);
        sorted = true;
    }
    
    Value read_value(size_t i) const {
        Value v;
        memcpy(&v, &reader[i*sizeof(Value)], sizeof(Value));
        return v;
    }

    // index
    void index(int num_threads) {
        close_writer();
        sort(num_threads);
        open_reader();
        n_records = record_count();
        indexed = true;
        close_reader();
        open_reader();
    }



    // todo we make a big change, pulling in all the stuff from mmmultimap
    
    
    
    void for_each_value(const std::function<void(const Value&)>& lambda) const {
        for (size_t i = 0; i < n_records; ++i) {
            lambda(read_value(i));
        }
    }

    void for_each_value_count(const std::function<void(const Value&, const uint64_t& count)>& lambda) const {
        assert(sorted);
        uint64_t curr_count = 0;
        bool first = true;
        Value last;
        for_each_value([&](const Value& v) {
                if (first) {
                    last = v;
                    first = false;
                } else if (v != last) {
                    lambda(last, curr_count);
                    curr_count = 0;
                    last = v;
                }
                ++curr_count;
            });
        lambda(last, curr_count);
    }

    void for_each_unique_value(const std::function<void(const Value&)>& lambda) const {
        assert(sorted);
        bool first = true;
        Value last;
        for_each_value([&](const Value& v) {
                if (first) {
                    last = v;
                    first = false;
                } else if (v != last) {
                    lambda(last);
                    last = v;
                }
            });
        lambda(last);
    }

    /// a local reimplementation of a pointer iterator
    class iterator {
    public:
        iterator(Value* ptr) : ptr(ptr) {}
        iterator() : ptr(nullptr) {}
        iterator(const set<Value>::iterator& other) : ptr(other.ptr) {}
        iterator& operator=(const set<Value>::iterator& other) {
            ptr = other.ptr;
        }
        
        Value& operator*() {
            return *ptr;
        }
        
        bool operator==(const set<Value>::iterator& other) {
            return ptr == other.ptr;
        }
        
        bool operator!=(const set<Value>::iterator& other) {
            return !(*this == other);
        }
        
        iterator& operator++() {
            ++ptr;
            return *this;
        }
        
        iterator operator++(int) {
            return iterator(ptr++);
        }
        
        iterator& operator--() {
            --ptr;
            return *this;
        }
        
        iterator operator--(int) {
            return iterator(ptr--);
        }
        
    private:
        Value* ptr;
        
        friend class const_iterator;
    };
    
    /// a local reimplementation of a const pointer iterator
    class const_iterator {
    public:
        const_iterator(const Value* ptr) : ptr(ptr) {}
        const_iterator() : ptr(nullptr) {}
        const_iterator(const set<Value>::const_iterator& other) : ptr(other.ptr) {}
        const_iterator& operator=(const set<Value>::const_iterator& other) {
            ptr = other.ptr;
        }
        const_iterator(const set<Value>::iterator& other) : ptr(other.ptr) {}
        const_iterator& operator=(const set<Value>::iterator& other) {
            ptr = other.ptr;
        }
        
        const Value& operator*() {
            return *ptr;
        }
        
        bool operator==(const set<Value>::const_iterator& other) {
            return ptr == other.ptr;
        }
        
        bool operator!=(const set<Value>::const_iterator& other) {
            return !(*this == other);
        }
        
        const_iterator& operator++() {
            ++ptr;
            return *this;
        }
        
        const_iterator operator++(int) {
            return const_iterator(ptr++);
        }
        
        const_iterator& operator--() {
            --ptr;
            return *this;
        }
        
        const_iterator operator--(int) {
            return const_iterator(ptr--);
        }
        
    private:
        const Value* ptr;
    };
};

}
