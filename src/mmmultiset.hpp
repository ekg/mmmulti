#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <functional>
#include <unordered_set>
#include <utility>

#include <fcntl.h>
#include <cstdio>
#include <cstdlib>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "sdsl/bit_vectors.hpp"
#include "ips4o.hpp"

namespace mmmulti {

/*
'mmmulti::set' is a disk-backed multiset values are stored
in a binary file. The key space is assumed to be numeric, but values
may be of arbitrary size.  To build the multiset we first append
values. We then sort the values. We can now iterate over the unique values
of the multiset, such as to obtain their counts.
*/

template <typename Value> class set {

private:
    
    // memory mapped buffer struct
    struct mmap_buffer_t {
        int fd;
        off_t size;
        void *data;
    };

    // utilities used by mmmultimap
    int open_mmap_buffer(const char* path, mmap_buffer_t* buffer) {
        buffer->data = nullptr;
        buffer->fd = open(path, O_RDWR);
        if (buffer->fd == -1) {
            goto error;
        }
        struct stat stats;
        if (-1 == fstat(buffer->fd, &stats)) {
            goto error;
        }
        if (!(buffer->data = mmap(nullptr,
                                  stats.st_size,
                                  PROT_READ | PROT_WRITE,
                                  MAP_SHARED,
                                  buffer->fd,
                                  0
                  ))) {
            goto error;
        }
        madvise(buffer, stats.st_size, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);
        buffer->size = stats.st_size;
        return 0;

    error:
        perror(path);
        if (buffer->data)
            munmap(buffer->data, stats.st_size);
        if (buffer->fd != -1)
            close(buffer->fd);
        buffer->data = 0;
        buffer->fd = 0;
        return -1;
    }

    void close_mmap_buffer(mmap_buffer_t* buffer) {
        if (buffer->data) {
            munmap(buffer->data, buffer->size);
            buffer->data = 0;
            buffer->size = 0;
        }

        if (buffer->fd) {
            close(buffer->fd);
            buffer->fd = 0;
        }
    }

    int get_thread_count(void) {
        int thread_count = 1;
#pragma omp parallel
        {
#pragma omp master
            thread_count = omp_get_num_threads();
        }
        return thread_count;
    }
    
    std::ofstream writer;
    std::vector<std::ofstream> writers;
    char* reader = nullptr;
    int reader_fd = 0;
    std::string filename;
    std::string index_filename;
    bool sorted = false;
    // key information
    uint64_t n_records = 0;
    bool indexed = false;
    uint32_t OUTPUT_VERSION = 1; // update as we change our format

public:

    // constructor
    set(void) { }

    set(const std::string& f) : filename(f) { open_writers(f); }

    ~set(void) { close_writers(); }

    void set_base_filename(const std::string& f) {
        filename = f;
    }

    // close/open backing file
    void open_main_writer(void) {
        if (writer.is_open()) {
            writer.seekp(0, std::ios_base::end); // seek to the end for appending
            return;
        }
        assert(!filename.empty());
        // open in binary append mode as that's how we write into the file
        //writer.open(filename.c_str(), std::ios::binary | std::ios::app);
        // remove the file; we only call this when making an index, and it's done once
        writer.open(filename.c_str(), std::ios::binary | std::ios::trunc);
        if (writer.fail()) {
            throw std::ios_base::failure(std::strerror(errno));
        }
    }

    // per-thread writers
    void open_writers(const std::string& f) {
        set_base_filename(f);
        open_writers();
    }

    void open_writers(void) {
        assert(!filename.empty());
        writers.clear();
        writers.resize(get_thread_count());
        for (size_t i = 0; i < writers.size(); ++i) {
            auto& writer = writers[i];
            writer.open(writer_filename(i), std::ios::binary | std::ios::app);
            if (writer.fail()) {
                throw std::ios_base::failure(std::strerror(errno));
            }
        }
    }

    std::string writer_filename(size_t i) {
        std::stringstream wf;
        wf << filename << ".tmp_write" << "." << i;
        return wf.str();
    }

    void open_reader(void) {
        if (reader_fd) return; //open
        assert(!filename.empty());
        // open in binary mode as we are reading from this interface
        reader_fd = open(filename.c_str(), O_RDWR);
        if (reader_fd == -1) {
            assert(false);
        }
        struct stat stats;
        if (-1 == fstat(reader_fd, &stats)) {
            assert(false);
        }
        if (!(reader =
              (char*) mmap(NULL,
                            stats.st_size,
                            PROT_READ | PROT_WRITE,
                            MAP_SHARED,
                            reader_fd,
                            0))) {
            assert(false);
        }
        madvise((void*)reader, stats.st_size, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);
    }

    std::ofstream& get_writer(void) {
        return writers[omp_get_thread_num()];
    }
    
    void sync_writers(void) {
        // close the temp writers and cat them onto the end of the main file
        open_main_writer();
        for (size_t i = 0; i < writers.size(); ++i) {
            writers[i].close();
            std::ifstream if_w(writer_filename(i), std::ios_base::binary);
            writer << if_w.rdbuf();
            if_w.close();
            std::remove(writer_filename(i).c_str());
        }
        writers.clear();
        writer.close();
    }

    void close_writers(void) {
        for (size_t i = 0; i < writers.size(); ++i) {
            std::remove(writer_filename(i).c_str());
        }
    }
    
    void close_reader(void) {
        if (reader) {
            size_t c = record_count();
            munmap(reader, c);
            reader = 0;
        }
        if (reader_fd) {
            close(reader_fd);
            reader_fd = 0;
        }
    }

    /// write the pair to end of backing file
    void append(const Value& v) {
        sorted = false; // assume we break the sort
        // write to the end of the file
        auto& writer = get_writer();
        writer.write((char*)&v, sizeof(Value));
    }

    /// return the number of records, which will only work after indexing
    size_t size(void) const {
        return n_records;
    }

    /// get the record count
    size_t record_count(void) {
        int fd = open(filename.c_str(), O_RDWR);
        if (fd == -1) {
            assert(false);
        }
        struct stat stats;
        if (-1 == fstat(fd, &stats)) {
            assert(false);
        }
        assert(stats.st_size % record_size == 0); // must be even records
        size_t count = stats.st_size / sizeof(Value);
        return count;
    }

    /// sort the record in the backing file by key
    void sort(void) {
        sync_writers();
        close_reader();
        if (sorted) return;
        //std::cerr << "sorting!" << std::endl;
        mmap_buffer_t buffer;
        open_mmap_buffer(filename.c_str(), &buffer);
        typename std::vector<Value>::iterator begin_ptr((Value*)buffer.data);
        uint64_t data_len = buffer.size/sizeof(Value);
        typename std::vector<Value>::iterator end_ptr((Value*)buffer.data+data_len);
        // sort in parallel (uses OpenMP if available, std::thread otherwise)
        ips4o::parallel::sort(begin_ptr, end_ptr);
        close_mmap_buffer(&buffer);
        sorted = true;
    }

    Value read_value(size_t i) const {
        Value v;
        memcpy(&v, &reader[i*sizeof(Value)], sizeof(Value));
        return v;
    }

    // index
    void index(void) {
        sort();
        open_reader();
        n_records = record_count();
        indexed = true;
        close_reader();
        open_reader();
    }

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
};

}
