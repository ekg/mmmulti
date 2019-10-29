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
'mmmulti::map' is a disk-backed multimap where keys and values are stored
in a binary file. The key space is assumed to be numeric, but values
may be of arbitrary size.  To build the multimap we first append
key/value pairs.  To query the multimap we must first index it.  We
first sort by key using bsort.  Then we pad the key space so that we
have one entry per integer in the range [0, max(keys)], sorting again
to put the padding pairs in the right positions.  We record the key
space by marking a bitvector of length equal to max(keys) with 1 at
those positions corresponding to the first record of each key in the
sorted array.  We compress this bitvector and build select supports on
it We are now able to traverse the sorted array using select queries
on this bitvector.
*/

template <typename Key, typename Value> class map {

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
    
    typedef struct { Key key; Value value; } Entry;
    std::ofstream writer;
    std::vector<std::ofstream> writers;
    char* reader = nullptr;
    int reader_fd = 0;
    std::string filename;
    std::string index_filename;
    bool sorted = false;
    bool padded = false;
    size_t record_size = 0;
    // key information
    Key max_key = 0;
    uint64_t n_records = 0;
    // null key and value
    Key nullkey;
    Value nullvalue;
    // compressed bitvector marked at key starts
    sdsl::sd_vector<> key_cbv;
    // select support for the key cbv
    sdsl::sd_vector<>::select_1_type key_cbv_select;
    bool indexed = false;
    uint32_t OUTPUT_VERSION = 1; // update as we change our format

    void init(void) {
        record_size = sizeof(Key) + sizeof(Value);
        nullkey = 0;
        for (size_t i = 0; i < sizeof(Value); ++i) {
            ((uint8_t*)&nullvalue)[i] = 0;
        }
    }

public:
    
    // forward declaration for iterator types
    class iterator;
    class const_iterator;

    // constructor
    map(void) { init(); }

    map(const std::string& f) : filename(f) { init(); open_writers(f); }

    ~map(void) { close_writers(); close_reader(); }

    void set_base_filename(const std::string& f) {
        filename = f;
        index_filename = filename+".idx";
    }

    // load from base file name
    void load(const std::string& f) {
        open_reader();
        set_base_filename(f);
        std::ifstream in(index_filename.c_str());
        std::string magic;
        in.read((char*)magic.c_str(), 9);
        uint32_t version;
        in.read((char*) &version, sizeof(version));
        assert(version == OUTPUT_VERSION);
        size_t record_size_in_bytes;
        sdsl::read_member(record_size_in_bytes, in);
        assert(record_size_in_bytes == record_size);
        sdsl::read_member(n_records, in);
        assert(n_records == record_count());
        sdsl::read_member(max_key, in);
        assert(max_key == nth_key(n_records));
        key_cbv.load(in);
        key_cbv_select.load(in);
    }

    // save indexes
    size_t save(sdsl::structure_tree_node* s = NULL, std::string name = "") {
        assert(max_key && indexed);
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
        // open the sdsl index
        std::ofstream out(index_filename.c_str());
        size_t written = 0;
        out << "mmmultimap"; written += 9;
        uint32_t version_buffer = OUTPUT_VERSION;
        out.write((char*) &version_buffer, sizeof(version_buffer));
        written += sdsl::write_member(record_size, out, child, "record_size");
        written += sdsl::write_member(record_count(), out, child, "record_count");
        written += sdsl::write_member(max_key, out, child, "max_key");
        written += key_cbv.serialize(out, child, "key_cbv");
        written += key_cbv_select.serialize(out, child, "key_cbv_select");
        out.close();
        return written;
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
        // check to see if we ran single-threaded
        uint64_t used_writers = 0;
        uint64_t writer_that_wrote = 0;
        for (size_t i = 0; i < writers.size(); ++i) {
            writers[i].close();
            if (filesize(writer_filename(i).c_str())) {
                ++used_writers;
                writer_that_wrote = i;
            }
        }
        bool single_threaded = used_writers == 1;
        // close the temp writers and cat them onto the end of the main file
        if (single_threaded) {
            std::rename(writer_filename(writer_that_wrote).c_str(), filename.c_str());
            for (size_t i = 0; i < writers.size(); ++i) {
                if (i != writer_that_wrote) {
                    std::remove(writer_filename(i).c_str());
                }
            }
        } else {
            open_main_writer();
            for (size_t i = 0; i < writers.size(); ++i) {
                std::ifstream if_w(writer_filename(i), std::ios_base::binary);
                writer << if_w.rdbuf();
                if_w.close();
                std::remove(writer_filename(i).c_str());
            }
        }
        writers.clear();
        writer.close();
        for (size_t i = 0; i < writers.size(); ++i) {
            std::remove(writer_filename(i).c_str());
        }
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
    void append(const Key& k, const Value& v) {
        sorted = false; // assume we break the sort
        // write to the end of the file
        auto& writer = get_writer();
        writer.write((char*)&k, sizeof(Key));
        writer.write((char*)&v, sizeof(Value));
    }

    /// return the number of records, which will only work after indexing
    size_t size(void) const {
        return n_records;
    }

    /// return the size of each combined record
    size_t get_record_size(void) const {
        return record_size;
    }

    /// return the backing buffer
    char* get_buffer(void) const {
        return reader;
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
        size_t count = stats.st_size / record_size;
        close(fd);
        return count;
    }

    std::ifstream::pos_type filesize(const char* filename) {
        std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
        return in.tellg();
    }

    /// sort the record in the backing file by key
    void sort(void) {
        if (sorted) return;
        //std::cerr << "sorting!" << std::endl;
        mmap_buffer_t buffer;
        open_mmap_buffer(filename.c_str(), &buffer);
        uint64_t data_len = buffer.size/record_size;
        // sort in parallel (uses OpenMP if available, std::thread otherwise)
        ips4o::parallel::sort((std::pair<Key, Value>*)buffer.data,
                              ((std::pair<Key, Value>*)buffer.data)+data_len);
        close_mmap_buffer(&buffer);
        sorted = true;
    }

    Entry read_entry(size_t i) const {
        Entry e;
        memcpy(&e, &reader[i*record_size], sizeof(Entry));
        return e;
    }

    // pad our key space with empty records so that we can query it directly with select operations
    void padsort(void) {
        close_reader();
        // blindly fill with a single key/value pair for each entity in the key space
        open_writers();
        // running this in parallel causes a strange race condition and segfaults at high thread counts
        // and it does not provide any performance benefit
        for (size_t i = 1; i <= max_key; ++i) {
            append(i, nullvalue);
        }
        sync_writers();
        sort();
        padded = true;
    }

    void simplesort(void) {
        close_reader();
        sync_writers();
        sort();
        padded = false;
    }

    // index
    void index(Key new_max = 0) {
        if (new_max) {
            max_key = new_max;
            padsort();
        } else {
            simplesort();
        }
        open_reader();
        n_records = record_count();
        if (padded) {
            sdsl::bit_vector key_bv(n_records+1);
            // record the key starts
            Key last = nullkey, curr = nullkey;
            Value val = nullvalue;
            Entry entry;
            //reader.read((char*)&last, sizeof(Key));
            for (size_t i = 0; i < n_records; ++i) {
                entry = read_entry(i);
                curr = entry.key;
                val = entry.value;
                if (curr != last) {
                    key_bv[i] = 1;
                }
                last = curr;
            }
            // the last key in the sort is our max key
            max_key = nth_key(n_records-1);
            key_bv[n_records] = 1; // sentinel
            // build the compressed bitvector
            sdsl::util::assign(key_cbv, sdsl::sd_vector<>(key_bv));
            key_bv.resize(0); // memory could be tight
            // build the select supports on the key bitvector
            sdsl::util::assign(key_cbv_select, sdsl::sd_vector<>::select_1_type(&key_cbv));
        }
        indexed = true;
        close_reader();
        open_reader();
    }

    Key nth_key(size_t n) const {
        Entry e = read_entry(n);
        return e.key;
    }

    Value nth_value(size_t n) const {
        Entry e = read_entry(n);
        return e.value;
    }

    void for_each_pair_parallel(const std::function<void(const Key&, const Value&)>& lambda) const {
#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < n_records; ++i) {
            Entry entry = read_entry(i);
            if (!padded || !is_null(entry.value)) {
                lambda(entry.key, entry.value);
            }
        }
    }

    void for_each_pair(const std::function<void(const Key&, const Value&)>& lambda) const {
        Entry entry;
        for (size_t i = 0; i < n_records; ++i) {
            entry = read_entry(i);
            if (!padded || !is_null(entry.value)) {
                lambda(entry.key, entry.value);
            }
        }
    }
    
    std::vector<Value> values(const Key& key) const {
        std::vector<Value> values;
        for_values_of(key, [&values](const Value& v) { values.push_back(v); });
        return values;
    }

    std::vector<Value> unique_values(const Key& key) const {
        std::vector<Value> values;
        for_unique_values_of(key, [&values](const Value& v) { values.push_back(v); });
        return values;
    }

    void for_unique_values_of(const Key& key, const std::function<void(const Value&)>& lambda) const {
        // quirk: if we've sorted by the whole binary record,
        // then we can do a simple 'uniq' operation to get the unique values
        Value last = nullvalue;
        for_values_of(key, [this,&lambda,&last](const Value& value) {
                if (value != last) {
                    lambda(value);
                    last = value;
                }
            });
    }

    bool is_null(const Value& value) const {
        for (size_t i = 0; i < sizeof(Value); ++i) {
            if (((uint8_t*)&value)[i] != 0) {
                return false;
            }
        }
        return true;
    }

    void for_values_of(const Key& key, const std::function<void(const Value&)>& lambda) const {
        if (!padded || key == 0 || key > max_key) {
            return;
        }
        size_t i = key_cbv_select(key);
        for ( ; i < n_records; ++i) {
            Entry entry = read_entry(i);
            if (entry.key != key) break;
            if (!is_null(entry.value)) {
                lambda(entry.value);
            }
        }
    }
};

}
