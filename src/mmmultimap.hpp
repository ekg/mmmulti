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

    // an entry is a POD struct of Key, Value
    typedef struct { Key key; Value value; } Entry;

    // the comparator used to sort the backing array
    struct EntryLess {
        bool operator()(const std::pair<Key,Value>& a, const std::pair<Key,Value>& b) const { return a < b; }
	};

    // memory mapped buffer
    mio::mmap_source reader;
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
    // a single writer thread reads from an atomic queue and writes to our (unsorted) backing file
    std::thread writer_thread;
    atomic_queue::AtomicQueue2<Entry, 2 << 16> entry_queue;
    std::atomic<bool> work_todo;

    void init(Value nullv) {
        record_size = sizeof(Key) + sizeof(Value);
        nullkey = 0;
        nullvalue = nullv;
        work_todo.store(false);
    }

public:
    
    // forward declaration for iterator types
    class iterator;
    class const_iterator;

    // make sure we and our friend the compiler don't try anything silly
    map(void) = delete;
    map(const map& m) = delete;
    map(map&& m) = delete;
    map& operator=(map&& m) = delete;

    // constructor
    map(Value nullv) { init(nullv); }

    map(const std::string& f, Value nullv) : filename(f) { init(nullv); }

    ~map(void) {
        close_writer();
        close_reader();
    }

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

    void writer_func(void) {
        assert(!filename.empty());
        // remove the file; we only call this when making an index, and it's done once
        std::ofstream writer(filename.c_str(), std::ios::binary | std::ios::trunc);
        if (writer.fail()) {
            throw std::ios_base::failure(std::strerror(errno));
        }
        Entry entry;
        while (work_todo.load() || !entry_queue.was_empty()) {
            if (entry_queue.try_pop(entry)) {
                do {
                    writer.write((char*)&entry, sizeof(Entry));
                } while (entry_queue.try_pop(entry));
            } else {
                std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            }
        }
        writer.close();
    }

    // close/open backing file
    void open_writer(void) {
        if (!work_todo.load()) {
            work_todo.store(true);
            writer_thread = std::thread(&map::writer_func, this);
        }
    }

    void close_writer(void) {
        if (work_todo.load()) {
            work_todo.store(false);
            while (!entry_queue.was_empty()) {
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

    /// write the pair to the backing file
    /// open_writer() must be called first to set up our buffer and writer
    void append(const Key& k, const Value& v) {
        entry_queue.push((Entry){k, v});
    }

    /// return the number of records, which will only work after indexing
    size_t size(void) const {
        return n_records;
    }

    /// return the size of each combined record
    size_t get_record_size(void) const {
        return record_size;
    }

    /// get the record count
    size_t record_count(void) {
        return reader.size() / sizeof(Entry);
    }

    std::ifstream::pos_type filesize(const char* filename) {
        std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
        return in.tellg();
    }

    /// sort the record in the backing file by key
    void sort(int num_threads) {
        if (sorted) return;
        //std::cerr << "sorting!" << std::endl;
        std::error_code error;
        mio::mmap_sink buffer = mio::make_mmap_sink(
            filename.c_str(), 0, mio::map_entire_file, error);
        if (error) { assert(false); }
        // sort in parallel (uses OpenMP if available, std::thread otherwise)
        ips4o::parallel::sort((std::pair<Key, Value>*)buffer.begin(),
                              (std::pair<Key, Value>*)buffer.end(),
                              EntryLess(),
                              num_threads);
        sorted = true;
    }

    Entry read_entry(size_t i) const {
        Entry e;
        memcpy(&e, &reader[i*record_size], sizeof(Entry));
        return e;
    }

    // pad our key space with empty records so that we can query it directly with select operations
    void padsort(int num_threads) {
        close_reader();
        // blindly fill with a single key/value pair for each entity in the key space
        // running this in parallel causes a strange race condition and segfaults at high thread counts
        // and it does not provide any performance benefit
        for (size_t i = 1; i <= max_key; ++i) {
            append(i, nullvalue);
        }
        close_writer();
        sort(num_threads);
        padded = true;
    }

    void simplesort(int num_threads) {
        close_reader();
        close_writer();
        sort(num_threads);
        padded = false;
    }

    // index
    void index(int num_threads, Key new_max = 0) {
        if (new_max) {
            max_key = new_max;
            padsort(num_threads);
        } else {
            simplesort(num_threads);
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
