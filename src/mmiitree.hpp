#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <functional>
#include <thread>
#include <mio/mmap.hpp>
#include "ips4o.hpp"
#include "atomic_queue.h"

/* Suppose there are N=2^(K+1)-1 sorted numbers in an array a[]. They
 * implicitly form a complete binary tree of height K+1. We consider leaves to
 * be at level 0. The binary tree has the following properties:
 *
 * 1. The lowest k-1 bits of nodes at level k are all 1. The k-th bit is 0.
 *    The first node at level k is indexed by 2^k-1. The root of the tree is
 *    indexed by 2^K-1.
 *
 * 2. For a node x at level k, its left child is x-2^(k-1) and the right child
 *    is x+2^(k-1).
 *
 * 3. For a node x at level k, it is a left child if its (k+1)-th bit is 0. Its
 *    parent node is x+2^k. Similarly, if the (k+1)-th bit is 1, x is a right
 *    child and its parent is x-2^k.
 *
 * 4. For a node x at level k, there are 2^(k+1)-1 nodes in the subtree
 *    descending from x, including x. The left-most leaf is x&~(2^k-1) (masking
 *    the lowest k bits to 0).
 *
 * When numbers can't fill a complete binary tree, the parent of a node may not
 * be present in the array. The implementation here still mimics a complete
 * tree, though getting the special casing right is a little complex. There may
 * be alternative solutions.
 *
 * As a sorted array can be considered as a binary search tree, we can
 * implement an interval tree on top of the idea. We only need to record, for
 * each node, the maximum value in the subtree descending from the node.
 *
 * This implementation allows the interval array to be stored in a memory
 * mapped file on disk. I've got a real lot of intervals to process. -EG
 */
namespace mmmulti {

template<typename S, typename T> // "S" is a scalar type; "T" is the type of data associated with each interval
class iitree {    

    struct StackCell {
		size_t x; // node
		int k, w; // k: level; w: 0 if left child hasn't been processed
		StackCell() {};
		StackCell(int k_, size_t x_, int w_) : x(x_), k(k_), w(w_) {};
	};

public:
	struct Interval {
		S st, en, max;
		T data;
	};
    // note that we have to set max to end initially
    Interval make_interval(const S &s, const S &e, const T &d) {
        return {s, e, e, d};
    }

private:
	struct IntervalLess {
		bool operator()(const Interval &a, const Interval &b) const { return a.st < b.st || (a.st == b.st && a.en < b.en); }
	};

    mio::mmap_source reader;
    std::string filename;
    std::string index_filename;
    bool sorted = false;
    // key information
    uint64_t n_records = 0;
    bool indexed = false;
    std::thread writer_thread;
    atomic_queue::AtomicQueue2<Interval, 2 << 16> interval_queue;
    std::atomic<bool> work_todo;
	uint64_t max_level = 0;

	uint64_t index_core(size_t a_size) {
        std::error_code error;
        mio::mmap_sink buf = mio::make_mmap_sink(
            filename.c_str(), 0, mio::map_entire_file, error);
        if (error) { assert(false); }
        Interval* a = (Interval*)buf.begin();
		size_t i, last_i; // last_i points to the rightmost node in the tree
		S last; // last is the max value at node last_i
		int64_t k;
		if (a_size == 0) return -1;
		for (i = 0; i < a_size; i += 2) last_i = i, last = a[i].max = a[i].en; // leaves (i.e. at level 0)
		for (k = 1; ((int64_t)1)<<k <= a_size; ++k) { // process internal nodes in the bottom-up order
			size_t x = ((int64_t)1)<<(k-1), i0 = (x<<1) - 1, step = x<<2; // i0 is the first node
			for (i = i0; i < a_size; i += step) { // traverse all nodes at level k
				S el = a[i - x].max;                          // max value of the left child
				S er = i + x < a_size? a[i + x].max : last; // of the right child
				S e = a[i].en;
				e = e > el? e : el;
				e = e > er? e : er;
				a[i].max = e; // set the max value for node i
			}
			last_i = last_i>>k&1? last_i - x : last_i + x; // last_i now points to the parent of the original last_i
			if (last_i < a_size && a[last_i].max > last) // update last accordingly
				last = a[last_i].max;
		}
		return k - 1;
	}
    
public:

    // forward declaration for iterator types
    class iterator;
    class const_iterator;
    
    iitree(const std::string& f) : filename(f) {
        work_todo.store(false);
    }

    ~iitree(void) {
        close_writer();
        close_reader();
    }

    // make sure we don't try anything silly
    iitree(void) = delete;
    iitree(const iitree& m) = delete;
    iitree(iitree&& m) = delete;
    iitree& operator=(iitree&& m) = delete;
    
    void set_base_filename(const std::string& f) {
        filename = f;
    }

    void writer_func(void) {
        std::ofstream writer(filename.c_str(), std::ios::binary | std::ios::trunc);
        if (writer.fail()) {
            throw std::ios_base::failure(std::strerror(errno));
        }
        Interval ival;
        while (work_todo.load() || !interval_queue.was_empty()) {
            if (interval_queue.try_pop(ival)) {
                do {
                    writer.write((char*)&ival, sizeof(Interval));
                } while (interval_queue.try_pop(ival));
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
            writer_thread = std::thread(&iitree::writer_func, this);
        }
    }

    void close_writer(void) {
        if (work_todo.load()) {
            work_todo.store(false);
            while (!interval_queue.was_empty()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            if (writer_thread.joinable()) {
                writer_thread.join();
            }
        }
    }

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

    /// return the number of records, which will only work after indexing
    size_t size(void) const {
        return n_records;
    }
    
    /// iterator to first value
    iterator begin(void) {
        return iterator((Interval*)reader.begin());
    }
    
    /// iterator to one past end
    iterator end(void) {
        return iterator((Interval*)reader.end());
    }

    /// const iterator to first value
    const_iterator begin(void) const {
        return const_iterator((Interval*)reader.begin());
    }

    /// const iterator to one past end
    const_iterator end(void) const {
        return const_iterator((Interval*)reader.end());
    }

    /// get the record count
    size_t record_count(void) {
        return reader.size() / sizeof(Interval);
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
        ips4o::parallel::sort((Interval*)buffer.begin(),
                              (Interval*)buffer.end(),
                              IntervalLess(),
                              num_threads);
        sorted = true;
    }

    Interval read_value(size_t i) const {
        Interval v;
        memcpy(&v, &reader[i*sizeof(Interval)], sizeof(Interval));
        return v;
    }

    void for_each_entry(const std::function<void(const Interval&)>& lambda) const {
        for (size_t i = 0; i < n_records; ++i) {
            lambda(read_value(i));
        }
    }

    void for_each_entry_count(const std::function<void(const Interval&, const uint64_t& count)>& lambda) const {
        assert(sorted);
        uint64_t curr_count = 0;
        bool first = true;
        Interval last;
        for_each_value([&](const Interval& v) {
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
    
    /// a local reimplementation of a pointer iterator
    class iterator {
    public:
        iterator(Interval* ptr) : ptr(ptr) {}
        iterator() : ptr(nullptr) {}
        iterator(const iitree<S, T>::iterator& other) : ptr(other.ptr) {}
        iterator& operator=(const iitree<S, T>::iterator& other) {
            ptr = other.ptr;
        }
        
        Interval& operator*() {
            return *ptr;
        }
        
        bool operator==(const iitree<S, T>::iterator& other) {
            return ptr == other.ptr;
        }
        
        bool operator!=(const iitree<S, T>::iterator& other) {
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
        Interval* ptr;
        
        friend class const_iterator;
    };
    
    /// a local reimplementation of a const pointer iterator
    class const_iterator {
    public:
        const_iterator(const Interval* ptr) : ptr(ptr) {}
        const_iterator() : ptr(nullptr) {}
        const_iterator(const iitree<S, T>::const_iterator& other) : ptr(other.ptr) {}
        const_iterator& operator=(const iitree<S, T>::const_iterator& other) {
            ptr = other.ptr;
        }
        const_iterator(const iitree<S, T>::iterator& other) : ptr(other.ptr) {}
        const_iterator& operator=(const iitree<S, T>::iterator& other) {
            ptr = other.ptr;
        }
        
        const Interval& operator*() {
            return *ptr;
        }
        
        bool operator==(const iitree<S, T>::const_iterator& other) {
            return ptr == other.ptr;
        }
        
        bool operator!=(const iitree<S, T>::const_iterator& other) {
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
        const Interval* ptr;
    };

    Interval* get_array(void) const {
        return (Interval*)reader.begin();
    }

public:
    /// write into our write buffer
    /// open_writer() must be called first to set up our buffer and writer
	void add(const S &s, const S &e, const T &d) {
        interval_queue.push(make_interval(s, e, d));
    }

	void index(int num_threads) {
        close_writer();
        sort(num_threads);
        open_reader();
        n_records = record_count();
        indexed = true;
        close_reader();
        max_level = index_core(n_records);
        open_reader();
	}

    /// get overlaps, callback takes index, start, end, data
	void overlap(const S &st, const S &en, const std::function<void(const size_t&, const S&, const S&, const T&)>& func) const {
		int64_t t = 0;
		StackCell stack[64];
        Interval* a = get_array();
        stack[t++] = StackCell(max_level, (1LL<<max_level) - 1, 0); // push the root; this is a top down traversal
		while (t) { // the following guarantees that numbers in out[] are always sorted
			StackCell z = stack[--t];
			if (z.k <= 3) { // we are in a small subtree; traverse every node in this subtree
				size_t i, i0 = z.x >> z.k << z.k, i1 = i0 + (1LL<<(z.k+1)) - 1;
				if (i1 >= n_records) i1 = n_records;
				for (i = i0; i < i1 && a[i].st < en; ++i)
					if (st < a[i].en) // if overlap, append to out[]
						func(i, a[i].st, a[i].en, a[i].data);
			} else if (z.w == 0) { // if left child not processed
				size_t y = z.x - (1LL<<(z.k-1)); // the left child of z.x; NB: y may be out of range (i.e. y>=n_records)
				stack[t++] = StackCell(z.k, z.x, 1); // re-add node z.x, but mark the left child having been processed
				if (y >= n_records || a[y].max > st) // push the left child if y is out of range or may overlap with the query
					stack[t++] = StackCell(z.k - 1, y, 0);
			} else if (z.x < n_records && a[z.x].st < en) { // need to push the right child
				if (st < a[z.x].en) func(z.x, a[z.x].st, a[z.x].en, a[z.x].data); // test if z.x overlaps the query; if yes, append to out[]
				stack[t++] = StackCell(z.k - 1, z.x + (1LL<<(z.k-1)), 0); // push the right child
			}
		}
	}

    /// callback takes only start, end, and data
    void overlap(const S &st, const S &en, const std::function<void(const S&, const S&, const T&)>& func) const {
        overlap(
            st, en,
            [&func](const size_t& idx, const S& start, const S& end, const T& data) {
                func(start, end, data);
            });
    }

	//size_t size(void) const { return a.size(); }
	const S &start(size_t i) const { return get_array()[i].st; }
	const S &end(size_t i) const { return get_array()[i].en; }
	const T &data(size_t i) const { return get_array()[i].data; }
};

}
