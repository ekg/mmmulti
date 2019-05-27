#include <iostream>
#include <random>
#include <limits>
#include <cassert>
#include <chrono>
#include "ips4o.hpp"
#include "args.hxx"
#include <endian.h>


#include <fcntl.h>
#include <cstdio>
#include <cstdlib>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>


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

//using namespace mmmultimap;

int main(int argc, char** argv) {

    args::ArgumentParser parser("memmapped multimap interface");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> in_file(parser, "FILE", "use this input file for a uint64_t sort", {'i', "in"});
    //args::Flag test_sort(parser, "", "test ips4o in memory", {'t', "test"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    }
    
    if (!args::get(in_file).empty()) {
        mmap_buffer_t buffer;
        open_mmap_buffer(args::get(in_file).c_str(), &buffer);
        std::vector<uint64_t>::iterator begin_ptr((uint64_t*)buffer.data);
        uint64_t data_len = buffer.size/sizeof(uint64_t);
        std::vector<uint64_t>::iterator end_ptr((uint64_t*)buffer.data+data_len);
        std::cerr << "data length " << data_len << std::endl;

        std::cerr << "starting sort" << std::endl;
        auto start = std::chrono::system_clock::now();
        // sort in parallel (uses OpenMP if available, std::thread otherwise)
        ips4o::parallel::sort(begin_ptr, end_ptr);
        // sort sequentially
        //ips4o::sort(x.begin(), x.end());
        auto end = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds = end-start;
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        std::cerr << "completed in " << elapsed_seconds.count() << "s" << std::endl;

        // check the sort
        /*
        std::cerr << "checking sort" << std::endl;
        uint64_t* x = (uint64_t*)buffer.data;
        for (int n=1; n<data_len; ++n) {
            assert(x[n-1] <= x[n]);
        }
        std::cerr << "ok" << std::endl;
        */

        close_mmap_buffer(&buffer);
    }

    return 0;

}
