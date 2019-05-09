#include "buffer.hpp"

namespace mmmultimap {

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

}
