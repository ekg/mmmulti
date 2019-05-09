#pragma once

#include <fcntl.h>
#include <cstdio>
#include <cstdlib>
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
//#include <time.h>
#include <unistd.h>
//#include <omp.h>

namespace mmmultimap {

struct mmap_buffer_t {
  int fd;
  off_t size;
  void *data;
};

int open_mmap_buffer(const char* path, mmap_buffer_t* buffer);
void close_mmap_buffer(mmap_buffer_t* buffer);

}
