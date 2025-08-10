#ifndef MEM_ALLOC_UTILS_H
#define MEM_ALLOC_UTILS_H
#include <cstddef>
#include <mm_malloc.h>
#include <iostream>
#include <cstring>
// decrease page faults and assert memory aligned
void alloc_aligned(void** mem_ptr, size_t alignment, size_t size);
void free_aligned(void* ptr);
struct aligned_buff_t {
  // 成员声明
  void* buff;
  size_t buff_size;
  aligned_buff_t() {
    buff = nullptr;
    buff_size = 1024;
  };
  ~aligned_buff_t() {
    free_aligned(buff);
  }
  void alloc_aligned(void** mem_ptr, size_t alignment, size_t size) {
    int ret;
    if (buff == nullptr) {
      buff_size = buff_size >= size * 2 ? buff_size : size * 2;
      ::alloc_aligned(&buff, alignment, buff_size);
      memset(buff, 0, buff_size);
    }
    if (size > buff_size) {
      free_aligned(buff);
      std::cerr << "*2\n";
      buff_size *= 2; // 倍增
      buff_size = buff_size >= size ? buff_size : size;
      ::alloc_aligned(&buff, alignment, buff_size);
      memset(buff, 0, buff_size);
    }
    *mem_ptr = buff;
    // 检查分配是否成功

  }
};

#endif