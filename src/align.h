#ifndef ALIGN_H
#define ALIGN_H
#include <vector>
#include "graph.h"
#include "parameter.h"
#include "result.h"
#include "mem_alloc_utils.h"
#include <stdio.h>
std::vector<res_t> POA(para_t* para, const graph& DAG, const std::string& seq);
std::vector<res_t> POA_SIMD(para_t* para, const graph& DAG, const std::string& _seq, aligned_buff_t* mpool = nullptr);
std::vector<res_t> POA_SIMD_ORIGIN(const para_t* para, const graph* DAG, const std::string& _seq, aligned_buff_t* mpool = nullptr);
std::vector<res_t> abPOA(const para_t* para, const graph* DAG, const std::string& _seq, aligned_buff_t* mpool = nullptr);
#endif