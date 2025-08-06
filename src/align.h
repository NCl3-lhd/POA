#ifndef ALIGN_H
#define ALIGN_H
#include <vector>
#include "graph.h"
#include "parameter.h"
#include "result.h"
std::vector<res_t> POA(para_t* para, const graph& DAG, const std::string& seq);
#endif