#ifndef HANDLE_INPUT_H
#define HANDLE_INPUT_H
#include <vector>
#include "sequence.h"
#include "parameter.h"

std::vector<seq_t> readFile(const char* path);
void initPara(para_t* para);
#endif