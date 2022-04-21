#ifndef _FLA_H
#define _FLA_H

#include "globalvar.h"
#include "src/src.h"
#include "ls/ls.h"

#include<functional>
#include<sstream>
#include<iterator>
#include<string>
#include<set>
#include<map>
#include<vector>
#include<algorithm>
#include<cmath>
#include "omp.h"


int random_choose_index_from(int lower, int upper);
int random_walk(int (*Route)[MAX_TASK_SEQ_LENGTH], int *Loads, const Task * inst_tasks);

void auto_corrlation(CARPInd sol, const Task * inst_tasks, int instance);
void auto_correlation_analysis(int instance, int tau);
double auto_correlation_analysis1(int instance, int tau);

void local_optimum_related(const Task * inst_tasks, int instance);
void lo_analysis(int instance, const Task * inst_tasks);
void lo_dis_analysis(int instance, const Task * inst_tasks);
size_t solution_hash_id(int *indi_seq, const Task* inst_tasks);
int hamming_distance(int *sol1_seq, int *sol2_seq, const Task* inst_tasks);

void return_probability(int instance, const Task *inst_tasks);

int global_optimum_related(int instance, int *seq, const Task *inst_tasks);

void check_solution_valid(CARPInd solution, const Task *inst_task);
void check_Route_valid(int (*Route)[MAX_TASK_SEQ_LENGTH], int *Loads, int cost, const Task *inst_task);

int cmp(const std::pair<size_t, int>& x, const std::pair<size_t, int>& y);
#endif