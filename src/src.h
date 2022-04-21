#ifndef  _SRC_H
#define _SRC_H

#include "../globalvar.h"
#include <fstream>

// function for read map data
void readMap(Task *inst_tasks, Arc *inst_arcs, const char *map1);
void save_instance_to_xml(const Task *inst_tasks, const Arc *inst_arcs, const Vehicles state, int instance_idx);
void read_instance_from_xml(Task *inst_tasks, Arc *inst_arcs, Vehicles *state, int instance_idx);


// calculate cost of map
void mod_dijkstra();
void update_cost(const Task *inst_tasks, const Arc *inst_arcs);
int get_additional_cost(Vehicles state);
int get_task_seq_total_cost(int *task_seq, const Task *inst_tasks);


// heuristic algorithms
void path_scanning(CARPInd *ps_indi, const Task *inst_tasks, const int *serve_mark);
void FredericksonHeuristic(int *FHRoute, int *Route, const Task *inst_tasks);
int split(int *split_task_seq, int *one_task_seq, int *split_route_loads, const Task *inst_tasks);
void augment_merge(CARPInd *am_indi, const Task *inst_tasks);
void rand_scanning(CARPInd *rs_indi, const Task *inst_tasks);


// module for processing virtual task related solution/individual
void indi_route_converter(CARPInd *dst, CARPInd *src, const Task *inst_tasks);
void construct_virtual_task(const Task *inst_tasks, Task *tasks_vt, const int *stop, const int *remain_capacity);
void inher_solution(CARPInd *inhrSolution, Vehicles state, const Task *inst_tasks_vt);
int repair_solution_greedy_insertion(CARPInd *solution, int *remain_seq, const int *stop, const Task *inst_tasks_vt);
void remain_solution(CARPInd *inhrSolution, Vehicles state, const Task *inst_tasks_vt);

// array operations functions
void delete_element(int *a, int k);
void find_ele_positions(int *positions, int *a, int e);
void add_element(int *a, int e, int k);
int rand_choose(int num);
void rand_perm(int *a, int num);
int max(int *Array);
void rand_selection(int *id1, int *id2, int popsize);
void AssignArray(int *Array1, int *Array2);
void AssignSubArray(int *Array1, int k1, int k2, int *Array2);
void JoinArray(int *JointArray, int *Array);
void ReverseDirection(int *Array, int k1, int k2);
int find_min(int *Array);


// simulator to generate the new DCARP instance
void nextScenario(CARPInd *Solution, Task *inst_tasks_vt, Task *inst_tasks, Arc *inst_arcs, Vehicles *state, unsigned int seed);

// other 'solution'/'individual' related functions
void clear_solution(CARPInd *solution);
void copy_individual(CARPInd *dest, CARPInd *src);
// void check_solution_valid(CARPInd solution, const Task *inst_task);
void check_seq_valid(CARPInd solution, const Task *inst_task);
int check_task_valid(int *seq);
int check_cost(CARPInd solution, const Task *inst_tasks);
int FindTask(int a, int b, const Task *inst_tasks, int NO_Task);
int get_total_vio_load(int *route_seg_load);
void saveResult(char *algorithm, int best);
void saveERT(char *algorithm, int state, double tpoint, int best);



#endif