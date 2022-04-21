#ifndef _HYLS_H
#define _HYLS_H

#include "../globalvar.h"
#include "../src/src.h"


#define LAMBDA 100000
#define AXVSIZE 400


typedef struct lns_route{
    int Route[101][MAX_TASK_SEQ_LENGTH];
    int total_cost;
    int loads[101];
    int total_vio_loads;
    double fitness;
} lns_route;



int SingleInsertion(CARPInd *BestSINeighbor, CARPInd *CurrSolution, const Task *inst_tasks);
int DoubleInsertion(CARPInd *BestDINeighbor, CARPInd *CurrSolution, const Task *inst_tasks);
int SWAP(CARPInd *BestDINeighbor, CARPInd *CurrSolution, const Task *inst_tasks);

int check_cost1(lns_route curr_solution, const Task *inst_tasks);


void hillclimbing(CARPInd *indi, CARPInd *LocalOptimal, int move_type, const Task *inst_tasks);

int expand(lns_route *curr_solution, int depth, const Task *inst_tasks);
double exautive_single_insertion(CARPInd indi, CARPInd *local, int *count, const Task *inst_tasks);
void update_global_best_solution(lns_route *route_solution, CARPInd *global_solution, const Task *inst_tasks, const char *type);


double HyLS(CARPInd InitSolution, CARPInd *bestSolution, int *count, const Task *inst_tasks);

void exhaustive_ls(CARPInd *InitSolution, CARPInd *IndiList, int *IndiListLength, const Task *inst_tasks);
int client_cross(CARPInd *indi, CARPInd *best_nb_solution, const Task *inst_tasks);
int client_reverse(CARPInd *indi, CARPInd *best_nb_solution, const Task *inst_tasks);
int client_swap(CARPInd *indi, CARPInd *best_nb_solution, const Task *inst_tasks);
int client_double_insertion(CARPInd *indi, CARPInd *best_nb_solution, const Task *inst_tasks);
int client_single_insertion(CARPInd *indi, CARPInd *best_nb_solution, const Task *inst_tasks);
int ls_single_insertion(lns_route *curr_solution, lns_route *next_solution, int u, int v, int trip_u, int trip_v, int pos_u, int pos_v, const Task *inst_tasks);
int ls_double_insertion(lns_route *curr_solution, lns_route *next_solution, int u, int x, int v, int trip_u, int trip_v, int pos_u, int pos_v, const Task *inst_tasks);
int ls_swap(lns_route *curr_solution, lns_route *next_solution, int u, int v, int trip_u, int trip_v, int pos_u, int pos_v, const Task *inst_tasks);
int ls_reverse(lns_route *curr_solution, lns_route *next_solution, int u, int v, int trip_u, int trip_v, int pos_u, int pos_v, const Task *inst_tasks);
int ls_cross(lns_route *curr_solution, lns_route *next_solution, int u, int v, int trip_u, int trip_v, int pos_u, int pos_v, const Task *inst_tasks);

#endif