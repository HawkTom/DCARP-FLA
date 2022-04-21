#ifndef _GLOBALVAR_H
#define _GLOBALVAR_H

#include<iostream>
#include<fstream>
#include<cstring>
#include<cmath>
#include<unordered_map>
#include<time.h>


#define MAX_TASKS_TAG_LENGTH 500
#define MAX_ARCS_TAG_LENGTH 1001
#define MAX_NODE_TAG_LENGTH 300
#define INF 1000000000

#define MAX_TASK_SEG_LENGTH 550
#define MAX_ROUTE_TAG_LENGTH 50

#define MAX_TASK_TAG_LENGTH 500
#define MAX_TASK_SEQ_LENGTH 500


extern int terminal_condition; //if the best has nevere been changed over 20 iterations, stop the algorithm
extern int terminal_duration;
extern int quick_response;

extern int req_arc_num;
extern int nonreq_arc_num;
extern int vertex_num;
extern int req_edge_num;
extern int nonreq_edge_num;
extern int vehicle_num;
extern int capacity;
extern int task_num;
extern int total_arc_num;
extern int DEPOT;
extern int trav_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
extern int serve_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];

extern int min_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
extern int shortest_path[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];

extern int cost_backup[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];

extern char map[20];


extern int costlb;
extern int costub;
extern int demandlb;
extern int demandub;
extern int edge_index[2];

extern std::unordered_map<std::string, unsigned int> edgesIDs;

typedef struct task
{
    int head_node;
    int tail_node;
    int head_node_r;
    int tail_node_r;
    int dead_cost;
    int serv_cost;
    int demand;
    int inverse;
    int vt;
} Task;

// format: head ----> tail
typedef struct arc
{
    int head_node;
    int tail_node;
    int trav_cost;
    unsigned int change;
    int link;
} Arc;

typedef struct individual
{
    int Sequence[MAX_TASK_SEQ_LENGTH]; // task id with depot inside
    int Assignment[250]; // when the representation is chromosome,this is used as the chromosome.
    int TotalCost;
    int Loads[80];
    int TotalVioLoad;
    double Fitness;
//    int start[101]; // for DCARP
} CARPInd;

typedef struct vehicle_state
{
    int remain_seqs[MAX_TASK_SEQ_LENGTH];
    int stop[101];
    int remain_capacity[101];
    int not_served_task_seq[MAX_TASK_SEQ_LENGTH];
} Vehicles;

typedef struct localoptimal
{
    int cost;
    int sequence[MAX_TASK_SEQ_LENGTH];
    int hitnum;
    float usedtime;
    int count;
} LOptimal;



#endif