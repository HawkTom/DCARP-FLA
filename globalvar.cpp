#include "globalvar.h"



int req_arc_num = 0; //NRA
int req_edge_num; //NRE
int nonreq_arc_num = 0;
int nonreq_edge_num;
int vertex_num;
int vehicle_num;
int capacity;
int task_num;
int total_arc_num;
int DEPOT;
int trav_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int serve_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int min_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int shortest_path[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];

int cost_backup[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
char map[20];


int costlb = INF;
int costub = 0;
int demandlb = INF;
int demandub = 0;

int edge_index[2];


std::unordered_map<std::string, unsigned int> edgesIDs;
std::unordered_map<std::string, unsigned int> SOLARV;