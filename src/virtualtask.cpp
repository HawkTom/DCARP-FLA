
// module for processing virtual task related solution/individual

#include "src.h"

void indi_route_converter(CARPInd *dst, CARPInd *src, const Task *inst_tasks)
{
    int i, load;
    load = 0;
    memset(dst->Sequence, 0, sizeof(int) * MAX_TASK_SEQ_LENGTH);
    memset(dst->Loads, 0, sizeof(int) * 50);
    dst->Sequence[0] = 1;
    dst->Sequence[1] = 0;
    for (i = 2; i <= src->Sequence[0]; i++)
    {
        if (src->Sequence[i] == 0)
        {
            dst->Sequence[0] ++;
            dst->Sequence[dst->Sequence[0]] = 0;
            dst->Loads[0] ++;
            dst->Loads[dst->Loads[0]] = load;
            load = 0;
            continue;
        }
        if (inst_tasks[src->Sequence[i]].vt > 0 && src->Sequence[i-1] != 0)
        {
            dst->Sequence[0] ++;
            dst->Sequence[dst->Sequence[0]] = 0;
            dst->Loads[0] ++;
            dst->Loads[dst->Loads[0]] = load;
            load = 0;
        }

        load += inst_tasks[src->Sequence[i]].demand;
        dst->Sequence[0] ++;
        dst->Sequence[dst->Sequence[0]] = src->Sequence[i];
    }
}

void construct_virtual_task(const Task *inst_tasks, Task *tasks_vt, const int *stop, const int *remain_capacity)
{

    int i;

    int actual_req_edge_num = req_edge_num;
    req_edge_num += stop[0];
    task_num = 2 * req_edge_num;

    
    for (i = 1; i <= actual_req_edge_num; i++)
    {
        tasks_vt[i].head_node = inst_tasks[i].head_node;
        tasks_vt[i].tail_node = inst_tasks[i].tail_node;
        tasks_vt[i].dead_cost = inst_tasks[i].dead_cost;
        tasks_vt[i].serv_cost = inst_tasks[i].serv_cost;
        tasks_vt[i].demand = inst_tasks[i].demand;
        tasks_vt[i].inverse = i + req_edge_num;
        tasks_vt[i].vt = 0;
        tasks_vt[i].head_node_r = 2*i-1;
        tasks_vt[i].tail_node_r = 2*i;

        tasks_vt[i + req_edge_num].head_node = tasks_vt[i].tail_node;
        tasks_vt[i + req_edge_num].tail_node = tasks_vt[i].head_node;
        tasks_vt[i + req_edge_num].head_node_r = tasks_vt[i].tail_node_r;
        tasks_vt[i + req_edge_num].tail_node_r = tasks_vt[i].head_node_r;

        tasks_vt[i + req_edge_num].dead_cost = tasks_vt[i].dead_cost;
        tasks_vt[i + req_edge_num].serv_cost = tasks_vt[i].serv_cost;
        tasks_vt[i + req_edge_num].demand = tasks_vt[i].demand;
        tasks_vt[i + req_edge_num].inverse = i;
        tasks_vt[i + req_edge_num].vt = 0;

    }
    if (stop[0] == 0)
    {
        return;
    }
    for (i = actual_req_edge_num+1; i <= req_edge_num; i++)
    {
        tasks_vt[i].head_node = DEPOT;
        tasks_vt[i].tail_node = stop[i-actual_req_edge_num];
        tasks_vt[i].dead_cost = 0;
        tasks_vt[i].serv_cost = min_cost[tasks_vt[i].head_node][tasks_vt[i].tail_node];
        tasks_vt[i].demand = capacity - remain_capacity[i-actual_req_edge_num];
        if (tasks_vt[i].demand == 0)
        {
            tasks_vt[i].demand = 1;
        }
        tasks_vt[i].inverse = i + req_edge_num;
        tasks_vt[i].vt = 1;
        tasks_vt[i].head_node_r = 2*i-1;
        tasks_vt[i].tail_node_r = 2*i;

        tasks_vt[i + req_edge_num].head_node = tasks_vt[i].head_node;
        tasks_vt[i + req_edge_num].tail_node = tasks_vt[i].tail_node;
        tasks_vt[i + req_edge_num].head_node_r = tasks_vt[i].tail_node_r;
        tasks_vt[i + req_edge_num].tail_node_r = tasks_vt[i].head_node_r;
        tasks_vt[i + req_edge_num].dead_cost = tasks_vt[i].dead_cost;
        tasks_vt[i + req_edge_num].serv_cost = tasks_vt[i].serv_cost;
        tasks_vt[i + req_edge_num].demand = tasks_vt[i].demand;
        tasks_vt[i + req_edge_num].inverse = i;
        tasks_vt[i + req_edge_num].vt = 1;
    }
    tasks_vt[0].head_node = DEPOT;
    tasks_vt[0].tail_node = DEPOT;
    tasks_vt[0].dead_cost = 0;
    tasks_vt[0].serv_cost = 0;
    tasks_vt[0].demand = 0;
    tasks_vt[0].inverse = 0;
}

// /**
//  * @brief inher solution: path-sacnning for new tasks
//  * 
//  * @param inhrSolution 
//  * @param state 
//  * @param inst_tasks_vt 
//  */
// void remain_solution(CARPInd *inhrSolution, Vehicles state, const Task *inst_tasks_vt)
// {
//     // check stop points with routes
//     int i, j;

//     int Route[101][MAX_TASK_SEQ_LENGTH];
//     memset(Route, 0, sizeof(Route));

//     int used[MAX_TASKS_TAG_LENGTH];
//     memset(used, 0, sizeof(used));

//     int vtIdx = req_edge_num - state.stop[0];
//     int actual_req_edge_num = vtIdx;
//     for (i = 1; i < state.remain_seqs[0]; i++)
//     {
//         if (state.remain_seqs[i] < 0)
//             continue;
        
//         if (state.remain_seqs[i] > actual_req_edge_num)
//         {
//             Route[Route[0][0]][0] ++;
//             Route[Route[0][0]][Route[Route[0][0]][0]] = state.remain_seqs[i] + state.stop[0];
//         } else if (state.remain_seqs[i] == 0)
//         {
//             vtIdx ++;
//             Route[0][0] ++;
//             Route[Route[0][0]][0] = 1;
//             Route[Route[0][0]][1] = vtIdx;
//         } else
//         {
//             Route[Route[0][0]][0] ++;
//             Route[Route[0][0]][Route[Route[0][0]][0]] = state.remain_seqs[i];
//         }

//         used[Route[Route[0][0]][Route[Route[0][0]][0]]] = 1;
//         used[inst_tasks_vt[Route[Route[0][0]][Route[Route[0][0]][0]]].inverse] = 1;
//         used[0]++;
        
//     }

//     for (i = 1; i < state.not_served_task_seq[0]; i++)
//     {
//         if (state.not_served_task_seq[i] < 0 )
//             continue;

//         if (state.not_served_task_seq[i] > actual_req_edge_num)
//         {
//             // state.not_served_task_seq[i] += state.stop[0];
//             Route[Route[0][0]][0] ++;
//             Route[Route[0][0]][Route[Route[0][0]][0]] = state.not_served_task_seq[i] + state.stop[0];
//         } else if (state.not_served_task_seq[i] == 0)
//         {
//             Route[0][0] ++;
//             Route[Route[0][0]][0] = 0;
//             continue;
//         } else
//         {
//             Route[Route[0][0]][0] ++;
//             Route[Route[0][0]][Route[Route[0][0]][0]] = state.not_served_task_seq[i];
//         }

//         used[Route[Route[0][0]][Route[Route[0][0]][0]]] = 1;
//         used[inst_tasks_vt[Route[Route[0][0]][Route[Route[0][0]][0]]].inverse] = 1;
//         used[0]++;
//     }

//     int addi_tasks[MAX_TASK_SEQ_LENGTH];
//     memset(addi_tasks, 0, sizeof(addi_tasks));

//     int load;
//     int cost[MAX_TASK_SEQ_LENGTH];
//     int aa = 0;
//     for (i = 1; i <= Route[0][0]; i++)
//     {
//         load = 0;
//         memset(cost, INF, sizeof(cost));
//         for (j = 1; j <= Route[i][0]; j++)
//         {
//             load += inst_tasks_vt[Route[i][j]].demand;
//             if (inst_tasks_vt[Route[i][j]].vt > 0)
//             {
//                 cost[j] = INF;
//             } else{
//                 cost[j] = min_cost[inst_tasks_vt[Route[i][j]].head_node][DEPOT] + min_cost[inst_tasks_vt[Route[i][j]].tail_node][DEPOT];
//             }
            
//         }
//         cost[0] = Route[i][0];
//         if (load <= capacity)
//             continue;

//         int min_pos;
//         while (load > capacity)
//         {
//             min_pos = find_min(cost); // serve the farther task firstly
//             load -= inst_tasks_vt[Route[i][min_pos]].demand;
//             cost[min_pos] = INF;
//             addi_tasks[0] ++;
//             addi_tasks[addi_tasks[0]] = Route[i][min_pos];
//             Route[i][min_pos] = -1;
//             aa ++;
//         }
//     }

//     int bb = 0;
//     for(i = 1; i <= req_edge_num; i++)
//     {
//         if (used[i])
//             continue;
        
//         addi_tasks[0] ++;
//         addi_tasks[addi_tasks[0]] = i;
//         bb ++;
//     }

//     int serv_mark[2*req_edge_num+1]; // the mark for served tasks
//     memset(serv_mark, 0, sizeof(serv_mark));

//     for (i=1; i <= addi_tasks[0]; i++)
//     {
//         serv_mark[addi_tasks[i]] = 1;
//         serv_mark[ inst_tasks_vt[addi_tasks[i]].inverse ] = 1;
//     }
//     CARPInd tmp_solution;
//     path_scanning(&tmp_solution, inst_tasks_vt, serv_mark); // new tasks and the tasks violate routes' capacity


//     for (i = 1; i <= Route[0][0]; i++)
//     {
//         inhrSolution->Sequence[0] ++;
//         inhrSolution->Sequence[inhrSolution->Sequence[0]] = 0;
//         load = 0;
//         for (j = 1; j <= Route[i][0]; j++)
//         {
//             if (Route[i][j] < 0)
//                 continue;
            
//             inhrSolution->Sequence[0] ++;
//             inhrSolution->Sequence[inhrSolution->Sequence[0]] = Route[i][j];
//             load += inst_tasks_vt[Route[i][j]].demand;
//         }
//         inhrSolution->Loads[0]++;
//         inhrSolution->Loads[inhrSolution->Loads[0]] = load;
//     }
//     for (i=1; i<=tmp_solution.Sequence[0]; i++)
//     {
//         inhrSolution->Sequence[0] ++;
//         inhrSolution->Sequence[inhrSolution->Sequence[0]] = tmp_solution.Sequence[i];
//     }

//     for (i=1; i<=tmp_solution.Loads[0]; i++)
//     {
//         inhrSolution->Loads[0] ++;
//         inhrSolution->Loads[inhrSolution->Loads[0]] = tmp_solution.Loads[i];
//     }

//     inhrSolution->TotalCost = get_task_seq_total_cost(inhrSolution->Sequence, inst_tasks_vt);
//     inhrSolution->TotalVioLoad = get_total_vio_load(inhrSolution->Loads);
//     check_solution_valid(*inhrSolution, inst_tasks_vt);
//     if (inhrSolution->TotalVioLoad > 0)
//     {
//         printf("inherted solution error \n");
//         exit(-1);
//     }

// }

// /**
//  * @brief inher_solution: split for previous routes, sequence transfer strategy
//  * 
//  * @param inhrSolution 
//  * @param state 
//  * @param inst_tasks_vt 
//  */

// void inher_solution(CARPInd *inhrSolution, Vehicles state, const Task *inst_tasks_vt)
// {
//     // inherience a solution
//     // Individual inhrSolution;
//     int remain_seq_num;
//     remain_seq_num = repair_solution_greedy_insertion(inhrSolution, state.remain_seqs, state.stop, inst_tasks_vt);
    
//     inhrSolution->TotalCost = split(inhrSolution->Sequence, inhrSolution->Assignment, inhrSolution->Loads, inst_tasks_vt);
//     memset(inhrSolution->Loads, 0, sizeof(inhrSolution->Loads));
//     if(!check_task_valid(inhrSolution->Sequence))
//     {
//         printf("Inhrsolution sequence not valid. \n");
//         exit(0);
//     }
//     if (inhrSolution->TotalCost != get_task_seq_total_cost(inhrSolution->Sequence, inst_tasks_vt))
//     {
//         printf("Inhrsolution cost not valid. \n");
//         exit(0);
//     }
//     int load = 0;
//     for (int i = 2; i <= inhrSolution->Sequence[0]; i++)
//     {   
//         if (inhrSolution->Sequence[i] == 0)
//         {
//             inhrSolution->Loads[0] ++;
//             inhrSolution->Loads[inhrSolution->Loads[0]] = load; 
//             load = 0;
//         } else
//         {
//             load += inst_tasks_vt[inhrSolution->Sequence[i]].demand;
//         }
//     }
// }

// int repair_solution_greedy_insertion(CARPInd *solution, int *remain_seq, const int *stop, const Task *inst_tasks_vt)
// {
//     int i, j;
//     int remain_seq_num = 0;

//     int used[MAX_TASKS_TAG_LENGTH];
//     memset(used, 0, sizeof(used));

//     int new_seq[MAX_TASK_SEQ_LENGTH];
//     memset(new_seq, 0, sizeof(new_seq));

//     // check stop points with routes
//     int vtIdx = req_edge_num - stop[0];
//     int actual_req_edge_num = vtIdx;
//     for (i = 1; i < remain_seq[0]; i++)
//     {
//         if (remain_seq[i] < 0)
//             continue;
        
//         if (remain_seq[i] > actual_req_edge_num)
//         {
//             remain_seq[i] += stop[0];
//         }
//         if (remain_seq[i] == 0)
//         {
//             vtIdx ++;
//             remain_seq[i] = vtIdx;
//         }

//         new_seq[0] ++;
//         new_seq[new_seq[0]] = remain_seq[i];
//         used[remain_seq[i]] = 1;
//         used[inst_tasks_vt[remain_seq[i]].inverse] = 1;
//         used[0]++;

//     }
//     remain_seq_num = new_seq[0] - stop[0];

//     int increase_cost, minimal_cost;
//     int minimal_index, flag;
//     for(i = 1; i <= req_edge_num; i++)
//     {
//         if (used[i])
//             continue;


//         increase_cost = min_cost[DEPOT][inst_tasks_vt[i].head_node]
//             + min_cost[inst_tasks_vt[i].tail_node][inst_tasks_vt[new_seq[1]].head_node];
//         minimal_cost = increase_cost;
//         minimal_index = 1;

        

//         for (j = 1; j <= new_seq[0]; j++)
//         {
//             increase_cost = min_cost[inst_tasks_vt[new_seq[j]].tail_node][inst_tasks_vt[i].head_node]
//             + min_cost[inst_tasks_vt[i].tail_node][inst_tasks_vt[new_seq[j+1]].head_node];

//             if (increase_cost < minimal_cost)
//             {
//                 minimal_cost = increase_cost;
//                 minimal_index = j+1;
//             }
//         }

//         flag = 0;
//         increase_cost = min_cost[DEPOT][inst_tasks_vt[i].tail_node]
//             + min_cost[inst_tasks_vt[i].head_node][inst_tasks_vt[new_seq[1]].head_node];
//         if (increase_cost < minimal_cost)
//         {
//             minimal_cost = increase_cost;
//             minimal_index = 1;
//             flag = 1;
//         }

//         for (j = 1; j <= new_seq[0]; j++)
//         {
//             increase_cost = min_cost[inst_tasks_vt[new_seq[j]].tail_node][inst_tasks_vt[i].tail_node]
//                             + min_cost[inst_tasks_vt[i].head_node][inst_tasks_vt[new_seq[j+1]].head_node];

//             if (increase_cost < minimal_cost)
//             {
//                 minimal_cost = increase_cost;
//                 minimal_index = j+1;
//                 flag = 1;
//             }
//         }
//         if (flag)
//         {
//             add_element(new_seq, inst_tasks_vt[i].inverse, minimal_index);
// //            printf("%d %d\n", inst_tasks_vt[i].inverse, new_seq[0]);
//         } else {
//             add_element(new_seq, i, minimal_index);
//         }
//         used[i] = 1;
//         used[inst_tasks_vt[i].inverse] = 1;
//         used[0]++;
//     }

//     memset(solution->Sequence, 0, sizeof(solution->Sequence));
//     memset(solution->Assignment, 0, sizeof(solution->Assignment));
//     for (i = 1; i <= new_seq[0]; i++)
//     {
//         if (new_seq[i] == 0)
//         {
//             continue;
//         }
//         solution->Assignment[0]++;
//         solution->Assignment[solution->Assignment[0]] = new_seq[i];
//     }
//     return remain_seq_num;
// }