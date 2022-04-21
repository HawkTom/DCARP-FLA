#include "ls.h"


int client_single_insertion(CARPInd *indi, CARPInd *best_nb_solution, const Task *inst_tasks)
{

    int i, j, k;
    lns_route curr_solution, next_solution;

    char type[30];

    strcpy(type, "Single Insertion");

    int Positions[101];
    find_ele_positions(Positions, indi->Sequence, 0);
    curr_solution.Route[0][0] = Positions[0] - 1;
    curr_solution.loads[0] = Positions[0] - 1;
    for (i = 1; i < Positions[0]; i++)
    {
        AssignSubArray(indi->Sequence, Positions[i], Positions[i+1], curr_solution.Route[i]); // Route[i]: 0 x x x x 0
        curr_solution.loads[i] = 0;
        for (j = Positions[i]+1; j < Positions[i+1]; j++)
        {
            curr_solution.loads[i] += inst_tasks[indi->Sequence[j]].demand;
        }
    }
    curr_solution.fitness = indi->TotalCost + LAMBDA * indi->TotalVioLoad;;
    curr_solution.total_vio_loads = indi->TotalVioLoad;
    curr_solution.total_cost = indi->TotalCost;

    int u, pos_u, trip_u, v, pos_v, trip_v, x, y, inv_u;
    int improve = 0;
    int global_improve = 0;
    // check_cost1(curr_solution, inst_tasks);
    
    while (1)
    {
        next_solution = curr_solution;
        for (trip_u = 1; trip_u <= curr_solution.Route[0][0]; trip_u ++)
        {
            for (pos_u = 2; pos_u < curr_solution.Route[trip_u][0]; pos_u++)
            {
                u = curr_solution.Route[trip_u][pos_u];

                for (trip_v = trip_u; trip_v <= curr_solution.Route[0][0]; trip_v ++)
                {
                    for (pos_v = pos_u + 1; pos_v < curr_solution.Route[trip_v][0]; pos_v ++)
                    {

                        v = curr_solution.Route[trip_v][pos_v];

                        improve = ls_single_insertion(&curr_solution, &next_solution, u, v, trip_u, trip_v, pos_u, pos_v, inst_tasks);

                        if (improve){
//                            printf("single insertion \n");
                            goto new_step;
                        }
                    }
                }
            }
        }


    new_step:
        if (improve)
        {
            // update to the server
            update_global_best_solution(&curr_solution, best_nb_solution, inst_tasks, type);
            global_improve = 1;
            improve = 0;
            continue;
        } else {
            break;
        }
    }



    update_global_best_solution(&curr_solution, best_nb_solution, inst_tasks, type);
    
    CARPInd mted_child;
    // route -> sequence
    mted_child.Sequence[0] = 1;
    k = 0;
    for (i = 1; i <= curr_solution.Route[0][0]; i++)
    {
        if (curr_solution.Route[i][0] > 2)
        {
            mted_child.Sequence[0] --;
            JoinArray(mted_child.Sequence, curr_solution.Route[i]);
            k ++;
            mted_child.Loads[k] = curr_solution.loads[i];
        }
    }
    mted_child.Loads[0] = k;
    mted_child.TotalCost = curr_solution.total_cost;
    mted_child.TotalVioLoad = curr_solution.total_vio_loads;

    if ( mted_child.TotalCost != get_task_seq_total_cost(mted_child.Sequence, inst_tasks))
    {
        printf("101: lmsals total cost error \n");
        // longjmp(buf, 2);
        exit(0);
    }
    return global_improve;

}


int client_double_insertion(CARPInd *indi, CARPInd *best_nb_solution, const Task *inst_tasks)
{

    int i, j, k;
    lns_route curr_solution, next_solution;

    char type[30];

    strcpy(type, "Double Insertion");

    int Positions[101];
    find_ele_positions(Positions, indi->Sequence, 0);
    curr_solution.Route[0][0] = Positions[0] - 1;
    curr_solution.loads[0] = Positions[0] - 1;
    for (i = 1; i < Positions[0]; i++)
    {
        AssignSubArray(indi->Sequence, Positions[i], Positions[i+1], curr_solution.Route[i]); // Route[i]: 0 x x x x 0
        curr_solution.loads[i] = 0;
        for (j = Positions[i]+1; j < Positions[i+1]; j++)
        {
            curr_solution.loads[i] += inst_tasks[indi->Sequence[j]].demand;
        }
    }
    curr_solution.fitness = indi->TotalCost + LAMBDA * indi->TotalVioLoad;;
    curr_solution.total_vio_loads = indi->TotalVioLoad;
    curr_solution.total_cost = indi->TotalCost;

    int u, pos_u, trip_u, v, pos_v, trip_v, x, y, inv_u;
    int improve = 0;
    int global_improve = 0;

    
    while (1)
    {
        next_solution = curr_solution;
        for (trip_u = 1; trip_u <= curr_solution.Route[0][0]; trip_u ++)
        {
            for (pos_u = 2; pos_u < curr_solution.Route[trip_u][0]; pos_u++)
            {
                u = curr_solution.Route[trip_u][pos_u];

                for (trip_v = trip_u; trip_v <= curr_solution.Route[0][0]; trip_v ++)
                {
                    for (pos_v = pos_u + 1; pos_v < curr_solution.Route[trip_v][0]; pos_v ++)
                    {

                        v = curr_solution.Route[trip_v][pos_v];

                        if (trip_u == trip_v && pos_v - pos_u > 1 || trip_v != trip_u)
                        {

                            x = curr_solution.Route[trip_u][pos_u+1];
                            if (trip_v != trip_u && x != 0)
                            {
                                improve = ls_double_insertion(&curr_solution, &next_solution, u, x, v, trip_u, trip_v, pos_u, pos_v, inst_tasks);
                                // check_cost1(curr_solution, inst_tasks);
                                if (improve){
//                                    printf("double insertion \n");
                                    goto new_step;
                                }
                            }
                        }



                    }
                }
            }
        }


    new_step:
        if (improve)
        {
            // update to the server
            update_global_best_solution(&curr_solution, best_nb_solution, inst_tasks, type);
            global_improve = 1;
            improve = 0;
            continue;
            // break;
        } else {
            break;
        }
    }


    update_global_best_solution(&curr_solution, best_nb_solution, inst_tasks, type);
    
    CARPInd mted_child;
    // route -> sequence
    mted_child.Sequence[0] = 1;
    k = 0;
    for (i = 1; i <= curr_solution.Route[0][0]; i++)
    {
        if (curr_solution.Route[i][0] > 2)
        {
            mted_child.Sequence[0] --;
            JoinArray(mted_child.Sequence, curr_solution.Route[i]);
            k ++;
            mted_child.Loads[k] = curr_solution.loads[i];
        }
    }
    mted_child.Loads[0] = k;
    mted_child.TotalCost = curr_solution.total_cost;
    mted_child.TotalVioLoad = curr_solution.total_vio_loads;

    if ( mted_child.TotalCost != get_task_seq_total_cost(mted_child.Sequence, inst_tasks))
    {
        printf("214: lmsals total cost error \n");
        // longjmp(buf, 2);
        exit(0);
    }
    return global_improve;

}


int client_swap(CARPInd *indi, CARPInd *best_nb_solution, const Task *inst_tasks)
{

    int i, j, k;
    lns_route curr_solution, next_solution;

    char type[30];

    strcpy(type, "SWAP");

    int Positions[101];
    find_ele_positions(Positions, indi->Sequence, 0);
    curr_solution.Route[0][0] = Positions[0] - 1;
    curr_solution.loads[0] = Positions[0] - 1;
    for (i = 1; i < Positions[0]; i++)
    {
        AssignSubArray(indi->Sequence, Positions[i], Positions[i+1], curr_solution.Route[i]); // Route[i]: 0 x x x x 0
        curr_solution.loads[i] = 0;
        for (j = Positions[i]+1; j < Positions[i+1]; j++)
        {
            curr_solution.loads[i] += inst_tasks[indi->Sequence[j]].demand;
        }
    }
    curr_solution.fitness = indi->TotalCost + LAMBDA * indi->TotalVioLoad;;
    curr_solution.total_vio_loads = indi->TotalVioLoad;
    curr_solution.total_cost = indi->TotalCost;

    int u, pos_u, trip_u, v, pos_v, trip_v, x, y, inv_u;
    int improve = 0;
    int global_improve = 0;

    
    while (1)
    {
        next_solution = curr_solution;
        for (trip_u = 1; trip_u <= curr_solution.Route[0][0]; trip_u ++)
        {
            for (pos_u = 2; pos_u < curr_solution.Route[trip_u][0]; pos_u++)
            {
                u = curr_solution.Route[trip_u][pos_u];

                for (trip_v = trip_u; trip_v <= curr_solution.Route[0][0]; trip_v ++)
                {
                    for (pos_v = pos_u + 1; pos_v < curr_solution.Route[trip_v][0]; pos_v ++)
                    {

                        v = curr_solution.Route[trip_v][pos_v];

                        if (trip_u == trip_v && pos_v - pos_u > 1 || trip_v != trip_u)
                        {
                            improve = ls_swap(&curr_solution, &next_solution, u, v, trip_u, trip_v, pos_u, pos_v, inst_tasks);
                            // check_cost1(curr_solution, inst_tasks);
                            if (improve){
//                                printf("swap \n");
                                goto new_step;
                            }
                        }

                    }
                }
            }
        }


    new_step:
        if (improve)
        {
            // update to the server
            update_global_best_solution(&curr_solution, best_nb_solution, inst_tasks, type);
            global_improve = 1;
            improve = 0;
            continue;
            // break;
        } else {
            break;
        }
    }


    update_global_best_solution(&curr_solution, best_nb_solution, inst_tasks, type);
    
    CARPInd mted_child;
    // route -> sequence
    mted_child.Sequence[0] = 1;
    k = 0;
    for (i = 1; i <= curr_solution.Route[0][0]; i++)
    {
        if (curr_solution.Route[i][0] > 2)
        {
            mted_child.Sequence[0] --;
            JoinArray(mted_child.Sequence, curr_solution.Route[i]);
            k ++;
            mted_child.Loads[k] = curr_solution.loads[i];
        }
    }
    mted_child.Loads[0] = k;
    mted_child.TotalCost = curr_solution.total_cost;
    mted_child.TotalVioLoad = curr_solution.total_vio_loads;

    if ( mted_child.TotalCost != get_task_seq_total_cost(mted_child.Sequence, inst_tasks))
    {
        printf("320: lmsals total cost error \n");
        // longjmp(buf, 2);
        exit(0);
    }

    return global_improve;
}

int client_reverse(CARPInd *indi, CARPInd *best_nb_solution, const Task *inst_tasks)
{

    int i, j, k;
    lns_route curr_solution, next_solution;

    char type[30];

    strcpy(type, "Reverse");

    int Positions[101];
    find_ele_positions(Positions, indi->Sequence, 0);
    curr_solution.Route[0][0] = Positions[0] - 1;
    curr_solution.loads[0] = Positions[0] - 1;
    for (i = 1; i < Positions[0]; i++)
    {
        AssignSubArray(indi->Sequence, Positions[i], Positions[i+1], curr_solution.Route[i]); // Route[i]: 0 x x x x 0
        curr_solution.loads[i] = 0;
        for (j = Positions[i]+1; j < Positions[i+1]; j++)
        {
            curr_solution.loads[i] += inst_tasks[indi->Sequence[j]].demand;
        }
    }
    curr_solution.fitness = indi->TotalCost + LAMBDA * indi->TotalVioLoad;;
    curr_solution.total_vio_loads = indi->TotalVioLoad;
    curr_solution.total_cost = indi->TotalCost;

    int u, pos_u, trip_u, v, pos_v, trip_v, x, y, inv_u;
    int improve = 0;
    int global_improve = 0;

    
    while (1)
    {
        next_solution = curr_solution;
        for (trip_u = 1; trip_u <= curr_solution.Route[0][0]; trip_u ++)
        {
            for (pos_u = 2; pos_u < curr_solution.Route[trip_u][0]; pos_u++)
            {
                u = curr_solution.Route[trip_u][pos_u];

                for (trip_v = trip_u; trip_v <= curr_solution.Route[0][0]; trip_v ++)
                {
                    // printf("%d, %d ||", curr_solution.Route[trip_u][0],  curr_solution.Route[trip_v][0]);
                    for (pos_v = 2; pos_v < curr_solution.Route[trip_v][0]; pos_v ++)
                    {

                        // printf("pos: %d %d\n", pos_u, pos_v);
                        v = curr_solution.Route[trip_v][pos_v];

                        x = curr_solution.Route[trip_u][pos_u+1];

                        y = curr_solution.Route[trip_v][pos_v+1];
                        // printf("u: %d x:%d v:%d y:%d\n", u,x,v,y);

                        // printf("trip: %d %d\n", trip_u, trip_v);
                        // if (pos_u == 2 && pos_v == 2 && trip_u == 1 && trip_v == 4)
                        // {
                        //     printf("stop\n");
                        // }
                        if (inst_tasks[curr_solution.Route[trip_u][pos_u]].vt > 0)
                            continue;

                        if (trip_u == trip_v && y != 0 && pos_v > pos_u)
                        {
                            // check_cost1(curr_solution, inst_tasks);
                            improve = ls_reverse(&curr_solution, &next_solution, u, v, trip_u, trip_v, pos_u, pos_v, inst_tasks);
                            // check_cost1(curr_solution, inst_tasks);
                            if (improve){
                                goto new_step;
                            }
                        }               
                    }
                }
            }
        }


    new_step:
        if (improve)
        {
            // update to the server
            update_global_best_solution(&curr_solution, best_nb_solution, inst_tasks, type);
            global_improve = 1;
            improve = 0;
            continue;
        } else {
            break;
        }
    }


    update_global_best_solution(&curr_solution, best_nb_solution, inst_tasks, type);
    
    CARPInd mted_child;
    // route -> sequence
    mted_child.Sequence[0] = 1;
    k = 0;
    for (i = 1; i <= curr_solution.Route[0][0]; i++)
    {
        if (curr_solution.Route[i][0] > 2)
        {
            mted_child.Sequence[0] --;
            JoinArray(mted_child.Sequence, curr_solution.Route[i]);
            k ++;
            mted_child.Loads[k] = curr_solution.loads[i];
        }
    }
    mted_child.Loads[0] = k;
    mted_child.TotalCost = curr_solution.total_cost;
    mted_child.TotalVioLoad = curr_solution.total_vio_loads;

    if ( mted_child.TotalCost != get_task_seq_total_cost(mted_child.Sequence, inst_tasks))
    {
        printf("436: lmsals total cost error \n");
        // longjmp(buf, 2);
        exit(0);
    }
    return global_improve;
}

int client_cross(CARPInd *indi, CARPInd *best_nb_solution, const Task *inst_tasks)
{
    int i, j, k;
    lns_route curr_solution, next_solution;

    char type[30];

    strcpy(type, "Cross");

    int Positions[101];
    find_ele_positions(Positions, indi->Sequence, 0);
    curr_solution.Route[0][0] = Positions[0] - 1;
    curr_solution.loads[0] = Positions[0] - 1;
    for (i = 1; i < Positions[0]; i++)
    {
        AssignSubArray(indi->Sequence, Positions[i], Positions[i+1], curr_solution.Route[i]); // Route[i]: 0 x x x x 0
        curr_solution.loads[i] = 0;
        for (j = Positions[i]+1; j < Positions[i+1]; j++)
        {
            curr_solution.loads[i] += inst_tasks[indi->Sequence[j]].demand;
        }
    }
    curr_solution.fitness = indi->TotalCost + LAMBDA * indi->TotalVioLoad;;
    curr_solution.total_vio_loads = indi->TotalVioLoad;
    curr_solution.total_cost = indi->TotalCost;

    int u, pos_u, trip_u, v, pos_v, trip_v, x, y, inv_u;
    int improve = 0;
    int global_improve = 0;

    
    while (1)
    {
        next_solution = curr_solution;
        for (trip_u = 1; trip_u <= curr_solution.Route[0][0]; trip_u ++)
        {
            for (pos_u = 2; pos_u < curr_solution.Route[trip_u][0]; pos_u++)
            {
                u = curr_solution.Route[trip_u][pos_u];

                for (trip_v = trip_u; trip_v <= curr_solution.Route[0][0]; trip_v ++)
                {
                    // printf("%d, %d ||", curr_solution.Route[trip_u][0],  curr_solution.Route[trip_v][0]);
                    for (pos_v = 2; pos_v < curr_solution.Route[trip_v][0]; pos_v ++)
                    {

                        // printf("pos: %d %d\n", pos_u, pos_v);
                        v = curr_solution.Route[trip_v][pos_v];

                        x = curr_solution.Route[trip_u][pos_u+1];

                        y = curr_solution.Route[trip_v][pos_v+1];
                        // printf("u: %d x:%d v:%d y:%d\n", u,x,v,y);

                        // printf("trip: %d %d\n", trip_u, trip_v);
                        // if (pos_u == 2 && pos_v == 2 && trip_u == 1 && trip_v == 4)
                        // {
                        //     printf("stop\n");
                        // }

                        if (trip_v != trip_u && x != 0 && y != 0)
                        {
                            improve = ls_cross(&curr_solution, &next_solution, u, v, trip_u, trip_v, pos_u, pos_v, inst_tasks);
                            if (improve){
                                goto new_step;
                            }
                        }

                        

                    }
                }
            }
        }


    new_step:
        if (improve)
        {
            // update to the server
            update_global_best_solution(&curr_solution, best_nb_solution, inst_tasks, type);
            global_improve = 1;
            improve = 0;
            continue;
        } else {
            break;
        }
    }

    // printf("xxxx\n");
    update_global_best_solution(&curr_solution, best_nb_solution, inst_tasks, type);
    
    CARPInd mted_child;
    // route -> sequence
    mted_child.Sequence[0] = 1;
    k = 0;
    for (i = 1; i <= curr_solution.Route[0][0]; i++)
    {
        if (curr_solution.Route[i][0] > 2)
        {
            mted_child.Sequence[0] --;
            JoinArray(mted_child.Sequence, curr_solution.Route[i]);
            k ++;
            mted_child.Loads[k] = curr_solution.loads[i];
        }
    }
    mted_child.Loads[0] = k;
    mted_child.TotalCost = curr_solution.total_cost;
    mted_child.TotalVioLoad = curr_solution.total_vio_loads;

    if ( mted_child.TotalCost != get_task_seq_total_cost(mted_child.Sequence, inst_tasks))
    {
        printf("553: lmsals total cost error \n");
        // longjmp(buf, 2);
        exit(0);
    }
    if ( best_nb_solution->Sequence[0] == 0)
    {
        printf("564: lmsals total cost error \n");
        // longjmp(buf, 2);
        exit(0);
    }
    return global_improve;

}


int ls_single_insertion(lns_route *curr_solution, lns_route *next_solution, int u, int v, int trip_u, int trip_v, int pos_u, int pos_v, const Task *inst_tasks)
{
    next_solution->loads[trip_u] -= inst_tasks[u].demand;
    next_solution->loads[trip_v] += inst_tasks[u].demand;


    if (curr_solution->loads[trip_u] > capacity)
        next_solution->total_vio_loads -= curr_solution->loads[trip_u] - capacity;
    if (curr_solution->loads[trip_v] > capacity)
        next_solution->total_vio_loads -= curr_solution->loads[trip_v] - capacity;

    if (next_solution->loads[trip_u] > capacity)
        next_solution->total_vio_loads += next_solution->loads[trip_u] - capacity;
    if (next_solution->loads[trip_v] > capacity)
        next_solution->total_vio_loads += next_solution->loads[trip_v] - capacity;



    int flag;
    if (pos_v == 2) // first task
    {
        flag = 1;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                    - min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[v].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[u].head_node]
                                    + min_cost[inst_tasks[u].tail_node][inst_tasks[v].head_node];

    }
    else
    {
        flag = 0;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                    - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                    + min_cost[inst_tasks[v].tail_node][inst_tasks[u].head_node]
                                    + min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];
    }

    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;
    if (next_solution->fitness < curr_solution->fitness)
    {
        // printf("-----\n");
        // check_cost1(*curr_solution, inst_tasks);
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        curr_solution->loads[trip_u] = next_solution->loads[trip_u];
        curr_solution->loads[trip_v] = next_solution->loads[trip_v];

        if (flag) {
            add_element(curr_solution->Route[trip_v], u, pos_v);
        } else {
            add_element(curr_solution->Route[trip_v], u, pos_v+1);
        }
        delete_element(curr_solution->Route[trip_u], pos_u);
        // check_cost1(*curr_solution, inst_tasks);
        return 1;
    }

    int inv_u = inst_tasks[u].inverse;
    if (pos_v == 2) // first task
    {
        flag = 1;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                    - min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[v].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[inv_u].head_node]
                                    + min_cost[inst_tasks[inv_u].tail_node][inst_tasks[v].head_node];

    }
    else
    {
        flag = 0;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                    - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                    + min_cost[inst_tasks[v].tail_node][inst_tasks[inv_u].head_node]
                                    + min_cost[inst_tasks[inv_u].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];
    }
    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;

    if (next_solution->fitness < curr_solution->fitness)
    {
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        curr_solution->loads[trip_u] = next_solution->loads[trip_u];
        curr_solution->loads[trip_v] = next_solution->loads[trip_v];

        if (flag) {
            add_element(curr_solution->Route[trip_v], inv_u, pos_v);
        } else {
            add_element(curr_solution->Route[trip_v], inv_u, pos_v + 1);
        }
        delete_element(curr_solution->Route[trip_u], pos_u);
        // check_cost1(*curr_solution, inst_tasks);
        return 1;
    }


    next_solution->total_cost = curr_solution->total_cost;
    next_solution->total_vio_loads = curr_solution->total_vio_loads;
    next_solution->fitness = curr_solution->fitness;
    next_solution->loads[trip_u] = curr_solution->loads[trip_u];
    next_solution->loads[trip_v] = curr_solution->loads[trip_v];
    return 0;
}

int ls_double_insertion(lns_route *curr_solution, lns_route *next_solution, int u, int x, int v, int trip_u, int trip_v, int pos_u, int pos_v, const Task *inst_tasks)
{

    next_solution->loads[trip_u] -= inst_tasks[u].demand + inst_tasks[x].demand;
    next_solution->loads[trip_v] += inst_tasks[u].demand + inst_tasks[x].demand;

    if (curr_solution->loads[trip_u] > capacity)
        next_solution->total_vio_loads -= curr_solution->loads[trip_u] - capacity;
    if (curr_solution->loads[trip_v] > capacity)
        next_solution->total_vio_loads -= curr_solution->loads[trip_v] - capacity;

    if (next_solution->loads[trip_u] > capacity)
        next_solution->total_vio_loads += next_solution->loads[trip_u] - capacity;
    if (next_solution->loads[trip_v] > capacity)
        next_solution->total_vio_loads += next_solution->loads[trip_v] - capacity;

    int flag;
    if (pos_v == 2)
    {
        flag = 1;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[x].head_node]
                                    - min_cost[inst_tasks[x].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    - min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[v].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[u].head_node]
                                    + min_cost[inst_tasks[u].tail_node][inst_tasks[x].head_node]
                                    + min_cost[inst_tasks[x].tail_node][inst_tasks[v].head_node];
    } else {
        flag = 0;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[x].head_node]
                                    - min_cost[inst_tasks[x].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                    + min_cost[inst_tasks[v].tail_node][inst_tasks[u].head_node]
                                    + min_cost[inst_tasks[u].tail_node][inst_tasks[x].head_node]
                                    + min_cost[inst_tasks[x].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];
    }
    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;
    if (next_solution->fitness < curr_solution->fitness)
    {
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        curr_solution->loads[trip_u] = next_solution->loads[trip_u];
        curr_solution->loads[trip_v] = next_solution->loads[trip_v];
        if (flag) {
            add_element(curr_solution->Route[trip_v], x, pos_v);
            add_element(curr_solution->Route[trip_v], u, pos_v);
        } else {
            add_element(curr_solution->Route[trip_v], u, pos_v + 1);
            add_element(curr_solution->Route[trip_v], x, pos_v + 2);
        }

        delete_element(curr_solution->Route[trip_u], pos_u+1);
        delete_element(curr_solution->Route[trip_u], pos_u);
        return 1;
    }

    int inv_u = inst_tasks[u].inverse;
    if (pos_v == 2)
    {
        flag = 1;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[x].head_node]
                                    - min_cost[inst_tasks[x].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    - min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[v].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[inv_u].head_node]
                                    + min_cost[inst_tasks[inv_u].tail_node][inst_tasks[x].head_node]
                                    + min_cost[inst_tasks[x].tail_node][inst_tasks[v].head_node];
    } else {
        flag = 0;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[x].head_node]
                                    - min_cost[inst_tasks[x].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                    + min_cost[inst_tasks[v].tail_node][inst_tasks[inv_u].head_node]
                                    + min_cost[inst_tasks[inv_u].tail_node][inst_tasks[x].head_node]
                                    + min_cost[inst_tasks[x].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];
    }
    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;
    if (next_solution->fitness < curr_solution->fitness)
    {
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        curr_solution->loads[trip_u] = next_solution->loads[trip_u];
        curr_solution->loads[trip_v] = next_solution->loads[trip_v];
        if (flag) {
            add_element(curr_solution->Route[trip_v], x, pos_v);
            add_element(curr_solution->Route[trip_v], inv_u, pos_v);
        } else {
            add_element(curr_solution->Route[trip_v], inv_u, pos_v + 1);
            add_element(curr_solution->Route[trip_v], x, pos_v + 2);
        }

        delete_element(curr_solution->Route[trip_u], pos_u+1);
        delete_element(curr_solution->Route[trip_u], pos_u);
        return 1;
    }


    int inv_x = inst_tasks[x].inverse;
    if (pos_v == 2)
    {
        flag = 1;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[x].head_node]
                                    - min_cost[inst_tasks[x].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    - min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[v].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[u].head_node]
                                    + min_cost[inst_tasks[u].tail_node][inst_tasks[inv_x].head_node]
                                    + min_cost[inst_tasks[inv_x].tail_node][inst_tasks[v].head_node];
    } else {
        flag = 0;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[x].head_node]
                                    - min_cost[inst_tasks[x].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                    + min_cost[inst_tasks[v].tail_node][inst_tasks[u].head_node]
                                    + min_cost[inst_tasks[u].tail_node][inst_tasks[inv_x].head_node]
                                    + min_cost[inst_tasks[inv_x].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];
    }
    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;
    if (next_solution->fitness < curr_solution->fitness)
    {
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        curr_solution->loads[trip_u] = next_solution->loads[trip_u];
        curr_solution->loads[trip_v] = next_solution->loads[trip_v];
        if (flag) {
            add_element(curr_solution->Route[trip_v], inv_x, pos_v);
            add_element(curr_solution->Route[trip_v], u, pos_v);
        } else {
            add_element(curr_solution->Route[trip_v], u, pos_v + 1);
            add_element(curr_solution->Route[trip_v], inv_x, pos_v + 2);
        }

        delete_element(curr_solution->Route[trip_u], pos_u+1);
        delete_element(curr_solution->Route[trip_u], pos_u);
        return 1;
    }

    if (pos_v == 2)
    {
        flag = 1;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[x].head_node]
                                    - min_cost[inst_tasks[x].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    - min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[v].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[inv_u].head_node]
                                    + min_cost[inst_tasks[inv_u].tail_node][inst_tasks[inv_x].head_node]
                                    + min_cost[inst_tasks[inv_x].tail_node][inst_tasks[v].head_node];
    } else {
        flag = 0;
        next_solution->total_cost = curr_solution->total_cost
                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[x].head_node]
                                    - min_cost[inst_tasks[x].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+2]].head_node]
                                    - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                    + min_cost[inst_tasks[v].tail_node][inst_tasks[inv_u].head_node]
                                    + min_cost[inst_tasks[inv_u].tail_node][inst_tasks[inv_x].head_node]
                                    + min_cost[inst_tasks[inv_x].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];
    }
    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;
    if (next_solution->fitness < curr_solution->fitness)
    {
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        curr_solution->loads[trip_u] = next_solution->loads[trip_u];
        curr_solution->loads[trip_v] = next_solution->loads[trip_v];
        if (flag) {
            add_element(curr_solution->Route[trip_v], inv_x, pos_v);
            add_element(curr_solution->Route[trip_v], inv_u, pos_v);
        } else {
            add_element(curr_solution->Route[trip_v], inv_u, pos_v + 1);
            add_element(curr_solution->Route[trip_v], inv_x, pos_v + 2);
        }

        delete_element(curr_solution->Route[trip_u], pos_u+1);
        delete_element(curr_solution->Route[trip_u], pos_u);
        return 1;
    }

    next_solution->total_cost = curr_solution->total_cost;
    next_solution->total_vio_loads = curr_solution->total_vio_loads;
    next_solution->fitness = curr_solution->fitness;
    next_solution->loads[trip_u] = curr_solution->loads[trip_u];
    next_solution->loads[trip_v] = curr_solution->loads[trip_v];
    return 0;
}

int ls_swap(lns_route *curr_solution, lns_route *next_solution, int u, int v, int trip_u, int trip_v, int pos_u, int pos_v, const Task *inst_tasks)
{

    next_solution->loads[trip_u] -= inst_tasks[u].demand - inst_tasks[v].demand;
    next_solution->loads[trip_v] += inst_tasks[u].demand - inst_tasks[v].demand;

    if (curr_solution->loads[trip_u] > capacity)
        next_solution->total_vio_loads -= curr_solution->loads[trip_u] - capacity;
    if (curr_solution->loads[trip_v] > capacity)
        next_solution->total_vio_loads -= curr_solution->loads[trip_v] - capacity;

    if (next_solution->loads[trip_u] > capacity)
        next_solution->total_vio_loads += next_solution->loads[trip_u] - capacity;
    if (next_solution->loads[trip_v] > capacity)
        next_solution->total_vio_loads += next_solution->loads[trip_v] - capacity;

    
    if (next_solution->total_vio_loads < 0)
    {
        printf("stop");
    }

    next_solution->total_cost = curr_solution->total_cost
                                - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[v].head_node]
                                + min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                - min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[v].head_node]
                                - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                + min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[u].head_node]
                                + min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];




    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;
    if (next_solution->fitness < curr_solution->fitness)
    {
        // check_cost1(*curr_solution, inst_tasks);
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        curr_solution->loads[trip_u] = next_solution->loads[trip_u];
        curr_solution->loads[trip_v] = next_solution->loads[trip_v];
        curr_solution->Route[trip_u][pos_u] = v;
        curr_solution->Route[trip_v][pos_v] = u;
        // check_cost1(*curr_solution, inst_tasks);
        if (curr_solution->loads[trip_u] > capacity || curr_solution->loads[trip_v] > capacity)
        {
            printf("stop");
        }
        return 1;
    }

    int inv_u = inst_tasks[u].inverse;
    next_solution->total_cost = curr_solution->total_cost
                                - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[v].head_node]
                                + min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                - min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[v].head_node]
                                - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                + min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[inv_u].head_node]
                                + min_cost[inst_tasks[inv_u].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];


    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;
    if (next_solution->fitness < curr_solution->fitness)
    {
        // check_cost1(*curr_solution, inst_tasks);
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        curr_solution->loads[trip_u] = next_solution->loads[trip_u];
        curr_solution->loads[trip_v] = next_solution->loads[trip_v];
        curr_solution->Route[trip_u][pos_u] = v;
        curr_solution->Route[trip_v][pos_v] = inv_u;
        // check_cost1(*curr_solution, inst_tasks);
        if (curr_solution->loads[trip_u] > capacity || curr_solution->loads[trip_v] > capacity)
        {
            printf("stop");
        }
        return 1;
    }


    int inv_v = inst_tasks[v].inverse;
    next_solution->total_cost = curr_solution->total_cost
                                - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[inv_v].head_node]
                                + min_cost[inst_tasks[inv_v].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                - min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[v].head_node]
                                - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                + min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[u].head_node]
                                + min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];

    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;
    if (next_solution->fitness < curr_solution->fitness)
    {
        // check_cost1(*curr_solution, inst_tasks);
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        curr_solution->loads[trip_u] = next_solution->loads[trip_u];
        curr_solution->loads[trip_v] = next_solution->loads[trip_v];
        curr_solution->Route[trip_u][pos_u] = inv_v;
        curr_solution->Route[trip_v][pos_v] = u;
        // check_cost1(*curr_solution, inst_tasks);
        if (curr_solution->loads[trip_u] > capacity || curr_solution->loads[trip_v] > capacity)
        {
            printf("stop");
        }
        return 1;
    }

    next_solution->total_cost = curr_solution->total_cost
                                - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[inv_v].head_node]
                                + min_cost[inst_tasks[inv_v].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                - min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[v].head_node]
                                - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                + min_cost[inst_tasks[curr_solution->Route[trip_v][pos_v-1]].tail_node][inst_tasks[inv_u].head_node]
                                + min_cost[inst_tasks[inv_u].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];

    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;
    if (next_solution->fitness < curr_solution->fitness)
    {
        // check_cost1(*curr_solution, inst_tasks);
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        curr_solution->loads[trip_u] = next_solution->loads[trip_u];
        curr_solution->loads[trip_v] = next_solution->loads[trip_v];
        curr_solution->Route[trip_u][pos_u] = inv_v;
        curr_solution->Route[trip_v][pos_v] = inv_u;
        // check_cost1(*curr_solution, inst_tasks);
        if (curr_solution->loads[trip_u] > capacity || curr_solution->loads[trip_v] > capacity)
        {
            printf("stop");
        }
        return 1;
    }

    next_solution->total_cost = curr_solution->total_cost;
    next_solution->total_vio_loads = curr_solution->total_vio_loads;
    next_solution->fitness = curr_solution->fitness;
    next_solution->loads[trip_u] = curr_solution->loads[trip_u];
    next_solution->loads[trip_v] = curr_solution->loads[trip_v];
    return 0;
}

int ls_reverse(lns_route *curr_solution, lns_route *next_solution, int u, int v, int trip_u, int trip_v, int pos_u, int pos_v, const Task *inst_tasks)
{
    next_solution->total_cost = curr_solution->total_cost
                                - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                + min_cost[inst_tasks[u].tail_node][inst_tasks[v].tail_node]
                                + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];
    // printf("(%d %d), (%d, %d)\n", inst_tasks[u].tail_node, inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node, inst_tasks[v].tail_node, inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node);
    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;

    if (next_solution->fitness < curr_solution->fitness)
    {
        check_cost1(*curr_solution, inst_tasks);
        // print_seq(*curr_solution, inst_tasks, trip_u);
        // printf("next:%d, curr:%d\n", next_solution->total_cost, curr_solution->total_cost);
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        for (int i = pos_u + 1; i <= pos_v; i++)
        {
            curr_solution->Route[trip_u][i] = inst_tasks[curr_solution->Route[trip_u][i]].inverse;
        }
        ReverseDirection(curr_solution->Route[trip_u], pos_u+1, pos_v);
        // print_seq(*curr_solution, inst_tasks, trip_u);
        check_cost1(*curr_solution, inst_tasks);


        return 1;
    }
    next_solution->total_cost = curr_solution->total_cost;
    next_solution->total_vio_loads = curr_solution->total_vio_loads;
    next_solution->fitness = curr_solution->fitness;
    next_solution->loads[trip_u] = curr_solution->loads[trip_u];
    next_solution->loads[trip_v] = curr_solution->loads[trip_v];
    return 0;
}


void print_seq(lns_route curr_solution, const Task *inst_tasks, int u)
{
    for (int i=1; i<= curr_solution.Route[u][0]; i++)
    {
        if (curr_solution.Route[u][i] == 0)
        {
            printf("0  ");
        } else {
            printf("(%d-->%d)  ", inst_tasks[curr_solution.Route[u][i]].head_node, inst_tasks[curr_solution.Route[u][i]].tail_node);
        }
    }
    printf("\n");
}

int ls_cross(lns_route *curr_solution, lns_route *next_solution, int u, int v, int trip_u, int trip_v, int pos_u, int pos_v, const Task *inst_tasks)
{
    int load_seg1 = 0, load_seg2 = 0;
    for (int i = 2; i <= pos_u; i++)
    {
        load_seg1 += inst_tasks[curr_solution->Route[trip_u][i]].demand;
    }
    for (int i = 2; i <= pos_v; i++)
    {
        load_seg2 += inst_tasks[curr_solution->Route[trip_v][i]].demand;
    }

    // case 1  // case 2 is not suitable for the directed CARP
    next_solution->loads[trip_u] = load_seg1 + curr_solution->loads[trip_v] - load_seg2;
    next_solution->loads[trip_v] = load_seg2 + curr_solution->loads[trip_u] - load_seg1;

    if (curr_solution->loads[trip_u] > capacity)
        next_solution->total_vio_loads -= curr_solution->loads[trip_u] - capacity;
    if (curr_solution->loads[trip_v] > capacity)
        next_solution->total_vio_loads -= curr_solution->loads[trip_v] - capacity;

    if (next_solution->loads[trip_u] > capacity)
        next_solution->total_vio_loads += next_solution->loads[trip_u] - capacity;
    if (next_solution->loads[trip_v] > capacity)
        next_solution->total_vio_loads += next_solution->loads[trip_v] - capacity;

    next_solution->total_cost = curr_solution->total_cost
                                - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                + min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                + min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node];

    next_solution->fitness = next_solution->total_cost + LAMBDA * next_solution->total_vio_loads;



    if (next_solution->fitness < curr_solution->fitness)
    {
//            printf("case1 \n");
        curr_solution->total_cost = next_solution->total_cost;
        curr_solution->total_vio_loads = next_solution->total_vio_loads;
        curr_solution->fitness = next_solution->fitness;
        curr_solution->loads[trip_u] = next_solution->loads[trip_u];
        curr_solution->loads[trip_v] = next_solution->loads[trip_v];

        int tmp_route[500];
        memcpy(tmp_route, curr_solution->Route[trip_u], sizeof(curr_solution->Route[trip_u]));
        curr_solution->Route[trip_u][0] = pos_u;
        for (int i = pos_v + 1; i <= curr_solution->Route[trip_v][0]; i++)
        {
            curr_solution->Route[trip_u][0] ++;
            curr_solution->Route[trip_u][curr_solution->Route[trip_u][0]] = curr_solution->Route[trip_v][i];
        }

        curr_solution->Route[trip_v][0] = pos_v;
        for (int i = pos_u + 1; i <= tmp_route[0]; i++)
        {
            curr_solution->Route[trip_v][0] ++;
            curr_solution->Route[trip_v][curr_solution->Route[trip_v][0]] = tmp_route[i];
        }
        // check_cost1(*curr_solution, inst_tasks);
        return 1;
    }
    next_solution->total_cost = curr_solution->total_cost;
    next_solution->total_vio_loads = curr_solution->total_vio_loads;
    next_solution->fitness = curr_solution->fitness;
    next_solution->loads[trip_u] = curr_solution->loads[trip_u];
    next_solution->loads[trip_v] = curr_solution->loads[trip_v];
    return 0;
}