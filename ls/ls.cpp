#include "ls.h"

void clear_solution(CARPInd *solution);
void copy_individual(CARPInd *dest, CARPInd *src);

void hillclimbing(CARPInd *indi, CARPInd *LocalOptimal, int move_type, const Task *inst_tasks)
{
    CARPInd BestNeighbor;
    int count;
    while (1)
    {
        switch (move_type)
        {
            case 1:
                /* single insertion */
                count = SingleInsertion(&BestNeighbor, indi, inst_tasks);
                break;

            case 2:
                /* single insertion */
                count = DoubleInsertion(&BestNeighbor, indi, inst_tasks);
                break;

            case 3:
                /* single insertion */
                count = SWAP(&BestNeighbor, indi, inst_tasks);
                break;
        }
        printf("%d %d\n", count, BestNeighbor.TotalCost);
        if (BestNeighbor.TotalCost < indi->TotalCost)
        {
            copy_individual(indi, &BestNeighbor);
            clear_solution(&BestNeighbor);
        } else {
            break;
        }
    }
    
    
}

void clear_solution(CARPInd *solution)
{
    memset(solution->Sequence, 0, sizeof(solution->Sequence));
    memset(solution->Assignment, 0, sizeof(solution->Assignment));
    memset(solution->Loads, 0, sizeof(solution->Loads));
    solution->TotalCost = 0;
    solution->TotalVioLoad = 0;
    solution->Fitness = 0;
}

void copy_individual(CARPInd *dest, CARPInd *src)
{
    memcpy(dest->Sequence, src->Sequence, sizeof(src->Sequence));
    memcpy(dest->Loads, src->Loads, sizeof(src->Loads));
    dest->TotalCost = src->TotalCost;
}

int SingleInsertion(CARPInd *BestSINeighbor, CARPInd *CurrSolution, const Task *inst_tasks)
{
    int i, j, k, u, v, z;
    double PenPrmt = LAMBDA;

    int count = 0;  // record neighbourhood moves (number of fitness evaluation)
    int Besti, Bestj, Bestu, Bestv, RID1, RID2;
    int Positions[101], Route[101][MAX_TASK_SEQ_LENGTH];
    find_ele_positions(Positions, CurrSolution->Sequence, 0);

    Route[0][0] = Positions[0] - 1; // PA*
    for(i = 1; i < Positions[0]; i++)
    {
        AssignSubArray(CurrSolution->Sequence, Positions[i], Positions[i+1], Route[i]);
    }

    CARPInd Neighbor;
    int MovedTasks[3];

    MovedTasks[0] = 1;
    BestSINeighbor->Fitness = INF;
    Bestu = -1;

    for (i = 1; i < Positions[0]; i++)
    {
        RID1 = i;
        for (j = 2; j < Route[i][0]; j++) // Route[i]: 0, t1, t2, t3, ..., tn, 0
        {
            for (u = 0; u < Positions[0]; u++)
            {
                if (u == i)
                    continue; // insertion happens in two different routes

                RID2 = u;
                if (u == 0 && Route[0][0]<vehicle_num && Route[i][0] > 3)
                {
                    // assign a new route, the route number < max vehicle, over one tasks in Route[i];
                    // AssignArray(CurrSolution->Loads, Neighbor.Loads);
                    memcpy(Neighbor.Loads, CurrSolution->Loads, sizeof(CurrSolution->Loads));
                    Neighbor.Loads[i] -= inst_tasks[Route[i][j]].demand;
                    Neighbor.TotalVioLoad = CurrSolution->TotalVioLoad;

                    if (Neighbor.Loads[i] > capacity)
                    {
                        Neighbor.TotalVioLoad -= inst_tasks[Route[i][j]].demand;
                    } else if ( Neighbor.Loads[i] > capacity - inst_tasks[Route[i][j]].demand)
                    {
                        Neighbor.TotalVioLoad -= Neighbor.Loads[i] + inst_tasks[Route[i][j]].demand - capacity;
                    }

                    Neighbor.Loads[0] ++;
                    Neighbor.Loads[Neighbor.Loads[0]] = inst_tasks[Route[i][j]].demand;

                    if (Neighbor.Loads[Neighbor.Loads[0]] > capacity)
                    {
                        Neighbor.TotalVioLoad += Neighbor.Loads[Neighbor.Loads[0]] - capacity;
                    }
                    int tmp_vio_load = get_total_vio_load(Neighbor.Loads);

                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                                - min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                                                + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                                                + min_cost[DEPOT][inst_tasks[Route[i][j]].head_node]
                                                + min_cost[inst_tasks[Route[i][j]].tail_node][DEPOT];
                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;

                    // check tabu list
                    if (Neighbor.Fitness < BestSINeighbor->Fitness)
                    {
                        MovedTasks[1] = Route[i][j];
                        Besti = i;
                        Bestj = j;
                        Bestu = u;
                        Bestv = 0;

                        BestSINeighbor->TotalCost = Neighbor.TotalCost;
                        // AssignArray(Neighbor.Loads, BestSINeighbor->Loads);
                        memcpy(BestSINeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                        BestSINeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                        BestSINeighbor->Fitness = Neighbor.Fitness;
                    }
                }
                if (u == 0)
                    continue;

                for (v = 2; v <= Route[u][0]; v++) // Route[i]: 0, t1, t2, t3, ..., tn, 0
                {
                    if ( u == i && v == j)
                        continue;

                    // AssignArray(CurrSolution->Loads, Neighbor.Loads);
                    memcpy(Neighbor.Loads, CurrSolution->Loads, sizeof(CurrSolution->Loads));
                    Neighbor.Loads[i] -= inst_tasks[Route[i][j]].demand;
                    Neighbor.Loads[u] += inst_tasks[Route[i][j]].demand;
                    Neighbor.TotalVioLoad = CurrSolution->TotalVioLoad;

                    if (Neighbor.Loads[i] > capacity)
                    {
                        Neighbor.TotalVioLoad -= inst_tasks[Route[i][j]].demand;
                    } else if (Neighbor.Loads[i] > capacity - inst_tasks[Route[i][j]].demand)
                    {
                        Neighbor.TotalVioLoad -= (Neighbor.Loads[i] + inst_tasks[Route[i][j]].demand - capacity);
                    }

                    if (Neighbor.Loads[u] > capacity + inst_tasks[Route[i][j]].demand)
                    {
                        Neighbor.TotalVioLoad += inst_tasks[Route[i][j]].demand;

                    } else if (Neighbor.Loads[u] > capacity)
                    {
                        Neighbor.TotalVioLoad += Neighbor.Loads[u] - capacity;
                    }

                    if (Neighbor.Loads[i] == 0)
                    {
                        RID1 = 0;
                        delete_element(Neighbor.Loads, i);
                    }
                    


                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                            - min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                            + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                            - min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                            + min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                            + min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[u][v]].head_node];

                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;

                    if (Neighbor.Fitness < BestSINeighbor->Fitness)
                    {
                        MovedTasks[1] = Route[i][j];
                        Besti = i;
                        Bestj = j;
                        Bestu = u;
                        Bestv = v;

                        BestSINeighbor->TotalCost = Neighbor.TotalCost;
                        // AssignArray(Neighbor.Loads, BestSINeighbor->Loads);
                        memcpy(BestSINeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                        BestSINeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                        BestSINeighbor->Fitness = Neighbor.Fitness;
                    }
                    // invert selected task
                    z = inst_tasks[Route[i][j]].inverse;
                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                         - min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                                         + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                                         - min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                                         + min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[z].head_node]
                                         + min_cost[inst_tasks[z].tail_node][inst_tasks[Route[u][v]].head_node];
                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;
                    if (Neighbor.Fitness < BestSINeighbor->Fitness)
                    {
                        MovedTasks[1] = z;
                        Besti = i;
                        Bestj = j;
                        Bestu = u;
                        Bestv = v;

                        BestSINeighbor->TotalCost = Neighbor.TotalCost;
                        // AssignArray(Neighbor.Loads, BestSINeighbor->Loads);
                        memcpy(BestSINeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                        BestSINeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                        BestSINeighbor->Fitness = Neighbor.Fitness;
                    }
                }
            }
        }
    }

    if (Bestu < 0)
    {
        return count;
    }

    // Assign Sequence
    if (Bestu == 0)
    {
        // new routes
        BestSINeighbor->Sequence[0] = 1;
        for(k = 1; k < Positions[0]; k++)
        {
            if (k == Besti)
            {
                delete_element(Route[k], Bestj);

                BestSINeighbor->Sequence[0] --; // Route[1] = 0
                JoinArray(BestSINeighbor->Sequence, Route[k]);
            } else {
                BestSINeighbor->Sequence[0] --; // Route[1] = 0
                JoinArray(BestSINeighbor->Sequence, Route[k]);
            }
        }
        BestSINeighbor->Sequence[0] ++;
        BestSINeighbor->Sequence[BestSINeighbor->Sequence[0]] = MovedTasks[1];
        BestSINeighbor->Sequence[0] ++;
        BestSINeighbor->Sequence[BestSINeighbor->Sequence[0]] = 0;

    } else{

        BestSINeighbor->Sequence[0] = 1;
        for(k = 1; k < Positions[0]; k++)
        {
            if (k == Besti)
            {
                delete_element(Route[k], Bestj);
                if (Route[k][0] == 2)
                    continue;

                BestSINeighbor->Sequence[0] --; // Route[1] = 0
                JoinArray(BestSINeighbor->Sequence, Route[k]);
            } else if (k == Bestu){
                add_element(Route[k], MovedTasks[1], Bestv);
                BestSINeighbor->Sequence[0] --; // Route[1] = 0
                JoinArray(BestSINeighbor->Sequence, Route[k]);
            } else {
                BestSINeighbor->Sequence[0] --; // Route[1] = 0
                JoinArray(BestSINeighbor->Sequence, Route[k]);
            }
        }
    }

    return count;

}


int DoubleInsertion(CARPInd *BestDINeighbor, CARPInd *CurrSolution, const Task *inst_tasks)
{
    int i, j, k, u, v, z, w;

    int count = 0;  // record neighbourhood moves (number of fitness evaluation)
    int Besti, Bestj, Bestu, Bestv, RID1, RID2;
    int Positions[101], Route[101][MAX_TASK_SEQ_LENGTH];
    find_ele_positions(Positions, CurrSolution->Sequence, 0);

    Route[0][0] = Positions[0] - 1; // PA*
    for(i = 1; i < Positions[0]; i++)
    {
        AssignSubArray(CurrSolution->Sequence, Positions[i], Positions[i+1], Route[i]);
    }

    CARPInd Neighbor;
    int MovedTasks[3];
    double PenPrmt = LAMBDA;

    MovedTasks[0] = 2;
    BestDINeighbor->Fitness = INF;

    Bestu = -1;


    for (i = 1; i < Positions[0]; i++)
    {
        if (Route[i][0] < 4)
            continue;
        RID1 = i;
        for (j = 2; j < Route[i][0]-1; j++) // Route[i]: 0, t1, t2, t3, ..., tn, 0
        {
            for (u = 0; u < Positions[0]; u++)
            {
                if (u == i)
                    continue;

                RID2 = u;

                if ( u==0 && Route[0][0] < vehicle_num && Route[i][0] > 4)
                {
                    // new routes
                    // AssignArray(CurrSolution->Loads, Neighbor.Loads);
                    memcpy(Neighbor.Loads, CurrSolution->Loads, sizeof(CurrSolution->Loads));
                    Neighbor.Loads[i] -= (inst_tasks[Route[i][j]].demand + inst_tasks[Route[i][j+1]].demand);
                    Neighbor.TotalVioLoad = CurrSolution->TotalVioLoad;

                    if (Neighbor.Loads[i] > capacity)
                    {
                        Neighbor.TotalVioLoad -= (inst_tasks[Route[i][j]].demand + inst_tasks[Route[i][j+1]].demand);
                    } else if ( Neighbor.Loads[i] > (capacity - inst_tasks[Route[i][j]].demand - inst_tasks[Route[i][j+1]].demand))
                    {
                        Neighbor.TotalVioLoad -= (Neighbor.Loads[i] - capacity + inst_tasks[Route[i][j]].demand + inst_tasks[Route[i][j+1]].demand);
                    }

                    Neighbor.Loads[0] ++;
                    Neighbor.Loads[Neighbor.Loads[0]] = inst_tasks[Route[i][j]].demand + inst_tasks[Route[i][j+1]].demand;

                    if (Neighbor.Loads[Neighbor.Loads[0]] > capacity)
                    {
                        Neighbor.TotalVioLoad += (Neighbor.Loads[Neighbor.Loads[0]] -capacity);
                    }

                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                                                 - min_cost[inst_tasks[Route[i][j+1]].tail_node][inst_tasks[Route[i][j+2]].head_node]
                                                                 + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j+2]].head_node]
                                                                 + min_cost[DEPOT][inst_tasks[Route[i][j]].head_node]
                                                                 + min_cost[inst_tasks[Route[i][j+1]].tail_node][DEPOT];
                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;

                    if (Neighbor.Fitness < BestDINeighbor->Fitness)
                    {
                        MovedTasks[1] = Route[i][j];
                        MovedTasks[2] = Route[i][j+1];
                        Besti = i;
                        Bestj = j;
                        Bestu = u;
                        Bestv = 0;

                        BestDINeighbor->TotalCost = Neighbor.TotalCost;
                        // AssignArray(Neighbor.Loads, BestDINeighbor->Loads);
                        memcpy(BestDINeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                        BestDINeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                        BestDINeighbor->Fitness = Neighbor.Fitness;
                    }

                }


                if (u == 0)
                    continue;

                for (v = 2; v <= Route[u][0]; v++)
                {
                    
                    if (u == i && v == j )
                        continue;
                    // AssignArray(CurrSolution->Loads, Neighbor.Loads);
                    memcpy(Neighbor.Loads, CurrSolution->Loads, sizeof(CurrSolution->Loads));
                    Neighbor.Loads[i] -= (inst_tasks[Route[i][j]].demand + inst_tasks[Route[i][j+1]].demand);
                    Neighbor.Loads[u] += (inst_tasks[Route[i][j]].demand + inst_tasks[Route[i][j+1]].demand);
                    Neighbor.TotalVioLoad = CurrSolution->TotalVioLoad;

                    if ( Neighbor.Loads[i] > capacity)
                    {
                        Neighbor.TotalVioLoad -= (inst_tasks[Route[i][j]].demand + inst_tasks[Route[i][j+1]].demand);
                    } else if (Neighbor.Loads[i] > capacity - (inst_tasks[Route[i][j]].demand + inst_tasks[Route[i][j+1]].demand))
                    {
                        Neighbor.TotalVioLoad -= (Neighbor.Loads[i] - capacity + inst_tasks[Route[i][j]].demand + inst_tasks[Route[i][j+1]].demand);
                    }

                    if (Neighbor.Loads[u] > capacity + (inst_tasks[Route[i][j]].demand + inst_tasks[Route[i][j+1]].demand))
                    {
                        Neighbor.TotalVioLoad += inst_tasks[Route[i][j]].demand + inst_tasks[Route[i][j+1]].demand;
                    } else if (Neighbor.Loads[u] > capacity)
                    {
                        Neighbor.TotalVioLoad += (Neighbor.Loads[u] - capacity);
                    }
                    if (Neighbor.Loads[i] == 0)
                    {
                        RID1 = 0;
                        delete_element(Neighbor.Loads, i);
                    }
                    
                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                         - min_cost[inst_tasks[Route[i][j+1]].tail_node][inst_tasks[Route[i][j+2]].head_node]
                                         + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j+2]].head_node]
                                         - min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                                         + min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                         + min_cost[inst_tasks[Route[i][j+1]].tail_node][inst_tasks[Route[u][v]].head_node];
                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;
                    if (Neighbor.Fitness < BestDINeighbor->Fitness)
                    {
                        MovedTasks[1] = Route[i][j];
                        MovedTasks[2] = Route[i][j+1];
                        Besti = i;
                        Bestj = j;
                        Bestu = u;
                        Bestv = v;

                        BestDINeighbor->TotalCost = Neighbor.TotalCost;
                        // AssignArray(Neighbor.Loads, BestDINeighbor->Loads);
                        memcpy(BestDINeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                        BestDINeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                        BestDINeighbor->Fitness = Neighbor.Fitness;
                    }

                    w = inst_tasks[Route[i][j]].inverse;
                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                         - min_cost[inst_tasks[Route[i][j+1]].tail_node][inst_tasks[Route[i][j+2]].head_node]
                                         - min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                                         + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j+2]].head_node]
                                         - min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                                         + min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[w].head_node]
                                         + min_cost[inst_tasks[w].tail_node][inst_tasks[Route[i][j+1]].head_node]
                                         + min_cost[inst_tasks[Route[i][j+1]].tail_node][inst_tasks[Route[u][v]].head_node];
                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;
                    if (Neighbor.Fitness < BestDINeighbor->Fitness)
                    {
                        MovedTasks[1] = w;
                        MovedTasks[2] = Route[i][j+1];
                        Besti = i;
                        Bestj = j;
                        Bestu = u;
                        Bestv = v;

                        BestDINeighbor->TotalCost = Neighbor.TotalCost;
                        // AssignArray(Neighbor.Loads, BestDINeighbor->Loads);
                        memcpy(BestDINeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                        BestDINeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                        BestDINeighbor->Fitness = Neighbor.Fitness;
                    }

                    z = inst_tasks[Route[i][j+1]].inverse;
                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                         - min_cost[inst_tasks[Route[i][j+1]].tail_node][inst_tasks[Route[i][j+2]].head_node]
                                         - min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                                         + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j+2]].head_node]
                                         - min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                                         + min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                         + min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[z].head_node]
                                         + min_cost[inst_tasks[z].tail_node][inst_tasks[Route[u][v]].head_node];
                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;

                    if (Neighbor.Fitness < BestDINeighbor->Fitness)
                    {
                        MovedTasks[1] = Route[i][j];
                        MovedTasks[2] = z;
                        Besti = i;
                        Bestj = j;
                        Bestu = u;
                        Bestv = v;

                        BestDINeighbor->TotalCost = Neighbor.TotalCost;
                        // AssignArray(Neighbor.Loads, BestDINeighbor->Loads);
                        memcpy(BestDINeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                        BestDINeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                        BestDINeighbor->Fitness = Neighbor.Fitness;
                    }

                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                         - min_cost[inst_tasks[Route[i][j+1]].tail_node][inst_tasks[Route[i][j+2]].head_node]
                                         - min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                                         + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j+2]].head_node]
                                         - min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                                         + min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[w].head_node]
                                         + min_cost[inst_tasks[w].tail_node][inst_tasks[z].head_node]
                                         + min_cost[inst_tasks[z].tail_node][inst_tasks[Route[u][v]].head_node];
                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;

                        if (Neighbor.Fitness < BestDINeighbor->Fitness)
                        {
                            MovedTasks[1] = w;
                            MovedTasks[2] = z;
                            Besti = i;
                            Bestj = j;
                            Bestu = u;
                            Bestv = v;

                            BestDINeighbor->TotalCost = Neighbor.TotalCost;
                            // AssignArray(Neighbor.Loads, BestDINeighbor->Loads);
                            memcpy(BestDINeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                            BestDINeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                            BestDINeighbor->Fitness = Neighbor.Fitness;
                        }
                }
            }
        }
    }

    if (Bestu < 0)
        return count;
    
    if (Bestu == 0)
    {
        BestDINeighbor->Sequence[0] = 1;
        for (k = 1; k < Positions[0]; k++)
        {
            if ( k == Besti )
            {
                delete_element(Route[k], Bestj+1);
                delete_element(Route[k], Bestj);
                BestDINeighbor->Sequence[0]--;
                JoinArray(BestDINeighbor->Sequence, Route[k]);
            } else {
                BestDINeighbor->Sequence[0]--;
                JoinArray(BestDINeighbor->Sequence, Route[k]);
            }
        }
        BestDINeighbor->Sequence[0] ++;
        BestDINeighbor->Sequence[BestDINeighbor->Sequence[0]] = MovedTasks[1];
        BestDINeighbor->Sequence[0] ++;
        BestDINeighbor->Sequence[BestDINeighbor->Sequence[0]] = MovedTasks[2];
        BestDINeighbor->Sequence[0] ++;
        BestDINeighbor->Sequence[BestDINeighbor->Sequence[0]] = 0;
    } else{

        BestDINeighbor->Sequence[0] = 1;
        for (k = 1; k < Positions[0]; k++)
        {
            if (k == Besti)
            {
                delete_element(Route[k], Bestj+1);
                delete_element(Route[k], Bestj);
                if (Route[k][0] == 2)
                    continue;

                BestDINeighbor->Sequence[0]--;
                JoinArray(BestDINeighbor->Sequence, Route[k]);
            } else if (k == Bestu)
            {
                add_element(Route[k], MovedTasks[2], Bestv);
                add_element(Route[k], MovedTasks[1], Bestv);
                BestDINeighbor->Sequence[0] --;
                JoinArray(BestDINeighbor->Sequence, Route[k]);
            } else{
                BestDINeighbor->Sequence[0] --;
                JoinArray(BestDINeighbor->Sequence, Route[k]);
            }
        }
    }

    return count;

}


int SWAP(CARPInd *BestSWAPNeighbor, CARPInd *CurrSolution, const Task *inst_tasks)
{
    int i, j, k, u, v, z, w;

    int count = 0;  // record neighbourhood moves (number of fitness evaluation)
    int Besti, Bestj, Bestu, Bestv, RID1, RID2;
    int Positions[101], Route[101][MAX_TASK_SEQ_LENGTH];
    find_ele_positions(Positions, CurrSolution->Sequence, 0);

    Route[0][0] = Positions[0] - 1; // PA*
    for(i = 1; i < Positions[0]; i++)
    {
        AssignSubArray(CurrSolution->Sequence, Positions[i], Positions[i+1], Route[i]);
    }

    CARPInd Neighbor;
    double PenPrmt=LAMBDA;
    int MovedTasks[3];

    MovedTasks[0] = 2;
    BestSWAPNeighbor->Fitness = INF;
    Bestu = -1;

    for (i = 1;  i < Route[0][0]; i++)
    {
        RID1 = i;
        for (j = 2; j < Route[i][0]; j++) // Route[i]: 0, t1, t2, t3, ..., tn, 0
        {
            for (u = i+1; u <= Route[0][0]; u++)
            {
                RID2 = u;
                for (v = 2; v < Route[u][0]; v++)
                {
                    // AssignArray(CurrSolution->Loads, Neighbor.Loads);
                    memcpy(Neighbor.Loads, CurrSolution->Loads, sizeof(CurrSolution->Loads));
                    Neighbor.Loads[i] -= inst_tasks[Route[i][j]].demand - inst_tasks[Route[u][v]].demand;
                    Neighbor.Loads[u] += inst_tasks[Route[i][j]].demand - inst_tasks[Route[u][v]].demand;
                    Neighbor.TotalVioLoad = CurrSolution->TotalVioLoad;

                    if (inst_tasks[Route[i][j]].demand > inst_tasks[Route[u][v]].demand)
                    {
                        if (Neighbor.Loads[i] > capacity)
                        {
                            Neighbor.TotalVioLoad -= inst_tasks[Route[i][j]].demand - inst_tasks[Route[u][v]].demand;
                        } else if (Neighbor.Loads[i] > capacity - inst_tasks[Route[i][j]].demand + inst_tasks[Route[u][v]].demand)
                        {
                            Neighbor.TotalVioLoad -= Neighbor.Loads[i] + inst_tasks[Route[i][j]].demand - inst_tasks[Route[u][v]].demand - capacity;
                        }

                        if (Neighbor.Loads[u] > capacity + inst_tasks[Route[i][j]].demand - inst_tasks[Route[u][v]].demand)
                        {
                            Neighbor.TotalVioLoad += inst_tasks[Route[i][j]].demand - inst_tasks[Route[u][v]].demand;
                        } else if (Neighbor.Loads[u] > capacity)
                        {
                            Neighbor.TotalVioLoad += Neighbor.Loads[u] - capacity;
                        }
                    } else{
                        if (Neighbor.Loads[u] > capacity)
                        {
                            Neighbor.TotalVioLoad -= inst_tasks[Route[u][v]].demand - inst_tasks[Route[i][j]].demand;
                        }
                        else if (Neighbor.Loads[u] > capacity - inst_tasks[Route[u][v]].demand + inst_tasks[Route[i][j]].demand)
                        {
                            Neighbor.TotalVioLoad -= Neighbor.Loads[u] + inst_tasks[Route[u][v]].demand - inst_tasks[Route[i][j]].demand - capacity;
                        }

                        if (Neighbor.Loads[i] > capacity + inst_tasks[Route[u][v]].demand - inst_tasks[Route[i][j]].demand)
                        {
                            Neighbor.TotalVioLoad += inst_tasks[Route[u][v]].demand - inst_tasks[Route[i][j]].demand;
                        }
                        else if (Neighbor.Loads[i] > capacity)
                        {
                            Neighbor.TotalVioLoad += Neighbor.Loads[i] - capacity;
                        }
                    }
                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                            - min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                            - min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                            - min_cost[inst_tasks[Route[u][v]].tail_node][inst_tasks[Route[u][v+1]].head_node]
                            + min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                            + min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[u][v+1]].head_node]
                            + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                            + min_cost[inst_tasks[Route[u][v]].tail_node][inst_tasks[Route[i][j+1]].head_node];

                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;
                    if (Neighbor.Fitness < BestSWAPNeighbor->Fitness)
                    {
                        MovedTasks[1] = Route[i][j];
                        MovedTasks[2] = Route[u][v];
                        Besti = i;
                        Bestj = j;
                        Bestu = u;
                        Bestv = v;

                        BestSWAPNeighbor->TotalCost = Neighbor.TotalCost;
                        // AssignArray(Neighbor.Loads, BestSWAPNeighbor->Loads);
                        memcpy(BestSWAPNeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                        BestSWAPNeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                        BestSWAPNeighbor->Fitness = Neighbor.Fitness;
                    }

                    w = inst_tasks[Route[i][j]].inverse;
                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                         - min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                                         - min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                                         - min_cost[inst_tasks[Route[u][v]].tail_node][inst_tasks[Route[u][v+1]].head_node]
                                         + min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[w].head_node]
                                         + min_cost[inst_tasks[w].tail_node][inst_tasks[Route[u][v+1]].head_node]
                                         + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                                         + min_cost[inst_tasks[Route[u][v]].tail_node][inst_tasks[Route[i][j+1]].head_node];

                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;
                    if (Neighbor.Fitness < BestSWAPNeighbor->Fitness)
                    {
                        MovedTasks[1] = w;
                        MovedTasks[2] = Route[u][v];
                        Besti = i;
                        Bestj = j;
                        Bestu = u;
                        Bestv = v;

                        BestSWAPNeighbor->TotalCost = Neighbor.TotalCost;
                        // AssignArray(Neighbor.Loads, BestSWAPNeighbor->Loads);
                        memcpy(BestSWAPNeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                        BestSWAPNeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                        BestSWAPNeighbor->Fitness = Neighbor.Fitness;
                    }

                    z = inst_tasks[Route[u][v]].inverse;
                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                         - min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                                         - min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                                         - min_cost[inst_tasks[Route[u][v]].tail_node][inst_tasks[Route[u][v+1]].head_node]
                                         + min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                         + min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[u][v+1]].head_node]
                                         + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[z].head_node]
                                         + min_cost[inst_tasks[z].tail_node][inst_tasks[Route[i][j+1]].head_node];

                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;
                    if (Neighbor.Fitness < BestSWAPNeighbor->Fitness)
                    {
                        MovedTasks[1] = Route[i][j];
                        MovedTasks[2] = z;
                        Besti = i;
                        Bestj = j;
                        Bestu = u;
                        Bestv = v;

                        BestSWAPNeighbor->TotalCost = Neighbor.TotalCost;
                        // AssignArray(Neighbor.Loads, BestSWAPNeighbor->Loads);
                        memcpy(BestSWAPNeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                        BestSWAPNeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                        BestSWAPNeighbor->Fitness = Neighbor.Fitness;
                    }

                    Neighbor.TotalCost = CurrSolution->TotalCost - min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[Route[i][j]].head_node]
                                         - min_cost[inst_tasks[Route[i][j]].tail_node][inst_tasks[Route[i][j+1]].head_node]
                                         - min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[Route[u][v]].head_node]
                                         - min_cost[inst_tasks[Route[u][v]].tail_node][inst_tasks[Route[u][v+1]].head_node]
                                         + min_cost[inst_tasks[Route[u][v-1]].tail_node][inst_tasks[w].head_node]
                                         + min_cost[inst_tasks[w].tail_node][inst_tasks[Route[u][v+1]].head_node]
                                         + min_cost[inst_tasks[Route[i][j-1]].tail_node][inst_tasks[z].head_node]
                                         + min_cost[inst_tasks[z].tail_node][inst_tasks[Route[i][j+1]].head_node];

                    count ++;
                    Neighbor.Fitness = Neighbor.TotalCost + PenPrmt * Neighbor.TotalVioLoad;

                    if (Neighbor.Fitness < BestSWAPNeighbor->Fitness)
                    {
                        MovedTasks[1] = w;
                        MovedTasks[2] = z;
                        Besti = i;
                        Bestj = j;
                        Bestu = u;
                        Bestv = v;

                        BestSWAPNeighbor->TotalCost = Neighbor.TotalCost;
                        // AssignArray(Neighbor.Loads, BestSWAPNeighbor->Loads);
                        memcpy(BestSWAPNeighbor->Loads, Neighbor.Loads, sizeof(CurrSolution->Loads));
                        BestSWAPNeighbor->TotalVioLoad = Neighbor.TotalVioLoad;
                        BestSWAPNeighbor->Fitness = Neighbor.Fitness;
                    }
                }
            }
        }
    }

    if(Bestu < 0)
    {
        return count;
    }
    BestSWAPNeighbor->Sequence[0] = 1;
    for (k = 1; k < Positions[0]; k++)
    {
        if (k == Besti)
        {
            delete_element(Route[k], Bestj);
            add_element(Route[k], MovedTasks[2], Bestj);
            BestSWAPNeighbor->Sequence[0] --;
            JoinArray(BestSWAPNeighbor->Sequence, Route[k]);
        } else if (k == Bestu)
        {
            delete_element(Route[k], Bestv);
            add_element(Route[k], MovedTasks[1], Bestv);
            BestSWAPNeighbor->Sequence[0] --;
            JoinArray(BestSWAPNeighbor->Sequence, Route[k]);
        } else{
            BestSWAPNeighbor->Sequence[0] --;
            JoinArray(BestSWAPNeighbor->Sequence, Route[k]);
        }
    }
    return count;

}



double exautive_single_insertion(CARPInd indi, CARPInd *local, int *count, const Task *inst_tasks)
{

    int i, j, k;
    lns_route curr_solution, next_solution;

    char type[30];

    strcpy(type, "Single Insertion");

    int Positions[101];
    find_ele_positions(Positions, indi.Sequence, 0);
    curr_solution.Route[0][0] = Positions[0] - 1;
    curr_solution.loads[0] = Positions[0] - 1;
    for (i = 1; i < Positions[0]; i++)
    {
        AssignSubArray(indi.Sequence, Positions[i], Positions[i+1], curr_solution.Route[i]); // Route[i]: 0 x x x x 0
        curr_solution.loads[i] = 0;
        for (j = Positions[i]+1; j < Positions[i+1]; j++)
        {
            curr_solution.loads[i] += inst_tasks[indi.Sequence[j]].demand;
        }
    }
    curr_solution.fitness = indi.TotalCost + LAMBDA * indi.TotalVioLoad;;
    curr_solution.total_vio_loads = indi.TotalVioLoad;
    curr_solution.total_cost = indi.TotalCost;

    int u, pos_u, trip_u, v, pos_v, trip_v, x, y, inv_u;
    int improve = 0;
    int global_improve = 0;
    // check_cost1(curr_solution, inst_tasks);
    
    int flag;
    *count = 0;

    double start_t, finish_t;
    start_t = clock();
    while (1)
    {
        /* code */
        
        improve = expand(&curr_solution, 0, inst_tasks);
        // std::cout << improve << ":" << curr_solution.total_cost << std::endl;
        
        if (improve < 0)
        {
            break;
        } else{
            (*count) ++;
        }
    }
    finish_t = clock();
    double duration = (double)(finish_t - start_t) / CLOCKS_PER_SEC;
    
    // record the totoal count, min, media, mean, max, skewness
    update_global_best_solution(&curr_solution, local, inst_tasks, type);
    
    return duration;

}

int expand(lns_route *curr_solution, int depth, const Task *inst_tasks)
{
    depth ++;
    if (depth > 5)
    {
        return -1;
    }

    int u, pos_u, trip_u, v, pos_v, trip_v, x, y, inv_u;

    int next_loads_u, next_loads_v, next_total_vio_loads, next_total_cost;

    int flag;
    int equal_move[100][8];
    equal_move[0][0] = 0;

    for (trip_u = 1; trip_u <= curr_solution->Route[0][0]; trip_u ++)
    {
        for (pos_u = 2; pos_u < curr_solution->Route[trip_u][0]; pos_u++)
        {
            u = curr_solution->Route[trip_u][pos_u];

            for (trip_v = trip_u; trip_v <= curr_solution->Route[0][0]; trip_v ++)
            {
                for (pos_v = pos_u + 1; pos_v < curr_solution->Route[trip_v][0]; pos_v ++)
                {

                    v = curr_solution->Route[trip_v][pos_v];

                    next_loads_u = curr_solution->loads[trip_u] - inst_tasks[u].demand;
                    next_loads_v = curr_solution->loads[trip_v] + inst_tasks[u].demand;

                    next_total_vio_loads = curr_solution->total_vio_loads;


                    if (curr_solution->loads[trip_u] > capacity)
                        next_total_vio_loads -= curr_solution->loads[trip_u] - capacity;
                    if (curr_solution->loads[trip_v] > capacity)
                        next_total_vio_loads -= curr_solution->loads[trip_v] - capacity;

                    if (next_loads_u > capacity)
                        next_total_vio_loads += next_loads_u - capacity;
                    if (next_loads_v > capacity)
                        next_total_vio_loads += next_loads_v - capacity;

                    if (next_total_vio_loads > 0)
                    {
                        continue;
                    }


                    
                    if (pos_v == 2) // first task
                    {
                        flag = 1;
                        next_total_cost = curr_solution->total_cost
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
                        next_total_cost = curr_solution->total_cost
                                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                                    - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                                    + min_cost[inst_tasks[v].tail_node][inst_tasks[u].head_node]
                                                    + min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];
                    }

                    // next_solution.fitness = next_solution.total_cost + LAMBDA * next_solution.total_vio_loads;
                    if (next_total_vio_loads == 0 && next_total_cost < curr_solution->total_cost)
                    {
                        // printf("-----\n");
                        // check_cost1(*curr_solution, inst_tasks);
                        curr_solution->total_cost = next_total_cost;
                        curr_solution->total_vio_loads = next_total_vio_loads;
                        curr_solution->loads[trip_u] = next_loads_u;
                        curr_solution->loads[trip_v] = next_loads_v;

                        if (flag) {
                            add_element(curr_solution->Route[trip_v], u, pos_v);
                        } else {
                            add_element(curr_solution->Route[trip_v], u, pos_v+1);
                        }
                        delete_element(curr_solution->Route[trip_u], pos_u);
                        // check_cost1(*curr_solution, inst_tasks);
                        return 1;
                    }
                    if (next_total_vio_loads == 0 && next_total_cost == curr_solution->total_cost)
                    {
                        equal_move[0][0] ++;
                        equal_move[equal_move[0][0]][0] = trip_u;
                        equal_move[equal_move[0][0]][1] = pos_u;
                        equal_move[equal_move[0][0]][2] = trip_v;
                        equal_move[equal_move[0][0]][3] = pos_v;
                        equal_move[equal_move[0][0]][4] = u;
                        equal_move[equal_move[0][0]][5] = flag;
                        equal_move[equal_move[0][0]][6] = inst_tasks[u].demand;
                        equal_move[equal_move[0][0]][7] = next_total_cost;
                        // printf("%d, %d\n", inst_tasks[inv_u].inverse, curr_solution->Route[trip_u][pos_u]);
                    }

                    int inv_u = inst_tasks[u].inverse;
                    if (pos_v == 2) // first task
                    {
                        flag = 1;
                        next_total_cost = curr_solution->total_cost
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
                        next_total_cost = curr_solution->total_cost
                                                    - min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[u].head_node]
                                                    - min_cost[inst_tasks[u].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                                    + min_cost[inst_tasks[curr_solution->Route[trip_u][pos_u-1]].tail_node][inst_tasks[curr_solution->Route[trip_u][pos_u+1]].head_node]
                                                    - min_cost[inst_tasks[v].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node]
                                                    + min_cost[inst_tasks[v].tail_node][inst_tasks[inv_u].head_node]
                                                    + min_cost[inst_tasks[inv_u].tail_node][inst_tasks[curr_solution->Route[trip_v][pos_v+1]].head_node];
                    }
                    // next_solution.fitness = next_solution.total_cost + LAMBDA * next_solution.total_vio_loads;

                    if (next_total_vio_loads == 0 && next_total_cost < curr_solution->total_cost)
                    {
                        // check_cost1(*curr_solution, inst_tasks);
                        curr_solution->total_cost = next_total_cost;
                        curr_solution->total_vio_loads = next_total_vio_loads;
                        curr_solution->loads[trip_u] = next_loads_u;
                        curr_solution->loads[trip_v] = next_loads_v;

                        if (flag) {
                            add_element(curr_solution->Route[trip_v], inv_u, pos_v);
                        } else {
                            add_element(curr_solution->Route[trip_v], inv_u, pos_v + 1);
                        }
                        delete_element(curr_solution->Route[trip_u], pos_u);
                        // check_cost1(*curr_solution, inst_tasks);
                        return 1;
                    }

                    if (next_total_vio_loads == 0 && next_total_cost == curr_solution->total_cost)
                    {
                        equal_move[0][0] ++;
                        equal_move[equal_move[0][0]][0] = trip_u;
                        equal_move[equal_move[0][0]][1] = pos_u;
                        equal_move[equal_move[0][0]][2] = trip_v;
                        equal_move[equal_move[0][0]][3] = pos_v;
                        equal_move[equal_move[0][0]][4] = inv_u;
                        equal_move[equal_move[0][0]][5] = flag;
                        equal_move[equal_move[0][0]][6] = inst_tasks[u].demand;
                        equal_move[equal_move[0][0]][7] = next_total_cost;
                        // printf("%d, %d\n", inst_tasks[inv_u].inverse, curr_solution->Route[trip_u][pos_u]);
                    }
                }
            }
        }
    }

    if (equal_move[0][0] > 0)
    {   
        if (equal_move[0][0] > 100)
        {
            printf("equal_move: %d\n", equal_move[0][0]);
            exit(0);
        }
        
        lns_route next_solution;
        lns_route best_solution;
        best_solution.total_cost = curr_solution->total_cost;

        int improve=0;
        for (int i = 1; i <= equal_move[0][0]; i++)
        {
            next_solution.total_cost = curr_solution->total_cost;
            next_solution.total_vio_loads = curr_solution->total_vio_loads;
            memcpy(next_solution.loads, curr_solution->loads, sizeof(curr_solution->loads));
            memcpy(next_solution.Route, curr_solution->Route, sizeof(curr_solution->Route));

            next_solution.loads[trip_u] -= equal_move[i][6];
            next_solution.loads[trip_v] += equal_move[i][6];

            if (equal_move[i][5]) {
                add_element(next_solution.Route[equal_move[i][2]], equal_move[i][4], equal_move[i][3]);
            } else {
                add_element(next_solution.Route[equal_move[i][2]], equal_move[i][4], equal_move[i][3] + 1);
            }
            delete_element(next_solution.Route[equal_move[i][0]], equal_move[i][1]);
            // check_cost1(*curr_solution, inst_tasks);
            // check_cost1(next_solution, inst_tasks);
            improve = expand(&next_solution, depth, inst_tasks);

            if (improve > 0 && next_solution.total_cost < best_solution.total_cost)
            {
                best_solution.total_cost = next_solution.total_cost;
                best_solution.total_vio_loads = next_solution.total_vio_loads;
                memcpy(best_solution.loads, next_solution.loads, sizeof(next_solution.loads));
                memcpy(best_solution.Route, next_solution.Route, sizeof(next_solution.Route));
            }
        }
        if (improve > 0)
        {
            curr_solution->total_cost = best_solution.total_cost;
            curr_solution->total_vio_loads = best_solution.total_vio_loads;
            memcpy(curr_solution->loads, best_solution.loads, sizeof(best_solution.loads));
            memcpy(curr_solution->Route, best_solution.Route, sizeof(best_solution.Route));
            return 1;
        } else{
            return -1;
        }
    } else{
        return -1;
    }
}

int check_cost1(lns_route curr_solution, const Task *inst_tasks)
{
    int i;
    int cost = 0, tmp_cost;
    for (i = 1; i <= curr_solution.Route[0][0]; i++)
    {
        tmp_cost = get_task_seq_total_cost(curr_solution.Route[i], inst_tasks);
        // printf("%d \t", tmp_cost);
        cost += tmp_cost;
    }
    // printf("\n");
    if (cost != curr_solution.total_cost)
    {
        printf("cost error. \n");
        return 0;
    }
    return 1;
}

void update_global_best_solution(lns_route *route_solution, CARPInd *global_solution, const Task *inst_tasks, const char *type)
{
    memset(global_solution->Sequence, 0, sizeof(int)*MAX_TASK_SEQ_LENGTH);
    memset(global_solution->Loads, 0, sizeof(int)*80);
    global_solution->Sequence[0] = 0;
    int i, j, load = 0;
    for (i = 1; i <= route_solution->Route[0][0]; i++)
    {
        if (route_solution->Route[i][0] <= 2) continue;
        for (j = 1; j < route_solution->Route[i][0]; j++)
        {
            if( inst_tasks[route_solution->Route[i][j]].vt > 0 && route_solution->Route[i][j-1] != 0)
            {
                global_solution->Sequence[0] ++;
                global_solution->Sequence[global_solution->Sequence[0]] = 0;
                global_solution->Sequence[0] ++;
                global_solution->Sequence[global_solution->Sequence[0]] = route_solution->Route[i][j];
                continue;
            }
            global_solution->Sequence[0] ++;
            global_solution->Sequence[global_solution->Sequence[0]] = route_solution->Route[i][j];
        }
    }
    global_solution->Sequence[0] ++;
    global_solution->Sequence[global_solution->Sequence[0]] = 0;


    global_solution->Loads[0] = -1;
    for (i = 1; i <= global_solution->Sequence[0]; i++)
    {
        if (global_solution->Sequence[i] == 0)
        {
            global_solution->Loads[0] ++;
            global_solution->Loads[global_solution->Loads[0]] = load;
            load = 0;            
            continue;
        }
        load += inst_tasks[global_solution->Sequence[i]].demand; 
    }


    global_solution->TotalCost = route_solution->total_cost;
    global_solution->TotalVioLoad = route_solution->total_vio_loads;
    int cost = get_task_seq_total_cost(global_solution->Sequence, inst_tasks);
    if (cost != global_solution->TotalCost)
    {
        for (int i=1; i <= route_solution->Route[0][0]; i++)
        {
            for (int j=1; j < route_solution->Route[i][0]; j++)
            {
                printf("%d ", route_solution->Route[i][j]);
            }
        }
        printf("\n");
        for (int i=1; i <= global_solution->Sequence[0]; i++)
        {
            printf("%d ", global_solution->Sequence[i]);
        }
        printf("\n");
        printf("616: cost error\n");
    }

    route_solution->Route[0][0] = 1;
    route_solution->Route[1][0] = 1;
    route_solution->Route[1][1] = 0;

    for (i = 2; i <= global_solution->Sequence[0]; i++)
    {
        if (global_solution->Sequence[i] == 0)
        {
            route_solution->Route[route_solution->Route[0][0]][0] ++;
            route_solution->Route[route_solution->Route[0][0]][route_solution->Route[route_solution->Route[0][0]][0]] = 0;
            route_solution->Route[0][0] ++;
            route_solution->Route[route_solution->Route[0][0]][0] = 1;
            route_solution->Route[route_solution->Route[0][0]][1] = 0;
            continue;
        }
        route_solution->Route[route_solution->Route[0][0]][0] ++;
        route_solution->Route[route_solution->Route[0][0]][route_solution->Route[route_solution->Route[0][0]][0]] = global_solution->Sequence[i];
    }
    route_solution->Route[0][0] --;
    memcpy(route_solution->loads, global_solution->Loads, sizeof(global_solution->Loads));
    
}