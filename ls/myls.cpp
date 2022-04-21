#include "ls.h"

double HyLS(CARPInd InitSolution, CARPInd *bestSolution, int *count, const Task *inst_tasks)
{

    CARPInd IndiList[AXVSIZE+1];
    int IndiListLength = 0;

    CARPInd solution;
    copy_individual(&solution, &InitSolution);

    double start_t, finish_t;
    start_t = clock();

    int best= INF;
    for (int i = 1; i < AXVSIZE; i++)
    {
        exhaustive_ls(&solution, IndiList, &IndiListLength, inst_tasks);
        if (solution.TotalCost < best)
        {
            best = solution.TotalCost;
            copy_individual(bestSolution, &solution);
            // printf("%d %d \n",i, best);
        }
        copy_individual(&solution, &IndiList[i]);
        if (i >= IndiListLength)
        {
            break;
        }
        // printf("%d %d \n",i, best);
    }
    // printf("%d \n", IndiListLength);
    
    finish_t = clock();
    *count = IndiListLength;
    double duration = (double)(finish_t - start_t) / CLOCKS_PER_SEC;
    // printf("HyLS: %d\n", best);
    return duration;
}

void exhaustive_ls(CARPInd *InitSolution, CARPInd *IndiList, int *IndiListLength, const Task *inst_tasks)
{

    CARPInd best_si_solution;
    CARPInd best_di_solution;
    CARPInd best_swap_solution;
    CARPInd best_rev_solution;
    CARPInd best_cross_solution;

    copy_individual(&best_si_solution, InitSolution);
    copy_individual(&best_di_solution, InitSolution);
    copy_individual(&best_swap_solution, InitSolution);
    copy_individual(&best_rev_solution, InitSolution);
    copy_individual(&best_cross_solution, InitSolution);

    int improves[6];
    int terminate_flag = 0;

    clock_t start_t, finish_t;
    start_t = clock();
    double duration = 0.0, duration1=0.0;
    
    while (1)
    {
        memset(improves, 0, sizeof(improves));
        // printf("Log: %d %d\n", InitSolution.TotalCost, additional_cost);

        improves[1] = client_single_insertion(InitSolution, &best_si_solution, inst_tasks);
        if (improves[1] && *IndiListLength <= AXVSIZE - 1)
        {
            (*IndiListLength) ++;
            copy_individual(&IndiList[(*IndiListLength)], &best_si_solution);
        }
        
            
        improves[2] = client_double_insertion(InitSolution, &best_di_solution, inst_tasks);
        if (improves[2]  && (*IndiListLength) <= AXVSIZE - 1)
        {
            (*IndiListLength) ++;
            copy_individual(&IndiList[(*IndiListLength)], &best_di_solution);
        }

        improves[3] = client_swap(InitSolution, &best_swap_solution, inst_tasks);
        if (improves[3]  && (*IndiListLength) <= AXVSIZE - 1)
        {
            (*IndiListLength) ++;
            copy_individual(&IndiList[(*IndiListLength)], &best_swap_solution);
        }


        improves[4] = client_reverse(InitSolution, &best_rev_solution, inst_tasks);
        if (improves[4]  && (*IndiListLength) <= AXVSIZE - 1)
        {
            (*IndiListLength) ++;
            copy_individual(&IndiList[(*IndiListLength)], &best_rev_solution);
        }

        improves[5] = client_cross(InitSolution, &best_cross_solution, inst_tasks);
        if (improves[5]  && (*IndiListLength) <= AXVSIZE - 1)
        {
            (*IndiListLength) ++;
            copy_individual(&IndiList[(*IndiListLength)], &best_cross_solution);
        }

        finish_t = clock();
        duration = (double)(finish_t - start_t) / CLOCKS_PER_SEC;
        terminate_flag = 0;
        for (int i = 1; i <= 5; i++)
        {
            // printf("%d\t", improves[i]);
            terminate_flag += improves[i];
        }
        // printf("\n");
        if (terminate_flag == 0 || duration > 10)
        {
            break;
        }


        int min_flag = 0, mincost=INF;
        if(best_si_solution.TotalCost < mincost)
        {
            mincost = best_si_solution.TotalCost;
            min_flag = 1;
        }
        if(best_di_solution.TotalCost < mincost)
        {
            mincost = best_di_solution.TotalCost;
            min_flag = 2;
        }
        if(best_swap_solution.TotalCost < mincost)
        {
            mincost = best_swap_solution.TotalCost;
            min_flag = 3;
        }
        if(best_rev_solution.TotalCost < mincost)
        {
            mincost = best_rev_solution.TotalCost;
            min_flag = 4;
        }
        if(best_cross_solution.TotalCost < mincost)
        {
            mincost = best_cross_solution.TotalCost;
            min_flag = 5;
        }

        switch (min_flag)
        {
        case 1:
            copy_individual(InitSolution, &best_si_solution);
            break;
        case 2:
            copy_individual(InitSolution, &best_di_solution);
            break;
        case 3:
            copy_individual(InitSolution, &best_swap_solution);
            break;
        case 4:
            copy_individual(InitSolution, &best_rev_solution);
            break;
        case 5:
            copy_individual(InitSolution, &best_cross_solution);
            break;
        default:
            break;
        }

        
        // printf("1 -> Log: %d %d\n", InitSolution->TotalCost, terminate_flag);
    }
}

void move_vt_to_first(CARPInd *dst, CARPInd *src, const Task *inst_tasks)
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
    dst->TotalCost = src->TotalCost;
}