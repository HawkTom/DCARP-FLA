#include "fla.h"

LOptimal indiList[100000];

void auto_corrlation(CARPInd sol, const Task * inst_tasks, int instance)
{

    int i;

    // Convert to Route[][]
    int Positions[101], Route[101][MAX_TASK_SEQ_LENGTH];
    find_ele_positions(Positions, sol.Sequence, 0);

    memset(Route, 0, sizeof(Route));

    Route[0][0] = Positions[0] - 1;
    for(i = 1; i < Positions[0]; i++)
    {
        AssignSubArray(sol.Sequence, Positions[i], Positions[i+1], Route[i]);
    }

    char path[50];
    memset(path, 0, sizeof(path));
    // sprintf(path, "result/autoc/instance%d.bin", instance);
    sprintf(path, "result/autoc/%s/instance%d.bin", map,  instance);


    int step;
    long int max_steps = 1000000000; //1000000000
    // clock_t start_t, finish_t;
    // start_t = clock();
    unsigned short int curr_cost = sol.TotalCost;
    int change_cost = 0;
    FILE *fp;
    fp = fopen(path, "wb");
    for (step=1; step <= max_steps; step++)
    {
        if (step % 1000000 == 0)
        {
            std::cout << step << std::endl;
        }
        change_cost = random_walk(Route, sol.Loads, inst_tasks);
        
        check_Route_valid(Route, sol.Loads, curr_cost+change_cost, inst_tasks);

        // record the cost of each step
        curr_cost = curr_cost + change_cost;
        
        fwrite(&curr_cost, sizeof(unsigned short int), 1, fp);
    }
    fclose(fp);
    // finish_t = clock();
    // double duration = (double)(finish_t - start_t) / CLOCKS_PER_SEC;
    // std::cout << duration << std::endl;
    auto_correlation_analysis1(instance, 100);
    
}

int random_walk(int (*Route)[MAX_TASK_SEQ_LENGTH], int *Loads, const Task * inst_tasks)
{   
    int i;
    int trip_u, trip_v, pos_u, pos_v;
    // random select a (u, v, pos_u, pos_v)
    while (1)
    {
        trip_u = random_choose_index_from(1, Route[0][0]);
        pos_u = random_choose_index_from(2, Route[trip_u][0]-1);

        if (Route[trip_u][0] == 3)
        {
            trip_v = random_choose_index_from(1, Route[0][0]);
            while (trip_u == trip_v)
            {
                trip_v = random_choose_index_from(1, Route[0][0]);
            }
        } else{
            trip_v = random_choose_index_from(1, Route[0][0]+1);
        }

        if (trip_u == trip_v)
        {
            pos_v = random_choose_index_from(2, Route[trip_v][0]-1);
            while(pos_u == pos_v)
            {
                pos_v = random_choose_index_from(2, Route[trip_v][0]-1);
            }
        } else if(trip_v > Route[0][0])
        {
            pos_v = 2;
        } else{
            pos_v = random_choose_index_from(2, Route[trip_v][0]);
        }

        if (inst_tasks[Route[trip_u][pos_u]].vt > 0)
        {
            continue;
        }

        // if not feasible, randomly select (u, v, pos_u, pos_v) again
        if ((Loads[trip_v] + inst_tasks[Route[trip_u][pos_u]].demand) < capacity)
        {
            break; // if feasible, move it
        }
    }

    if (pos_v == 2 && inst_tasks[Route[trip_v][2]].vt > 0)
    {
        trip_v = Route[0][0] + 1;
    }

    int cost0 = get_task_seq_total_cost(Route[trip_u], inst_tasks) + get_task_seq_total_cost(Route[trip_v], inst_tasks);
    
    int move_task = Route[trip_u][pos_u];

    Loads[trip_v] += inst_tasks[move_task].demand;
    Loads[trip_u] -= inst_tasks[move_task].demand;
    if (trip_u == trip_v)
    {
        delete_element(Route[trip_u], pos_u);
        add_element(Route[trip_v], move_task, pos_v);
    } else {
        if (trip_v > Route[0][0])
        {
            Loads[0]++;
            Loads[Loads[0]] = inst_tasks[move_task].demand;
            Route[0][0] ++;
            Route[Route[0][0]][0] = 3;
            Route[Route[0][0]][1] = 0;
            Route[Route[0][0]][2] = move_task;
            Route[Route[0][0]][3] = 0;
        } else{
            add_element(Route[trip_v], move_task, pos_v);
        }
        delete_element(Route[trip_u], pos_u);       
    }

    // if (trip_v > Route[0][0] || (pos_v == 2 && inst_tasks[Route[trip_v][2]].vt > 0))
    // {
    //     Loads[0]++;
    //     Loads[Loads[0]] = inst_tasks[move_task].demand;
    //     Route[0][0] ++;
    //     Route[Route[0][0]][0] = 3;
    //     Route[Route[0][0]][1] = 0;
    //     Route[Route[0][0]][2] = move_task;
    //     Route[Route[0][0]][3] = 0;
    //     delete_element(Route[trip_u], pos_u);
    //     Loads[trip_u] -= inst_tasks[move_task].demand;
    // } else if(trip_u == trip_v)
    // {
    //     delete_element(Route[trip_u], pos_u);
    //     add_element(Route[trip_v], move_task, pos_v);
    // } else {
    //     add_element(Route[trip_v], move_task, pos_v);
    //     delete_element(Route[trip_u], pos_u);
    //     Loads[trip_v] += inst_tasks[move_task].demand;
    //     Loads[trip_u] -= inst_tasks[move_task].demand;
    // }

    int cost1 = get_task_seq_total_cost(Route[trip_u], inst_tasks) + get_task_seq_total_cost(Route[trip_v], inst_tasks);

    if (trip_u == trip_v)
    {
        cost1/=2;
        cost0/=2;
    }

    // delete route with task num equal to 0
    int j, k;
    for (i = 1; i <= Route[0][0]; i++)
    {
        if (Route[i][0] < 3)
        {
            for (j=i+1; j <= Route[0][0]; j++)
            {   
                memset(Route[j-1], 0, sizeof(MAX_TASK_SEQ_LENGTH));
                for (k=0; k<= Route[j][0];k++)
                {
                    Route[j-1][k] = Route[j][k];
                }
                Loads[j-1] = Loads[j];
            }
            memset(Route[j-1], 0, sizeof(MAX_TASK_SEQ_LENGTH));
            Loads[j-1] = 0;
            Loads[0] --;
            Route[0][0]--;
        }
    }


    return (cost1 - cost0);
    // printf("trip_u:%d, pos_u:%d, trip_v:%d, pos_v:%d\n", trip_u, pos_u, trip_v, pos_v);
}

int random_choose_index_from(int lower, int upper)
{
    int x = rand();
    int k = x%(upper+1-lower);

    // printf("%d %d %d\n", x, num, k);

    return (k+lower); 
}

double auto_correlation_analysis1(int instance, int max_tau)
{
    char path[50];
    memset(path, 0, sizeof(path));
    sprintf(path, "result/autoc/%s/instance%d.bin", map,  instance);

    std::cout << path << std::endl;
    clock_t start_t, finish_t;
    start_t = clock();

    FILE *fp;
    fp = fopen(path, "rb");

    unsigned short int curr_cost = 0;
    double sum=0; 
    long int i=0;
    while ( fread(&curr_cost, sizeof(unsigned short int), 1, fp) == 1 )
    {
        sum += curr_cost;
        i++;
        if (curr_cost > 50000)
        {
            fclose(fp);
            std::cout << "error " << i << "," << curr_cost << std::endl;
            exit(0);
        }
    }
    double avr = sum/i;
    rewind(fp);
    i = 0;

    int arr_size = 10000;

    std::cout << "No problem " << sum << "," << avr << std::endl;

    max_tau += 1;
    double square_sum=0, y_sum[max_tau];
    
    unsigned short int curr_cost_arr[arr_size+max_tau];
    unsigned short int *p1, *p2;

    memset(y_sum, 0, sizeof(y_sum));

    int j, k;
    while ( fread(&curr_cost, sizeof(unsigned short int), 1, fp) == 1 )
    {
        // record the cost of each step
        curr_cost_arr[i] = curr_cost;
        i++;
        if(i <= arr_size)
        {
            square_sum += (curr_cost - avr) * (curr_cost - avr); 
        }
        if (i == arr_size+max_tau)
        {
            p1 = curr_cost_arr;
            for (k=0; k < max_tau; k++)
            {
                p2 = curr_cost_arr + k;
                for (j=0; j < arr_size; j++)
                {
                    y_sum[k] += (*(p1+j) - avr) * (*(p2+j) - avr);
                    // std::cout << *(p1+j) << "," << *(p2+j) << std::endl;
                }
            }
            fseek(fp, -max_tau*sizeof(unsigned short int), 1);
            i = 0;
        }
    }
    p1 = curr_cost_arr;
    for (k=0; k < max_tau; k++)
    {
        p2 = curr_cost_arr+k;
        for (j=0; j < i-k; j++)
        {
            y_sum[k] += (*(p1+j) - avr) * (*(p2+j) - avr);
            // std::cout << *(p1+j) << "," << *(p2+j) << std::endl;
        }
    }
    fclose(fp);

    memset(path, 0, sizeof(path));
    sprintf(path, "analysis/autoc/%s/instance%d.txt", map, instance);
    fp = fopen(path, "w");


    double autoc[max_tau];
    for (k=0; k<max_tau;k++)
    {
        autoc[k] = y_sum[k]/square_sum;
        // std::cout << "auto correlation:" << k << "," << y_sum[k] << "," << square_sum << ","<< autoc[k] << std::endl;
        fprintf(fp, "%.4f;", autoc[k]);
    }
    double corr_length = -1 / log(autoc[1]);
    // fseek(fp, -1*sizeof(char), 1);
    fprintf(fp, "%.4f\n", corr_length);
    fclose(fp);
    std::cout << "auto correlation:" << y_sum << "," << square_sum << ","<<autoc << std::endl;
    
    memset(path, 0, sizeof(path));
    sprintf(path, "result/autoc/%s/instance%d.bin", map,  instance);
    std::remove(path);


    finish_t = clock();
    double duration = (double)(finish_t - start_t) / CLOCKS_PER_SEC;
    std::cout << duration << "," <<i<< std::endl;
    
    return corr_length;

}


void local_optimum_related(const Task * inst_tasks, int instance)
{
    int runs = 100000;
    double start_t, finish_t;
    start_t = omp_get_wtime();    

    int i, dis, flag;

    #pragma omp parallel for num_threads(110)
    for (i=0; i<runs; i++)
    {
        CARPInd indi, local;
        double t;
        int count;
        rand_scanning(&indi, inst_tasks);
        // t = exautive_single_insertion(indi, &local, &count, inst_tasks);
        t = HyLS(indi, &local, &count, inst_tasks);
        check_solution_valid(local, inst_tasks);
        indiList[i].cost = local.TotalCost;
        memcpy(indiList[i].sequence, local.Sequence, sizeof(local.Sequence));
        indiList[i].hitnum = 1;
        indiList[i].count = count;
        indiList[i].usedtime = t;
        // std::cout << i << std::endl;
    }
    char path[50];
    memset(path, 0, sizeof(path));
    sprintf(path, "result/lor/%s/instance%d.bin", map, instance);
    printf("%s\n", path);
    FILE *fp;
    fp = fopen(path, "wb");

    for (i=0; i<runs; i++)
    {
        fwrite(&indiList[i], sizeof(LOptimal), 1, fp);
    }
    fclose(fp);

    finish_t = omp_get_wtime(); ;
    double duration = finish_t - start_t;
    std::cout << duration << std::endl;
}

void lo_analysis(int instance, const Task * inst_tasks)
{
    char path[50];
    memset(path, 0, sizeof(path));
    sprintf(path, "/home/hao/DCARP/data/lorHy/%s/instance%d.bin", map, instance);
    // sprintf(path, "result/lor/%s/instance%d.bin", map, instance);
    FILE *fp;
    fp = fopen(path, "rb");
    std::map<size_t, LOptimal> LOS;
    size_t id;
    LOptimal local;
    int flag;
    while( fread(&local, sizeof(LOptimal), 1, fp) == 1 )
    {
        size_t id = solution_hash_id(local.sequence, inst_tasks);
        auto lo = LOS.find(id);
        if (lo != LOS.end())
        {
            // already have
            int dis = hamming_distance(local.sequence, lo->second.sequence, inst_tasks);
            if (dis == 0)
            {
                flag = 1; // already exist
            } else{
                flag = 0; // not exist
            }
        } else{
            flag = 0; // not exist
        }
        if (flag)
        {
            lo->second.hitnum ++; // Record the number of times that each local optimal is hit
            lo->second.usedtime += local.usedtime;
        } else{
            local.hitnum = 1;
            LOS.emplace(id, local);
        }
    }
    fclose(fp);

    

    int total_num_lo = LOS.size();

    std::map<size_t, LOptimal>::iterator iter;
    int min_value=INF, max_value = 0;
    for (iter = LOS.begin(); iter!=LOS.end(); iter++)
    {
        if (iter->second.cost < min_value)
        {
            min_value = iter->second.cost;
        }
        if (iter->second.cost > max_value)
        {
            max_value = iter->second.cost;
        }
    }
    int bins = 100;
    double interval = (max_value - min_value)/bins;
    int offset = 0;
    while(true)
    {   
        offset ++;
        if (min_value + (bins+offset)*interval > max_value)
        {
            break;
        }
    }
    printf("instance: %d, offset: %d\n", instance, offset);

    int hist[bins+offset], hist_cost[bins+offset];
    int hit_hist[bins+offset];
    float time_hist[bins+offset];

    memset(hist_cost, 0, sizeof(hist_cost));
    for(int i = 0; i < bins+offset; i++)
    {
        hist_cost[i] = min_value + i*interval;
        if (hist_cost[i] > max_value)
        {
            break;
        }
    }

    memset(hist, 0, sizeof(hist));
    memset(hit_hist,0, sizeof(hit_hist));
    memset(time_hist,0, sizeof(time_hist));
    int idx;
    float avr_time = 0;
    for (iter = LOS.begin(); iter!=LOS.end(); iter++)
    {
        idx = (int)(iter->second.cost - min_value)/interval;
        hist[idx] ++; // local_optimal_num ~ cost : proportion

        hit_hist[idx] += iter->second.hitnum; // hit_num ~ cost : probability

        avr_time += iter->second.usedtime; 

        time_hist[idx] += (iter->second.usedtime/iter->second.hitnum);
    }
    avr_time /= 100000; // average time 
    for(int i = 0; i < bins+offset; i++)
    {
        if (hist[i] == 0)
        {
            continue;
        }
        time_hist[i] /= hist[i]; // cost ~ used time
    }

    // validation
    int sum1 = 0, sum2 = 0;
    for (int i=0; i<bins+offset; i++)
    {
        sum1+=hist[i];
        sum2+=hit_hist[i];
    }
    if (sum1 != total_num_lo)
    {
        
        std::cout << "hist error" << sum1 << "," << total_num_lo << std::endl;
    }
    if (sum2 != 100000)
    {
        std::cout << "hit_hist error" << sum1 << "," << total_num_lo << std::endl;
    }

    memset(path, 0, sizeof(path));
    sprintf(path, "result/best/%s/instance%d.txt", map, instance);
    
    fp = fopen(path, "r");
    int best;
    fscanf(fp, "%d\n", &best);
    int best_sequence[MAX_TASK_SEQ_LENGTH];
    int i = -1, value;
    fscanf(fp, "%d,", &best_sequence[0]);
    for (i=1; i <= best_sequence[0]; i++)
    {
        fscanf(fp, "%d,", &value);
        best_sequence[i] = value;
        // std::cout << best_sequence[i] << std::endl;
    }
    fclose(fp);

    memset(path, 0, sizeof(path));
    sprintf(path, "analysis/lor/%s/instance%d.bin", map, instance);
    fp = fopen(path, "wb");
    fwrite(&total_num_lo, sizeof(total_num_lo), 1, fp);
    
    unsigned short int F[100000];
    unsigned short int D[100000];
    i=0;
    int dis_tmp;
    // calculate fitness distance coefficients
    float avr_f = 0, avr_d=0;
    for (iter = LOS.begin(); iter!=LOS.end(); iter++)
    {
        dis_tmp = hamming_distance(iter->second.sequence, best_sequence, inst_tasks);
        iter->second.cost;
        F[i] = iter->second.cost;
        D[i] = dis_tmp;

        avr_f += F[i];
        avr_d += D[i];
        i++;
        for (int j=0; j < iter->second.hitnum; j++)
        {
            fwrite(&dis_tmp, sizeof(dis_tmp), 1, fp);
            fwrite(&iter->second.cost, sizeof(iter->second.cost), 1, fp);
            fwrite(&iter->second.usedtime, sizeof(iter->second.usedtime), 1, fp);
        }
    }
    
    avr_f /= i;
    avr_d /= i; 

    float CFD = 0, sigma_F = 0, sigma_D = 0;
    for (int j=0; j < i; j++)
    {
        CFD += (F[j]-avr_f) * (D[j] - avr_d);
        sigma_F += (F[j]-avr_f) * (F[j]-avr_f);
        sigma_D += (D[j]-avr_d) * (D[j]-avr_d);
    }
    CFD = CFD / i;
    sigma_F = sqrt(sigma_F/i);
    sigma_D = sqrt(sigma_D/i);
    float fdc = CFD / (sigma_F * sigma_D);

    fwrite(&fdc, sizeof(fdc), 1, fp);

    fclose(fp);


    memset(path, 0, sizeof(path));
    sprintf(path, "analysis/lor/%s/instance%d.txt", map, instance);
    fp = fopen(path, "w");

    fprintf(fp, "hist,");
    for(int i = 0; i < bins+offset; i++)
    {
        fprintf(fp, "%d,", hist[i]);
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "hithist,");
    for(int i = 0; i < bins+offset; i++)
    {
        fprintf(fp, "%d,", hit_hist[i]);
    }
    fprintf(fp, "\n");

    fprintf(fp, "costhist,");
    for(int i = 0; i < bins+offset; i++)
    {
        fprintf(fp, "%d,", hist_cost[i]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "timehist");
    for(int i = 0; i < bins+offset; i++)
    {
        fprintf(fp, "%.2f,", time_hist[i]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "avrtime,%.2f\nfdc,%.4f\n", avr_time, fdc);


    fclose(fp);
    std::cout << "end" << std::endl;

}

void lo_dis_analysis(int instance, const Task * inst_tasks)
{
    char path[50];
    memset(path, 0, sizeof(path));
    sprintf(path, "/home/hao/DCARP/data/lorHy/%s/instance%d.bin", map, instance);
    // sprintf(path, "result/lor/instance0.bin");
    FILE *fp;
    fp = fopen(path, "rb");
    std::map<size_t, LOptimal> LOS;
    size_t id;
    LOptimal local;
    int flag, dis;
    while( fread(&local, sizeof(LOptimal), 1, fp) == 1 )
    {
        size_t id = solution_hash_id(local.sequence, inst_tasks);
        auto lo = LOS.find(id);
        if (lo != LOS.end())
        {
            // already have
            dis = hamming_distance(local.sequence, lo->second.sequence, inst_tasks);
            if (dis == 0)
            {
                flag = 1; // already exist
            } else{
                flag = 0; // not exist
            }
        } else{
            flag = 0; // not exist
        }
        if (flag)
        {
            lo->second.hitnum ++; // Record the number of times that each local optimal is hit
            lo->second.usedtime += local.usedtime;
        } else{
            local.hitnum = 1;
            LOS.emplace(id, local);
        }
    }
    fclose(fp);

    

    int samples = 2000;
    CARPInd solution1, solution2;
    double sum_dis = 0;
    int x = 0;
    for (int i=1; i<=samples;i++)
    {
        rand_scanning(&solution1, inst_tasks);
        rand_scanning(&solution2, inst_tasks);
        x = hamming_distance(solution1.Sequence, solution2.Sequence, inst_tasks);
        sum_dis += x;
    }
    int avr_dis1 = int(sum_dis/samples);

    

    std::map<size_t, LOptimal>::iterator iter1, iter2;
    int max_num = LOS.size();
    int idx1, idx2;

    sum_dis = 0;
    for (int i=1; i <= samples; i++){
        iter1 = LOS.begin();
        iter2 = LOS.begin();
        idx1 = rand_choose(max_num-1);
        idx2 = rand_choose(max_num-1);
        while (idx2 == idx1)
        {
            idx2 = rand_choose(max_num-1);
        }
        std::advance(iter1, idx1);
        std::advance(iter2, idx2);
        x = hamming_distance(iter1->second.sequence, iter2->second.sequence, inst_tasks);
        sum_dis += x;
    }
    int avr_dis2 = int(sum_dis/samples);

    // printf("%d %d \n", avr_dis1, avr_dis2);

    std::vector<std::pair<size_t, int>> LOVector;
    for (iter1 = LOS.begin(); iter1!=LOS.end(); iter1++)
    {
        LOVector.push_back(std::make_pair(iter1->first, iter1->second.cost));
    }
    std::sort(LOVector.begin(), LOVector.end(), cmp);

    int k = 0;
    int avr_dis[10];
    int avr_cost[10];
    for (k=0; k <= 9; k++)
    {   
        sum_dis = 0;
        int sum_cost = 0;
        for (int i=1; i<=1000; i++)
        {
            idx1 = int(max_num*0.1)*k + rand_choose(int(max_num*0.1));
            idx2 = int(max_num*0.1)*k + rand_choose(int(max_num*0.1));
            while (idx2 == idx1)
            {
                idx2 = int(max_num*0.1)*k + rand_choose(int(max_num*0.1));
            }
            x = hamming_distance(LOS[LOVector[idx1].first].sequence, LOS[LOVector[idx2].first].sequence, inst_tasks);
            sum_dis += x;
            sum_cost += (LOS[LOVector[idx1].first].cost + LOS[LOVector[idx2].first].cost);
        }
        int tmp_dis = int(sum_dis/1000);
        int tmp_cost = int(sum_cost/2000);
        // printf("%d %d\n", k, tmp_dis);
        avr_dis[k] = tmp_dis;
        avr_cost[k] = tmp_cost;
    }

    memset(path, 0, sizeof(path));
    sprintf(path, "analysis/lor/%s/instance%d.txt", map, instance);
    // sprintf(path, "analysis/lor/instance0.txt");
    fp = fopen(path, "a");
    fprintf(fp, "dis,random,%d,optimum,%d\noptimumdis,", avr_dis1, avr_dis2);
    for (int i=0; i<=9; i++)
    {
        fprintf(fp, "%d,", avr_dis[i]);
    }
    fprintf(fp, "\noptimucost,");
    for (int i=0; i<=9; i++)
    {
        fprintf(fp, "%d,", avr_cost[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

int cmp(const std::pair<size_t, int>& x, const std::pair<size_t, int>& y)  
{  
    return x.second < y.second;  
}  

// save all local optimum and its information
size_t solution_hash_id(int *indi_seq, const Task* inst_tasks)
{
    int i, j;
    std::string edge;
    int head, tail;
    std::vector<int> arr;

    for (i=1; i<indi_seq[0]; i++)
    {
        head = inst_tasks[indi_seq[i]].tail_node_r;
        tail = inst_tasks[indi_seq[i+1]].head_node_r;

        edge = std::to_string(head) + "," + std::to_string(tail);
        arr.push_back(edgesIDs[edge]);
    }
    std::sort(arr.begin(), arr.end());
    
    std::ostringstream vts;
    if(!arr.empty())
    {
        std::copy(arr.begin(), arr.end()-1, std::ostream_iterator<int>(vts, ","));
        vts << arr.back();
    }
    
    size_t id;
    std::hash<std::string> h;
    id = h(vts.str());

    // std::cout << vts.str() << std::endl;
    // std::cout << id << std::endl;
    return id;
}

int hamming_distance(int *sol1_seq, int *sol2_seq, const Task* inst_tasks)
{
    int i, j;
    std::string edge;
    int head, tail;
    std::set<int> set1;
    std::set<int> set2;

    for (i=1; i<sol1_seq[0]; i++)
    {
        head = inst_tasks[sol1_seq[i]].tail_node_r;
        tail = inst_tasks[sol1_seq[i+1]].head_node_r;

        edge = std::to_string(head) + "," + std::to_string(tail);
        set1.emplace(edgesIDs[edge]);
    }
    for (i=1; i<sol2_seq[0]; i++)
    {
        head = inst_tasks[sol2_seq[i]].tail_node_r;
        tail = inst_tasks[sol2_seq[i+1]].head_node_r;

        edge = std::to_string(head) + "," + std::to_string(tail);
        set2.emplace(edgesIDs[edge]);
    }

    std::set<int> set_inter = {};
    std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(set_inter, set_inter.begin()));
    int dis = set1.size()+set2.size() - 2*set_inter.size();
    return dis;
}

void return_probability(int instance, const Task *inst_tasks)
{
    int samples = 1000;
    int distance = 3; // neighbourhood moves 3~50;
    
    char path[50];
    memset(path, 0, sizeof(path));
    // sprintf(path, "result/best/instance%d.txt", instance);
    sprintf(path, "result/best/%s/instance%d.txt", map, instance);
    
    FILE *fp = fopen(path, "r");
    int best;
    fscanf(fp, "%d\n", &best);
    int best_sequence[MAX_TASK_SEQ_LENGTH];
    int i = -1, value;
    fscanf(fp, "%d,", &best_sequence[0]);
    for (i=1; i <= best_sequence[0]; i++)
    {
        fscanf(fp, "%d,", &value);
        best_sequence[i] = value;
        // std::cout << best_sequence[i] << std::endl;
    }
    fclose(fp);

    CARPInd best_indi, local;
    memcpy(best_indi.Sequence, best_sequence, sizeof(best_sequence));
    best_indi.TotalCost = get_task_seq_total_cost(best_indi.Sequence, inst_tasks);
    check_solution_valid(best_indi, inst_tasks);
    int count;
    HyLS(best_indi, &local, &count, inst_tasks);
    if(local.TotalCost < best_indi.TotalCost)
    {
        memcpy(best_indi.Sequence, local.Sequence, sizeof(local.Sequence));
    }

    // Convert to Route[][]
    int Positions[101], bestRoute[101][MAX_TASK_SEQ_LENGTH];
    find_ele_positions(Positions, best_indi.Sequence, 0);

    memset(bestRoute, 0, sizeof(bestRoute));

    bestRoute[0][0] = Positions[0] - 1;
    for(i = 1; i < Positions[0]; i++)
    {
        AssignSubArray(best_indi.Sequence, Positions[i], Positions[i+1], bestRoute[i]);
    }

    int j, k;
    int best_loads[101], tmp_load;
    memset(best_loads, 0, sizeof(best_loads));
    for (i=1; i<= bestRoute[0][0]; i++)
    {
        tmp_load = 0;
        for (j=1; j<=bestRoute[i][0]; j++)
        {
            tmp_load += inst_tasks[bestRoute[i][j]].demand;
        }
        best_loads[0]++;
        best_loads[best_loads[0]] = tmp_load;
    }

    printf("%s-%d return probability. \n", map, instance);
    CARPInd indi;
    int Route[101][MAX_TASK_SEQ_LENGTH], loads[101];
    int return_samples[10];
    memset(return_samples, 0, sizeof(return_samples));
    for (i=0; i<10; i++)
    {   
        for(k=1; k <= samples; k++)
        {   
            // change another sample
            memcpy(Route, bestRoute, sizeof(bestRoute));
            memcpy(loads, best_loads, sizeof(best_loads));
            for (j=0; j<=i; j++)
            {
                random_walk(Route, loads, inst_tasks);
                // check_Route_valid(Route, loads, inst_tasks);
            }
            // Route -> Solution
            memset(indi.Sequence, 0, sizeof(indi.Sequence));
            memset(indi.Loads, 0, sizeof(indi.Loads));
            for (int ii=1; ii<=Route[0][0]; ii++)
            {
                tmp_load = 0;
                for (int jj=1; jj<Route[ii][0]; jj++)
                {
                    indi.Sequence[0]++;
                    indi.Sequence[indi.Sequence[0]] = Route[ii][jj];
                    tmp_load += inst_tasks[Route[ii][jj]].demand;
                }
                indi.Loads[0]++;
                indi.Loads[indi.Loads[0]] = tmp_load;
            }
            indi.Sequence[0]++;
            indi.Sequence[indi.Sequence[0]] = 0;
            indi.TotalCost = get_task_seq_total_cost(indi.Sequence, inst_tasks);
            // check_solution_valid(indi, inst_tasks);
            
            HyLS(indi, &local, &count, inst_tasks);
            int dis = hamming_distance(local.Sequence, best_indi.Sequence, inst_tasks);
            if (dis == 0)
            {
                return_samples[i] ++;
            }
        }
        printf("iteration: %d, %d\n", i, return_samples[i]);
        if (1.0*return_samples[i]/samples < 0.001)
        {
            break;
        }
    }
    memset(path, 0, sizeof(path));
    // sprintf(path, "result/returnp/instance%d.txt", instance);
    sprintf(path, "result/returnp/%s/instance%d.txt", map, instance);
    fp = fopen(path, "w");
    for (k=0; k<=i; k++)
    {
        fprintf(fp, "%d,", return_samples[k]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

void check_Route_valid(int (*Route)[MAX_TASK_SEQ_LENGTH], int *Loads, int checkCost, const Task *inst_task)
{
    int i, j;
    int flag = 0;
    int used[MAX_TASK_SEQ_LENGTH];
    memset(used, 0, sizeof(used));

    int load = 0;
    int cost = 0;
    for(i = 1; i <= Route[0][0]; i++)
    {
        cost += get_task_seq_total_cost(Route[i], inst_task);
        for (j = 2; j < Route[i][0]; j++)
        {
            if (Route[i][j] <= req_edge_num)
            {
                used[Route[i][j]] = 1;
            } else
            {
                used[inst_task[Route[i][j]].inverse] = 1;
            }
            load += inst_task[Route[i][j]].demand;
        }
        if (load != Loads[i])
        {
            printf("load error %d \t", i);
            flag = 1;
        }

        if (load > capacity)
        {
            printf("move is not feasible %d \t", i);
            flag = 1;
        }
        load = 0;

    }

    if (cost!=checkCost)
    {
        printf("cost error %d %d \t", cost, checkCost);
        flag = 1;
    }
    
    // printf("lack of: ");
    for (i = 1; i <= req_edge_num; i++)
    {
        if(used[i] == 0)
        {
            printf("not used %d \t", i);
            flag = 1;
        }
    }

    /** 
     * @brief check cost
     * 
     * */
    // if ( solution.TotalCost != get_task_seq_total_cost(solution.Sequence, inst_task))
    // {
    //     flag = 1;
    //     printf("solution's cost error\n");
    // }
    if (flag)
    {
        printf("\n");
        exit(0);
    }
}

void check_solution_valid(CARPInd solution, const Task *inst_task)
{
    int i, j;
    int used[MAX_TASK_SEQ_LENGTH];
    memset(used, 0, sizeof(used));
    int route_num = 0;
    for(i = 1; i <= solution.Sequence[0]; i++)
    {
        if (solution.Sequence[i] == 0)
        {
            route_num ++;
            continue;
        }
        if (solution.Sequence[i] <= req_edge_num)
        {
            used[solution.Sequence[i]] = 1;
        } else
        {
            used[inst_task[solution.Sequence[i]].inverse] = 1;
        }
    }
    int flag = 0;
    // if((route_num-1) != solution.Loads[0])
    // {
    //     printf("route num error\n");
    //     flag = 1;
    // }
    // printf("lack of: \n");
    for (i = 1; i <= req_edge_num; i++)
    {
        if(used[i] == 0)
        {
            printf("not used %d \t", i);
            flag = 1;
        }
    }
    // printf("\n");
    // get_each_load(solution.Sequence, inst_task);
    if ( solution.TotalCost != get_task_seq_total_cost(solution.Sequence, inst_task))
    {
        flag = 1;
        printf("solution's cost error\n");
    }
    if (flag)
    {
        printf("\n");
        exit(0);
    }
}