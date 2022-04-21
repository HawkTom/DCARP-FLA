#include "globalvar.h"
#include "src/src.h"
#include "ls/ls.h"
#include "libs/cmdline.h"
#include "fla.h"

// #include <limits> std::cout << sizeof(unsigned short) << "," << sizeof(unsigned short) << "," << std::numeric_limits<unsigned short>::max();



void available_edges_map(const Task* inst_tasks);



void dcarp_instance(int instance);
void scarp_instance();



// void instance_test(int instance)
int main(int argc, char *argv[])
{

    
    
    cmdline::parser inst_info;
    inst_info.add<std::string>("map", 'm', "map name", true, "");
    inst_info.add<int>("ins", 'i', "instance index", true, 0, cmdline::range(0, 65535));
    inst_info.parse_check(argc, argv);
    int instance = inst_info.get<int>("ins");
    strcpy(map, inst_info.get<std::string>("map").c_str());

    // strcpy(map, "egl-e2-A-1");
    // int instance = 2;
    

    DEPOT = 1;
    dcarp_instance(instance);

    // auto_correlation_analysis1(instance, 100);
    // auto_correlation_analysis(instance, 1);

    // #pragma omp parallel for num_threads(4)
    // for (instance = 1; instance <= 1; instance++)
    // {
    //     std::cout << map << ":" << instance << " start." << std::endl;
    //     if (instance == 0)
    //     {
    //         scarp_instance();
    //     } else{
    //         dcarp_instance(instance);
    //     }
    //     std::cout << map << ":" << instance << " complete." << std::endl;
    //     // break;
    // }



    
    

    
    
    return 0;

}

void scarp_instance()
{
    Task inst_tasks[MAX_TASKS_TAG_LENGTH];
    Arc inst_arcs[MAX_ARCS_TAG_LENGTH];

    readMap(inst_tasks, inst_arcs, map);
    update_cost(inst_tasks, inst_arcs);
    mod_dijkstra();

    available_edges_map(inst_tasks);

    

    // CARPInd solution1, solution2;
    // double avr_dis = 0;
    // int x = 0;
    // for (int i=1; i<=5000;i++)
    // {
    //     rand_scanning(&solution1, inst_tasks);
    //     avr_dis += solution1.TotalCost;
    //     // rand_scanning(&solution2, inst_tasks_vt);
    //     // x = hamming_distance(solution1.Sequence, solution2.Sequence, inst_tasks_vt);
    //     // avr_dis += x;
    // }
    // int avr_cost = int(avr_dis/5000);
    // printf("%d\n", avr_cost);


    // CARPInd solution;
    // int serve_mark[MAX_TASK_TAG_LENGTH];
    // memset(serve_mark, 0, sizeof(serve_mark));
    // for (int i = 1; i <= task_num; i++)
    // {
    //     serve_mark[i] = 1;
    // }
    // path_scanning(&solution, inst_tasks, serve_mark);
    
    // auto_corrlation(solution, inst_tasks, 0);
    // local_optimum_related(inst_tasks, 0);
    // lo_analysis(0, inst_tasks);
    // lo_dis_analysis(0, inst_tasks);
    // return_probability(0, inst_tasks);
    // return_probability(0, inst_tasks);
}

void dcarp_instance(int instance)
{
    Task inst_tasks[MAX_TASKS_TAG_LENGTH];
    Task inst_tasks_vt[MAX_TASKS_TAG_LENGTH];
    Arc inst_arcs[MAX_ARCS_TAG_LENGTH];
    Vehicles state;

    memset(state.stop, 0, sizeof(state.stop));
    memset(state.remain_capacity, 0, sizeof(state.remain_capacity));
    memset(state.remain_seqs, 0, sizeof(state.remain_seqs));
    memset(state.not_served_task_seq, 0, sizeof(state.not_served_task_seq));
    read_instance_from_xml(inst_tasks, inst_arcs, &state, instance);

    update_cost(inst_tasks, inst_arcs);
    mod_dijkstra();
    construct_virtual_task(inst_tasks, inst_tasks_vt, state.stop, state.remain_capacity);
    int additional_cost = get_additional_cost(state);

    available_edges_map(inst_tasks_vt);

    // CARPInd solution, converted_solution;
    // int serve_mark[MAX_TASK_TAG_LENGTH];
    // memset(serve_mark, 0, sizeof(serve_mark));
    // for (int i = 1; i <= task_num; i++)
    // {
    //     serve_mark[i] = 1;
    // }
    // path_scanning(&solution, inst_tasks_vt, serve_mark);
    // indi_route_converter(&converted_solution, &solution, inst_tasks_vt);

    // auto_corrlation(solution, inst_tasks_vt, instance);

    // // srand(2002);
    // CARPInd solution1, solution2;
    // double avr_dis = 0;
    // int x = 0;
    // for (int i=1; i<=5000;i++)
    // {
    //     rand_scanning(&solution1, inst_tasks_vt);
    //     avr_dis += solution1.TotalCost;
    //     // rand_scanning(&solution2, inst_tasks_vt);
    //     // x = hamming_distance(solution1.Sequence, solution2.Sequence, inst_tasks_vt);
    //     // avr_dis += x;
    // }
    // int avr_cost = int(avr_dis/5000);
    // printf("%d\n", avr_cost);
    // char path[50];
    // memset(path, 0, sizeof(path));
    // sprintf(path, "result/best/%s/instance%d.txt", map, instance);
    // FILE *fp;
    // fp = fopen(path, "a");
    // fprintf(fp, "expected:%d\n", avr_cost);
    // fclose(fp);

    // local_optimum_related(inst_tasks_vt, instance);
    // lo_analysis(instance, inst_tasks_vt);
    lo_dis_analysis(instance, inst_tasks_vt);
    // return_probability(instance, inst_tasks_vt);
}




void available_edges_map(const Task* inst_tasks)
{
    int VRtemp[MAX_NODE_TAG_LENGTH];
    memset(VRtemp, 0, sizeof(VRtemp));

    int i, j;
    int RNum = inst_tasks[1].inverse-1;

    int idx = 1;
    std::string edge;
    for (i=0; i < 2*RNum; i++)
    {
        for (j = i+1; j <= 2*RNum; j++)
        {   
            if (i % 2 == 1 && j == i+1)
                continue;
            edge = std::to_string(i) + "," + std::to_string(j);
            edgesIDs.insert({edge, idx});
            edge = std::to_string(j) + "," + std::to_string(i);
            edgesIDs.insert({edge, idx});
            // std::cout << idx << ":" << edge << std::endl;
            idx ++;
        }
    }
}

