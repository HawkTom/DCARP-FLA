#include "src.h"
#include <cstring>
#include "../libs/tinyxml2.h"
using namespace tinyxml2;

void readMap(Task *inst_tasks, Arc *inst_arcs, const char *map1)
{
    FILE *fp;

    char dummy[101];

// //    fp = fopen("../instances/example.dat", "r");
//     char path[101];
//     // strcpy(path, "../map/");.dat
//     strcpy(path, "maps/");
//     strcat(path, map1);

    char path[50];
    memset(path, 0, sizeof(path));
    sprintf(path, "maps/%s.dat", map);

    fp = fopen(path, "r");
    if (fp == NULL)
    {
        printf("The file <%s> can't be open\n", path);
        exit(0);
    }

    while (fscanf(fp, "%s", dummy) != EOF)
    {

        if (strcmp(dummy, "VERTICES")==0)
        {
            fscanf(fp, "%s", dummy);
            fscanf(fp, "%d", &vertex_num);
        }
        else if (strcmp(dummy, "ARISTAS_REQ") == 0)
        {
            fscanf(fp, "%s", dummy);
            fscanf(fp, "%d", &req_edge_num);
        }
        else if (strcmp(dummy, "ARISTAS_NOREQ")==0)
        {
            fscanf(fp, "%s", dummy);
            fscanf(fp, "%d", &nonreq_edge_num);
        }
        else if (strcmp(dummy, "VEHICULOS")==0)
        {
            fscanf(fp, "%s", dummy);
            fscanf(fp, "%d", &vehicle_num);
        }
        else if (strcmp(dummy, "CAPACIDAD")==0)
        {
            fscanf(fp, "%s", dummy);
            fscanf(fp, "%d", &capacity);
        }
        else if (strcmp(dummy, "LISTA_ARISTAS_REQ")==0) {

            fscanf(fp, "%s", dummy);
            task_num = 2 * req_edge_num + req_arc_num;
            total_arc_num = task_num + 2 * nonreq_edge_num + nonreq_arc_num;
            for (int i = 1; i <= req_edge_num; i++) {
                fscanf(fp, "%s", dummy);
                fscanf(fp, "%d,", &inst_tasks[i].head_node);
                fscanf(fp, "%d)", &inst_tasks[i].tail_node);
                fscanf(fp, "%s", dummy);
                fscanf(fp, "%d", &inst_tasks[i].serv_cost);
                fscanf(fp, "%s", dummy);
                fscanf(fp, "%d", &inst_tasks[i].demand);

                inst_tasks[i].dead_cost = inst_tasks[i].serv_cost;
                inst_tasks[i].inverse = i + req_edge_num;
                inst_tasks[i].vt = 0;
                inst_tasks[i].head_node_r = 2*i-1;
                inst_tasks[i].tail_node_r = 2*i;

                inst_tasks[i + req_edge_num].head_node = inst_tasks[i].tail_node;
                inst_tasks[i + req_edge_num].tail_node = inst_tasks[i].head_node;
                inst_tasks[i + req_edge_num].dead_cost = inst_tasks[i].dead_cost;
                inst_tasks[i + req_edge_num].serv_cost = inst_tasks[i].serv_cost;
                inst_tasks[i + req_edge_num].demand = inst_tasks[i].demand;
                inst_tasks[i + req_edge_num].inverse = i;
                inst_tasks[i + req_edge_num].vt = 0;
                inst_tasks[i + req_edge_num].head_node_r = inst_tasks[i].tail_node_r;
                inst_tasks[i + req_edge_num].tail_node_r = inst_tasks[i].head_node_r;


                inst_arcs[i].head_node = inst_tasks[i].head_node;
                inst_arcs[i].tail_node = inst_tasks[i].tail_node;
                inst_arcs[i].trav_cost = inst_tasks[i].dead_cost;
                inst_arcs[i].change = 0;
                inst_arcs[i].link = 0;
                inst_arcs[i + req_edge_num].head_node = inst_arcs[i].tail_node;
                inst_arcs[i + req_edge_num].tail_node = inst_arcs[i].head_node;
                inst_arcs[i + req_edge_num].trav_cost = inst_arcs[i].trav_cost;
                inst_arcs[i + req_edge_num].change = 0;
                inst_arcs[i+ req_edge_num].link = 0;

                if (costlb > inst_tasks[i].dead_cost)
                    costlb = inst_tasks[i].dead_cost;

                if (costub < inst_tasks[i].dead_cost)
                    costub = inst_tasks[i].dead_cost;

                if (demandlb > inst_tasks[i].demand)
                    demandlb = inst_tasks[i].demand;

                if (demandub < inst_tasks[i].demand)
                    demandub = inst_tasks[i].demand;

            }
        }
        else if (strcmp(dummy, "LISTA_ARISTAS_NOREQ")==0)
        {
            fscanf(fp, "%s", dummy);
            for (int i=task_num+1; i<=task_num+nonreq_edge_num;i++)
            {
                fscanf(fp, "%s", dummy);
                fscanf(fp, "%d,", &inst_arcs[i].head_node);
                fscanf(fp, "%d)", &inst_arcs[i].tail_node);
                fscanf(fp, "%s", dummy);
                fscanf(fp, "%d", &inst_arcs[i].trav_cost);
                inst_arcs[i].change = 0;
                inst_arcs[i].link = 0;

                inst_arcs[i + nonreq_edge_num].head_node = inst_arcs[i].tail_node;
                inst_arcs[i + nonreq_edge_num].tail_node = inst_arcs[i].head_node;
                inst_arcs[i + nonreq_edge_num].trav_cost = inst_arcs[i].trav_cost;
                inst_arcs[i + nonreq_edge_num].change = 0;
                inst_arcs[i + nonreq_edge_num].link = 0;

                if (costlb > inst_arcs[i].trav_cost)
                    costlb = inst_arcs[i].trav_cost;

                if (costub < inst_arcs[i].trav_cost)
                    costub = inst_arcs[i].trav_cost;

            }
        }
        else if (strcmp(dummy, "DEPOSITO")==0)
        {
            fscanf(fp, "%s", dummy);
            fscanf(fp, "%d", &DEPOT);
        }

    }

    fclose(fp);

    inst_tasks[0].head_node = DEPOT;
    inst_tasks[0].tail_node = DEPOT;
    inst_tasks[0].head_node_r = 0;
    inst_tasks[0].tail_node_r = 0;
    inst_tasks[0].dead_cost = 0;
    inst_tasks[0].serv_cost = 0;
    inst_tasks[0].demand = 0;
    inst_tasks[0].inverse = 0;
    inst_arcs[0].head_node = DEPOT;
    inst_arcs[0].tail_node = DEPOT;
    inst_arcs[0].trav_cost = 0;

    for (int i=1; i<=total_arc_num; i++)
    {
        cost_backup[inst_arcs[i].head_node][inst_arcs[i].tail_node] = inst_arcs[i].trav_cost;
    }

    memset(edge_index, 0, sizeof(edge_index));
    edge_index[0] = req_edge_num;
    edge_index[1] = nonreq_edge_num;

}

void read_instance_from_xml(Task *inst_tasks, Arc *inst_arcs, Vehicles *state, int instance_idx)
{
    char path[50];
    memset(path, 0, sizeof(path));
    sprintf(path, "instance/%s/instance%d.xml", map, instance_idx);
    // sprintf(path, "xml/%s-instance%d.xml", map, instance_idx);

    // std::cout << path << std::endl;

    XMLDocument doc;
    doc.LoadFile(path);
    XMLElement* root = doc.FirstChildElement( "Instance" );

    XMLElement* vehicles;
    vehicles = root->FirstChildElement("State");

    char rem_seq_str[2000];

    state->stop[0] = 0;
    state->remain_capacity[0] = 0;
    state->remain_seqs[0] = 1;
    state->remain_seqs[1] = 0;
    state->not_served_task_seq[0] = 1;
    state->not_served_task_seq[1] = 0;

    XMLElement* curr_veh;
    curr_veh = vehicles->FirstChildElement("vehicle");
    while (true)
    {
        if (curr_veh == 0) break;

        int stop_point = atoi(curr_veh->Attribute("stop"));

        if (stop_point == 1)
        {
            strcpy(rem_seq_str, curr_veh->Attribute("seq"));
            char *token = strtok(rem_seq_str, " ");
            while (token != NULL) {
                state->not_served_task_seq[0]++;
                state->not_served_task_seq[state->not_served_task_seq[0]] = atoi(token);
                token = strtok(NULL, " ");
            }
            state->not_served_task_seq[0]++;
            state->not_served_task_seq[state->not_served_task_seq[0]] = 0;
        } else{
            state->stop[0] ++;
            state->stop[state->stop[0]] = stop_point;

            state->remain_capacity[0] ++;
            state->remain_capacity[state->remain_capacity[0]] = atoi(curr_veh->Attribute("rcapacity"));

            strcpy(rem_seq_str, curr_veh->Attribute("seq"));
            char *token = strtok(rem_seq_str, " ");
            while (token != NULL) {
                state->remain_seqs[0]++;
                state->remain_seqs[state->remain_seqs[0]] = atoi(token);
                token = strtok(NULL, " ");
            }
            state->remain_seqs[0]++;
            state->remain_seqs[state->remain_seqs[0]] = 0;
        }
        
        XMLElement* next_veh = curr_veh->NextSiblingElement();
        curr_veh = next_veh;
    }

    XMLElement* tasks;
    tasks = root->FirstChildElement("Tasks");

    task_num = atoi(tasks->Attribute("num"));

    XMLElement* curr_task;
    curr_task = tasks->FirstChildElement("task");
    int k = -1;
    while(true)
    {
        if (curr_task == 0) break;
        k ++;
        inst_tasks[k].head_node = atoi(curr_task->Attribute("head_node"));
        inst_tasks[k].tail_node = atoi(curr_task->Attribute("tail_node"));
        inst_tasks[k].demand = atoi(curr_task->Attribute("demand"));
        inst_tasks[k].dead_cost = atoi(curr_task->Attribute("dead_cost"));
        inst_tasks[k].serv_cost = atoi(curr_task->Attribute("serv_cost"));
        inst_tasks[k].inverse = atoi(curr_task->Attribute("inverse"));
        // inst_tasks[k].head_node_r = 2*k-1;
        // inst_tasks[k].tail_node_r = 2*k;

        // inst_tasks[k].vt = atoi(curr_task->Attribute("vt"));

        XMLElement* next_task = curr_task->NextSiblingElement();
        curr_task = next_task;
    }

    XMLElement* arcs;
    arcs = root->FirstChildElement("Arcs");
    XMLElement* curr_arc;
    curr_arc = arcs->FirstChildElement("arc");
    k = -1;
    while(true)
    {
        if (curr_arc == 0) break;
        k ++;
        inst_arcs[k].head_node = atoi(curr_arc->Attribute("head_node"));
        inst_arcs[k].tail_node = atoi(curr_arc->Attribute("tail_node"));
        inst_arcs[k].trav_cost = atoi(curr_arc->Attribute("trav_cost"));
        // inst_arcs[k].change = atoi(curr_arc->Attribute("change"));
        // inst_arcs[k].link = atoi(curr_arc->Attribute("link"));

        XMLElement* next_arc = curr_arc->NextSiblingElement();
        curr_arc = next_arc;
    }

    XMLElement* info;
    info = root->FirstChildElement("Info");
    XMLElement* tmp;
    tmp = info->FirstChildElement("req_edge");
    req_edge_num = atoi(tmp->Attribute("num"));
    // std::cout << req_edge_num << std::endl;

    tmp = info->FirstChildElement("non_req_edge");
    nonreq_edge_num = atoi(tmp->Attribute("num"));

    tmp = info->FirstChildElement("vertex");
    vertex_num = atoi(tmp->Attribute("num"));

    tmp = info->FirstChildElement("total_arc");
    total_arc_num = atoi(tmp->Attribute("num"));

    tmp = info->FirstChildElement("capacity");
    capacity = atoi(tmp->Attribute("num"));

    tmp = info->FirstChildElement("vehicle_num");
    vehicle_num = atoi(tmp->Attribute("num"));

    edge_index[0] = req_edge_num;
    edge_index[1] = nonreq_edge_num;
}


void save_instance_to_xml(const Task *inst_tasks, const Arc *inst_arcs, const Vehicles state, int instance_idx)
{
    int i;

    XMLDocument doc;
    XMLElement* root = doc.NewElement("Instance");
    root->SetAttribute("map", map);
    root->SetAttribute("instance", instance_idx);
    doc.InsertFirstChild(root);
    
    XMLElement* vehicles;
    vehicles = root->InsertNewChildElement("State");

    for (i = 1; i <= state.stop[0]; i++)
    {
        XMLElement* veh = vehicles->InsertNewChildElement("vehicle");
        veh->SetAttribute("stop", state.stop[i]);
        veh->SetAttribute("rcapacity", state.remain_capacity[i]);
    }

    char str[2000];
    char stm[10];
    memset(str, 0, sizeof(str));
    for (i=1; i <= state.remain_seqs[0]; i++)
    {
        sprintf(stm, "%d", state.remain_seqs[i]);
        strcat(str, stm);
        strcat(str, " ");
    }
    vehicles->SetAttribute("seq", str);


    XMLElement* tasks;
    tasks = root->InsertNewChildElement("Tasks");
    tasks->SetAttribute("num", task_num);
    for (i = 0; i <= task_num; i++)
    {   
        XMLElement* task;
        task = tasks->InsertNewChildElement("task");
        task->SetAttribute("head_node", inst_tasks[i].head_node);
        task->SetAttribute("tail_node", inst_tasks[i].tail_node);
        task->SetAttribute("demand", inst_tasks[i].demand);
        task->SetAttribute("dead_cost", inst_tasks[i].dead_cost);
        task->SetAttribute("serv_cost", inst_tasks[i].serv_cost);
        task->SetAttribute("inverse", inst_tasks[i].inverse);
        task->SetAttribute("vt", inst_tasks[i].vt);
    }

    XMLElement* arcs;
    arcs = root->InsertNewChildElement("Arcs");
    for (i = 0; i <= total_arc_num; i++)
    {   
        XMLElement* arc;
        arc = arcs->InsertNewChildElement("arc");
        arc->SetAttribute("head_node", inst_arcs[i].head_node);
        arc->SetAttribute("tail_node", inst_arcs[i].tail_node);
        arc->SetAttribute("trav_cost", inst_arcs[i].trav_cost);
        arc->SetAttribute("change", inst_arcs[i].change);
        arc->SetAttribute("link", inst_arcs[i].link);
    }

    XMLElement* info;
    info = root->InsertNewChildElement("Info");

    XMLElement* tmp;
    tmp = info->InsertNewChildElement("req_edge");
    tmp->SetAttribute("num", req_edge_num);

    tmp = info->InsertNewChildElement("non_req_edge");
    tmp->SetAttribute("num", nonreq_edge_num);

    tmp = info->InsertNewChildElement("vertex");
    tmp->SetAttribute("num", vertex_num);

    tmp = info->InsertNewChildElement("total_arc");
    tmp->SetAttribute("num", total_arc_num);


    char path[50];
    memset(path, 0, sizeof(path));

    sprintf(path, "xml/%s-instance%d.xml", map, instance_idx);
    doc.SaveFile(path);
}