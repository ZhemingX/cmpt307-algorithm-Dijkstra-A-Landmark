//
// Created by Zheming Xu on 2019/10/30.
// ID:301414922
//This .cpp will perform the Dijkstra's Algorithm, A*-search Algorithm and ,landmark algorithm each for the given nodes
// and the adjacent lists, and return the number of nodes visited for each algorithm for comparison.
//The input form is:
//Enter the information about nodes:
//1 : 49.2622862226, -122.920346717
//2 : 49.2628916127, -122.914482463
//...
//Enter the infomation about adjacent lists:
//1 : 863, 864, 
//2 : 90, 
//3 : 949, 
//...
//Here are the results of 20 pairs of nodes generated randomly:
//(final results)

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>

#define MAX 2000
#define INF 10000000

using namespace std;

// construct graph ADT
typedef struct node linknode; //define node
typedef struct node_head header; //define header of each adjacent list
typedef struct graph Graph; // define the graph

struct node{
    int num;
    double info[2]; //record info[0]=latitude and info[1]=longitude
    double dist; //record dist[u];
    double Adist;//record dist[u]+to_t[u];
    linknode* next;
};

struct node_head{
    linknode* first;
};

struct graph{
    int n;
    header adjaTable[MAX];
};

linknode nv[MAX];
Graph G;

void createGraphNode(int num_node){ //determine the info of each node
    char c1,c2;
    for(int i = 1;i <= num_node;i++){
        cin>>nv[i].num>>c1>>nv[i].info[0]>>c2>>nv[i].info[1];
        nv[i].next = NULL;
        nv[i].dist = 0;
        nv[i].Adist = 0;
    }
};

void Graph_init(int num_node){ //initial the graph
    G.n = num_node;
    for(int i = 1;i <= G.n;i++){
        G.adjaTable[i].first = &nv[i];
    }
}

void createGraphEdge(int num_node){//determine the edge of the graph

    char buffer[100];
    char s[10];
    linknode* temp;
    for(int i = 1;i <= num_node+1;i++){
        cin.getline(buffer,100);
        int tableheader;
        int start = 0;
        int end = start;
        int len = strlen(buffer);
        while(start<len && end<len){
            end = start;
            if(isdigit(buffer[start])){
                while(isdigit(buffer[end+1])){
                    end++;
                } //get a data type string

                strncpy(s,buffer+start,end-start+1);
                s[end-start+1] = '\0';
                //string to int
                int data = atoi(s);
                if(start==0){
                    tableheader = data;
                    start = end + 1;
                    temp = G.adjaTable[tableheader].first;
                }
                else{
                    linknode *p = (linknode*)malloc(sizeof(linknode));
                    p->num = data;
                    p->info[0] = nv[data].info[0];
                    p->info[1] = nv[data].info[1];
                    p->next = NULL;
                    p->dist = 0;
                    p->Adist = 0;
                    temp->next = p;
                    temp = p;
                    start = end + 1;
                }
            }
            else
                start++;
        }
    }
}
//

//construct distance function
double distance(linknode n1,linknode n2){
    double p1_0 = n1.info[0], p1_1 = n1.info[1], p2_0 = n2.info[0], p2_1 = n2.info[1];
    double dlat = 2*M_PI*(p2_0-p1_0)/360;
    double mlat = 2*M_PI*(p1_0+p2_0)/2/360;
    double dlon = 2*M_PI*(p2_1-p1_1)/360;
    return 6371009*pow(pow(dlat,2)+pow(cos(mlat)*dlon,2),0.5);
}
//

//construct heap ADT
#define Mindata -1
struct HeapStruct;
typedef struct HeapStruct *PriorityQueue;

PriorityQueue Initialize(int Maxelements);
int isEmpty(PriorityQueue H);
void PrecolateDown1(PriorityQueue H, int position);
void PrecolateDown2(PriorityQueue H, int position);
void CreateHeap1(PriorityQueue H, linknode* nv, int num);
void CreateHeap2(PriorityQueue H, linknode* nv, int num);
linknode DeleteMin1(PriorityQueue H);
linknode DeleteMin2(PriorityQueue H);
void DecreaseKey1(PriorityQueue H,int id, double new_dist);
void DecreaseKey2(PriorityQueue H,int id, double new_dist, double new_Adist);

struct HeapStruct{
    int Capacity;
    int Size;
    linknode Elements[MAX];
};

PriorityQueue Initialize(){
    PriorityQueue H;
    H = (PriorityQueue)malloc(sizeof(struct HeapStruct));
    H->Capacity = MAX;
    H->Size = 0;
    (H->Elements[0]).dist = Mindata;
    (H->Elements[0]).Adist = Mindata;

    return H;
}

int isEmpty(PriorityQueue H){
    return H->Size == 0;
}

void PrecolateDown1(PriorityQueue H, int position){ //for dist as key value
    int i,Child;
    linknode temp = H->Elements[position];
    for(i=position;i*2<=H->Size;i=Child){
        Child=i*2;
        if(Child!=H->Size&&(H->Elements[Child+1]).dist<(H->Elements[Child]).dist)
            Child++;
        if(temp.dist>(H->Elements[Child]).dist)
            H->Elements[i] = H->Elements[Child];
        else
            break;
    }
    H->Elements[i] = temp;
};

void PrecolateDown2(PriorityQueue H, int position){ //for dist[u]+to_t[u] as key value
    int i,Child;
    linknode temp = H->Elements[position];
    for(i=position;i*2<=H->Size;i=Child){
        Child=i*2;
        if(Child!=H->Size&&(H->Elements[Child+1]).Adist<(H->Elements[Child]).Adist)
            Child++;
        if(temp.Adist>(H->Elements[Child]).Adist)
            H->Elements[i] = H->Elements[Child];
        else
            break;
    }
    H->Elements[i] = temp;
};

void CreateHeap1(PriorityQueue H, linknode* nv, int num){
    for(int i=1;i<=num;i++){
        H->Elements[i] = nv[i];
        H->Size++;
    }
    for(int i=num/2;i>0;i--)
        PrecolateDown1(H,i);
}

void CreateHeap2(PriorityQueue H, linknode* nv, int num){
    for(int i=1;i<=num;i++){
        H->Elements[i] = nv[i];
        H->Size++;
    }
    for(int i=num/2;i>0;i--)
        PrecolateDown2(H,i);
}

linknode DeleteMin1(PriorityQueue H){
    linknode MinElement;

    if(isEmpty(H))
        return H->Elements[0];

    MinElement = H->Elements[1];
    H->Elements[1] = H->Elements[H->Size--];
    PrecolateDown1(H,1);

    return MinElement;
};

linknode DeleteMin2(PriorityQueue H){
    linknode MinElement;

    if(isEmpty(H))
        return H->Elements[0];

    MinElement = H->Elements[1];
    H->Elements[1] = H->Elements[H->Size--];
    PrecolateDown2(H,1);

    return MinElement;
};

void DecreaseKey1(PriorityQueue H,int id, double new_dist){ // for Dijkstra's algorithm
    for(int i=1;i<=H->Size;i++){
        if((H->Elements[i]).num==id){
            (&(H->Elements[i]))->dist = new_dist;
        }
    }
    for(int i=(H->Size)/2;i>0;i--)
        PrecolateDown1(H,i);
};

void DecreaseKey2(PriorityQueue H,int id, double new_dist, double new_Adist){ // A*-search
    for(int i=1;i<=H->Size;i++){
        if((H->Elements[i]).num==id){
            (&(H->Elements[i]))->dist = new_dist;
            (&(H->Elements[i]))->Adist = new_Adist;
        }
    }
    for(int i=(H->Size)/2;i>0;i--)
        PrecolateDown2(H,i);
};
//
int parent[MAX];//record parent of visited node

//Dijkstra's Algorithm: output count, and dis[t] if you want
int Dijkstra(Graph G, int s, int t) {
    //record if node has been visited
    int visited[MAX] = {0};
    //count record the number of nodes visited
    int count = 0;
    //initialization
    for (int i = 1; i <= G.n; i++) {
        if (i == s) {
            nv[i].dist = 0;
        }
        else {
            nv[i].dist = INF;
        }
    }
    //initialize priority queue Q to contain all nodes
    // using dist values as keys
    PriorityQueue H = Initialize();
    CreateHeap1(H, nv, G.n);
    while (!isEmpty(H)) {
        linknode temp = DeleteMin1(H);
        if (visited[temp.num] == 0){
            count++;
            visited[temp.num] = 1;
        }

        if (temp.num == t) //end flag
            return count;
        int temp_index = temp.num;
        linknode *v = G.adjaTable[temp_index].first->next;
        while (v != NULL) {
            int v_nv = v->num;
            if (nv[v_nv].dist > temp.dist + distance(temp, nv[v_nv])) {
                //relax
                nv[v_nv].dist = temp.dist + distance(temp, nv[v_nv]);
                DecreaseKey1(H, v_nv, nv[v_nv].dist);
                parent[v_nv] = temp.num;
            }
            v = v->next;
        }
    }

    return count;
};

//A* search Algorithm: output count, and dis[t] if you want
int A_search(Graph G, int s, int t){
    //record if node has been visited
    int visited[MAX] = {0};
    //count record the number of nodes visited
    int count = 0;
    //initialization
    for (int i = 1; i <= G.n; i++) {
        if (i == s) {
            nv[i].dist = 0;
            nv[i].Adist = nv[i].dist + distance(nv[i],nv[t]); //calculate Adist[u] = dist[u] + to_t[u];
        }
        else {
            nv[i].dist = INF;
            nv[i].Adist = nv[i].dist + distance(nv[i],nv[t]); //calculate Adist[u] = dist[u] + to_t[u];
        }
    }
    //initialize priority queue Q to contain all nodes
    // using dist[u]+to_t[u] values as keys
    PriorityQueue H = Initialize();
    CreateHeap2(H, nv, G.n);
    while (!isEmpty(H)) {
        linknode temp = DeleteMin2(H);
        if (visited[temp.num] == 0){
            count++;
            visited[temp.num] = 1;
        }

        if (temp.num == t) //end flag
            return count;
        int temp_index = temp.num;
        linknode *v = G.adjaTable[temp_index].first->next;
        while (v != NULL) {
            int v_nv = v->num;
            if (nv[v_nv].dist > temp.dist + distance(temp, nv[v_nv])) {
                //relax
                nv[v_nv].dist = temp.dist + distance(temp, nv[v_nv]);
                nv[v_nv].Adist = nv[v_nv].dist + distance(nv[v_nv], nv[t]);
                DecreaseKey2(H, v_nv, nv[v_nv].dist, nv[v_nv].Adist);
                parent[v_nv] = temp.num;
            }
            v = v->next;
        }
    }

    return count;
};

//landmark, this time randomly choose three landmarks nv[250], nv[500], nv[750];
double dist_250[MAX], dist_500[MAX], dist_750[MAX]; //distz[u]
//determine to_t[u]
double to_t(int u, int t){
    double to_t1 = abs(dist_250[u]-dist_250[t]);
    double to_t2 = abs(dist_500[u]-dist_500[t]);
    double to_t3 = abs(dist_750[u]-dist_750[t]);

    double submax = (to_t1>to_t2)?to_t1:to_t2;
    double max = (submax>to_t3)?submax:to_t3;

    return max;
}
//landmark Algorithm: output count, and dis[t] if you want
int landmark(Graph G, int s, int t){
    //record if node has been visited
    int visited[MAX] = {0};
    //count record the number of nodes visited
    int count = 0;
    //initialization
    for (int i = 1; i <= G.n; i++) {
        if (i == s) {
            nv[i].dist = 0;
            nv[i].Adist = nv[i].dist + to_t(i,t); //calculate Adist[u] = dist[u] + to_t[u];
        }
        else {
            nv[i].dist = INF;
            nv[i].Adist = nv[i].dist + to_t(i,t); //calculate Adist[u] = dist[u] + to_t[u];
        }
    }
    //initialize priority queue Q to contain all nodes
    // using dist[u]+to_t[u] values as keys
    PriorityQueue H = Initialize();
    CreateHeap2(H, nv, G.n);
    while (!isEmpty(H)) {
        linknode temp = DeleteMin2(H);
        if (visited[temp.num] == 0){
            count++;
            visited[temp.num] = 1;
        }

        if (temp.num == t) //end flag
            return count;
        int temp_index = temp.num;
        linknode *v = G.adjaTable[temp_index].first->next;
        while (v != NULL) {
            int v_nv = v->num;
            if (nv[v_nv].dist > temp.dist + distance(temp, nv[v_nv])) {
                //relax
                nv[v_nv].dist = temp.dist + distance(temp, nv[v_nv]);
                nv[v_nv].Adist = nv[v_nv].dist + to_t(v_nv,t);
                DecreaseKey2(H, v_nv, nv[v_nv].dist, nv[v_nv].Adist);
                parent[v_nv] = temp.num;
            }
            v = v->next;
        }
    }

    return count;
}

int main(){
    //generate the graph
    int nver;
    cin>>nver;
    cout<<"Enter the information about nodes:"<<endl;
    createGraphNode(nver);
    Graph_init(nver);
    cout<<"Enter the infomation about adjacent lists:"<<endl;
    createGraphEdge(nver);

    //generate dist_250[], dist_500[], dist_750[] for landmark algorithm
    for(int i=1;i<=1000;i++){
        A_search(G,250,i);
        dist_250[i] = G.adjaTable[i].first->dist;
    }
    for(int i=1;i<=1000;i++){
        A_search(G,500,i);
        dist_500[i] = G.adjaTable[i].first->dist;
    }
    for(int i=1;i<=1000;i++){
        A_search(G,750,i);
        dist_750[i] = G.adjaTable[i].first->dist;
    }

    cout<<"Here are the results of 20 pairs of nodes generated randomly:"<<endl;

    //generate source node and destination node
    int s,t;
    //count_1, count_2 count_3 are used to record number of nodes visited for Dijkstra, A* search and landmark Algorithm
    int count_1;
    int count_2;
    int count_3;

    int total = 1;
    while(total<=20){ //randomly generate 20 pairs
        s = rand()%1000+1;
        t = rand()%1000+1;
        count_1 = Dijkstra(G,s,t);
        count_2 = A_search(G,s,t);
        count_3 = landmark(G,s,t);
        cout<<"random choose pair("<<s<<","<<t<<")"<<endl;
        cout<<"the number of vertices Dijkstra algorithm has visited is "<<count_1<<endl;
        cout<<"the number of vertices A* algorithm has visited is "<<count_2<<endl;
        cout<<"the number of vertices Landmark algorithm has visited is "<<count_3<<endl;

        total++;
    }
    return 0;
}