/*******************
initialize:04/07/2024 
update: 04/02/2026
Lara D. Hissa
********************/


#ifndef GRAPH_H
#define GRAPH_H

#include <unordered_map>
#include <set>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <random>
#include "node.h"


/*********************************************************************************/
/*          CLASS GRAPH -- Manages the spatial environment, node positioning,                
                           and connectivity based on radius S
/*********************************************************************************/   

class graph{
    private:
        map<pair<int,int>,int> Density; //Stores the number of individuals per site;
        vector<node> NodeList; //Stores the characteristics about node;
        vector<pair<int,int>> RegionS; //Neighboring region;
        
       
    public:
        int get_numnodes(); //Return the number of the nodes in the graph;
        int get_degree(int v); //Returns the degree of the node v;
        pair<int,int> get_position(int v); //Returns the position of the node;
        void generate_random_positions(const int &NumNodes, const int &SizeX, const int &SizeY, std::mt19937& gen); //Function that places individuals in random positions on an X by Y grid;
        int get_nodeid(int v); //Returns the id of the node ;
        void print_density_nodes_position(); //Prints the density of nodes at each location in the graph;
        void print_node_position();  //Prints the node's position;
        void region_S(const double &S, const int &SizeX, const int &SizeY); //Function to create the region S of radius;
        void add_neighbors(const int &SizeX, const int &SizeY); //Graph creation, adds the neighbors of each node;
        void print_list_neighbors();
        void clear_node_list(); //Clean the NodeList;
        void add_node(const node& newNode); //Add a node in the graph;
        const vector<pair<int, int>>&get_regionS() const; //Returns the RegionS;
        const vector<node>& get_node_list() const; //Returns NodeList;
        vector<int> get_neighbor_ids(int nodeId); //Returns the neighboring ids according to the id of a node; 
        void save_node_neighbors(const string &filename);
        void print_region_S() const;
        int available_positions(int v, const int &SizeX, const int &SizeY);
        float mean_degree();
        double local_density(int v, const int &SizeX, const int &SizeY);
             
};


inline int graph::get_numnodes(){return NodeList.size();}
inline int graph::get_degree(int v){return NodeList[v].degree();} 
inline pair<int,int> graph::get_position(int v){return NodeList[v].get_position();}
inline int graph::get_nodeid(int v){return NodeList[v].get_nodeid();}
inline void graph::clear_node_list(){NodeList.clear();}
inline void graph::add_node(const node& newNode){NodeList.push_back(newNode);}
inline const vector<node>&graph::get_node_list() const{return NodeList;}
inline const vector<pair<int,int>>&graph::get_regionS() const{return RegionS;}


float graph::mean_degree(){
    int sum = 0;
    for(int i = 0; i < get_numnodes(); i++){
        int NodeId = NodeList[i].degree();
        sum +=NodeId;
    }
    float mean = sum/NodeList.size();
return mean;
}

vector<int> graph::get_neighbor_ids(int nodeId) {
    vector<int> neighborIds;
    //Find the nodeId
    auto it = find_if(NodeList.begin(), NodeList.end(), [&](const node& n) {
        return n.get_nodeid() == nodeId;
    });

    //Check if the node was found, if found, get the IDs of its neighbors
    if (it != NodeList.end()) {
        const vector<int>& neighbors = it->get_neighbors_id();
        neighborIds.insert(neighborIds.end(), neighbors.begin(), neighbors.end());
    }

    return neighborIds;
}        





void graph:: generate_random_positions(const int &NumNodes, const int &SizeX, const int &SizeY, std::mt19937& gen){
  
  std::uniform_int_distribution<int>distX(0,SizeX-1);
  std::uniform_int_distribution<int>distY(0,SizeY-1);
    
    for(int i = 0; i<NumNodes;i++){
        int CoordX;
        int CoordY;
        
        CoordX = distX(gen);
        CoordY = distY(gen);        

        Density[{CoordX, CoordY}]++;
        NodeList.push_back(node(i,{CoordX, CoordY}));            
    }

}
    

void graph::region_S(const double &S, const int &SizeX, const int &SizeY) {
    int p = static_cast<int>(S + 1); 

    //Iterate over all positions in the current node's neighborhood square
    for (int i = -p; i <= p; i++) {
        for (int j = -p; j <= p; j++){
            if (i * i + j * j <= S * S) {
                RegionS.push_back({i, j});
            }
        }
    }
}


void graph::print_region_S() const {
    cout << "RegionS:" << endl;
    for (const auto& pos : RegionS) {
        cout << "(" << pos.first << ", " << pos.second << ")" << endl;
    }

}

int graph::available_positions(int v, const int &SizeX, const int &SizeY){
    int X = get_position(v).first;
    int Y = get_position(v).second;
    vector<pair<int,int>> positions;

    for(const auto &regionPos: RegionS){
        int PossibleX = X + regionPos.first;
        int PossibleY = Y + regionPos.second;
        if(PossibleX >= 0 && PossibleX < SizeX &&  PossibleY >= 0 && PossibleY < SizeY){
            positions.push_back({PossibleX,PossibleY});
        }
    }

   

    return positions.size();
}


inline void graph::add_neighbors(const int &SizeX, const int &SizeY){
    for(auto &node : NodeList){
        node.clear_neighbor_id();

        int positionX = node.get_position().first;
        int positionY = node.get_position().second;
       
        for(const auto &regionPos: RegionS){
            //Without periodic conditions
            int neighborX = positionX + regionPos.first;
            int neighborY = positionY + regionPos.second;

            //With periodic conditions
            //int neighborX = (positionX + regionPos.first + SizeX) % SizeX;
            //int neighborY = (positionY + regionPos.second + SizeY) % SizeY;


            if(neighborX >= 0 && neighborX < SizeX && neighborY >= 0 && neighborY < SizeY){ 
                for( auto& otherNode : NodeList){
                    if(otherNode.get_position() == make_pair(neighborX, neighborY) && otherNode.get_nodeid() != node.get_nodeid()){
                        node.add_neighbor_position({neighborX, neighborY});
                        node.add_neighbor_id(otherNode.get_nodeid());
                    }
                }
            }
        }
    }
}







inline void graph::print_density_nodes_position(){
    for(auto it = Density.begin(); it!= Density.end();it++){
        cout << "Position( " << it-> first.first << "," << it -> first.second <<"): " << it ->second << "individals" << endl; 
    }
}

inline void graph::print_node_position(){
    for(int i = 0; i < NodeList.size();i++){
        cout << "Node " <<NodeList[i].get_nodeid() << " at position(" << NodeList[i].get_position().first << "," << NodeList[i].get_position().second<<")" << endl;
    }
}

void graph::save_node_neighbors(const string &filename){
    ofstream outFile(filename);

    if (!outFile.is_open()) {
        cerr << "Error opening file for writing." << endl;
        return;
    }

    for (const auto& node : NodeList) {
        outFile << node.get_nodeid() << " " << node.degree() << " ";
        auto neighbors = node.get_neighbors_id();
        for (const auto& neighbor_id : neighbors) {
            outFile << neighbor_id << " ";
        }
        outFile << endl;
        }

    cout << "SAVE NODE NEIGHBORS " << filename << endl;
}

inline void graph::print_list_neighbors() {
    for (const auto& node : NodeList) {
        //node.print_neighbors();  //Uncomment this loop to print the connections for each node
        node.print_neighbor_ids();
       
    }
     cout << NodeList.size() << endl;

}


double graph::local_density(int v, const int &SizeX, const int &SizeY){
    int X = get_position(v).first;
    int Y = get_position(v).second;
    double total_density = 0;
    int count = 0;

    //Get the positions around the 'v' node
    for (const auto &regionPos : RegionS) {
        int PossibleX = X + regionPos.first;
        int PossibleY = Y + regionPos.second;

        if (PossibleX >= 0 && PossibleX < SizeX && PossibleY >= 0 && PossibleY < SizeY) {
            pair<int, int> pos = {PossibleX, PossibleY};

            if (Density.find(pos) != Density.end()) {
                total_density += Density[pos];  //Add position density
            } else {
                total_density += 0;  
            }

            count++;
        }
    }

    return (count > 0) ? (total_density / count) : 0;
}


#endif
