/*******************
initialize:04/07/2024 
update: 04/02/2026
********************/


#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <vector>
#include <utility>

using namespace std;


/************************************************************************************/
/*	CLASSE NODE -- Represents an individual in the simulation, storing its identity, *
 *                 spatial coordinates, and local connectivity	                     *
/***********************************************************************************/

class node{
	private:
        int NodeId; //Node's index
        pair<int,int> NodePosition; //Node's position;
		vector<pair<int,int>> NeighPositions; //Stores the positions that are neighbors of the node;
        vector<int> NeighId; //Stores the index of neighbors;
    public:
        int get_nodeid() const; //Return the id of the node;
		pair<int,int> get_position() const ; //Return the position x and y of the node;
        int degree() const; // Returns the node degree (number of neighbors)
		void add_neighbor_position(pair<int,int> v); // Adds a neighbor index position 'v' to the list of neighbors positions;
        void add_neighbor_id(int v); // Adds a neighbor index 'v'to the list of neighbors id;
        const vector<pair<int,int>> get_neighbors() const; //Returns the neighboring position of the node;
        const vector<int> get_neighbors_id() const; //Return the index of neighbors;
        void print_neighbors() const ; //Prints the neighboring nodes (positions);
        void clear_neighbor_id(); //Clean the vector;
        void print_neighbor_ids() const; //Prints the neighboring nodes (Id);
        vector<int> get_neighbor_ids(int nodeId) const;
        node(int _NodeId, pair<int,int> _NodePosition): NodeId(_NodeId), NodePosition(_NodePosition){} //Definition of the node which depends of the id, position x and y. (Constructor)
        void clear(); //clean all properties of the node;
};


void node::clear() {
    NodeId = 0;
    NodePosition = {0, 0};
    NeighPositions.clear();
    NeighId.clear();
}

inline int node::get_nodeid() const {return NodeId;}
inline pair<int,int> node::get_position() const {return NodePosition;}
inline int node::degree() const { return NeighId.size(); }
inline void node::add_neighbor_position(pair<int,int>v){NeighPositions.push_back(v); }
inline void node::add_neighbor_id(int v){NeighId.push_back(v);}
inline const vector<pair<int,int>> node::get_neighbors() const {return NeighPositions;}
inline const vector<int>node::get_neighbors_id() const {return NeighId;}
inline void node::print_neighbors() const{
    cout << "Node " << NodeId << " neighbors:" << endl;
    for (const auto& neighbor : NeighPositions) {
        cout << "    (" << neighbor.first << ", " << neighbor.second << ")" << endl;
    }
};

inline void node::clear_neighbor_id(){
    NeighId.clear();
}

inline void node::print_neighbor_ids() const{
    cout << "Neighbor IDs of Node " << NodeId << ":" << endl;
    for (const auto& id : NeighId) {
        cout << "    " << id << endl;
    }
}


#endif
