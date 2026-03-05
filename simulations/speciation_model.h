/*******************
initialize:04/07/2024 
update: 04/02/2026
Lara D. Hissa
********************/




#ifndef SPECIATION_MODEL
#define SPECIATION_MODEL

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <vector>
#include <map>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <fstream>
#include <ctime>
#include <cmath>
#include <climits>
#include<numeric>
#include "graph.h"


/*******************************************************************************************************************\
                        CLASS SPECIATION_MODEL -- Determines the speciation dynamics 
********************************************************************************************************************/

class speciation_model{
    private:   
        vector<vector<int>> RepGenome; // Population reproductive genomes
        vector<vector<int>> RepGenAux; // Auxiliary storage for offspring reproductive genomes
        vector<vector<int>> PheGenome; // Population phenotypic genomes
        vector<vector<int>> PheGenAux; // Auxiliary storage for offspring phenotypic genomes
        vector<double> Phenotype; //Expressed traits calculated from the sum of the loci by PheGenome
        vector<double> Fitness; // Environmental compatibility (selection pressure)
        
        
        map<pair<int,int>,int> DensityOffspring; //Stores the number of Offspring per site;
        vector<pair<int, int>>  OffspringPosition; 
        vector<node> OffspringList; //List of Offspring like NodeList 
        vector<vector<int>> GeneticNeigh; //Genetic Neighbors
        map<int, int> SpeciesId;
        vector<tuple<int, int, int>> OffspringParents; 
        map<int, int> SpeciesCount;
        map<int, vector<int>> SpeciesMembers;
        std::mt19937& gen_;
      




        public:
        
        vector<double> DIFF;      // Phenotypic differences at the environment boundary
        vector<double> DIFFTOTAL; // Global phenotypic variance
    
        speciation_model(graph &net, const int &SizeRepGenome, const int &SizePheGenome, std::mt19937& gen_ref);
    
        void dynamics(graph &net, const float &qmating, const float &MutationRate, const double &S, const int &SizeX, const int &SizeY,
                      const double &G, const int &NumNodes, const int &SizeRepGenome, const int &SizePheGenome, const float &Phenotype0, 
                      float &sigma, int currentGeneration, const int &maxGenerations, float &kmed);

        bool try_reproduce(graph &net, int i, const float &qmating, const float &MutationRate, const double &S, const int &SizeX, const int &SizeY, 
                           double rho0, const int &NumNodes, float &kmed);
       
        void reproduce(int father, int mother, const float &MutationRate, graph &net, const int &SizeX, const int &S);
        void offspring_region(graph &net, int i, const double &S, const int &SizeX, const int &SizeY); // Offspring spatial dispersal
        void new_graph(graph &net, const int &SizeX, const int &SizeY); // Transitions population to the next generation
        void update_genomes_from_offspring(graph &net, const int &SizeRepGenome, const int &SizePheGenome);
        bool find_new_partner(graph &net, int &i, int &k, int &mother, const float &qmating, float &sigma, int currentGeneration, const int &maxGenerations);

       
        // --- Biological and environmental calculations ---
        void phenotype(); // Translates genotype to phenotype (additive model)
        void compatibility_environmental(graph &net, const double &S, const int &SizeY, const int &SizeX, float &sigma, const float &Phenotype0);         
        double calculate_dissimilarity(int v, int u); // Genetic distance calculation
        void calculate_genetic_neighbors(graph &net, double G);
        pair<map<int,int>,int> check_species(graph &net, const double &G); // Clustering
    
        // --- auxiliaries---
        long double fitness_mean(int j, const int &SizeX, graph &net);
        int diff_phenotype_boundary(graph&net, int i, int j);
    
        // --- Getters ---
        vector<vector<int>> get_repgenome(); 
        vector<vector<int>> get_phegenome();
        vector<vector<int>> get_geneticneigh();
        map<pair<int,int>,int> get_density_offspring();
        double get_phenotype(int i);
        double get_fitness(int i);
        int get_speciesid(int i);
        tuple<int, int, int> get_parents(int i);
        vector<int> get_speciesmembers(int i);
        int get_numspecies_parents();
    };




speciation_model::speciation_model(graph &net, const int &SizeRepGenome, const int &SizePheGenome, std::mt19937& gen_ref): gen_(gen_ref){
   std::uniform_real_distribution<> dis(0.3, 0.7);

    int NumNodes = net.get_numnodes();
    RepGenome.assign(NumNodes, vector<int>(SizeRepGenome, 0));  // Initialize RepGenome
    PheGenome.assign(NumNodes, vector<int>(SizePheGenome, 0));  // Initialize PheGenome

    for(auto &vec : PheGenome){
        float randFactor = dis(gen_ref);

        int onesCount = randFactor*SizePheGenome;

        fill(vec.begin(), vec.begin() + onesCount, 1);
        shuffle(vec.begin(), vec.end(),gen_ref);
    
    }
    
    for (int i = 0; i < NumNodes; i++) {
        int nodeid = net.get_nodeid(i);
        SpeciesId[nodeid] = 0;
    }
    
  

}  


inline vector<vector<int>> speciation_model::get_repgenome(){return RepGenome;}
inline vector<vector<int>> speciation_model::get_phegenome(){return PheGenome;}
inline map<pair<int,int>,int> speciation_model::get_density_offspring(){return DensityOffspring;}
inline vector<vector<int>> speciation_model::get_geneticneigh(){return GeneticNeigh;}
inline double speciation_model::get_phenotype(int i){return Phenotype[i];}
inline double speciation_model::get_fitness(int i){return Fitness[i];}
inline tuple<int, int, int> speciation_model::get_parents(int i){return OffspringParents[i];}
inline int speciation_model::get_numspecies_parents(){

    set<int> SpeciesSet;
    for(auto &entry:OffspringParents){
        int species = get<2>(entry);
        SpeciesSet.insert(species);
    }
    return SpeciesSet.size();
}



inline int speciation_model::get_speciesid(int i){
    auto it = SpeciesId.find(i);
    return it -> second;
}

inline vector<int> speciation_model::get_speciesmembers(int i){
    auto it = SpeciesMembers.find(i);
    return it -> second;
}


void speciation_model::phenotype(){
    Phenotype.clear();
    Phenotype.reserve(PheGenome.size());
 
    for (auto& genome : PheGenome) {
        float fit = 0;
        for (int value : genome) {
            fit += value;
        }
        fit /= genome.size();
        Phenotype.push_back(fit);
        
      
    }
}


void speciation_model::compatibility_environmental(graph &net, const double &S, const int &SizeY, const int &SizeX, float &sigma, const float &Phenotype0){

    // =========================================================================
    // SCENARIO 1: HETEROGENEOUS ENVIRONMENT (Two Ecological Niches: E1 and E2)
    // To use this scenario, keep the block below active.
    // =========================================================================
     
  double Phenotype00 = 0.0; 

    if (sigma != 0) {
        Fitness.clear(); 

        for (int i = 0; i < net.get_numnodes(); i++) {
            int X = net.get_position(i).first;

            if (X < (SizeX / 2)) {
                Phenotype00 = Phenotype0 - 0.10; // Niche E1
            } else {
                Phenotype00 = Phenotype0 + 0.10; // Niche E2
            }

            double a = pow(double(Phenotype[i] - Phenotype00), 2);
            double d = 2 * pow(sigma, 2);
            double dd = -a / d;
            double compatibility = exp(dd);
            Fitness.push_back(compatibility);      
        }
    } else {
       
        for (int i = 0; i < net.get_numnodes(); i++) {
            Fitness.push_back(0);
        }
    }

   /* // =========================================================================
    // SCENARIO 2: HOMOGENEOUS ENVIRONMENT (Single Ecological Niche)
    // To use this scenario, uncomment this block and comment the one above.
    // =========================================================================

    
    /*if (sigma != 0) {

        Fitness.clear();

        for (int i = 0; i < net.get_numnodes(); i++) {

            double a = pow(double(Phenotype[i] - Phenotype0), 2);
            double d = 2 * pow(sigma, 2);
            double dd = -a / d;
            double compatibility = exp(dd);
            Fitness.push_back(compatibility);

        }

        } else {

            for (int i = 0; i < net.get_numnodes(); i++) {
                Fitness.push_back(0);
    
            } 
        } */ 

}




long double speciation_model::fitness_mean(int j, const int &SizeX, graph&net){
    double sumx1 = 0;
    double sumx2 = 0;
    int countx1 = 0;
    int countx2 = 0;
    for (int i = 0; i < net.get_numnodes(); i++) {
        int X = net.get_position(i).first;
        if ((X - 1) < (SizeX / 2)) {
            sumx1 += Fitness[i];
            countx1++;
        } else {
            sumx2 += Fitness[i];
            countx2++;
        }
    }
    
    long double FitnessMean1 = sumx1/countx1;
    long double FitnessMean2 = sumx2/countx2;
    int X = net.get_position(j).first;
    
    
    long double FitnessMean = (X - 1 < (SizeX / 2)) ? FitnessMean1 : FitnessMean2;

return FitnessMean;
}




inline double speciation_model::calculate_dissimilarity(int v, int u){
    double dissimilarity = 0.0;
    const vector<int>& genome1 = RepGenome[v];
    const vector<int>& genome2 = RepGenome[u];
    for(int i = 0; i < genome1.size(); ++i) {
        dissimilarity += abs(genome1[i] - genome2[i]); 
    }
   
    dissimilarity = dissimilarity/genome1.size();
   

    return dissimilarity;
}

void speciation_model::calculate_genetic_neighbors(graph &net, double G) {
    int NumNodes = net.get_numnodes();
    GeneticNeigh.clear();
    GeneticNeigh.resize(NumNodes); 

    for (int i = 0; i < NumNodes; ++i) {
        vector<int> Neighbors;
        vector<int> NeighborsIds = net.get_neighbor_ids(i);

        for (int NeighborId : NeighborsIds) {
            double dissimilarity = calculate_dissimilarity(i, NeighborId);
            if (dissimilarity <= G) {
                Neighbors.push_back(NeighborId);
            }
        }

        GeneticNeigh[i] = Neighbors;
    }
}


void speciation_model::reproduce(int father, int mother, const float&MutationRate, graph &net, const int &SizeX, const int &S){
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    vector<int> OffspringGenome(RepGenome[father].size());
    vector<int> OffspringPhen(PheGenome[father].size());
    
    for (int kc = 0; kc < RepGenome[father].size(); kc++) {
        double aux = dist(gen_);
        // Gene inheritance
        if (aux < 0.5) {
            OffspringGenome[kc] = RepGenome[father][kc];
        } else {
            OffspringGenome[kc] = RepGenome[mother][kc];
        }

        aux = dist(gen_);

        //Genetic mutation
        if (aux < MutationRate) {
            OffspringGenome[kc] = 1 - OffspringGenome[kc]; // Inverts the bit value
        }
    
    }

    RepGenAux.push_back(OffspringGenome);

    for (int kc = 0; kc < PheGenome[father].size(); kc++) {
        double aux = dist(gen_);
        if (aux < 0.5) {
            OffspringPhen[kc] = PheGenome[father][kc];
        } else {
            OffspringPhen[kc] = PheGenome[mother][kc];
        }

        aux = dist(gen_);  

        //Genetic mutation
        if (aux < MutationRate) {
            OffspringPhen[kc] = 1 - OffspringPhen[kc]; 
        }
    }
    
    PheGenAux.push_back(OffspringPhen);

    OffspringParents.push_back(make_tuple(net.get_nodeid(father), net.get_nodeid(mother),get_speciesid(father))); 

    double diff_val = 0;
    double diff_total=0;
    int xfather = net.get_position(father).first;
    int xmother = net.get_position(mother).first;
    diff_total = Phenotype[father] - Phenotype[mother];
    DIFFTOTAL.push_back(diff_total);

    if((SizeX/2 - S < xfather) && (xfather < SizeX/2 + S)){
        if ((SizeX/2 - S < xmother) && (xmother < SizeX/2 + S)){
            diff_val = Phenotype[father] - Phenotype[mother];
            DIFF.push_back(diff_val);
        }
    }


   

}
    

//Without periodic conditions
inline void speciation_model::offspring_region(graph& net, int i, const double& S, const int& SizeX, const int& SizeY) {
    pair<int, int> FatherPosition = net.get_position(i);
    const vector<pair<int, int>>& RegionS = net.get_regionS();
   

    vector<vector<bool>> grid(SizeX, vector<bool>(SizeY, false));
    vector<pair<int, int>> PossiblePositions;

    for (const auto& regionPos : RegionS) {

        int OffspringX = FatherPosition.first + regionPos.first; 
        int OffspringY = FatherPosition.second + regionPos.second;


        if (OffspringX >= 0 && OffspringX < SizeX && OffspringY >= 0 && OffspringY < SizeY) {
            grid[OffspringX][OffspringY] = true;
            PossiblePositions.push_back({OffspringX, OffspringY});
        }
    }

    if (!PossiblePositions.empty()) {
        uniform_int_distribution<int> pos_dist(0,PossiblePositions.size() - 1);
        size_t idx = pos_dist(gen_);
        const auto& pos = PossiblePositions[idx];
        DensityOffspring[pos]++;
        OffspringPosition.push_back(pos);

        return; 
    
    }

}



void speciation_model::new_graph(graph &net, const int &SizeX, const int &SizeY) {
    net.clear_node_list();  //Cleans the old graph

   
    if (OffspringParents.size() != OffspringPosition.size()) {
        cerr << "Erro" << endl;
        return;
    }

    for( int i = 0; i < OffspringParents.size(); i++){
        const auto& position = OffspringPosition[i];
        OffspringList.push_back(node(i, position));  //Adds all nodes to a List OffspringList
    }
      
    //Create a new graph
    for (const auto& Offspring : OffspringList) {
        net.add_node(Offspring);  
    }
    net.add_neighbors(SizeX, SizeY); 

}


inline void speciation_model::update_genomes_from_offspring(graph &net, const int &SizeRepGenome, const int &SizePheGenome){
    int NumNodes = net.get_numnodes();  
    RepGenome.assign(NumNodes, vector<int>(SizeRepGenome, 0));  
    PheGenome.assign(NumNodes, vector<int>(SizePheGenome, 0));

    for (int i = 0; i < NumNodes; ++i) {
        RepGenome[i] = RepGenAux[i];
        PheGenome[i] = PheGenAux[i];
    }
    
    RepGenAux.clear(); 
    PheGenAux.clear();

}
bool speciation_model::try_reproduce(graph &net, int i, const float &qmating, const float &MutationRate, const double &S,
    const int &SizeX, const int &SizeY, double rho0, const int &NumNodes, float &kmed) {
        
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    // Population constants for recruitment control
    double Nt = static_cast<double>(net.get_numnodes()); 
    double N0 = static_cast<double>(NumNodes);     
    
    // Mating probability based on individual fitness and mating cost (qmating)
    double qnorm = 1.0 - Fitness[i] + qmating * Fitness[i];
    
    vector<int> neighborsI = GeneticNeigh[i];
    if (neighborsI.empty()) return false; 
    
    // Check if reproduction occurs for this individual
    if (dist(gen_) < qnorm) return false;
    
    std::uniform_int_distribution<int> int_dist(0, neighborsI.size() - 1);
    int k = neighborsI[int_dist(gen_)];
    reproduce(i, k, MutationRate, net, SizeX, S);
    offspring_region(net, i, S, SizeX, SizeY);
    
    // Extra Reproduction (Population Recruitment)
    // Calculates how many extra offspring are needed to maintain population stability
    int Factor = static_cast<int>(std::ceil(N0 / Nt));
    int Radd = std::max(0, Factor); 
    
    if (net.get_degree(i) < kmed) {
        double Pextra = 1.0 - qnorm; 
         
        for (int r = 0; r < Radd; r++) {
            if (dist(gen_) < Pextra) { 
                reproduce(i, k, MutationRate, net, SizeX, S);
                offspring_region(net, i, S, SizeX, SizeY);
            }
        }
    }
    
    return true; 
}


void speciation_model::dynamics(graph &net, const float &qmating, const float &MutationRate, 
                                const double &S, const int &SizeX, const int &SizeY, 
                                const double &G, const int &NumNodes, const int &SizeRepGenome, const int &SizePheGenome, const float &Phenotype0, float &sigma,
                                 int currentGeneration, const int &maxGenerations, float &kmed) {
   
    OffspringPosition.clear();
    OffspringList.clear();
    DensityOffspring.clear();
    OffspringParents.clear();
    DIFF.clear();
    
    double rho0 = static_cast<double>(net.get_numnodes()) / (SizeX * SizeY);
    
    for (int i = 0; i < net.get_numnodes(); i++) {
        if (net.get_degree(i) == 0) continue;
        
        if (try_reproduce(net, i, qmating, MutationRate, S, SizeX, SizeY, rho0, NumNodes, kmed)) {
            continue;
        }

    }

    new_graph(net, SizeX, SizeY);
    update_genomes_from_offspring(net, SizeRepGenome, SizePheGenome);
    phenotype();
    compatibility_environmental(net, S, SizeY, SizeX, sigma, Phenotype0);
    calculate_genetic_neighbors(net, G);
}


pair<map<int, int>, int> speciation_model::check_species(graph &net, const double &G) {
    int NumNodes = net.get_numnodes();
    vector<bool> visited(NumNodes, false);
    SpeciesId.clear(); 
    SpeciesMembers.clear();
    int NumSpecies = 0;
    SpeciesMembers.clear();
    

    for (int i = 0; i < NumNodes; ++i) {
        if (visited[i]) continue;

        //Initialize BFS
        queue<int> q;
        q.push(i);
        visited[i] = true;
        SpeciesId[i] = NumSpecies + 1;
        SpeciesMembers[NumSpecies].push_back(i); 
        
        //Expand specie to the first node
        while (!q.empty()) {
            int current = q.front();
            q.pop();

            for(int j = 0; j < NumNodes; ++j) {
                if (!visited[j] && calculate_dissimilarity(current, j) <= G) {
                    visited[j] = true;
                    SpeciesId[net.get_nodeid(j)] = NumSpecies + 1;
                    q.push(j);
                    SpeciesMembers[NumSpecies].push_back(j);

                }
            }
        }

        NumSpecies++;
    }

   
    return {SpeciesId, NumSpecies};
}







#endif
