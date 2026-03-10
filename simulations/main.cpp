/*******************
initialize: 04/07/2024
update: 28/02/2026
********************/

/*********************************************************************************************************
 * SIMULATION MAIN --  Have the speciation model, handling parameters, initialization, 
 *                    and the generational loop.
 *********************************************************************************************************/

 #include <iostream>
 #include <cstdlib>
 #include <fstream>
 #include <string>
 #include <filesystem>
 #include <random>
 #include <functional>
 #include "graph.h"
 #include "speciation_model.h"
 #include "auxiliar_functions.h"
 
 using namespace std;
 namespace fs = std::filesystem;
 
 int main(int argc, char *argv[]) {
    
     if (argc < 5 || argc > 6) {
         cerr << "Usage: " << argv[0] << " <NumNodes> <G> <sigma> <directory> [run_id]" << endl;
         return 1;
     }
 
     string run_id = (argc == 6) ? argv[5] : "";
     std::mt19937 gen;
     
     unsigned int seed_value = (!run_id.empty()) ? 
         std::hash<std::string>{}(run_id + argv[2] + argv[3]) : std::random_device{}();
     gen.seed(seed_value);
     
     // --- 1. PARAMETERS ---
     int NumNodes     = atoi(argv[1]);
     double G         = atof(argv[2]);
     float sigma      = atof(argv[3]);
     string directory = argv[4];
 
     // Model Constants
     const float MutationRate = 0.00025f;
     const float qmating      = 0.37f;
     const int t      = 5000;
     const int SizeX = 100, SizeY = 100;
     const double S           = 6.0;
     const int SizeRepGenome  = 1500;
     const int SizePheGenome  = 100;
     const float Phenotype0   = 0.50f;
 
     if (directory.back() != '/' && directory.back() != '\\') directory += '/';
 
     // --- 2. DIRECTORIES ---
     try {
         fs::create_directories(directory + "All");
         fs::create_directories(directory + "diff");
     } catch (const fs::filesystem_error& e) {
         cerr << "Error creating directories: " << e.what() << endl;
         return 1;
     }
 
     string global_file = directory + "evolution_global" + run_id + ".dat";
     ofstream outFile(global_file);
     outFile << "Generation\tSpeciesCount\tPopulationSize\n";
 
     // --- 3. INITIALIZATION ---
     graph net;
     net.generate_random_positions(NumNodes, SizeX, SizeY, gen); 
     net.region_S(S, SizeX, SizeY);
     net.add_neighbors(SizeX, SizeY);
 
     float rho = NumNodes / (float)(SizeX * SizeY);   
     float kmed = rho * 3.14159f * S * S; 
     
     speciation_model speciation(net, SizeRepGenome, SizePheGenome, gen); 
     auxiliar_functions auxiliar;
 
     speciation.phenotype();
     speciation.compatibility_environmental(net, S, SizeY, SizeX, sigma, Phenotype0);
     speciation.calculate_genetic_neighbors(net, G);
 
     auto [SpeciesId, NumSpecies] = speciation.check_species(net, G); 
    
     // --- 4. MAIN SIMULATION LOOP ---
    for (int i = 1; i <= t; i++) {
       speciation.dynamics(net, qmating, MutationRate, S, SizeX, SizeY, G, NumNodes,
                           SizeRepGenome, SizePheGenome, Phenotype0,
                           sigma, i, t, kmed);
    
       tie(SpeciesId, NumSpecies) = speciation.check_species(net, G);
       
       string allFile = directory + "All/data" + to_string(i) + ".dat";
       string diffPath = directory + "diff/"; 

       auxiliar.save_node_positions_and_species(allFile, net, SpeciesId, speciation);
       //auxiliar.save_diff_phenotype_analysis(speciation, diffPath, i);
       
       cerr << "t: " << i << " | Species: " << NumSpecies << " | Pop: " << net.get_numnodes() << endl;
       
       outFile << i << "\t" << NumSpecies << "\t" << net.get_numnodes() << endl; 
    }

    string PheFile = directory + "phenotype/phegenome" + ".dat";
    string diffPath = directory + "diff/"; 

    auxiliar.save_phe_genomes(speciation, PheFile);
    auxiliar.save_diff_phenotype_analysis(speciation, diffPath, 5000);
}
