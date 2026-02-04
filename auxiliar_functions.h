 /*******************
initialize:04/07/2024
update: 04/02/2026
Lara D. Hissa

********************/ 



/*********************************************************************************************************
 * CLASS AUXILIAR_FUNCTIONS
 * Purpose: Provides utility methods for data export, visualization, and statistical analysis.
 *********************************************************************************************************/

 class auxiliar_functions {
    public:
        // --- Printing & Console Debug ---
        void print_genomes(speciation_model &spec_model); // Print the RepGenomes of all individuals;
        void print_density_offspring_position(speciation_model &spec_model);
        void print_offspring_position(speciation_model &spec_model); // Prints the positions of all the offspring;

        // --- Data Export (File I/O) ---
        void save_rep_genomes(speciation_model &spec_model, const string& filename); 
        void save_node_positions_and_species(const string& filename, graph& net, map<int, int>& SpeciesId, speciation_model &spec_model); // Saves x,y, species, degree, phenotype and fitness;
        void save_density(speciation_model &spec_model, const string &filename, graph &net); // Records offspring density per grid site;
        void save_genetic_neighbors(speciation_model &spec_model, const string &filename); // Exports the list of genetic neighbors;
        void save_tudo(const string &filename, graph &net, map<int, int>& SpeciesId, speciation_model &spec_model); // Comprehensive export for network analysis;
        
        // --- Statistical & Evolutionary Analysis ---
        map<int, int> get_all_node_degrees(speciation_model &spec_model);
        void save_diff_phenotype_analysis(speciation_model &spec_model, const string& pasta, int iteration); // Exports boundary vs total phenotype variance;
};

// =================================================================================================
//                                     IMPLEMENTATIONS
// =================================================================================================



inline void auxiliar_functions::save_diff_phenotype_analysis(speciation_model &spec_model, const string& pasta, int iteration) {
    if (spec_model.DIFF.empty() && spec_model.DIFFTOTAL.empty()) return;

    // 1. Save DIFF (Boundary)
    string pathBoundary = pasta + "diff_boundary_" + to_string(iteration) + ".txt";
    ofstream fileBoundary(pathBoundary);
    if (fileBoundary.is_open()) {
        for (const auto& value : spec_model.DIFF) fileBoundary << value << "\n";
        fileBoundary.close();
        spec_model.DIFF.clear(); 
    }

    // 2. Save DIFFTOTAL (Full Population)
    string pathTotal = pasta + "diff_total_" + to_string(iteration) + ".txt";
    ofstream fileTotal(pathTotal);
    if (fileTotal.is_open()) {
        for (const auto& value : spec_model.DIFFTOTAL) fileTotal << value << "\n";
        fileTotal.close();
        spec_model.DIFFTOTAL.clear(); 
    }
}

inline void auxiliar_functions::print_genomes(speciation_model &spec_model) {
    vector<vector<int>> RepGenome = spec_model.get_repgenome();
    for(int i = 0; i < RepGenome.size(); i++){
        cout << "Node ID " << i << " Genome: ";
        for(int j = 0; j < RepGenome[i].size(); j++) cout << RepGenome[i][j];
        cout << "\n************************" << endl;
    }
}

void auxiliar_functions::save_rep_genomes(speciation_model &spec_model, const string& filename) {
    vector<vector<int>> RepGenome = spec_model.get_repgenome();
    ofstream outFile(filename);
    if (!outFile.is_open()) return;

    for (const auto& genome : RepGenome) {
        for (int bit : genome) outFile << bit;
        outFile << endl;
    }
    outFile.close();
}

void auxiliar_functions::save_node_positions_and_species(const string& filename, graph& net, map<int, int>& SpeciesId, speciation_model &spec_model) {
    ofstream outFile(filename);
    if (!outFile.is_open()) return;

    for (const auto& node_entry : SpeciesId) {
        int NodeId = node_entry.first;
        int SpId = node_entry.second;
        auto pos = net.get_position(NodeId);
        outFile << NodeId << " " << pos.first << " " << pos.second << " " 
                << SpId << " " << net.get_degree(NodeId) << " " 
                << spec_model.get_phenotype(NodeId) << " " << spec_model.get_fitness(NodeId) << endl;
    }
    outFile.close();
}

void auxiliar_functions::save_density(speciation_model &spec_model, const string &filename, graph &net) {
    map<pair<int,int>,int> DensityOffspring = spec_model.get_density_offspring();
    ofstream outFile(filename);
    if (!outFile.is_open()) return;

    for(auto const& [pos, count] : DensityOffspring) {
        outFile << pos.first << " " << pos.second << " " << count << endl;
    }
    outFile.close();
}

map<int, int> auxiliar_functions::get_all_node_degrees(speciation_model &spec_model) {
    vector<vector<int>> GeneticNeigh = spec_model.get_geneticneigh();
    map<int, int> NodeDegrees;
    for (size_t i = 0; i < GeneticNeigh.size(); i++) {
        NodeDegrees[i] = GeneticNeigh[i].size();
    }
    return NodeDegrees;
}

void auxiliar_functions::save_genetic_neighbors(speciation_model &spec_model, const string &filename) {
    vector<vector<int>> GeneticNeigh = spec_model.get_geneticneigh();
    ofstream outFile(filename);
    if (!outFile.is_open()) return;

    for (size_t i = 0; i < GeneticNeigh.size(); i++) {
        outFile << i << " " << GeneticNeigh[i].size() << " -- ";
        for (int neighbor : GeneticNeigh[i]) outFile << neighbor << " ";
        outFile << endl;
    }
    outFile.close();
}

inline void auxiliar_functions::print_density_offspring_position(speciation_model &spec_model) {
    map<pair<int,int>,int> DensityOffspring = spec_model.get_density_offspring();
    for(auto const& [pos, count] : DensityOffspring) {
        cout << "Position(" << pos.first << "," << pos.second << "): " << count << " offspring" << endl; 
    }
}

void auxiliar_functions::save_tudo(const string &filename, graph &net, map<int, int>& SpeciesId, speciation_model &spec_model) {
    vector<vector<int>> GeneticNeigh = spec_model.get_geneticneigh();
    ofstream outFile(filename);
    if (!outFile.is_open()) return;

    for (const auto& node_entry : SpeciesId) {
        int i = node_entry.first;
        if (GeneticNeigh[i].empty()) continue;

        auto pos_i = net.get_position(i);
        double pheno_i = spec_model.get_phenotype(i);

        for (int neighbor : GeneticNeigh[i]) {
            auto pos_n = net.get_position(neighbor);
            outFile << i << " " << neighbor << " " << pos_i.first << " " << pos_i.second << " "
                    << pos_n.first << " " << pos_n.second << " " << pheno_i << " "
                    << spec_model.get_phenotype(neighbor) << " " << node_entry.second << " " 
                    << SpeciesId[neighbor] << endl;
        }
    }
    outFile.close();
}