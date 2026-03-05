#ifndef SPECIES_TREE
#define SPECIES_TREE

#include <iostream>
#include <map>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <numeric>  
#include "speciation_model.h"
#include "graph.h"

using namespace std;

class species_tree {
    private:
        
        vector<int> ExtinctionTime;
        vector<int> ExtinctSize;
        vector<int> SisterSp;
        vector<int> AllParentsSp;
        
      
        vector<int> OldId;     
        vector<int> OldSize;   
        
        int fullcount = 0;    
        int Oldfullcount = 0;
          
        vector<vector<int>> TimeMatrix;     
        vector<vector<int>> OldTimeMatrix;  

    public:
        void MRCAT(graph &net, speciation_model &spec,
                   int &NumSpeciesOld, int &NumSpeciesActual,
                   string &directory, int CurrentTime, int &TotalTime);
        void initialize_MRCAT();
        void save_full_phylogeny(const string &directory);
};

void species_tree::initialize_MRCAT(){
    int maxSp = 15000;
    ExtinctionTime.resize(maxSp,0);
    ExtinctSize.resize(maxSp,0);
    SisterSp.resize(maxSp,0);
    AllParentsSp.resize(maxSp,0);

    OldId.resize(maxSp);
    iota(OldId.begin(), OldId.end(),0);

    TimeMatrix.resize(maxSp, vector<int>(maxSp,0));
    OldTimeMatrix.resize(maxSp, vector<int>(maxSp,0));

    OldSize.resize(maxSp,0);
}

void species_tree::save_full_phylogeny(const string &directory) {
    int igtfull = Oldfullcount;

    // === matrix-fullphy.dat ===
    {
        ofstream out(directory + "matrix-fullphy.dat", ios::trunc);
        for (int i = 1; i <= igtfull; ++i) {
            for (int j = 1; j <= igtfull; ++j) {
                out << TimeMatrix[i][j] << " ";
            }
            out << "\n";
        }
    }

    // === ext-times.dat ===
    {
        ofstream out(directory + "ext-times.dat", ios::trunc);
        for (int i = 1; i <= igtfull; ++i) {
            out << ExtinctionTime[i] << "\n";
        }
    }

    // === ext-sizes.dat ===
    {
        ofstream out(directory + "ext-sizes.dat", ios::trunc);
        for (int i = 1; i <= igtfull; ++i) {
            out << ExtinctSize[i] << " " << ExtinctionTime[i] << " " << SisterSp[i] << "\n";
        }
    }

    // === hybridizations.dat ===
    {
        ofstream out(directory + "hybridizations.dat", ios::trunc);
        for (int ii = 1; ii <= igtfull; ++ii) {
            int ij = igtfull - ii + 1;
            if (SisterSp[ij] > 0) {
                out << ExtinctionTime[ij] << "\n"; // ou outro critério se precisar
            }
        }
    }
}


/*void species_tree::save_matrix_time(const string &path){
    ofstream out(path, ios::trunc); 
    int size = Oldfullcount;        // total (vivas + extintas) no passo atual
    for (int i = 1; i <= size; ++i){
        for (int j = 1; j <= size; ++j){
            out << TimeMatrix[i][j] << " ";
        }
        out << "\n";
    }
}*/

void species_tree::MRCAT(graph &net, speciation_model &spec, int &NumSpeciesOld, int &NumSpeciesActual, string &directory, int CurrentTime, int &TotalTime){
    /****************************************************************************
     * 1) Ancestral mais provável + irmãs e contagem de ramificações
     ****************************************************************************/
    vector<int> ParentSp(NumSpeciesActual + 1, 0);
    vector<int> sister1(NumSpeciesActual + 1, 0);
    vector<int> sister2(NumSpeciesActual + 1, 0);
    vector<int> sister3(NumSpeciesActual + 1, 0);
    vector<int> multi(NumSpeciesOld + 1, 0);

    for (int i = 1; i <= NumSpeciesActual; ++i) {
        map<int,int> hist;
        for (int individual : spec.get_speciesmembers(i-1)) {
            auto [father, mother, Fathersp] = spec.get_parents(individual);
            hist[Fathersp]++;
        }
        vector<pair<int,int>> vec(hist.begin(), hist.end());
        sort(vec.begin(), vec.end(), [](auto &a, auto &b){
            if (a.second != b.second) return a.second > b.second;
            return a.first < b.first;
        });
        if (!vec.empty()) {
            ParentSp[i] = vec[0].first;
            multi[vec[0].first]++;
            if (vec.size() > 1) sister1[i] = vec[1].first;
            if (vec.size() > 2) sister2[i] = vec[2].first;
            if (vec.size() > 3) sister3[i] = vec[3].first;
        }
    }

    /****************************************************************************
     * 2) Marca extinções verdadeiras (no passo anterior)
     ****************************************************************************/
    vector<int> IsExtinct(NumSpeciesOld + 1, 0);
    for (int i = 1; i <= NumSpeciesOld; ++i) {
        bool found = false;
        for (int j = 1; j <= NumSpeciesActual; ++j) {
            if (ParentSp[j] == i || sister1[j] == i || sister2[j] == i || sister3[j] == i) {
                found = true; break;
            }
        }
        if (!found) IsExtinct[i] = 1;
    }

    /****************************************************************************
     * 3) Atribui IDs históricos consistentes
     ****************************************************************************/
    int MaxSpecies = 15000;
    vector<int> TrueSizes(MaxSpecies+2, 0);     // por ID histórico
    vector<int> TrueParentId(MaxSpecies+2, 0);  // por ID histórico
    vector<int> TrueSpeciesId(MaxSpecies+2, 0); // índice atual -> ID histórico
    OldId.resize(MaxSpecies+2, 0);

    if (fullcount == 0) fullcount = NumSpeciesOld;

    for (int i = 1; i <= NumSpeciesActual; ++i) {
        int p = ParentSp[i];
        int CurrentSize = (int)spec.get_speciesmembers(i-1).size();
        if (multi[p] == 1) {
            TrueSpeciesId[i] = OldId[p];
            TrueSizes[TrueSpeciesId[i]] = CurrentSize;
            TrueParentId[TrueSpeciesId[i]] = TrueSpeciesId[i];
        } else if (multi[p] > 1) {
            fullcount++;
            TrueSpeciesId[i] = fullcount;
            TrueSizes[TrueSpeciesId[i]] = CurrentSize;
            TrueParentId[TrueSpeciesId[i]] = OldId[p];
        }
    }

    // grava histórico (por ID histórico)
    {
        ofstream out(directory + "global_abundances.dat", ios::app);
        for (int k = 1; k <= fullcount; ++k) out << TrueSizes[k] << " ";
        out << "\n";
    }
    {
        ofstream out(directory + "global_parents.dat", ios::app);
        for (int k = 1; k <= fullcount; ++k) out << TrueParentId[k] << " ";
        out << "\n";
    }

    /****************************************************************************
     * 4) Preparação para montar a nova lista (extantes + extintas copiadas)
     *    - snapshot por ÍNDICE do passo anterior (igual ao Fortran: *_old)
     *    - mapeia tamanho antigo por índice: PrevExtantSizes[i] = OldSize[OldId[i]]
     *    - zera arrays atuais (extantes devem ter ext_time=0)
     ****************************************************************************/

    vector<int> ext_time_old = ExtinctionTime;  
    vector<int> size_ext_old = ExtinctSize;     
    vector<int> sister_sp_old = SisterSp;      

    // tamanho antigo por ÍNDICE de espécie viva do passo anterior (1..NumSpeciesOld)
    vector<int> PrevExtantSizes(MaxSpecies+2, 0);
    for (int i = 1; i <= NumSpeciesOld; ++i) {
        int hid = (i < (int)OldId.size() ? OldId[i] : 0);
        if (hid >= 1 && hid < (int)OldSize.size())
            PrevExtantSizes[i] = OldSize[hid];
    }

    // ZERA arrays atuais (como no Fortran: ext_time = 0; sister_sp = 0)
    if ((int)ExtinctionTime.size() < MaxSpecies+2) ExtinctionTime.resize(MaxSpecies+2, 0);
    if ((int)ExtinctSize.size()    < MaxSpecies+2) ExtinctSize.resize(MaxSpecies+2, 0);
    if ((int)SisterSp.size()       < MaxSpecies+2) SisterSp.resize(MaxSpecies+2, 0);
    std::fill(ExtinctionTime.begin(), ExtinctionTime.end(), 0);
    std::fill(ExtinctSize.begin(),    ExtinctSize.end(),    0);
    std::fill(SisterSp.begin(),       SisterSp.end(),       0);

    // parent_sp por índice atual (para matriz de tempos)
    AllParentsSp.assign(MaxSpecies+2, 0);
    for (int i = 1; i <= NumSpeciesActual; ++i) AllParentsSp[i] = ParentSp[i];

    int Newfullcount = NumSpeciesActual;

    /****************************************************************************
     * 5) Copia espécies extintas do passo anterior
     ****************************************************************************/
    for (int i = 1; i <= Oldfullcount; ++i) {
        bool found = false;
        for (int j = 1; j <= NumSpeciesActual; ++j) {
            if (i == ParentSp[j]) { found = true; break; }
        }
        if (!found) {
            Newfullcount++;
            if ((int)ExtinctionTime.size() <= Newfullcount) ExtinctionTime.resize(Newfullcount+2, 0);
            if ((int)ExtinctSize.size()    <= Newfullcount) ExtinctSize.resize(Newfullcount+2, 0);
            if ((int)SisterSp.size()       <= Newfullcount) SisterSp.resize(Newfullcount+2, 0);

            AllParentsSp[Newfullcount] = i;

            // herda valores antigos
            int et_old = (i < (int)ext_time_old.size() ? ext_time_old[i] : 0);
            int sz_old = (i < (int)size_ext_old.size() ? size_ext_old[i] : 0);
            int ss_old = (i < (int)sister_sp_old.size()? sister_sp_old[i]: 0);

            ExtinctionTime[Newfullcount] = et_old;
            ExtinctSize[Newfullcount]    = sz_old;

            if (ss_old > 0) SisterSp[Newfullcount] = Newfullcount;
            else            SisterSp[Newfullcount] = ss_old;

            // se era "nova extinção" (ainda sem tempo registrado), define agora
            if (ExtinctionTime[Newfullcount] == 0)
                ExtinctionTime[Newfullcount] = TotalTime - CurrentTime + 1;

            // se não havia tamanho salvo, usa o tamanho do passo anterior
            if (ExtinctSize[Newfullcount] == 0 && i <= NumSpeciesOld)
                ExtinctSize[Newfullcount] = PrevExtantSizes[i];

            // ajuste de sister_sp para espécies que eram realmente extintas nesse passo
            if (i <= NumSpeciesOld) {
                if (IsExtinct[i] == 1) SisterSp[Newfullcount] = -1;
                else                    SisterSp[Newfullcount] = Newfullcount;
            }
        }
    }

    // prepara estado para o próximo passo
    OldSize      = TrueSizes;     
    Oldfullcount = Newfullcount;

    /****************************************************************************
     * 6) Atualiza matriz de tempos entre espécies (extantes + extintas)
     ****************************************************************************/
    for (int i = 1; i <= Oldfullcount; ++i)
        for (int j = 1; j <= Oldfullcount; ++j)
            TimeMatrix[i][j] = 0;

    for (int i1 = 1; i1 <= Oldfullcount; ++i1) {
        int k1 = AllParentsSp[i1]; // índice antigo
        for (int i2 = i1+1; i2 <= Oldfullcount; ++i2) {
            int k2 = AllParentsSp[i2]; // índice antigo
            int val = 1;
            if (k1 >= 1 && k1 < (int)OldTimeMatrix.size() &&
                k2 >= 1 && k2 < (int)OldTimeMatrix.size()) {
                val = OldTimeMatrix[k1][k2] + 1;
            }
            TimeMatrix[i1][i2] = val;
            TimeMatrix[i2][i1] = val;
        }
    }

    OldTimeMatrix = TimeMatrix;   
    OldId         = TrueSpeciesId; 
}

#endif
