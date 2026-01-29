#include "Protein.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

const std::vector<std::string> PROTEIN_NAMES = {"Phe","Phe","Leu","Leu","Ser","Ser","Ser","Ser","Tyr","Tyr","Stop","Stop","Cys","Cys","Stop","Trp","Leu","Leu","Leu","Leu","Pro","Pro","Pro","Pro","His","His","Gln","Gln","Arg","Arg","Arg","Arg","Ile","Ile","Ile","Met","Thr","Thr","Thr","Thr","Asn","Asn","Lys","Lys","Ser","Ser","Arg","Arg","Val","Val","Val","Val","Ala","Ala","Ala","Ala","Asp","Asp","Glu","Glu","Gly","Gly","Gly","Gly"};
const std::vector<std::string> CODONS = {"UUU","UUC","UUA","UUG","UCU","UCC","UCA","UCG","UAU","UAC","UAA","UAG","UGU","UGC","UGA","UGG","CUU","CUC","CUA","CUG","CCU","CCC","CCA","CCG","CAU","CAC","CAA","CAG","CGU","CGC","CGA","CGG","AUU","AUC","AUA","AUG","ACU","ACC","ACA","ACG","AAU","AAC","AAA","AAG","AGU","AGC","AGA","AGG","GUU","GUC","GUA","GUG","GCU","GCC","GCA","GCG","GAU","GAC","GAA","GAG","GGU","GGC","GGA","GGG"};


/**
 * @brief Construct a default Protein
 * 
 */
Protein::Protein() {
    _sourceSpecies = "Unknown";
    _codon = "";
}

/**
 * @brief Construct a protein from inputted species name and codon
 * 
 */
Protein::Protein(string speciesName, string codon) {
    _sourceSpecies = speciesName;
    _codon = codon;
}

/**
 * @brief return the protein name from the codon
 * 
 * @return std::string 
 */
string Protein::findProtein() const {
    for (size_t i = 0; i < CODONS.size(); i++) {
        if (_codon == CODONS.at(i)) {
            return PROTEIN_NAMES.at(i);
        }
    }

    return "?";
}

/**
 * @brief Get the Source Species object
 * 
 * @return std::string 
 */
string Protein::getSourceSpecies() const {
    return _sourceSpecies;
}

/**
 * @brief Get the Codon string
 * 
 * @return std::string 
 */
string Protein::getCodon() const {
    return _codon;
}