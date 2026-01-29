#ifndef PROTEIN_H
#define PROTEIN_H

#include <string>
#include <vector>

class Protein {
    public:
        /**
         * @brief Construct a default Protein
         * 
         */
        Protein();
        
        /**
         * @brief Construct a protein from inputted species name and codon
         * 
         */
        Protein(std::string, std::string);

        /**
         * @brief return the protein name from the codon
         * 
         * @return std::string 
         */
        std::string findProtein() const;

        /**
         * @brief Get the Source Species object
         * 
         * @return std::string 
         */
        std::string getSourceSpecies() const;

        /**
         * @brief Get the Codon string
         * 
         * @return std::string 
         */
        std::string getCodon() const;

    private:
        std::string _sourceSpecies;
        std::string _codon;
};

#endif