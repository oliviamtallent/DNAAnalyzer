#ifndef DNASTRAND_H
#define DNASTRAND_H

#include <string>
#include <vector>
#include "Protein.h"
#include <SFML/Graphics.hpp>

class DNAStrand {   
    public:
        /**
         * @brief Construct a new DNAStrand object
         * 
         */
        DNAStrand();

        /**
         * @brief Construct new DNAStrand object from species name and sequence
         * 
         */
        DNAStrand(std::string, std::string, int);

        /**
         * @brief Copy constructor
         * 
         * @param copy 
         */
        DNAStrand(const DNAStrand& copy);

        /**
         * @brief Destroy the DNAStrand object
         * 
         */
        ~DNAStrand();

        /**
         * @brief copy assignment operator
         * 
         * @param other 
         * @return DNAStrand& 
         */
        DNAStrand& operator=(const DNAStrand& other);

        /**
         * @brief helper to deallocate memory of dna strand
         * 
         */
        void deallocate();

        /**
         * @brief helper to deep copy memory of dna strand
         * 
         */
        void deepCopy(const DNAStrand&);

        /**
         * @brief Helper function to run setup code for creating the sequence variables
         * 
         */
        void setupData();

        /**
         * @brief Maps the current _sequence to its associated nucleotide pair sequence
         * 
         */
        void createPairSequence();

        /**
         * @brief Maps the current _sequence to a vector of its associated codons
         * 
         */
        void createCodonSequence();

        /**
         * @brief Maps the current _codonSequence to a vector of its associated proteins
         * 
         */
        void createProteinSequence();

        /**
         * @brief Get the Sequence object
         * 
         * @return std::string const DNA Sequence
         */
        std::string getSequence() const;

        /**
         * @brief Set the Sequence object
         * 
         */
        void setSequence(std::string);

        /**
         * @brief Get the Pair Sequence object
         * 
         * @return std::string pair sequence
         */
        std::string getPairSequence() const;

        /**
         * @brief Get the Protein Sequence object
         * 
         * @return vector<Protein> 
         */
        std::vector<Protein*> getProteinSequence() const;

        /**
         * @brief Modify the nucelotide at given index
         * 
         */
        void modifyNucleotide(int, char);

        /**
         * @brief Modify the codon at the given index
         * 
         */
        void modifyCodon(int, std::string);

        /**
         * @brief  Find the percentage of the two DNA strands that share similarity
         * 
         */
        double compareDNA(DNAStrand&) const;

        /**
         * @brief Find similarity clusters of nucleotides
         * 
         * @return std::vector<int> clusters
         */
        std::vector<int> findClusters(DNAStrand&) const;

        /**
         * @brief Find the percentage of the protein sequences of the two DNA strands that share similarity
         * 
         * @return double 
         */
        double compareProteins(DNAStrand&) const;

        /**
         * @brief Find similarity clusters of proteins
         * 
         * @return std::vector<int> clusters
         */
        std::vector<int> findProteinClusters(DNAStrand&) const;

        /**
         * @brief Uses SMFL Library to display the given DNA Strand
         * 
         */
        void drawNucleotides(sf::RenderWindow&, sf::Vector2f, int);

        /**
         * @brief Uses SMFL Library to display the given Protein Sequence
         * 
         */
        void drawProteins(sf::RenderWindow&, sf::Vector2f, int);

        /**
         * @brief highlights the nucleotide clusters found
         * 
         */
        void highlightNucleotideClusters(sf::RenderWindow&, sf::Vector2f, int, std::vector<int>&);

        /**
         * @brief highlights the protein clusters found
         * 
         */
        void highlightProteinClusters(sf::RenderWindow&, sf::Vector2f, int, std::vector<int>&);
    private:
        std::string _sourceSpecies;
        int _class;
        std::string _sequence;
        std::string _pairSequence;
        std::vector<std::string> _codonSequence;
        std::vector<Protein*> _proteinSequence;
};

std::ostream& operator<<(std::ostream&, const DNAStrand&);

#endif