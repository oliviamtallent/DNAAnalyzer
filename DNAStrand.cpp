#include "DNAStrand.h"
#include "dna_functions.h"
#include <string>
#include <vector>
#include <iostream>
#include <SFML/Graphics.hpp>
#include <cmath>

using namespace std;

DNAStrand::DNAStrand() {
    _sourceSpecies = "Unknown";
    _class = 0;
    _sequence = "";
    _pairSequence = "";
    _codonSequence = {};
    _proteinSequence = {};
}

DNAStrand::DNAStrand(string speciesName, string dnaSequence, int classNum) {
    if (speciesName == "") {
        _sourceSpecies = "Unknown";
    }

    _sequence = dnaSequence;
    _class = classNum;
    setupData();
}

DNAStrand::DNAStrand(const DNAStrand& copy) {
    deepCopy(copy);
}

DNAStrand& DNAStrand::operator=(const DNAStrand& other) {
    if (&other == this) {
        return *this;
    }

    deallocate();
    deepCopy(other);

    return *this;
}

DNAStrand::~DNAStrand() {
    deallocate();
}

void DNAStrand::deallocate() {
    for (size_t i = 0; i < _proteinSequence.size(); i++) {
        delete _proteinSequence.at(i);
    }
    _proteinSequence.clear();
}

void DNAStrand::deepCopy(const DNAStrand& copy) {
    _sourceSpecies = copy._sourceSpecies;
    _sequence = copy._sequence;
    _class = copy._class;
    _pairSequence = copy._pairSequence;
    _codonSequence = copy._codonSequence;
    _proteinSequence.clear();

    for (size_t i = 0; i < copy.getProteinSequence().size(); i++) {
        _proteinSequence.push_back(new Protein(copy.getProteinSequence().at(i)->getSourceSpecies(), copy.getProteinSequence().at(i)->getCodon()));
    }
}

void DNAStrand::setupData() {
    createPairSequence();
    createCodonSequence();
    createProteinSequence();
}

/**
 * @brief Maps the current _sequence to its associated nucleotide pair sequence
 * 
 */
void DNAStrand::createPairSequence() {
    _pairSequence.clear();
    for (size_t i = 0; i < _sequence.length(); i++) {
        _pairSequence.push_back(getPair(_sequence.at(i))); 
    }
}

/**
 * @brief Maps the current _sequence to a vector of its associated codons
 * 
 */
void DNAStrand::createCodonSequence() {
    _codonSequence.clear();
    for (size_t i = 0; i < _pairSequence.length(); i += 3) {
        _codonSequence.push_back(_pairSequence.substr(i, 3));
    }
}

/**
 * @brief Maps the current _codonSequence to a vector of its associated proteins
 * 
 */
void DNAStrand::createProteinSequence() {
    _proteinSequence.clear();
    for (size_t i = 0; i < _codonSequence.size(); i++) {
        _proteinSequence.push_back(new Protein(_sourceSpecies, _codonSequence.at(i)));
    }
}

/**
 * @brief Get the Sequence object
 * 
 * @return std::string const DNA Sequence
 */
string DNAStrand::getSequence() const {
    return _sequence;
}

/**
 * @brief Set the Sequence object
 * 
 */
void DNAStrand::setSequence(string dnaSequence) {
    _sequence = dnaSequence;
}

/**
 * @brief Get the Pair Sequence object
 * 
 * @return std::string pair sequence
 */
string DNAStrand::getPairSequence() const {
    return _pairSequence;
}

/**
 * @brief Get the Protein Sequence object
 * 
 * @return vector<Protein> 
 */
vector<Protein*> DNAStrand::getProteinSequence() const {
    return _proteinSequence;
}

/**
 * @brief Modify the nucelotide at given index
 * 
 */
void DNAStrand::modifyNucleotide(int index, char nucleotide) {
    _sequence.at(index) = nucleotide;
    _pairSequence.at(index) = getPair(nucleotide);

    // find codon
    int codonIndex = index / 3;
    int codonPosition = index % 3;

    // edit codon
    string modifiedCodon = _codonSequence.at(codonIndex);
    modifiedCodon.at(codonPosition) = nucleotide;
    _codonSequence.at(codonIndex) = modifiedCodon;

    // edit protein
    delete _proteinSequence.at(codonIndex);
    _proteinSequence.at(codonIndex) = new Protein(_sourceSpecies, modifiedCodon);
}

/**
 * @brief Modify the codon at the given index
 * 
 */
void DNAStrand::modifyCodon(int index, string codon) {
    for (int i = index; i < 3; i++) {
        _sequence.at(i) = codon.at(i-index);
        _pairSequence.at(i) = getPair(codon.at(i-index));
    }
    
    // edit protein
    delete _proteinSequence.at(index);
    _proteinSequence.at(index) = new Protein(_sourceSpecies, codon);
}

/**
 * @brief  Find the percentage of the two DNA strands that share similarity
 * 
 */
double DNAStrand::compareDNA(DNAStrand& other) const {
    double matches = 0;

    // find strand end
    size_t end = this->_sequence.size();
    if (other._sequence.size() < end) {
        end = other._sequence.size();
    }

    // find number of matches in strand
    for (size_t i = 0; i < end; i++) {
        if (other.getSequence().at(i) == _sequence.at(i)) {
            matches += 1;
        }
    }

    // calculate percentage similarity
    return matches/(int)end * 100;
}

/**
 * @brief Find similarity clusters of nucleotides
 * 
 * @return std::vector<int> clusters
 */
vector<int> DNAStrand::findClusters(DNAStrand& other) const {
    // find end of strand
    size_t end = _sequence.size();
    if (other._sequence.size() < end) {
        end = other._sequence.size();
    }

    // create vector holding the similarity percentages for subsections of 5 nucleotides
    vector<double> similarityVals;
    for (size_t i = 0; i <= end - 5; i++) {
        double matches = 0;
        for (size_t j = 0; j < 5; j++) {
            if (_sequence[i+j] == other.getSequence()[i+j]) {
                matches += 1;
            }
        }
        similarityVals.push_back(matches/5);
    }

    // check for top 5
    vector<int> topClusters;
    if (similarityVals.size() < 5) {
        for (size_t i = 0; i < similarityVals.size(); i++) {
            topClusters.push_back((int)i);
        }
    } else {
        topClusters = {0, 5, 10, 15, 20};
        size_t currMin = getMin(similarityVals, topClusters);
        for (size_t i = 5; i < similarityVals.size(); i++) {
            if (similarityVals.at(i) > similarityVals.at(topClusters.at(currMin))) {
                // only add unique clusters
                // TODO could override if near cluster has higher similarity
                bool isUnique = true;
                for (size_t j = 0; j < topClusters.size(); j++) {
                    if (abs((int)i - (int)topClusters.at(j)) < 5) {
                        isUnique = false;
                        break;
                    }
                }
                if (isUnique) {
                    topClusters.at(currMin) = (int)i;
                    currMin = getMin(similarityVals, topClusters);
                }
            }
        }
    }

    return topClusters;
}

/**
 * @brief Find the percentage of the protein sequences of the two DNA strands that share similarity
 * 
 * @return double 
 */
double DNAStrand::compareProteins(DNAStrand& other) const {
    double matches = 0;

    // find end of strand
    size_t end = this->_proteinSequence.size();
    if (other.getProteinSequence().size() < end) {
        end = other.getProteinSequence().size();
    }

    // find number of matches
    for (size_t i = 0; i < end; i++) {
        if (other.getProteinSequence().at(i)->findProtein() == _proteinSequence.at(i)->findProtein() && _proteinSequence[i]->findProtein() != "?") {
            matches += 1;
        }
    }

    return matches/(int)end * 100;
}

/**
 * @brief Find similarity clusters of proteins
 * 
 * @return std::vector<int> clusters
 */
vector<int> DNAStrand::findProteinClusters(DNAStrand& other) const {
    // find end of strand
    size_t end = _proteinSequence.size();
    if (other._proteinSequence.size() < end) {
        end = other._proteinSequence.size();
    }

    // create vector of similarity percentages in clusters of 5 proteins.
    vector<double> similarityVals;
    for (size_t i = 0; i < end - 5; i++) {
        double matches = 0;
        for (size_t j = 0; j < 5; j++) {
            if (_proteinSequence[i + j]->findProtein() == other._proteinSequence[i + j]->findProtein() && _proteinSequence[i + j]->findProtein() != "?") {
                matches+=1;
            }
        }
        similarityVals.push_back(matches/5);
    }

    // check for top 5
    vector<int> topClusters;
    if (similarityVals.size() < 5) {
        for (size_t i = 0; i < similarityVals.size(); i++) {
            topClusters.push_back((int)i);
        }
    } else {
        topClusters = {0, 1, 2, 3, 4};
        size_t currMin = getMin(similarityVals, topClusters);
        for (size_t i = 5; i < similarityVals.size(); i++) {
            if (similarityVals.at(i) > similarityVals.at(topClusters.at(currMin))) {
                // only add unique clusters
                // TODO could override if near cluster has higher similarity
                bool isUnique = true;
                for (size_t j = 0; j < topClusters.size(); j++) {
                    if (abs((int)i - (int)topClusters.at(j)) < 5) {
                        isUnique = false;
                        break;
                    }
                }
                if (isUnique) {
                    topClusters.at(currMin) = (int)i;
                    currMin = getMin(similarityVals, topClusters);
                }
            }
        }
    }

    return topClusters;
}

std::ostream& operator<<(ostream& os, const DNAStrand& WH) {
    os << "DNA Strand: " << WH.getSequence() << endl;
    return os;
}

void DNAStrand::drawNucleotides(sf::RenderWindow& rw, sf::Vector2f startPosition, int scrollPos) {
    size_t end = rw.getSize().x / 12;
    if (end > _sequence.length() - scrollPos) {
        end = _sequence.length() - scrollPos;
    }
    for (size_t i = scrollPos; i < end + scrollPos; i++) {
        sf::ConvexShape shape;
        shape.setPointCount(5);
        shape.setPoint(0, sf::Vector2f(0, 0));
        shape.setPoint(1, sf::Vector2f(10, 0));
        shape.setPoint(2, sf::Vector2f(10, 20));
        if (_sequence.at(i) == 'A' || _sequence.at(i) == 'G') {
            shape.setPoint(3, sf::Vector2f(5, 14));
        } else {
            shape.setPoint(3, sf::Vector2f(5, 24));
        }
        shape.setPoint(4, sf::Vector2f(0, 20));
        
        // set color based on nucleotide
        if (_sequence.at(i) == 'A') {
            shape.setFillColor(sf::Color::Red);
        } else if (_sequence.at(i) == 'T') {
            shape.setFillColor(sf::Color::Blue);
        } else if (_sequence.at(i) == 'C') {
            shape.setFillColor(sf::Color::Green);
        } else {
            shape.setFillColor(sf::Color::Yellow);
        }

        shape.setPosition(sf::Vector2f(startPosition.x + (float)(i - scrollPos) * 12, startPosition.y));
        rw.draw(shape);
    }
}

void DNAStrand::drawProteins(sf::RenderWindow& rw, sf::Vector2f startPosition, int scrollPos) {
    int proteinScroll = scrollPos/3;
    size_t end = rw.getSize().x / 36 + 1;
    if (end > _proteinSequence.size() - proteinScroll) {
        end = _proteinSequence.size() - proteinScroll;
    }
    for (size_t i = proteinScroll; i < end + proteinScroll; i++) {
        sf::ConvexShape shape;
        shape.setPointCount(5);
        shape.setPoint(0, sf::Vector2f(0, 0));
        shape.setPoint(1, sf::Vector2f(30, 0));
        shape.setPoint(2, sf::Vector2f(30, 20));
        shape.setPoint(3, sf::Vector2f(15, 14));
        shape.setPoint(4, sf::Vector2f(0, 20));
        
        // set color based on protein
        string protein = _proteinSequence.at(i)->findProtein();
        if (protein == "Phe") {
            shape.setFillColor(sf::Color::Green);
        } else if (protein == "Leu") {
            shape.setFillColor(sf::Color::Blue);
        } else if (protein == "Stop") {
            shape.setFillColor(sf::Color::Red);
        } else if (protein == "Ser") {
            shape.setFillColor(sf::Color::Yellow);
        } else if (protein == "Tyr") {
            shape.setFillColor(sf::Color::Magenta);
        } else if (protein == "Trp") {
            shape.setFillColor(sf::Color(141, 186, 224));
        } else if (protein == "Pro") {
            shape.setFillColor(sf::Color::Cyan);
        } else if (protein == "His") {
            shape.setFillColor(sf::Color(210, 250, 211));
        } else if (protein == "Gln") {
            shape.setFillColor(sf::Color(247, 233, 151));
        } else if (protein == "Arg") {
            shape.setFillColor(sf::Color(207, 133, 6));
        } else if (protein == "Met") {
            shape.setFillColor(sf::Color(122, 40, 57));
        } else if (protein == "Ile") {
            shape.setFillColor(sf::Color(168, 63, 176));
        } else if (protein == "Thr") {
            shape.setFillColor(sf::Color(119, 71, 161));
        } else if (protein == "Cys") {
            shape.setFillColor(sf::Color(35, 84, 17));
        } else if (protein == "Asn") {
            shape.setFillColor(sf::Color(69, 135, 111));
        } else if (protein == "Val") {
            shape.setFillColor(sf::Color(140, 106, 11));
        } else if (protein == "Ala") {
            shape.setFillColor(sf::Color(77, 16, 29));
        } else if (protein == "Asp") {
            shape.setFillColor(sf::Color(94, 138, 135));
        } else if (protein == "Gly") {
            shape.setFillColor(sf::Color(205, 255, 97));
        } else if (protein == "Lys") {
            shape.setFillColor(sf::Color(41, 0, 92));
        } else {
            shape.setFillColor(sf::Color::White);
        }

        // find offset for partial scroll
        float fractional = (float)(scrollPos % 3) * 12;
        float xOffset = startPosition.x - fractional;

        shape.setPosition(sf::Vector2f(2 + xOffset + startPosition.x + (float)(i - proteinScroll) * 36, startPosition.y));
        rw.draw(shape);
    }
}

void DNAStrand::highlightNucleotideClusters(sf::RenderWindow& rw, sf::Vector2f startPosition, int scrollPos, vector<int>& clusterIndexes) {
    size_t end = rw.getSize().x / 12;
    if (end > _sequence.length() - scrollPos) {
        end = _sequence.length() - scrollPos;
    }

    for (size_t i = 0; i < clusterIndexes.size(); i++) {
        int clusterStart = clusterIndexes.at(i);
        int clusterEnd = clusterIndexes.at(i) + 4;
        if (clusterEnd >= scrollPos && clusterStart <= (int)end + scrollPos) {
            for (size_t j = 0; j < 5; j++) {
                sf::ConvexShape shape;
                shape.setPointCount(4);
                shape.setPoint(0, sf::Vector2f(0, 0));
                shape.setPoint(1, sf::Vector2f(10, 0));
                shape.setPoint(2, sf::Vector2f(10, 5));
                shape.setPoint(3, sf::Vector2f(0, 5));

                shape.setFillColor(sf::Color::Yellow);

                shape.setPosition(sf::Vector2f(startPosition.x + (float)(clusterIndexes.at(i) - scrollPos + j) * 12, startPosition.y));
                rw.draw(shape);
            }
        }
    }
}

void DNAStrand::highlightProteinClusters(sf::RenderWindow& rw, sf::Vector2f startPosition, int scrollPos, vector<int>& clusterIndexes) {
    int proteinScroll = scrollPos/3;
    size_t end = rw.getSize().x / 36;
    if (end > _sequence.length() - proteinScroll) {
        end = _sequence.length() - proteinScroll;
    }

    for (size_t i = 0; i < clusterIndexes.size(); i++) {
        int clusterStart = clusterIndexes.at(i);
        int clusterEnd = clusterIndexes.at(i) + 4;
        if (clusterEnd >= proteinScroll && clusterStart <= (int)end + proteinScroll) {
            for (size_t j = 0; j < 5; j++) {
                sf::ConvexShape shape;
                shape.setPointCount(4);
                shape.setPoint(0, sf::Vector2f(0, 0));
                shape.setPoint(1, sf::Vector2f(30, 0));
                shape.setPoint(2, sf::Vector2f(30, 5));
                shape.setPoint(3, sf::Vector2f(0, 5));

                shape.setFillColor(sf::Color::Yellow);

                // find offset for partial scroll
                float fractional = (float)(scrollPos % 3) * 12;
                float xOffset = startPosition.x - fractional;

                shape.setPosition(sf::Vector2f(xOffset + startPosition.x + (float)(clusterIndexes.at(i) - proteinScroll + j) * 36, startPosition.y));
                rw.draw(shape);
            }
        }
    }
}