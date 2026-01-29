#include "dna_functions.h"

using namespace std;

char getPair(char nucleotide) {
    if (nucleotide == 'A') {
        return 'U';
    } else if (nucleotide == 'T') {
        return 'A';
    } else if (nucleotide == 'G') {
        return 'C';
    } else if (nucleotide == 'C') {
        return 'G';
    }
    return ' ';
}

size_t getMin(vector<double> vals, vector<int> indexes) {
    size_t min = 0;
    for (size_t i = 1; i < indexes.size(); i++) {
        if (vals.at(indexes.at(i)) < vals.at(indexes.at(min))) {
            min = i;
        }
    }
    return min;
}