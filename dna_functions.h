#ifndef DNA_FUNCTIONS_H
#define DNA_FUNCTIONS_H

#include <vector>

/**
 * @brief Get the DNA nucleotide pair of inputted nucleotide
 * 
 * @return char 
 */
char getPair(char);

/**
 * @brief Get the index of the minimum of the list of indexes give a vector of double weights
 * 
 * @return size_t 
 */
size_t getMin(std::vector<double>, std::vector<int>);

#endif