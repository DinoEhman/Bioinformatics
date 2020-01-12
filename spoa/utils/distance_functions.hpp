/*!
 * @file distance_functions.hpp
 *
 * @brief functions for calculating distance between two strings
 */
#pragma once

#include <string>

//int min(int a, int b, int c);

int hamming_distance(std::string s1, std::string s2);

//int levenstein_distance(std::string s1, std::string s2);

int weighted_distance(std::string s1, std::string s2);

