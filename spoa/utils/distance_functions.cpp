/*!
 * @file distance_functions.cpp
 *
 * @brief functions for calculating distance between two strings
 */

#include <string>

//int min(int a, int b, int c);

/**
 * @author Luka Justic
 * Calculate hamming distance between two strings
 * 
 * @param two strings
 * @cond same length
 * @return distance 
 */
int hamming_distance(std::string s1, std::string s2){

    if(s1.length() != s2.length()){
        return -1;
    }
    
    int i = 0, count = 0; 
    
    for(i = 0; i < s1.size(); ++i){
        if(s1[i] != s2[i]){
            count += 1;
        }
    }
    return count; 

}

//int levenstein_distance(std::string s1, std::string s2);

/**
 * @author Lovre Budimir
 * Calculate weighted distance between two strings
 * 
 * @param two strings
 * @cond same length
 * @return distance 
 */
float weighted_distance(std::string s1, std::string s2){

    if(s1.length() != s2.length()){
        return -1;
    }

    int i = 0;
    float count = 0; 
    
    for(i = 0; i < s1.size(); ++i){

        if(s1[i] == s2[i] and s1[i]=='-') {
            count += 0.75;
        }

        if((s1[i] != s2[i] and s1[i]=='-') || (s1[i] != s2[i] and s2[i]=='-')) {
            count += 0.75;
        }

        if(s1[i] != s2[i]){
            count += 1;
        }
    }
    return count; 

}

