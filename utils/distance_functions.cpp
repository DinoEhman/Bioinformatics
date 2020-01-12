#include "distance_functions.h"

/**
 * @author Luka Justic
 * Calculate hamming distance between two strings
 * 
 * @param two strings
 * @return distance 
 */
int hamming_distance(std::string s1, std::string s2){

    int i = 0, count = 0; 
    
    for(i = 0; i < s1.size(); ++i){
        if(s1[i] != s2[i]){
            count += 1;
        }
    }
    return count; 

}

/**
 * @author Lovre Budimir
 * Calculate same distance between two strings
 * 
 * @param two strings
 * @return distance 
 */
int match_distance(std::string s1, std::string s2){

    if(s1.length() != s2.length()){
         fprintf(stderr, "DIFFRENT STRING SIZES");
        return -1;
    }

    int i = 0, count = 0; 
    
    for(i = 0; i < s1.size(); ++i){
        if(s1[i] == s2[i]){
            count += 1;
        }
    }
    return count; 

}