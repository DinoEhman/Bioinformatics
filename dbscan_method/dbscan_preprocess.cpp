#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include "spoa/spoa.hpp"


// metoda pronalzi u potpunosti J29B-1_M13F-pUC i J29B-3_M13F-pUC, ali ne J29B-6_M13F-pUC za J29_B_CE_IonXpress_005
// metoda za jelenref01 napravi jednu zamjenu, jelenref02 u potpunosti, a kod jelenref04 jedno umetanje viska, 2 nepotrebna brisanja i 2 zamjene za J30_B_CE_IonXpress_006

/**
 * @author Dino Ehman i Luka Justic
 * Reads all sequences from fastq file and stores them in vector
 * 
 * @param file path
 * @return vector with string elements
 */
std::vector<std::string> readFastQFile(std::string file)
{
    std::ifstream fastq(file);
    std::vector<std::string> sequences;
    std::string line;

    int i = 0;
    if (fastq.is_open())
    {
        while (getline(fastq, line))
        {
            switch (i % 4)
            {
            case 1:
                sequences.push_back(line);
                break;

            default:
                break;
            }
            i++;
        }
    }
    return sequences;
}

/**
 * @author Lovre Budimir
 * Find sequences with most common length and all sequences with lenght in close range of most common length
 * 
 * @param sequences and range/2
 * @return filtered sequences
 */
std::vector<std::string> find_sequences_with_most_common_length_plus_minus_n(std::vector<std::string> allSequences, int n){

    std::map<int, std::vector<std::string>> commonLengthMap;

    for(std::vector<std::string>::iterator it = allSequences.begin(); it != allSequences.end(); ++it){
        int key = (*it).size();
        commonLengthMap[key].push_back(*it);
    }

    std::vector<std::string> result = commonLengthMap.begin()->second;
    int maxLentgh = commonLengthMap.begin()->first;
    int max = result.size();



    for(std::map<int, std::vector<std::string>>::iterator it = commonLengthMap.begin(); it != commonLengthMap.end(); it ++) {
        int currentSize = (it->second).size();
        if(currentSize > max){
            result = it->second;
            maxLentgh = it->first;
            max = currentSize;
        } else if(currentSize == max && it->first > maxLentgh){
            result = it->second;
            maxLentgh = it->first;
            max = currentSize;  
        }
    }

    std::vector<std::string> result_sequences;

    int first_index = maxLentgh - n;
    int last_index = maxLentgh + n;

    for(int i = first_index; i <= last_index; ++i){
        if(commonLengthMap.count(i) == 1) {
            result_sequences.reserve(result_sequences.size() + commonLengthMap[i].size());
            result_sequences.insert(result_sequences.end(), commonLengthMap[i].begin(), commonLengthMap[i].end());
        }
    }


    return result_sequences;

}


/**
 * @author Luka Justić
 * Read sequences from fastq file.
 * Filter sequences by length.
 * Write sequences to a file.
 */
int main(int argc, char **argv){

    std::vector<std::string> allSequences = readFastQFile(argv[1]);

    std::vector<std::string> sequences = find_sequences_with_most_common_length_plus_minus_n(allSequences, 5);

    std::ofstream outfile (argv[2]);

    for(std::vector<std::string>::iterator it = sequences.begin(); it != sequences.end(); ++it){
        std::string current = (*it);
        outfile << current << std::endl;
    }

    outfile.close();

}
