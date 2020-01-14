#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include "spoa/spoa.hpp"


// metoda pronalzi u potpunosti J29B-1_M13F-pUC i J29B-3_M13F-pUC, ali ne J29B-6_M13F-pUC za J29_B_CE_IonXpress_005
// metoda za jelenref01 napravi jednu zamjenu, jelenref02 u potpunosti, a kod jelenref04 jedno umetanje viska, 2 nepotrebna brisanja i 2 zamjene za J30_B_CE_IonXpress_006

/**
 * @author Dino Ehman
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
 * @author Luka Justic
 * Calculate hamming distance between two strings
 * 
 * @param two strings
 * @return distance 
 */
int calculate_distance(std::string s1, std::string s2){

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
 * Find elements from multiple sequence alignment with k distance and make them centroids for clusters
 * 
 * @param multiple sequence alignment and distance
 * @return centorids
 */
std::map<int, std::string> init_clusters(std::vector<std::string> msa, int k){

    std::map<int, std::string> centroids;

    int group_id_counter = 0;
    centroids[group_id_counter] = *msa.begin();

    for(std::vector<std::string>::iterator it = msa.begin() + 1; it != msa.end(); ++it){
        
        int flag = 1;
        std::string current = (*it);

        for(std::map<int, std::string>::iterator it2 = centroids.begin(); it2 != centroids.end(); it2++){

            std::string centroid = (it2->second);
            
            int distance = calculate_distance(current, centroid);

            if(distance < k) {
                flag = 0;
            }
        }

        if(flag == 1){
            group_id_counter += 1;
            centroids[group_id_counter] = current;
        }

    }

    return centroids;

}

int main(int argc, char **argv){

    std::vector<std::string> allSequences = readFastQFile(argv[6]);

    std::vector<std::string> sequences = find_sequences_with_most_common_length_plus_minus_n(allSequences, 5);

    auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(atoi(argv[1])),
                                                        atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    auto graph = spoa::createGraph();

    for (const auto &it : sequences)
    {
        auto alignment = alignment_engine->align(it, graph);
        graph->add_alignment(alignment, it);
    }

    std::vector<std::string> msa;
    graph->generate_multiple_sequence_alignment(msa);

    //pronadi sekvence koje su udaljene jedne od drugih za barem 30 i postavi ih kao pocetne centroide
    std::map<int, std::string> centroids = init_clusters(msa, 30);

    std::ofstream outfile (argv[7]);

    std::string number_of_clusters = std::to_string(centroids.size());
    outfile << number_of_clusters << std::endl;

    for(std::vector<std::string>::iterator it = msa.begin() + 1; it != msa.end(); ++it){
        std::string current = (*it);
        outfile << current << std::endl;
    }

    outfile.close();

}
