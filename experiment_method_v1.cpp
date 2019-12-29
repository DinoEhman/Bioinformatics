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
 * Remove all clusters with less then s elements
 * 
 * @param clusters and min number of elements in cluster
 * @return clusters
 */ 
std::map<int, std::vector<std::string>> filter_clusters(std::map<int, std::vector<std::string>> clusters, int s) {


    for(std::map<int, std::vector<std::string>>::iterator it = clusters.begin(); it != clusters.end(); ){

        if((it->second).size() < s){
            clusters.erase(it++);
        } else {
            ++it;
        }

    }

    return clusters;

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

/**
 * @author Lovre Budimir
 * Iterate over all elements in multiple sequence alignment and if the distance between centroid and element is less then k
 * add element to the cluster
 * 
 * @param multiple sequence alignment, centorid and distance
 * @return clusters
 */
std::map<int, std::vector<std::string>> create_clusters(std::vector<std::string> msa, std::map<int, std::string> centroids, int k){

    std::map<int, std::vector<std::string>> clusters;

    for(std::map<int, std::string>::iterator it = centroids.begin(); it != centroids.end(); it++){

        int cluster_id = it->first;
        std::string centroid = it->second;

        for(std::vector<std::string>::iterator it = msa.begin() + 1; it != msa.end(); ++it){

            std::string current = (*it);

            if(calculate_distance(current, centroid) < k){
                clusters[cluster_id].push_back(current);
            }
        }
    }

    return clusters;
}

/**
 * @author Lovre Budimir
 * Remove all '-' from strings
 * 
 * @param clusters
 * @return clusters with elements without '-' in string
 */
std::map<int, std::vector<std::string>> clean_clusters(std::map<int, std::vector<std::string>> clusters){

    for(std::map<int, std::vector<std::string>>::iterator it = clusters.begin(); it != clusters.end(); it++){

        for(std::vector<std::string>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); it2++){

            (*it2).erase(std::remove((*it2).begin(), (*it2).end(), '-'), (*it2).end());
        
        }

    }

    return clusters;
}

int main(int argc, char **argv){

    std::vector<std::string> allSequences = readFastQFile("./fastq/J30_B_CE_IonXpress_006.fastq");

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

    std::map<int, std::vector<std::string>> clusters = create_clusters(msa, centroids, 15);

    std::map<int, std::vector<std::string>> filtered_clusters = filter_clusters(clusters, 10);
    fprintf(stderr, "Number of filtered clusters (%zu)\n", filtered_clusters.size());

    std::map<int, std::vector<std::string>> final_clusters = clean_clusters(filtered_clusters);

    
    
    for(std::map<int, std::vector<std::string>>::iterator it = final_clusters.begin(); it != final_clusters.end(); it++){

        std::vector<std::string> cluster_sequences = it->second;

        auto alignment_engine2 = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(atoi(argv[1])),
                                                        atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

        auto graph2 = spoa::createGraph();

        for (const auto &it : cluster_sequences)
        {
            auto alignment2 = alignment_engine2->align(it, graph2);
            graph2->add_alignment(alignment2, it);
        }

        std::string consensus2 = graph2->generate_consensus();
        fprintf(stderr, "Alel (%zu)\n", consensus2.size());
        fprintf(stderr, "%s\n\n", consensus2.c_str());

    }


}