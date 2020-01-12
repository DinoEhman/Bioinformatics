#include "utils.h"
#include <algorithm>
#include "distance_functions.h"
#include "spoa/spoa.hpp"


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
            
            int distance = hamming_distance(current, centroid);

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

            if(hamming_distance(current, centroid) < k){
                clusters[cluster_id].push_back(current);
            }
        }
    }

    return clusters;
}

std::vector<std::string> find_allels(std::map<int, std::vector<std::string>> clusters,  int alg, int m, int n, int g, int e){

    std::vector<std::string> consensuses;

    for(std::map<int, std::vector<std::string>>::iterator it = clusters.begin(); it != clusters.end(); it++){

        std::vector<std::string> cluster_sequences = it->second;

        // parametri: AlignmentType (0 - lokalno, 1 - globalno, 2 - poluglobalno), m (podudaranje), n (zamjena), g (umetanje), e (brisanje)
        auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(alg), m, n, g, e);

        auto graph = spoa::createGraph();

        for (const auto &it : cluster_sequences)
        {
            auto alignment = alignment_engine->align(it, graph);
            graph->add_alignment(alignment, it);
        }

        std::string consensus = graph->generate_consensus();


        // fprintf(stderr, "Alel (%zu)\n", consensus.size());
        // fprintf(stderr, "%s\n\n", consensus.c_str());

        consensuses.push_back(consensus);

    }

    return consensuses;


}

std::map<int, std::vector<std::string>> merge_clusters(std::map<int, std::vector<std::string>> merged_centroids, std::map<std::string, std::vector<std::string>> cluster_map) {

    std::map<int, std::vector<std::string>> merged_clusters;

    for(int i = 0; i < merged_centroids.size(); ++i){

        std::vector<std::string> centroids = merged_centroids.find(i)->second;

        for(int j = 0; j < centroids.size(); ++j){

            std::vector<std::string> cluster = cluster_map.find(centroids[j])->second;

            for(int k = 0; k < cluster.size(); ++k){

                merged_clusters[i].push_back(cluster[k]);

            }

        }

    }

    return merged_clusters;
}
