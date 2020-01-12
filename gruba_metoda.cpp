#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include "spoa/spoa.hpp"

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
                line.erase(remove_if(line.begin(), line.end(), isspace), line.end());  //remove whitespace              
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
 * Reads all sequences from fasta file and stores them in vector
 * 
 * @param file path
 * @return vector with string elements
 */
std::vector<std::string> readFastaFile(std::string file)
{
    std::ifstream fastq(file);
    std::vector<std::string> alels;
    std::string line;

    int i = 0;
    if (fastq.is_open())
    {
        while (getline(fastq, line))
        {
            switch (i % 2)
            {
            case 1:
                line.erase(remove_if(line.begin(), line.end(), isspace), line.end());  //remove whitespace        
                alels.push_back(line);
                break;

            default:
                break;
            }
            i++;
        }
    }
    return alels;
}

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


std::map<int, std::vector<std::string>> group_sequences(std::vector<std::string> msa, int k) {

    int group_id_counter = 0;
    std::map<int, std::vector<std::string>> clusters;


    clusters[group_id_counter].push_back(*msa.begin());

    for(std::vector<std::string>::iterator it = msa.begin() + 1; it != msa.end(); ++it){

        int flag_found_cluster = 0;
        std::string current = (*it);
        
        for(std::map<int, std::vector<std::string>>::iterator it2 = clusters.begin(); it2 != clusters.end(); it2++){

            int max_k = 0;

            for(std::vector<std::string>::iterator it3 = (it2->second).begin(); it3 != (it2->second).end(); it3++){

                int k_current = hamming_distance(*it3, current);

                if(k_current > max_k) max_k = k_current;

            }

            if(max_k < k) {
                (it2->second).push_back(current);
                flag_found_cluster = 1;
                break;
            } 

        }

        if(flag_found_cluster == 0){
            group_id_counter += 1;
            clusters[group_id_counter].push_back(current);
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


        consensuses.push_back(consensus);

    }

    return consensuses;

}

/**
 * @author Lovre Budimir
 * find corresponding consensus for every expected allele if possible
 * 
 * @param consensuses
 */
void test_results(std::vector<std::string> expected, std::vector<std::string> predicted){

    printf("TESTING: find corresponding consensus for every expected allele if possible\n");

    int expected_id = 1;

    for (const auto &true_alel : expected){

        int predicted_id = 1;
        
        int true_length = true_alel.length();

        std::string best_predicted_msa;
        std::string best_expected_msa;
        int best_predicted_id = -1;
        double best_acc = 0;
        int best_matches = -1;


        for (const auto &pred_alel : predicted){

            // globalno poravnanje za testiranje
            auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(2), 4, -1, -2, -2);

            auto graph = spoa::createGraph();

            auto alignment1 = alignment_engine->align(true_alel, graph);
            graph->add_alignment(alignment1, true_alel);

            auto alignment2 = alignment_engine->align(pred_alel, graph);
            graph->add_alignment(alignment2, pred_alel);

            std::vector<std::string> msa;
            graph->generate_multiple_sequence_alignment(msa);

            int matches = match_distance(msa[0], msa[1]);
            double acc = (double)matches / (double) true_length;

            if(acc > best_acc){
                best_predicted_id = predicted_id;
                best_expected_msa = msa[0];
                best_predicted_msa = msa[1];
                best_acc = acc;
                best_matches = matches;
            }

            predicted_id++;
        }

        // nesto nije u redu, duljina expecteda i predicteda je ista, ali u consoli se ne ispise visak na kraju do kraja
        fprintf(stderr, "Expected(%d) and Predicted(%d)\n", expected_id, best_predicted_id);fflush(stdout);
        fprintf(stderr, "%s\n\n", best_expected_msa.c_str()); fflush(stderr);
        fprintf(stderr, "%s\n\n", best_predicted_msa.c_str()); fflush(stderr);

        double perc = best_acc * 100;
        fprintf(stderr, "Acc: %.2f%% (%d/%d)\n\n\n", perc, best_matches, true_length); fflush(stderr);        

        expected_id++;
    
    }

}


int main(int argc, char **argv)
{
    std::vector<std::string> allSequences = readFastQFile("./fastq/J29_B_CE_IonXpress_005.fastq");

    std::vector<std::string> sequences = find_sequences_with_most_common_length_plus_minus_n(allSequences, 5);

    fprintf(stderr, "Number of input sequences (%zu)\n", sequences.size());

    // 2. korak
    // parametri: AlignmentType (0 - lokalno, 1 - globalno, 2 - poluglobalno), m (podudaranje), n (zamjena), g (umetanje), e (brisanje)
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


    //3. korak
    std::map<int, std::vector<std::string>> clusters = group_sequences(msa, 15);
    fprintf(stderr, "Number of clusters (%zu)\n", clusters.size());

    for(std::map<int, std::vector<std::string>>::iterator it = clusters.begin(); it != clusters.end(); ++it){
        fprintf(stderr, "(%zu) ", (it->second).size());
    }

    clusters = clean_clusters(clusters);

    std::map<int, std::vector<std::string>> filtered_clusters = filter_clusters(clusters, 10);
    fprintf(stderr, "\nNumber of filtered clusters (%zu)\n", filtered_clusters.size());



    std::vector<std::string> consensuses = find_allels(filtered_clusters, atoi(argv[1]),atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    int test = 1;

    if(test){

        std::vector<std::string> expected = readFastaFile("./J29B_expected.fasta");

        test_results(expected, consensuses);

    }else{

        for (const auto &it : consensuses){
            fprintf(stderr, "Alel (%zu)\n", it.size());
            fprintf(stderr, "%s\n\n", it.c_str());
        }

    }

    return 0;
}