#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include "spoa/spoa.hpp"
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

std::vector<std::string> find_sequences_with_most_common_length(std::vector<std::string> allSequences){

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

    // fprintf(stderr, "NUmber of 296 sequences (%zu)\n", commonLengthMap[296].size());
    // fprintf(stderr, "NUmber of 295 sequences (%zu)\n", commonLengthMap[295].size());
    // fprintf(stderr, "NUmber of 294 sequences (%zu)\n", commonLengthMap[294].size());
    // fprintf(stderr, "NUmber of 293 sequences (%zu)\n", commonLengthMap[293].size());
    // fprintf(stderr, "NUmber of 292 sequences (%zu)\n", commonLengthMap[292].size());
    // fprintf(stderr, "NUmber of 291 sequences (%zu)\n", commonLengthMap[291].size());
    // fprintf(stderr, "NUmber of 297 sequences (%zu)\n", commonLengthMap[297].size());
    // fprintf(stderr, "NUmber of 298 sequences (%zu)\n", commonLengthMap[298].size());
    // fprintf(stderr, "NUmber of 299 sequences (%zu)\n", commonLengthMap[299].size());
    // fprintf(stderr, "NUmber of 300 sequences (%zu)\n", commonLengthMap[300].size());
    // fprintf(stderr, "NUmber of 301 sequences (%zu)\n", commonLengthMap[301].size());

    return result;

}

std::map<int, std::vector<std::string>> filter_clusters(std::map<int, std::vector<std::string>> clusters, int s) {

    //std::map<int, std::vector<std::string>> filtered;

    for(std::map<int, std::vector<std::string>>::iterator it = clusters.begin(); it != clusters.end(); ){

        if((it->second).size() < s){
            clusters.erase(it++);
        } else {
            ++it;
        }

    }

    return clusters;

}   


int calculate_distance(std::string s1, std::string s2){

    int i = 0, count = 0; 
    
    for(i = 0; i < s1.size(); ++i){
        if(s1[i] != s2[i]){
            count += 1;
        }
    }
    return count; 

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

                int k_current = calculate_distance((*it3), current);

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

int main(int argc, char **argv)
{
    std::vector<std::string> allSequences = readFastQFile("./fastq/J30_B_CE_IonXpress_006.fastq");

    // 1. korak
    std::vector<std::string> sequences = find_sequences_with_most_common_length(allSequences); //1

    fprintf(stderr, "Number of input sequences (%zu)\n", sequences.size());

    fprintf(stderr, "Most common size (%zu)\n", (*sequences.begin()).size());

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
    std::map<int, std::vector<std::string>> clusters = group_sequences(msa, 8);
    fprintf(stderr, "Number of clusters (%zu)\n", clusters.size());

    std::map<int, std::vector<std::string>> filtered_clusters = filter_clusters(clusters, 10);
    fprintf(stderr, "Number of filtered clusters (%zu)\n", filtered_clusters.size());

    //4. korak
    for(std::map<int, std::vector<std::string>>::iterator it = filtered_clusters.begin(); it != filtered_clusters.end(); it++){

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




    //std::string consensus = graph->generate_consensus();
    //fprintf(stderr, "Consensus (%zu)\n", consensus.size());
    //fprintf(stderr, "%s\n\n", consensus.c_str());




    //int counter = 1;
    //fprintf(stderr, "Multiple sequence alignment\n");
    //for (const auto& it: msa) {
    //    fprintf(stderr, "%d %s\n",counter, it.c_str());
    //    counter++;
    //}

    return 0;
}