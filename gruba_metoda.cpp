#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include "spoa/spoa.hpp"

struct seq_msa {
    std::string original;
    std::string msa;
};

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

    int first_index = maxLentgh - 5;
    int last_index = maxLentgh + 5;

    for(int i = first_index; i <= last_index; ++i){
        if(commonLengthMap.count(i) == 1) {
            result_sequences.reserve(result_sequences.size() + commonLengthMap[i].size());
            result_sequences.insert(result_sequences.end(), commonLengthMap[i].begin(), commonLengthMap[i].end());
        }
    }


    return result_sequences;

}

std::map<int, std::vector<seq_msa>> filter_clusters(std::map<int, std::vector<seq_msa>> clusters, int s) {


    for(std::map<int, std::vector<seq_msa>>::iterator it = clusters.begin(); it != clusters.end(); ){

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

std::map<int, std::vector<seq_msa>> group_sequences(std::vector<seq_msa> seq_msa_vector, int k) {

    int group_id_counter = 0;
    std::map<int, std::vector<seq_msa>> clusters;


    clusters[group_id_counter].push_back(*seq_msa_vector.begin());

    for(std::vector<seq_msa>::iterator it = seq_msa_vector.begin() + 1; it != seq_msa_vector.end(); ++it){

        int flag_found_cluster = 0;
        seq_msa current = (*it);
        
        for(std::map<int, std::vector<seq_msa>>::iterator it2 = clusters.begin(); it2 != clusters.end(); it2++){

            int max_k = 0;

            for(std::vector<seq_msa>::iterator it3 = (it2->second).begin(); it3 != (it2->second).end(); it3++){

                int k_current = calculate_distance(it3->msa, current.msa);

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
    std::vector<std::string> sequences = find_sequences_with_most_common_length_plus_minus_n(allSequences, 5); //1

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

    std::vector<seq_msa> seq_msa_vector;
    for(int i = 0; i < sequences.size(); ++i){
        seq_msa new_seq_msa;

        new_seq_msa.original = sequences[i];
        new_seq_msa.msa = msa[i];

        seq_msa_vector.push_back(new_seq_msa);
    }  

    //3. korak
    std::map<int, std::vector<seq_msa>> clusters = group_sequences(seq_msa_vector, 12);
    fprintf(stderr, "Number of clusters (%zu)\n", clusters.size());

    for(std::map<int, std::vector<seq_msa>>::iterator it = clusters.begin(); it != clusters.end(); ++it){
        fprintf(stderr, "(%zu) ", (*it).second.size());
    }

    std::map<int, std::vector<seq_msa>> filtered_clusters = filter_clusters(clusters, 100);
    fprintf(stderr, "\nNumber of filtered clusters (%zu)\n", filtered_clusters.size());


    //4. korak
    std::vector<std::string> result_consensuses;

    for(std::map<int, std::vector<seq_msa>>::iterator it = filtered_clusters.begin(); it != filtered_clusters.end(); it++){

        std::vector<seq_msa> cluster_sequences = it->second;

        auto alignment_engine2 = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(atoi(argv[1])),
                                                        atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

        auto graph2 = spoa::createGraph();

        for (const auto &it : cluster_sequences)
        {
            auto alignment2 = alignment_engine2->align(it.original, graph2);
            graph2->add_alignment(alignment2, it.original);
        }

        std::string consensus2 = graph2->generate_consensus();
        result_consensuses.push_back(consensus2);
        fprintf(stderr, "Alel (%zu)\n", consensus2.size());
        fprintf(stderr, "%s\n\n", consensus2.c_str());

    }



    //test
    std::vector<std::string> expected = readFastaFile("./J30B_expected.fasta");

    for(std::vector<std::string>::iterator it = result_consensuses.begin(); it != result_consensuses.end(); it++){
        std::string current_result = (*it);

        int rbr = 0;
        for(std::vector<std::string>::iterator it2 = expected.begin(); it2 != expected.end(); it2++){
            std::string current_expected = (*it2);

            auto alignment_engine2 = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(atoi(argv[1])),
                                                        atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

            auto graph2 = spoa::createGraph();

            auto alignment_result = alignment_engine2->align(current_result, graph2);
            graph2->add_alignment(alignment_result, current_result);

            auto alignment_expected = alignment_engine2->align(current_expected, graph2);
            graph2->add_alignment(alignment_expected, current_result);

            std::vector<std::string> msa2;
            graph2->generate_multiple_sequence_alignment(msa2);

            int dist = calculate_distance(msa[0], msa[1]);

            if(dist == 0){
                fprintf(stderr, "Result:\n%s\n", current_result.c_str());
                fprintf(stderr, "Rbr (%d)\n", rbr);
                fprintf(stderr, "Expected:\n%s\n", current_expected.c_str());
                fprintf(stderr, "Result msa:\n%s\n", msa[0].c_str());
                fprintf(stderr, "Expected msa:\n%s\n\n", msa[1].c_str());
            }

            rbr += 1;
        }
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