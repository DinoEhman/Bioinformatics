#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include "spoa/spoa.hpp"
#include "utils/readers.h"
#include "utils/distance_functions.h"
#include "utils/utils.h"
#include "utils/test.h"

std::map<int, std::vector<std::string>> group_sequences(std::vector<std::string> msa, int k)
{

    int group_id_counter = 0;
    std::map<int, std::vector<std::string>> clusters;

    clusters[group_id_counter].push_back(*msa.begin());

    for (std::vector<std::string>::iterator it = msa.begin() + 1; it != msa.end(); ++it)
    {

        int flag_found_cluster = 0;
        std::string current = (*it);

        for (std::map<int, std::vector<std::string>>::iterator it2 = clusters.begin(); it2 != clusters.end(); it2++)
        {

            int max_k = 0;

            for (std::vector<std::string>::iterator it3 = (it2->second).begin(); it3 != (it2->second).end(); it3++)
            {

                int k_current = hamming_distance(*it3, current);

                if (k_current > max_k)
                    max_k = k_current;
            }

            if (max_k < k)
            {
                (it2->second).push_back(current);
                flag_found_cluster = 1;
                break;
            }
        }

        if (flag_found_cluster == 0)
        {
            group_id_counter += 1;
            clusters[group_id_counter].push_back(current);
        }
    }

    return clusters;
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

    for (std::map<int, std::vector<std::string>>::iterator it = clusters.begin(); it != clusters.end(); ++it)
    {
        fprintf(stderr, "(%zu) ", (it->second).size());
    }

    clusters = clean_clusters(clusters);

    std::map<int, std::vector<std::string>> filtered_clusters = filter_clusters(clusters, 10);
    fprintf(stderr, "\nNumber of filtered clusters (%zu)\n", filtered_clusters.size());

    std::vector<std::string> consensuses = find_allels(filtered_clusters, atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    int test = 1;

    if (test)
    {

        std::vector<std::string> expected = readFastaFile("./J29B_expected.fasta");

        test_results(expected, consensuses);
    }
    else
    {

        for (const auto &it : consensuses)
        {
            fprintf(stderr, "Alel (%zu)\n", it.size());
            fprintf(stderr, "%s\n\n", it.c_str());
        }
    }

    return 0;
}