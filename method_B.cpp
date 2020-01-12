#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include "spoa/spoa.hpp"
#include "utils/readers.h"
#include "utils/distance_functions.h"
#include "utils/utils.h"

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

                int k_current = hamming_distance((*it3), current);

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

    std::vector<std::string> allSequences = readFastQFile(argv[6]);

    std::vector<std::string> sequences = find_sequences_with_most_common_length_plus_minus_n(allSequences, atoi(argv[7]));

    auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(atoi(argv[1])), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    auto graph = spoa::createGraph();

    for (const auto &it : sequences)
    {
        auto alignment = alignment_engine->align(it, graph);
        graph->add_alignment(alignment, it);
    }

    std::vector<std::string> msa;
    graph->generate_multiple_sequence_alignment(msa);

    //pronadi sekvence koje su udaljene jedne od drugih za barem 30 i postavi ih kao pocetne centroide
    std::map<int, std::string> centroids = init_clusters(msa, atoi(argv[8]));

    //u cluster stavi poravnatu sekvencu ako se nalazi unutar k od centroida
    std::map<int, std::vector<std::string>> clusters = create_clusters(msa, centroids, atoi(argv[9]));

    //izracunaj consensuse (nove centroide) iz clustera
    std::vector<std::string> new_centroids = find_allels(clusters, atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    // mapa u kojoj su kljucevi novi centroidi, a vrijednosti clusteri
    std::map<std::string, std::vector<std::string>> cluster_map;
    for (int i = 0; i < new_centroids.size(); ++i)
    {
        cluster_map[new_centroids[i]] = clusters.find(i)->second;
    }

    //grupiraj nove centroide ako su slicni
    std::map<int, std::vector<std::string>> merged_centroids = group_sequences(new_centroids, atoi(argv[10]));

    clusters = merge_clusters(merged_centroids, cluster_map);

    std::map<int, std::vector<std::string>> filtered_clusters = filter_clusters(clusters, atoi(argv[11]));
    fprintf(stderr, "Number of filtered clusters (%zu)\n", filtered_clusters.size());

    std::map<int, std::vector<std::string>> final_clusters = clean_clusters(filtered_clusters);

    std::vector<std::string> consensuses = find_allels(final_clusters, atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    std::ofstream outfile (argv[12]);
    for (const auto &it : consensuses)
    {
        fprintf(stderr, "Alel (%zu)\n", it.size());
        fprintf(stderr, "%s\n\n", it.c_str());
        outfile << "Alel size: " << it.size() << std::endl;
        outfile << it.c_str() << std::endl;
    }
}