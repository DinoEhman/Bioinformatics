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

// metoda pronalzi u potpunosti J29B-1_M13F-pUC i J29B-3_M13F-pUC, ali ne J29B-6_M13F-pUC za J29_B_CE_IonXpress_005
// metoda za jelenref01 napravi jednu zamjenu, jelenref02 u potpunosti, a kod jelenref04 jedno umetanje viska, 2 nepotrebna brisanja i 2 zamjene za J30_B_CE_IonXpress_006

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

    std::map<int, std::vector<std::string>> clusters = create_clusters(msa, centroids, atoi(argv[9]));

    std::map<int, std::vector<std::string>> filtered_clusters = filter_clusters(clusters, atoi(argv[10]));
    fprintf(stderr, "Number of filtered clusters (%zu)\n", filtered_clusters.size());

    std::map<int, std::vector<std::string>> final_clusters = clean_clusters(filtered_clusters);

    std::vector<std::string> consensuses = find_allels(final_clusters, atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    std::ofstream outfile(argv[11]);
    for (const auto &it : consensuses)
    {
        fprintf(stderr, "Alel (%zu)\n", it.size());
        fprintf(stderr, "%s\n\n", it.c_str());
        outfile << "Alel size: " << it.size() << std::endl;
        outfile << it.c_str() << std::endl;
    }
}