#include "test.h"
#include "distance_functions.h"
#include "spoa/spoa.hpp"
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