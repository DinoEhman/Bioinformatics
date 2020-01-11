#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include "spoa/spoa.hpp"

std::map<int, std::vector<std::string>> readKmeansOutput(std::string file){

    std::ifstream output(file);
    std::map<int, std::vector<std::string>> clusters;
    std::string line;

    int i = 0;
    int cluster_id;
    if (output.is_open())
    {
        while (getline(output, line))
        {
            switch (i % 2)
            {
            case 0:
                cluster_id = atoi(line.c_str());
                break;
            case 1:
                clusters[cluster_id].push_back(line);
                break;

            default:
                break;
            }
            i++;
        }
    }
    return clusters;

}
    
int main(int argc, char **argv){

    std::map<int, std::vector<std::string>> clusters = readKmeansOutput(argv[1]);

    std::cout << clusters.size() << std::endl;
    
    for(std::map<int, std::vector<std::string>>::iterator it = clusters.begin(); it != clusters.end(); it++){

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