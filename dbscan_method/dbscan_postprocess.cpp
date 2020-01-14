#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include "spoa/spoa.hpp"

/**
 * @author Luka Justić and Dino Ehman
 * Reads clustered sequences created by dbscan method and stores them in clusters (set of strings)
 * 
 * @param file path
 * @return map with vectors of strings elements
 */
std::map<int, std::vector<std::string>> readKmeansOutput(std::string file)
{

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

/**
 * @author Luka Justić
 * Read output of kmean method.
 * For every cluster generate concenzus and write it to a fasta file.
 */
int main(int argc, char **argv)
{

    std::map<int, std::vector<std::string>> clusters = readKmeansOutput(argv[6]);

    std::cout << clusters.size() << std::endl;
    std::vector<std::string> consensuses;

    for (std::map<int, std::vector<std::string>>::iterator it = clusters.begin(); it != clusters.end(); it++)
    {

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
        consensuses.push_back(consensus2);
    }
    
    std::ofstream outfile(argv[7]);
    for (const auto &it : consensuses)
    {
        fprintf(stderr, "Alel (%zu)\n", it.size());
        fprintf(stderr, "%s\n\n", it.c_str());
        outfile << "Alel size: " << it.size() << std::endl;
        outfile << it.c_str() << std::endl;
    }
}