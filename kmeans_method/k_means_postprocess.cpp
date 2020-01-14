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
    
int main(int argc, char **argv){

    std::map<int, std::vector<std::string>> clusters = readKmeansOutput(argv[6]);

    std::map<int, std::vector<std::string>> final_clusters = clean_clusters(clusters);
    
    std::vector<std::string> consensuses; 

    for(std::map<int, std::vector<std::string>>::iterator it = final_clusters.begin(); it != final_clusters.end(); it++){

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
        //fprintf(stderr, "Alel (%zu)\n", consensus2.size());
        //fprintf(stderr, "%s\n\n", consensus2.c_str());
	consensuses.push_back(consensus2);
    }

    std::ofstream outfile (argv[7]);
    for (const auto &it : consensuses)
    {
        fprintf(stderr, "Alel (%zu)\n", it.size());
        fprintf(stderr, "%s\n\n", it.c_str());
        outfile << "Alel size: " << it.size() << std::endl;
        outfile << it.c_str() << std::endl;
    }
}
