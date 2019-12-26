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

int main(int argc, char **argv)
{
    std::vector<std::string> allSequences = readFastQFile("./fastq/J29_B_CE_IonXpress_005.fastq");

    std::vector<std::string> sequences = find_sequences_with_most_common_length(allSequences);

    fprintf(stderr, "Number of input sequences (%zu)\n", sequences.size());

    fprintf(stderr, "Most common size (%zu)\n", (*sequences.begin()).size());

    // parametri: AlignmentType (0 - lokalno, 1 - globalno, 2 - poluglobalno), m (podudaranje), n (zamjena), g (umetanje), e (brisanje)
    auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(atoi(argv[1])),
                                                        atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    auto graph = spoa::createGraph();

    for (const auto &it : sequences)
    {
        auto alignment = alignment_engine->align(it, graph);
        graph->add_alignment(alignment, it);
    }

    std::string consensus = graph->generate_consensus();

    fprintf(stderr, "Consensus (%zu)\n", consensus.size());
    fprintf(stderr, "%s\n", consensus.c_str());

    std::vector<std::string> msa;
    graph->generate_multiple_sequence_alignment(msa);

    int counter = 1;
    fprintf(stderr, "Multiple sequence alignment\n");
    for (const auto& it: msa) {
        fprintf(stderr, "%d %s\n",counter, it.c_str());
        counter++;
    }

    return 0;
}