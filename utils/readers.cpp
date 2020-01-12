#include"readers.h"
#include<fstream>
#include<algorithm>

/**
 * @author Dino Ehman
 * Reads all sequences from fastq file and stores them in vector
 * 
 * @param file path
 * @return vector with string elements
 */
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
                line.erase(remove_if(line.begin(), line.end(), isspace), line.end());  //remove whitespace              
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

/**
 * @author Lovre Budimir
 * Reads all sequences from fast file and stores them in vector
 * 
 * @param file path
 * @return vector with string elements
 */
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
                line.erase(remove_if(line.begin(), line.end(), isspace), line.end());  //remove whitespace        
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