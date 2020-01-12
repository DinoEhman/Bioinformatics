#include<iostream>
#include "utils/readers.h"
#include "utils/test.h"
int main(int argc, char **argv)
{
    std::vector<std::string> expected = readFastaFile(argv[1]);
    std::vector<std::string> predicted = readFastaFile(argv[2]);

    test_results(expected,predicted);
}