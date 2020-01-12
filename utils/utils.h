#include<vector>
#include<string>
std::vector<std::string> find_sequences_with_most_common_length_plus_minus_n(std::vector<std::string> allSequences, int n);
std::map<int, std::vector<std::string>> clean_clusters(std::map<int, std::vector<std::string>> clusters);
std::map<int, std::vector<std::string>> filter_clusters(std::map<int, std::vector<std::string>> clusters, int s);
std::map<int, std::string> init_clusters(std::vector<std::string> msa, int k);
std::map<int, std::vector<std::string>> create_clusters(std::vector<std::string> msa, std::map<int, std::string> centroids, int k);
std::vector<std::string> find_allels(std::map<int, std::vector<std::string>> clusters,  int alg, int m, int n, int g, int e);
std::map<int, std::vector<std::string>> merge_clusters(std::map<int, std::vector<std::string>> merged_centroids, std::map<std::string, std::vector<std::string>> cluster_map);