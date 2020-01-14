rm -rf build
mkdir build
g++ -c -std=c++11 utils/distance_functions.cpp -o build/distance_functions.o
g++ -c -std=c++11 utils/readers.cpp -o build/readers.o
g++ -c -std=c++11 utils/utils.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/utils.o
g++ -c -std=c++11 utils/test.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/test.o
g++ -c -std=c++11 method_B.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/method_B.o
g++ -std=c++11 build/method_B.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o method_B
g++ -c -std=c++11 method_A.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/method_A.o
g++ -std=c++11 build/method_A.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o method_A
g++ -c -std=c++11 rough_method.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/rough_method.o
g++ -std=c++11 build/rough_method.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o rough_method
g++ -c -std=c++11 big_punishment.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/big_punishment.o
g++ -std=c++11 build/big_punishment.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o big_punishment
g++ -c -std=c++11 result_compare.cpp -o build/result_compare.o
g++ -std=c++11 build/result_compare.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o result_compare
g++ -std=c++11 kmeans_method/k_means_preprocess.cpp -Ispoa/include/ -Lspoa/build/lib -lspoa -o kmeans_method/k_means_preprocess
g++ -std=c++11 kmeans_method/k_means_postprocess.cpp -Ispoa/include/ -Lspoa/build/lib -lspoa -o kmeans_method/k_means_postprocess
g++ -std=c++11 dbscan_method/dbscan_preprocess.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o dbscan_method/dbscan_preprocess
g++ -std=c++11 dbscan_method/dbscan_postprocess.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o dbscan_method/dbscan_postprocess
