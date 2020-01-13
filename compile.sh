rm -rf build
mkdir build
g++ -c utils/distance_functions.cpp -o build/distance_functions.o
g++ -c utils/readers.cpp -o build/readers.o
g++ -c utils/utils.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/utils.o
g++ -c utils/test.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/test.o
g++ -c method_B.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/method_B.o
g++ build/method_B.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o method_B
g++ -c method_A.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/method_A.o
g++ build/method_A.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o method_A
g++ -c rough_method.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/rough_method.o
g++ build/rough_method.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o rough_method
g++ -c big_punishment.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/big_punishment.o
g++ build/big_punishment.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o big_punishment
g++ -c result_compare.cpp -o build/result_compare.o
g++ build/result_compare.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o result_compare