# Bioinformatics
Clustering sequenced genes

### Methods used:
1. Big punishment method
2. Rough method
3. Method A
4. Method B
5. Kmeans method
6. DBSCAN method

### Results on test set:

##### 1. Big punishment method
- J29 - finds 2 expected alleles completely but not third one
- J30 - finds first allele completely, second one with 1 change but doesn't find third expected allele

##### 2. Rough method
- J29 - finds 2 expected alleles completely but not third one
- J30 - finds first allele completely, second one with 1 change and third expected allele with 1 change, 2 deletion and 1 insertion

##### 3. Method A
- J29 - finds 2 expected alleles completely but not third one
- J30 - finds first allele completely, second one with 1 change and third expected allele with 1 change, 2 deletion and 1 insertion

##### 4. Method B
- J29 - finds first and third alleles completely and second one with 1 insertion
- J30 - finds first and second allele completely, and third expected allele with 1 change, 2 deletion and 1 insertion

##### 5. Kmeans method
- J29 - finds first allele with 1 insertion, second one with 1 insertion but doesn't find third expected allele
- J30 - finds first allele completely, second one with 1 change but doesn't find third expected allele

##### 6. DBSCAN method
- J29 - finds first allele completely, second one with 1 change and third allele with 1 deletion
- J30 - finds first allele completely, second one with 1 change and third allele with 1 insertion, 1 change and 2 deletion

### Compiling methods

After you clone the repository you will need to install spoa from https://github.com/rvaser/spoa in Bioinformatics directory.

All compilations are put together in script __compile.sh__ which does all the work or you can run it manually one by one.

``g++ -c utils/distance_functions.cpp -o build/distance_functions.o``

``g++ -c utils/readers.cpp -o build/readers.o``

``g++ -c utils/utils.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/utils.o``

``g++ -c utils/test.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/test.o``

``g++ -c method_B.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/method_B.o``

``g++ build/method_B.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o method_B``

``g++ -c method_A.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/method_A.o``

``g++ build/method_A.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o method_A``

``g++ -c rough_method.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/rough_method.o``

``g++ build/rough_method.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o rough_method``

``g++ -c big_punishment.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o build/big_punishment.o``

``g++ build/big_punishment.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o big_punishment``

``g++ -c result_compare.cpp -o build/result_compare.o``

``g++ build/result_compare.o build/distance_functions.o build/readers.o build/test.o build/utils.o -Ispoa/include -Lspoa/build/lib -lspoa -o result_compare``

``g++ kmeans_method/k_means_preprocess.cpp -Ispoa/include/ -Lspoa/build/lib -lspoa -o kmeans_method/k_means_preprocess``

``g++ kmeans_method/k_means_postprocess.cpp -Ispoa/include/ -Lspoa/build/lib -lspoa -o kmeans_method/k_means_postprocess``

``g++ dbscan_method/dbscan_preprocess.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o dbscan_method/dbscan_preprocess``

``g++ dbscan_method/dbscan_postprocess.cpp -Ispoa/include -Lspoa/build/lib -lspoa -o dbscan_method/dbscan_postprocess``


### Running programs

##### 1. Big punishment method

``./big_punishment 1 0 -1 -100 -100 files/fastq/J29_B_CE_IonXpress_005.fastq 10 5 output.txt``
``./big_punishment 1 0 -1 -100 -100 files/fastq/J30_B_CE_IonXpress_006.fastq 10 10 output.txt``

##### 2. Rough method

``./rough_method 1 0 -1 -1 -1 files/fastq/J29_B_CE_IonXpress_005.fastq 5 15 80 output.txt ``

``./rough_method 1 0 -1 -1 -1 files/fastq/J30_B_CE_IonXpress_006.fastq 5 15 10 output.txt`` 

##### 3. Method A

``./method_A 1 0 -1 -1 -1 files/fastq/J29_B_CE_IonXpress_005.fastq 5 30 20 6 output.txt``

``./method_A 1 0 -1 -1 -1 files/fastq/J30_B_CE_IonXpress_006.fastq 5 30 20 6 output.txt``

##### 4. Method B

``./method_B 2 1 -1 -1 -1 files/fastq/J29_B_CE_IonXpress_005.fastq 5 30 29 17 6 output.txt``

``./method_B 1 0 -1 -1 -1 files/fastq/J30_B_CE_IonXpress_006.fastq 5 30 10 16 6 output.txt``

##### 5. Kmeans method

``./kmeans_method/k_means_preprocess 1 0 -1 -1 -1 files/fastq/J30_B_CE_IonXpress_006.fastq kmeans_method/kmeans_input_30.txt``

``python kmeans_method/kmeans.py kmeans_method/kmeans_input_30.txt kmeans_method/kmeans_output_30.txt``

``./kmeans_method/k_means_postprocess 1 0 -1 -1 -1 kmeans_method/kmeans_output_30.txt kmeans_method/kmeans_results_30.txt``


#####  6. DBSCAN method

``./dbscan_method/dbscan_preprocess fastq/J29_B_CE_IonXpress_005.fastq dbscan_method/dbscan_data.txt``

``python dbscan_method/dbscan.py dbscan_method/dbscan_data.txt 5 5 dbscan_method/output.txt``

``./dbscan_method/dbscan_postprocess 0 1 -1 -1 -1 dbscan_method/output.txt``
