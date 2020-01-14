./kmeans_method/k_means_preprocess 0 5 -4 -8 -6 files/fastq/J30_B_CE_IonXpress_006.fastq kmeans_method/kmeans_input_30.txt

python kmeans_method/kmeans.py kmeans_method/kmeans_input_30.txt kmeans_method/kmeans_output_30.txt


./kmeans_method/k_means_postprocess 0 5 -4 -8 -6 kmeans_method/kmeans_output_30.txt kmeans_method/kmeans_results_30.txt