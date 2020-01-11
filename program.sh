g++ dbscan_method/dbscan_preprocess.cpp -Ispoa/include/ -Lspoa/build/lib -lspoa -o dbscan_method/dbscan_preprocess
./dbscan_method/dbscan_preprocess fastq/J29_B_CE_IonXpress_005.fastq dbscan_method/dbscan_data.txt
python dbscan_method/dbscan.py dbscan_method/dbscan_data.txt 5 5 dbscan_method/output.txt
g++ dbscan_method/dbscan_postprocess.cpp -Ispoa/include/ -Lspoa/build/lib -lspoa -o dbscan_method/dbscan_postprocess
./dbscan_method/dbscan_postprocess 0 5 -4 -8 -6 dbscan_method/output.txt