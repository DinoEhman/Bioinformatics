from leven import levenshtein       
import numpy as np
from sklearn.cluster import dbscan
import sys
import time


def load_data(path):
    with open(path, 'r') as file:  
        lines = file.readlines()

    sequences = []
    for i, line in enumerate(lines):
        sequences.append(line.strip())

    return sequences


path_to_data = "dbscan_data_j30.txt"
data = load_data(path_to_data)

print("dbscan started with",len(data),"alels")

start_time = time.time()
def lev_metric(x, y):
     i, j = int(x[0]), int(y[0])     # extract indices
     return levenshtein(data[i], data[j])

X = np.arange(len(data)).reshape(-1, 1)

EPS = [3,4,5]
MIN = [8,9,10,11,12]

for eps in EPS:
    for min_sample in MIN:
        items, labels = dbscan(X, metric=lev_metric, eps=eps, min_samples=min_sample, algorithm='brute')

        #print("Items: ", items)
        #print("Labels: ", labels)

        clusters = {}

        for index, label in enumerate(labels):
            if label == -1:
                continue
            
            if label not in clusters:
                clusters[label] = []

            #print(index)
            clusters[label].append(data[index])

        # Write clusters to a file
        with open("outputs_{}_eps_{}_min_{}.txt".format(30,eps, min_sample), 'w') as file:  # Use file to refer to the file object
            for key in clusters:
                values = clusters[key]
                for value in values:
                    file.write(str(key)+"\n")
                    file.write(value+"\n")

        print("Done with: ",eps," - ",min_sample)
        print("--- %s seconds ---" % (time.time() - start_time))


