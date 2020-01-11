from leven import levenshtein       
import numpy as np
from sklearn.cluster import dbscan
import sys

EPS = [3,4,5,6,7]

def load_data(path):
    with open(path, 'r') as file:  
        lines = file.readlines()

    sequences = []
    for i, line in enumerate(lines):
        sequences.append(line.strip())

    return sequences

#data = ["ACCTCCTAGAAG", "ACCTACTAGAAGTT", "GAATATTAGGCCGA"]

path_to_data = sys.argv[1]
data = load_data(path_to_data)
print(len(data))

def lev_metric(x, y):
     i, j = int(x[0]), int(y[0])     # extract indices
     return levenshtein(data[i], data[j])

X = np.arange(len(data)).reshape(-1, 1)

for eps in EPS:
    items, labels = dbscan(X, metric=lev_metric, eps=eps, min_samples=3, algorithm='brute')

    #print("Items: ", items)
    #print("Labels: ", labels)

    clusters = {}

    for index, label in enumerate(labels):
        if label == -1:
            continue
        
        if label not in clusters:
            clusters[label] = []

        print(index)
        clusters[label].append(data[index])

    #print(clusters)
    print("Done with eps: ",eps)

    # Write clusters to a file
    with open('output_30_eps_{}.txt'.format(eps), 'x') as file:  # Use file to refer to the file object
        for key in clusters:
            values = clusters[key]
            for value in values:
                file.write(str(key)+"\n")
                file.write(value+"\n")

