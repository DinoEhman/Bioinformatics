from leven import levenshtein       
import numpy as np
from sklearn.cluster import dbscan
import sys

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
print("dbscan started with",len(data),"alels")

def lev_metric(x, y):
     i, j = int(x[0]), int(y[0])     # extract indices
     return levenshtein(data[i], data[j])

X = np.arange(len(data)).reshape(-1, 1)

items, labels = dbscan(X, metric=lev_metric, eps=int(sys.argv[2]), min_samples=int(sys.argv[3]), algorithm='brute')

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

# Write clusters to a file
with open(sys.argv[4], 'w') as file:  # Use file to refer to the file object
    for key in clusters:
        values = clusters[key]
        for value in values:
            file.write(str(key)+"\n")
            file.write(value+"\n")

