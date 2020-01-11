from leven import levenshtein       
import numpy as np
from sklearn.cluster import dbscan

def load_data(path):
    with open(path, 'r') as file:  
        lines = file.readlines()

    sequences = []
    for i, line in enumerate(lines):
        if i % 4 == 1: 
            sequences.append(line.strip())

    return sequences

#data = ["ACCTCCTAGAAG", "ACCTACTAGAAGTT", "GAATATTAGGCCGA"]

data = load_data('J30_B_CE_IonXpress_006.fastq')
data = data[0:1000]
print(len(data))

def lev_metric(x, y):
     i, j = int(x[0]), int(y[0])     # extract indices
     return levenshtein(data[i], data[j])

X = np.arange(len(data)).reshape(-1, 1)

items, labels = dbscan(X, metric=lev_metric, eps=5, min_samples=50, algorithm='brute')

#print("Items: ", items)
print("Labels: ", labels)

clusters = {}

for index, label in enumerate(labels):
    if label == -1:
        continue
    
    if label not in clusters:
        clusters[label] = []

    print(index)
    clusters[label].append(data[index])

print(clusters)

