def numerize(line):
    numerized = []
    
    for char in line:
        if char == 'A':
            numerized.append(1)
        elif char == 'T':
            numerized.append(4)
        elif char == 'C':
            numerized.append(2)
        elif char == 'G':
            numerized.append(3)
        elif char == '-':
            numerized.append(0)
            
    return numerized

def denumerize(line):
    denumerized = ""
    
    for char in line:
        if char == 1:
            denumerized += "A"
        elif char == 4:
            denumerized += "T"
        elif char ==  2:
            denumerized += "C"
        elif char == 3:
            denumerized += "G"
        elif char == 0:
            denumerized += "-"
            
    return denumerized

with open('kmeans_data.txt', 'r') as file:  # Use file to refer to the file object
    lines = file.readlines()
    

k = int(lines[0].split("=")[1].strip())
print(k)

msa = []

for i in range(1,len(lines)):
    line = lines[i]
    numerized = numerize(line)
    msa.append(numerized)
    
print(msa)

for m in msa:
    print(denumerize(m))

from sklearn.cluster import KMeans
import numpy as np
X = np.array(msa)
kmeans = KMeans(n_clusters=k, random_state=0).fit(X)
labels = kmeans.labels_
print(labels)

clusters = {}

for label in labels:
    clusters[label] = []
    
for i in range(len(msa)):
    label = labels[i]
    m = msa[i]
    
    clusters[label].append(m)
    
print(clusters)

with open('output.txt', 'x') as file:  # Use file to refer to the file object
    for key in clusters:
        values = clusters[key]
        
        for value in values:
            de_value = denumerize(value)
            file.write(str(key)+"\n")
            file.write(de_value+"\n")