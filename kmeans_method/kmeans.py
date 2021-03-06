import sys

'''
@ Luka Justić
Numerize method turns a sequence of letters to sequence of numbers
'''
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

'''
@ Luka Justić
Denumerize method turns a sequence od numbers to a string of letters
'''
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

'''
@ Luka Justić
'''

# Open the file with aligned sequences
with open(sys.argv[1], 'r') as file:  
    lines = file.readlines()
    
# First line in document contains the number of classes
k = int(lines[0])
print("K:",k)

msa = []

# Turn aligned sequences to vectors of numbers
for i in range(1,len(lines)):
    line = lines[i]
    numerized = numerize(line)
    msa.append(numerized)
    
#print(msa)


from sklearn.cluster import KMeans
import numpy as np

X = np.array(msa)

# We use KMeans algorithm from sklearn library
kmeans = KMeans(n_clusters=k, random_state=0).fit(X)
labels = kmeans.labels_
#print(labels)

clusters = {}

# Initialize clusters
for label in labels:
    clusters[label] = []
    
# Add sequences to clusters using labels
for i in range(len(msa)):
    label = labels[i]
    m = msa[i]
    
    clusters[label].append(m)
    
#print(clusters)

# Write clusters to a file
with open(sys.argv[2], 'w') as file:  # Use file to refer to the file object
    for key in clusters:
        values = clusters[key]
        
        for value in values:
            de_value = denumerize(value)
            file.write(str(key)+"\n")
            file.write(de_value+"\n")
