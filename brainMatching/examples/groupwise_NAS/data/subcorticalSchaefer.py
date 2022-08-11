import os
import numpy as np

OLD = "./connectomes_schaefer"
# DIR = "./connectomes_new_schaefer"
DIR = "./connectomes2"

files = []
with open("schaefer_filtered.txt", "r") as f:
    files = f.read().splitlines()

# Filter 220x220 Schaefer connectomes to 200x200
for filename in os.listdir(OLD):
    connectome = []
    with open(f"{OLD}/{filename}", "r") as f:
        temp = f.read().splitlines()
        for i, line in enumerate(temp):
            connectome.append(line.split())

    with open(f"{DIR}/{filename}", "w") as f:
        for i in range(200):
            temp = connectome[i][:200] + connectome[i][218:]
            f.write(" ".join(temp) + "\n")
        
        for i in range(218, 220):
            temp = connectome[i][:200] + connectome[i][218:]
            f.write(" ".join(temp) + "\n")



""" # RUN ONLY ONCE - deletes extra connectome files
newFiles = []
for file in files:
    newFiles.append(f"{file}.txt")

print(len(os.listdir(OLD)))
count = 0
for filename in os.listdir(OLD):
    if filename not in newFiles:
        count += 1
        os.remove(f"{OLD}/{filename}")
print(count)
print(len(os.listdir(OLD))) """