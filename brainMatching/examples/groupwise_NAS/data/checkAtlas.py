# Script to check if subject names are same across Schaefer vs Desikan
# also contains misc code for checking if Schaefer220 --> Schaefer200/202

import os

# Sanity check to see if Schaefer connectomes were filtered properly
""" for file in os.listdir("./connectomes_schaefer_200"):
    with open(f"./connectomes_new_schaefer/{file}") as f:
        temp = f.read().splitlines()
        if len(temp) != 200 or len(temp[0].split()) != 200:
            print(file)
            # print(len(temp), len(temp[0]))
            break """

for file in os.listdir("./connectomes_schaefer/norm_schaefer_200_2"):
    with open(f"./connectomes_schaefer/norm_schaefer_200_2/{file}") as f:
        temp = f.read().splitlines()
        if len(temp) != 200 or len(temp[0].split()) != 200:
            print(file)
            print(len(temp), len(temp[0].split()))
            break
exit()

sFiles = []
with open("schaefer_filtered.txt", "r") as f:
    sFiles = f.read().splitlines()

dFiles = []
with open("desikan_filtered.txt", "r") as f:
    dFiles = f.read().splitlines()

# print(len(sFiles), len(dFiles))

count = 0
# Check if Schaefer & Desikan subject names are same
for i, _ in enumerate(sFiles):
    s = "_".join(sFiles[i].split("_", 2)[:2])
    d = "_".join(dFiles[i].split("_", 2)[:2])
    if s != d:
        count += 1
# print(count)