import os

sFiles = []
with open("schaefer_filtered.txt", "r") as f:
    sFiles = f.read().splitlines()

dFiles = []
with open("desikan_filtered.txt", "r") as f:
    dFiles = f.read().splitlines()

print(len(sFiles), len(dFiles))

count = 0
for i, _ in enumerate(sFiles):
    s = "_".join(sFiles[i].split("_", 2)[:2])
    d = "_".join(dFiles[i].split("_", 2)[:2])
    if s != d:
        count += 1
print(count)