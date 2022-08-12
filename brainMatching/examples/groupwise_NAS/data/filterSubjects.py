# Script to take intersection of subject names from csv & connectomes directory
    # Creates file of MISSING subject names
# Creates Schaefer & Desikan FILTERED subject names

import os

files = []
with open("all_subjects.txt", "r") as f:
    files = f.read().splitlines()

DIR = './connectomes_schaefer'
# print(len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]))
# print(len([name for name in os.listdir(DIR)]))

# check which demographics subjects are not in Schaefer
count = 0
with open("missing.txt", "w") as file:
    for f in files:
        if not os.path.isfile(f"connectomes_schaefer/{f}.txt"):
            count += 1
            print(f)
            file.write(f +"\n")

print("Missing:", count)

# intersection of Schaefer connectomes & tobacco demographics
filtered = []
for f in files:
    if os.path.isfile(f"connectomes_schaefer/{f}.txt"):
        filtered.append(f)

# print(len(files), len(filtered))

# store filtered connectome names to file
with open("schaefer_filtered.txt", "w") as f:
    for val in filtered:
        f.write(val +"\n")

# create & store Desikan subject names
desikan = []
for i, _ in enumerate(filtered):
    temp = "_".join(filtered[i].split("_", 2)[:2])
    temp = f"{temp}_desikan_connmat_sift2"
    desikan.append(temp)

DIR2 = './connectomes_desikan'
count = 0
for f in desikan:
    if not os.path.isfile(f"connectomes_desikan/{f}.txt"):
        count += 1
        print(f)
print(count)

with open("desikan_filtered.txt", "w") as f:
    for val in desikan:
        f.write(val +"\n")

# print(os.path.isfile(f"connectomes_schaefer/{files[1]}.txt"))
# print(all([os.path.isfile(f"connectomes_schaefer/{f}.txt") for f in files]))