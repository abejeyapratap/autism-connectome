import os

files = []
with open("all_subjects.txt", "r") as f:
    files = f.read().splitlines()

# DIR = './connectomes2'
# print(len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]))

# check which demographics subjects are not in Schaefer
with open("missing.txt", "w") as file:
    for f in files:
        if not os.path.isfile(f"connectomes2/{f}.txt"):
            print(f)
            file.write(f +"\n")

# intersection of Schaefer connectomes & tobacco demographics
filtered = []
for f in files:
    if os.path.isfile(f"connectomes2/{f}.txt"):
        filtered.append(f)

# print(len(files), len(filtered))

# store filtered connectome names to file
""" with open("subjects_filtered.txt", "w") as f:
    for val in filtered:
        f.write(val +"\n") """

# print(os.path.isfile(f"connectomes2/{files[1]}.txt"))
# print(all([os.path.isfile(f"connectomes2/{f}.txt") for f in files]))