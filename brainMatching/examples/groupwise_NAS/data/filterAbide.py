# Script to take intersection of subject names from csv & connectomes directory
    # Creates file of MISSING subject names

import os

files = []
with open("abide_all_subjects.txt", "r") as f:
    files = f.read().splitlines()

# print(files[:3])
DIR = './abide_nyu'
# print(len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]))
# print(len([name for name in os.listdir(DIR)]))

# check which demographics subjects are not in ABIDE NYU
count = 0
with open("abide_nyu_missing.txt", "w") as file:
    for f in files:
        if not os.path.isfile(f"{DIR}/{f}.1D"):
            count += 1
            # print(f"{DIR}/{f}.1D")
            file.write(f +"\n")

print("Missing:", count)
quit()

""" 
# check if 200x200 connectome (it's actually 176x200)
with open(f"{DIR}/{files[102]}.1D") as f:
    temp = f.read().splitlines()
    print(len(temp), len(temp[0].split()))

quit() """

# intersection of ABIDE connectomes & tobacco demographics
filtered = []
for f in files:
    if os.path.isfile(f"{DIR}/{f}.1D"):
        filtered.append(f)

# print(len(files), len(filtered))

# store filtered connectome names to file
with open("abide_nyu_filtered.txt", "w") as f:
    for val in filtered:
        f.write(val +"\n")

# print(os.path.isfile(f"connectomes_schaefer/{files[1]}.txt"))
# print(all([os.path.isfile(f"connectomes_schaefer/{f}.txt") for f in files]))