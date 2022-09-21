import os

CONNECTOME_LEN = 220

# Load connectome into matrix
connectome = []
with open(f"R0006_V0021_DTI_Schaefer2018_200_7Networks_connmat_sift2.txt", "r") as f:
    lines = f.read().splitlines()
    if len(lines) != CONNECTOME_LEN:
        print("fileName")
    for i, stringLine in enumerate(lines):
        line = stringLine.split()
        if len(line) != CONNECTOME_LEN:
            print("noo")
        connectome.append([float(num) for num in stringLine.split()])

# set diagonals to 0
for i in range(CONNECTOME_LEN):
    connectome[i][i] = 0

edgesSummation = sum(map(sum, connectome)) / 20

# Normalize connectome by summation
for i in range(CONNECTOME_LEN):
    for j in range(CONNECTOME_LEN):
        connectome[i][j] = connectome[i][j] / edgesSummation
print(max(map(max, connectome)))
quit()
# Save normalized connectome for use in brainMatch
with open("duplicate.txt", "w") as f:
    for i in range(CONNECTOME_LEN):
        f.write(" ".join(map(str, connectome[i])) + "\n")


""" count = 0
for i in range(220):
    for j in range(220):
        if i == j:
            continue
        if connectome[i][j] != connectome[j][i]:
            count += 1
print(count) """