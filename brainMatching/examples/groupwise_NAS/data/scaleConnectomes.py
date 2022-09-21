import os

CONNECTOME_LEN = 220

OLD = "./connectomes_schaefer"
DIR = "./connectomes_norm_schaefer"

os.mkdir(DIR)

for filename in os.listdir(OLD):
    # Load connectome into matrix
    connectome = []
    sysPath = f"{OLD}/{filename}"
    with open(sysPath, "r") as f:
        lines = f.read().splitlines()
        if len(lines) != CONNECTOME_LEN:
            print(sysPath)
        for i, stringLine in enumerate(lines):
            line = stringLine.split()
            if len(line) != CONNECTOME_LEN:
                print(sysPath)
            connectome.append([float(num) for num in stringLine.split()])

    # set diagonals to 0
    for i in range(CONNECTOME_LEN):
        connectome[i][i] = 0

    edgesSummation = sum(map(sum, connectome)) / 2

    # Normalize connectome by summation
    for i in range(CONNECTOME_LEN):
        for j in range(CONNECTOME_LEN):
            connectome[i][j] = connectome[i][j] / edgesSummation
    
    # print(max(map(max, connectome)))

    # Save normalized connectome for use in brainMatch
    with open(f"{DIR}/{filename}", "w") as f:
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