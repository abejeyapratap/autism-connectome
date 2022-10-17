import os

# Schaefer
# OLD = "./connectomes_schaefer/schaefer_200"
# DIR = "./connectomes_schaefer/norm_schaefer_200_2"

# Desikan
OLD = "./connectomes_desikan/desikan"
DIR = "./connectomes_desikan/norm_desikan_2"

CONNECTOME_LEN = 86

if not os.path.exists(DIR):
    os.mkdir(DIR)

for filename in os.listdir(OLD):
    # Load connectome into matrix
    connectome = []
    sysPath = f"{OLD}/{filename}"
    with open(sysPath, "r") as f:
        lines = f.read().splitlines()
        if len(lines) != CONNECTOME_LEN:
            print(sysPath)
            exit(1)
        for i, stringLine in enumerate(lines):
            line = stringLine.split()
            if len(line) != CONNECTOME_LEN:
                print(sysPath)
                exit(2)
            connectome.append([float(num) for num in stringLine.split()])

    # set diagonals to 0
    for i in range(CONNECTOME_LEN):
        connectome[i][i] = 0

    edgesSummation = sum(map(sum, connectome)) / 2
    # edgesSummation = edgesSummation / 50 # scale to avoid super small numbers

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