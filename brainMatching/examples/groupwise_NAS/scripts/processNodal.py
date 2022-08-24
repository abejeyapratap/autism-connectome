# Script to identify ROIs that consistently mismatch more in patients vs healthy
# We calculate node mismatches 1 subject relative to ALL healthy controls

# define healthy matrix (220 x 150)
# start at control 0
    # for all healthy "matching pairs" of control 0 (should be 150-1), count # of mismatches of node i
    # normalize by 150-1

# define healthy matrix (220 x 163)
# start at patient 0
    # for all healthy "matching pairs" of patient 0 (should be 163), count # of mismatches of node i
    # normalize by 150

# for each of the 220 nodes, calculate group difference b/w healthy & patient mismatches

# report p-value & fdr corrected