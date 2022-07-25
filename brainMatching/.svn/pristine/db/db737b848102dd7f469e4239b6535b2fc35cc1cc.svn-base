Sample scripts to run structure-function matching on the cluster.

Note that, in "Osmanlıoğlu, Y., Tunç, B., Parker, D., Elliott, M.A., Baum, G.L., Ciric, R., Satterthwaite, T.D., Gur, R.E., Gur, R.C. and Verma, R., 2019. System-level matching of structural and functional connectomes in the human brain. NeuroImage, 199, pp.93-104.", we used the hardest version of matching algorithm: edgesIncludeDiagRandDiag

This version of the problem randomly generates numbers on the diagonals of the connectomes and then calculates matching. The purpose in here was to remove the bias that would help identify the correct matching node due to the only zero entry in funcitonal connectomes being in the diagonals. In order to balance this bias after removing the diagonals and randomly setting the values, this experiment is needed to be repeated several times (1000 is what we did in NIMG paper) and then the matching results needs to be averaged for the rest of the evaluation of results. Thus, this is a job that is needed to be run on the cluster.

Additionally, we ran a permutation testing for the experiment in NIMG paper (which might not have  been necessary in the first place, but we did). This test requires shuffling of the connectomes and then doing the matching again. We repeated this experiemnt 1000 times as well, thus requires to be run on the cluster.

I'm adding the scripts and codes that I used to evaluate the results of these experiment, just in case someone would like to run structure-funciton coupling one day.
