# BrainMatch: A Program for Graph Matching of Brain Networks

**BrainMatch** is a program, that has two major uses: calculating brain network similarity of a subject relative to a group, and structure-function coupling of a subject. For the details of how to generate the executable and potential use cases, see information below.

## Folder Structure
- `groupwise_NAS`: contains experiments and scripts for groupwise graph matching of brain networks
- `subjectwise_strFunc`: contains prototype experiments and scripts for structure-function coupling of brain networks

## Setup
1. Generating the executable brain matching program:
	1. Extract the DLIB library files: `cd lib`
  	2. Extract the tar ball: `tar -zxvf dlib-19.4.tar.gz`

2. Compiling the brainMatch code to generate the executable
	1. Compile the brain matching code: `make -f makefile`
      	- This will generate Build and bin folders. Generated executable will be located under `bin/`
3. If you are making changes to the code, you might want to remove previous builds: `make -f makefile clean`

*Note*: It is absolutely necessary to have the `brainMatch` executable. Otherwise, all other scripts will fail.


### Information Archive
3) Example codes
----------------
Although the brain matching program does many things, two of its major uses are exemplified under ./examples/ folder
In order to generate results all at once, type:
	$> ./runAll.sh
under each experiment folder.

a) SUBJECTWISE: structure-function matching at connectome and system level
This use case of the code matches two connectomes belonging to the same subject (such as matching structure and function of an individual) across the dataset. Thus, this is a subjectwise matching use of the program.
This use is previously presented in "Osmanlıoğlu, Y., Tunç, B., Parker, D., Elliott, M.A., Baum, G.L., Ciric, R., Satterthwaite, T.D., Gur, R.E., Gur, R.C. and Verma, R., 2019. System-level matching of structural and functional connectomes in the human brain. NeuroImage, 199, pp.93-104."
 - In the example folder, I copied 10 structural and functional connectomes from the PNC dataset.
 - I use edgesIncludeDiagRandDiag as the assignment cost in this project.
 - I also did permutation testing in this project to calculate significance of matchings (although was not necessary)
 - Example code calculates 10 iterations of matching, 10 permutations testinggs, then calculates average matching matrices oput of these iterations
 - finally, we generate some heatmaps showing significant matchings at connectome and system level.

Please cite the aforementioned paper for this use case of the code

b) GROUPWISE: Network anomaly score (NAS)
This use case of the code matches connectome of a single subject with everybody else and in its evaluation, we consider average similarity of the individual relative to a group (such as healthycontrols). Thus, this use case of the code is groupwise matching.

In this experiment I use TBI data of Hoon. We use structural connectomes, that are generated using 500 seeds per voxel, probabilistic tractography. Details of the dataset and of all experiments can be found in "Osmanlioglu, Y., Parker, D., Alappatt, J.A., Gugger, J.J., Diaz-Arrastia, R.R., Whyte, J., Kim, J.J. and Verma, R., 2021. Connectomic Assessment of Injury Burden and Longitudinal Structural Network Alterations in Moderate-to-severe Traumatic Brain Injury. bioRxiv."

 - I first calculate matching scores among every two subject in the sample
 - I then calculate NAS score of patients relative to the healthy.
 - Then I plot box plots of patients at 3/6/12 months and healthy controls.
 - Finally I do mixed linear effects model (Using R) to calculate and plot change of NAS with days since injury. Note that, if this stage fails, check the error messages under the results folder to see which packages are missing in your R installation.

Please cite the aforementioned paper for this use case of the code


3) Netbeans project
nbproject folder and "Makefile" (with capital M) can be used to work on the project using netbeans.

4) the code uses DLIB library, which is provided with the code under ./lib/ folder. Compiler will take care of incorporating this library into the main code.


