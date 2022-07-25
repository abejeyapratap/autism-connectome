# -*- coding: utf-8 -*-
"""
This program takes a matrix that keeps the maping numbers of structural
regions to function regions and gerenates a heat map of the matrix with
region names written in the axis
sample: python heatMapStrFunc.py -i ./significantMappingsReduced.txt -o ./asd.png --display
"""
import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
import argparse

##### get command line parameters
parser = argparse.ArgumentParser(description='Plot heat map of given matrix')
parser.add_argument('-i','--inputPath', help='source path to the matrix to be plotted', required=True, type=str)
parser.add_argument('-o','--outputPath', help='destination path to save the png file', required=False, type=str)
parser.add_argument('--permutation', help='Reoreder the matrix', required=False, type=str)
parser.add_argument('--log', help='log scale', required=False,action='store_true',default=False)
parser.add_argument('--normalize', help='normalize the data to fit in [0,1] range', required=False,action='store_true',default=False)
parser.add_argument('--dataSign', help='Matrix consists of positive/mixed/negative values', required=False,type=str,choices=['pos','mixed','neg'],default='pos')
parser.add_argument('--threshold', help='lower/upper threshold for the cell values', required=False, type=float)
parser.add_argument('--scale', help='scale matrix', required=False, default = 1.0, type=float)
parser.add_argument('--bcolor', help='background color', required=False, default='k', type=str)
parser.add_argument('--colorMap', help='color map', required=False, default='jet', type=str)
parser.add_argument('--interpolation', help='interpolation', required=False, default='nearest', type=str, choices=['none','nearest'])
parser.add_argument('--display', help='display plot', required=False,action='store_true')
parser.add_argument('--roiNames', help='write ROI names in axes', required=False,action='store_true',default=False)


args = vars(parser.parse_args())
inputPath=args['inputPath']
threshold=args['threshold']
backgroundColor=args['bcolor']
colorMap=args['colorMap']
displayPlotFlag=args['display']
logScaleFlag=args['log']
normalizeFlag=args['normalize']
dataSign=args['dataSign']
interpolation=args['interpolation']
scale=args['scale']
roiNamesFlag=args['roiNames']

#if output path is not provided, then save the output image under the same folder 
#with input image by only changing the file extension
if parser.parse_args().outputPath:
    outputPath = args['outputPath'] 
else:
	outputPath = inputPath.rsplit(".",1)[0] + ".png"

print("read from " + inputPath)
matrix = np.loadtxt(inputPath)
matrix *= scale

numRoi=len(matrix)

if normalizeFlag==True:
    matrix /= matrix.max()

if logScaleFlag==True:
    matrix1 = np.log(matrix)
    matrix1[matrix<=0] = 0
    matrix = matrix1
    
if parser.parse_args().permutation:
    mapping = np.loadtxt(args['permutation'],dtype=int)
    permutation = np.array(mapping[:,1]-1)
    matrix = matrix[:,permutation][permutation,:]

cmap_=mpl.cm.get_cmap(colorMap) # nipy_spectral
cmap_.set_under(backgroundColor)
cmap_.set_over(backgroundColor) #'k' for black, 'w' for white, #c8b7b7 for pinkish

#roiNamesYeo=['visual','somatomotor','dorsal','ventral','limbic','frontoparietal','default','subcortical']
roiNamesYeo=['visual','somato.','dorsal','ventral','limbic','fronto','DMN','subcor.'][:numRoi]
roiNamesLausanne=['lateralorbitofrontal-1','lateralorbitofrontal-2','lateralorbitofrontal-3','lateralorbitofrontal-4','parsorbitalis','frontalpole','medialorbitofrontal-1','medialorbitofrontal-2','medialorbitofrontal-3','parstriangularis-1','parstriangularis-2','parsopercularis-1','parsopercularis-2','rostralmiddlefrontal-1','rostralmiddlefrontal-2','rostralmiddlefrontal-3','rostralmiddlefrontal-4','rostralmiddlefrontal-5','rostralmiddlefrontal-6','superiorfrontal-1','superiorfrontal-2','superiorfrontal-3','superiorfrontal-4','superiorfrontal-5','superiorfrontal-6','superiorfrontal-7','superiorfrontal-8','caudalmiddlefrontal-1','caudalmiddlefrontal-2','caudalmiddlefrontal-3','precentral-1','precentral-2','precentral-3','precentral-4','precentral-5','precentral-6','paracentral-1','paracentral-2','paracentral-3','rostralanteriorcingulate-1','caudalanteriorcingulate-1','posteriorcingulate-1','posteriorcingulate-2','isthmuscingulate','postcentral-1','postcentral-2','postcentral-3','postcentral-4','postcentral-5','supramarginal-1','supramarginal-2','supramarginal-3','supramarginal-4','superiorparietal-1','superiorparietal-2','superiorparietal-3','superiorparietal-4','superiorparietal-5','superiorparietal-6','superiorparietal-7','inferiorparietal-1','inferiorparietal-2','inferiorparietal-3','inferiorparietal-4','inferiorparietal-5','inferiorparietal-6','precuneus-1','precuneus-2','precuneus-3','precuneus-4','precuneus-5','cuneus-1','cuneus-2','pericalcarine-1','pericalcarine-2','lateraloccipital-1','lateraloccipital-2','lateraloccipital-3','lateraloccipital-4','lateraloccipital-5','lingual-1','lingual-2','lingual-3','fusiform-1','fusiform-2','fusiform-3','fusiform-4','parahippocampal','entorhinal','temporalpole','inferiortemporal-1','inferiortemporal-2','inferiortemporal-3','inferiortemporal-4','middletemporal-1','middletemporal-2','middletemporal-3','middletemporal-4','bankssts','superiortemporal-1','superiortemporal-2','superiortemporal-3','superiortemporal-4','superiortemporal-5','transversetemporal','insula-1','insula-2','insula-3','thalamusproper','caudate','putamen','pallidum','accumbensarea','hippocampus','amygdala','lateralorbitofrontal-1','lateralorbitofrontal-2','lateralorbitofrontal-3','lateralorbitofrontal-4','parsorbitalis','frontalpole','medialorbitofrontal-1','medialorbitofrontal-2','parstriangularis-1','parsopercularis-1','parsopercularis-2','rostralmiddlefrontal-1','rostralmiddlefrontal-2','rostralmiddlefrontal-3','rostralmiddlefrontal-4','rostralmiddlefrontal-5','rostralmiddlefrontal-6','superiorfrontal-1','superiorfrontal-2','superiorfrontal-3','superiorfrontal-4','superiorfrontal-5','superiorfrontal-6','superiorfrontal-7','superiorfrontal-8','superiorfrontal-9','caudalmiddlefrontal-1','caudalmiddlefrontal-2','caudalmiddlefrontal-3','precentral-1','precentral-2','precentral-3','precentral-4','precentral-5','precentral-6','precentral-7','precentral-8','paracentral-1','paracentral-2','rostralanteriorcingulate','caudalanteriorcingulate','posteriorcingulate-1','posteriorcingulate-2','isthmuscingulate','postcentral-1','postcentral-2','postcentral-3','postcentral-4','postcentral-5','postcentral-6','postcentral-7','supramarginal-1','supramarginal-2','supramarginal-3','supramarginal-4','supramarginal-5','superiorparietal-1','superiorparietal-2','superiorparietal-3','superiorparietal-4','superiorparietal-5','superiorparietal-6','superiorparietal-7','inferiorparietal-1','inferiorparietal-2','inferiorparietal-3','inferiorparietal-4','inferiorparietal-5','precuneus-1','precuneus-2','precuneus-3','precuneus-4','precuneus-5','cuneus','pericalcarine','lateraloccipital-1','lateraloccipital-2','lateraloccipital-3','lateraloccipital-4','lateraloccipital-5','lingual-1','lingual-2','lingual-3','lingual-4','fusiform-1','fusiform-2','fusiform-3','fusiform-4','parahippocampal','entorhinal','temporalpole','inferiortemporal-1','inferiortemporal-2','inferiortemporal-3','inferiortemporal-4','middletemporal-1','middletemporal-2','middletemporal-3','middletemporal-4','bankssts-1','bankssts-2','superiortemporal-1','superiortemporal-2','superiortemporal-3','superiortemporal-4','superiortemporal-5','transversetemporal','insula-1','insula-2','insula-3','insula-4','thalamusproper','caudate','putamen','pallidum','accumbensarea','hippocampus','amygdala','brainstem']
roiNamesDesikan=['lateralorbitofrontal-left','parsorbitalis-left','frontalpole-left','medialorbitofrontal-left','parstriangularis-left','parsopercularis-left','rostralmiddlefrontal-left','superiorfrontal-left','caudalmiddlefrontal-left','precentral-left','paracentral-left','rostralanteriorcingulate-left','posteriorcingulate-left','isthmuscingulate-left','postcentral-left','supramarginal-left','superiorparietal-left','inferiorparietal-left','precuneus-left','cuneus-left','pericalcarine-left','lateraloccipital-left','lingual-left','fusiform-left','parahippocampal-left','entorhinal-left','temporalpole-left','inferiortemporal-left','middletemporal-left','bankssts-left','superiortemporal-left','transversetemporal-left','insula-left','thalamusproper-left','caudate-left','putamen-left','pallidum-left','accumbensarea-left','hippocampus-left','amygdala-left','lateralorbitofrontal-right','parsorbitalis-right','frontalpole-right','medialorbitofrontal-right','parstriangularis-right','Parsopercularis-right','rostralmiddlefrontal-right','superiorfrontal-right','caudalmiddlefrontal-right','precentral-right','paracentral-right','rostralanteriorcingulate-right','posteriorcingulate-right','isthmuscingulate-right','postcentral-right','supramarginal-right','superiorparietal-right','inferiorparietal-right','precuneus-right','cuneus-right','pericalcarine-right','lateraloccipital-right','lingual-right','fusiform-right','parahippocampal-right','entorhinal-right','temporalpole-right','inferiortemporal-right','middletemporal-right','bankssts-right','superiortemporal-right','transversetemporal-right','insula-right','thalamusproper-right','caudate-right','putamen-right','pallidum-right','accumbensarea-right','hippocampus-right','amygdala-right','brainstem']

size=matrix.shape[0]
if size==len(roiNamesYeo):
	roiNames=roiNamesYeo
elif size==len(roiNamesLausanne):
	roiNames=roiNamesLausanne
elif size==len(roiNamesDesikan):
	roiNames=roiNamesDesikan

fig, ax = plt.subplots(1, 1)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') #show labels of the xAxis on top
ax.tick_params(axis='both',labelsize=14)

if dataSign=='mixed':
    plt.imshow(matrix, interpolation=interpolation,cmap=cmap_)    
elif dataSign=='pos':
    if parser.parse_args().threshold==False:
        plt.imshow(matrix, interpolation=interpolation,vmin=matrix.min(),cmap=cmap_)    
    else:
        plt.imshow(matrix, interpolation=interpolation,vmin=threshold,cmap=cmap_)    
elif dataSign=='neg':
    if parser.parse_args().threshold==False:
        plt.imshow(matrix, interpolation=interpolation,vmax=matrix.max(),cmap=cmap_)    
    else:
        plt.imshow(matrix, interpolation=interpolation,vmax=threshold,cmap=cmap_)

if roiNamesFlag==True:
    plt.yticks(range(size), roiNames, size=16, rotation=45)
    plt.xticks(range(size), roiNames, size=16, rotation=45)# 'vertical' or 'horizontal' are other options
cbar=plt.colorbar(spacing='proportional',orientation='vertical', shrink=1.0, format="%.0f")
cbar.ax.tick_params(labelsize=16)

print("save to " + outputPath)
plt.savefig(outputPath, transparent=False, dpi=300, bbox_inches='tight')

#display the plot
if displayPlotFlag:
    plt.show()
