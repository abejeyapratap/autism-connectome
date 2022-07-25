# -*- coding: utf-8 -*-
"""
Given a matrix containing numerical values,
this program plots the heatmap of the matrix into a png file.

sample: python heatMap.py -i ./actualMatchingMatrixFull.txt -o ./heatmap.png --threshold 5 --color w --display
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
parser.add_argument('--interpolation', help='interpolation', required=False, default='none', type=str, choices=['none','nearest'])
parser.add_argument('--display', help='display plot', required=False,action='store_true')


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

#if output path is not provided, then save the output image under the same folder 
#with input image by only changing the file extension
if parser.parse_args().outputPath:
    outputPath = args['outputPath'] 
else:
	outputPath = inputPath.rsplit(".",1)[0] + ".png"

print("read from " + inputPath)
matrix = np.loadtxt(inputPath)
matrix *= scale

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
cmap_.set_over(backgroundColor)

fig, ax = plt.subplots(1, 1)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') #show labels of the xAxis on top

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

plt.colorbar()
print("save to " + outputPath)
plt.savefig(outputPath, transparent=False, dpi=300, bbox_inches='tight')

#display the plot
if displayPlotFlag:
    plt.show()
