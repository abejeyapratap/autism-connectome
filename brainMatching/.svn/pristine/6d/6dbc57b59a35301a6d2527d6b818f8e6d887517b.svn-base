# -*- coding: utf-8 -*-
"""
This program takes two matrices and draws their heatmap side by side
sample: python heatMapDouble.py -i ./actualMatchingMatrixReduced.txt ./permutedMatchingMatrixReduced.txt --title actual permuted -o ./asd.png --doubleColorBars

"""
import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
import argparse

##### get command line parameters
parser = argparse.ArgumentParser(description='Plot heat map of given matrix')
parser.add_argument('-i','--inputPath', help='source path to the matrices to be plotted', required=True,type=str,nargs=2)
parser.add_argument('-o','--outputPath', help='destination path to save the png file', required=False,type=str)
parser.add_argument('--threshold', help='lower threshold for the cell values', required=False, default=0.00001, type=float)
parser.add_argument('--scale', help='scale matrix', required=False, default = 1.0, type=float)
parser.add_argument('--bcolor', help='background color', required=False, default='k', type=str)
parser.add_argument('--colorMap', help='color map', required=False, default='jet', type=str)
parser.add_argument('--display', help='display plot', required=False,action='store_true')
parser.add_argument('--doubleColorBars', help='draw double colorbars instead of a common color bar', required=False,action='store_true')
parser.add_argument('--title', help='titles for plots', required=False, type=str,nargs=2) #You can escape white space in the title from command line with "\ " as in <two\ words>

args = vars(parser.parse_args())
filePath1 = args['inputPath'][0]
filePath2 = args['inputPath'][1]
threshold=args['threshold']
displayPlotFlag=args['display']
doubleColorBarFlag=args['doubleColorBars']
backgroundColor=args['bcolor']
colorMap=args['colorMap']
scale=args['scale']

#if output path is not provided, then save the output image under the same folder 
#with first input file by merging the filenames of the input files
if parser.parse_args().outputPath:
    outputPath = args['outputPath'] 
else:
    folderPath = filePath1.rsplit("/",1)[0]
    filename1 = filePath1.rsplit("/",1)[1].rsplit(".",1)[0]
    filename2 = filePath2.rsplit("/",1)[1].rsplit(".",1)[0]
    outputPath = folderPath + "/" + filename1 + "_" + filename2 + ".png"

#get/set tiles
if parser.parse_args().title:
    title1=args['title'][0]
    title2=args['title'][1]
else:
    title1=filename1
    title2=filename2

#load matrices
print("read from " + filePath1 + " and " + filePath2)
matrix1 = np.loadtxt(filePath1)
matrix2 = np.loadtxt(filePath2)
matrix1 *= scale
matrix2 *= scale


#set color scale and the background color
cmap_=mpl.cm.get_cmap(colorMap) # nipy_spectral
cmap_.set_under(backgroundColor)  #'k' for black, 'w' for white, #c8b7b7 for pinkish


#initialize plots and set titles
fig, xarr = plt.subplots(1,2)
xarr[0].set_title(title1)
xarr[1].set_title(title2)

#adjustments for colorbars
if doubleColorBarFlag==True:
    vmax1=matrix1.max()
    vmax2=matrix2.max()
    fig.subplots_adjust(wspace = 0.3)
else:
    vmax1=max(matrix1.max(), matrix2.max())
    vmax2=vmax1
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.29, 0.03, 0.41])

im = xarr[0].imshow(matrix1, interpolation='none',vmin=threshold, vmax=vmax1, cmap=cmap_)
if doubleColorBarFlag==True:
    fig.colorbar(im, ax=xarr[0],fraction=0.046, pad=0.04) #individual colorbar

im = xarr[1].imshow(matrix2, interpolation='none',vmin=threshold, vmax=vmax2, cmap=cmap_)
if doubleColorBarFlag==True:
	fig.colorbar(im, ax=xarr[1],fraction=0.046, pad=0.04) #individual colorbar
    
#single colorbar for both heatmaps
if doubleColorBarFlag==False:
    fig.colorbar(im, cax=cbar_ax)

print("write to " + outputPath)
plt.savefig(outputPath, transparent=False, dpi=300, bbox_inches='tight')

#display the plot
if displayPlotFlag:
    plt.show()
