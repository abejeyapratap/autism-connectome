#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 16:18:27 2017

@author: yusuf
"""

import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
import scipy.stats as stt
import statsmodels.api as sm


def fdr(pvalues, correction_type='Benjamini-Hochberg'):
    """
    fdr(pvalues, correction_type='Benjamini-Hochberg'):
    
    Returns an adjusted pvalue array same shape as the input.
    
    http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python
    """
    import numpy as np
    #retain a copy of the data for reshaping to original format
    pvalues_orig=np.array(pvalues)
    pvalues_out=np.ones((np.size(pvalues_orig),))  #1d array
    
    #grab only those p-values that are not exactly equal to 1
    # (pvalues are 1 when the model was not fit and the element was skipped over)
    w=np.where(pvalues_orig.flatten()<1)[0]
    pvalues=pvalues_orig[pvalues_orig<1]  #flattens
    
    n=float(pvalues.shape[0])
    new_pvalues=np.empty(int(n))
    if correction_type=='Bonferroni':
        new_pvalues=n*pvalues
        new_pvalues[new_pvalues>1]=1
    elif correction_type=='Bonferroni-Holm':
        values=[(pvalue,i) for i,pvalue in enumerate(pvalues)]
        values.sort()
        for rank,vals in enumerate(values):
            new_pvalues[vals[1]]=(n-rank)*vals[0]
    elif correction_type=='Benjamini-Hochberg':
        values=[(pvalue,i) for i,pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values=[]
        for i,vals in enumerate(values):
            rank=n-i
            pvalue,index=vals
            new_values.append((n/rank)*pvalue)
        for i in range(0,int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1]=new_values[i]
        for i,vals in enumerate(values):
            new_pvalues[vals[1]]=new_values[i]
    
    #fill in values
    pvalues_out[w]=new_pvalues
    
    return pvalues_out.reshape(pvalues_orig.shape)

def savetxt_compact(fname, mat, fmt="%.6f", delimiter=' ', fileAccessMode='a'):
    with open(fname, fileAccessMode) as fh:
        for row in mat:
            line = delimiter.join("0" if value == 0 else fmt % value for value in row)
            fh.write(line + '\n')

def calculateZScore(data,controlGroupOrder):
    newData=np.zeros(data.shape)
    for i in range(len(data)):
        if(i<len(controlGroupOrder)):
            controls=np.delete(controlGroupOrder,i)
        meanControl = np.mean(data[controls])
        stdControl = np.std(data[controls])
        newData[i] = (data[i]-meanControl)/stdControl
    return newData


def heatMap(matrix,displayPlotFlag,outputPath,backgroundColor,threshold):
    cmap_=mpl.cm.get_cmap('jet')
    cmap_.set_under(backgroundColor)
    
    fig=plt.figure()
    plt.imshow(matrix, interpolation='none',vmin=threshold,cmap=cmap_)
    plt.colorbar()
    if(outputPath!=''):
        print("save to " + outputPath)
        plt.savefig(outputPath, transparent=True, dpi=300, bbox_inches='tight')
    
    #display the plot
    if displayPlotFlag:
        plt.show()

def calculateGroupDifference(data1,data2,parametric=True,paired=False):
    import numpy as np
    import scipy.stats as stt
          
    mean = np.array([data1.mean(),data2.mean()])
    var = np.array([data1.var(),data2.var()])
    std = np.array([data1.std(),data2.std()])
    size = np.array([len(data1),len(data2)])
    
    if(var[0]==0 and var[1]==0):
        return 0,1
    
    if(parametric==True):
        if(paired==True):#Student's t test for repeated measures for normal distributions
            # calculate parametric paired t-test, with the assumption of distribution being normal
            statistic, pValue=stt.ttest_rel(data1,data2)
            #Calculate effect sie as explained in: http://www.real-statistics.com/students-t-distribution/paired-sample-t-test/cohens-d-paired-samples/
            r=stt.pearsonr(data1,data2)[0]
            s_z=np.sqrt(std[0]**2+std[1]**2-2*r*std[0]*std[1])
            s_rm=s_z/np.sqrt(2*(1-r))
            effectSize = (mean[0]-mean[1])/s_rm ##effect size for repeated measures
        elif(paired==False):#student's t test for independent variables having normal distribution
            # do F-test to see if the variances of two distributions are different
            if(var[1]!=0):
                F=var[0]/var[1]
                pVal=1-stt.f.cdf(F,len(data1)-1,len(data2)-1)#if, pVal<0.05, two distributions have different variance
            else:### if variance of one of the data is zero while the other is not, then consider this as the two dataset has different variance
                pVal=0
            if(pVal>0.05):
                statistic, pValue = stt.ttest_ind(data1,data2,equal_var=True)
            else:
                statistic, pValue = stt.ttest_ind(data1,data2,equal_var=False)
            
            #calculate effectSize using Cohen's D for t-test : https://en.wikipedia.org/wiki/Effect_size#Cohen's_d
            s_pooled=np.sqrt(((size[0]-1)*var[0]+(size[1]-1)*var[1])/float(size[0]+size[1]-2))
            effectSize = (mean[0]-mean[1])/s_pooled 
    elif(parametric==False):
        if(paired==True): #Wilcoxon Signed-rank test for repeated measures with non normal distribution
            # calculate non parametric t-test, when the data is not normally distributed
            statistic, pValue=stt.wilcoxon(data1,data2)
            # calculate effect size:  http://yatani.jp/teaching/doku.php?id=hcistats:wilcoxonsigned
            effectSize = statistic/np.sqrt(size[0])
        elif(paired==False):#Mann-Whitney U test
            statistic, pValue = stt.mannwhitneyu(data1,data2)
            # effectSize for U-test: https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Effect_sizes
            effectSize = 1-(2*statistic/float(size[0]*size[1])) 
        
    return effectSize,pValue


def drawCorrelationPlot(data1,data2,r,p,data1Label,data2Label,plotTitle,outputPath,pointNames=None):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_xlabel(data1Label,fontsize=28)
    ax.set_ylabel(data2Label,fontsize=28)
    ax.tick_params(axis='both',labelsize=18)
    mpl.rcParams['font.sans-serif'] = "Times New Roman"
    mpl.rcParams['font.family'] = "serif"
    ax.text(0.77,0.9, 'r='+str("{0:.3f}".format(r))+'\np='+("{0:.3f}".format(p)),
            horizontalalignment='left',verticalalignment='center',
            transform = ax.transAxes, color='grey',fontsize=18,
            bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    coeff=np.polyfit(data1,data2,1)
    data2Regressed=np.polyval(coeff,data1)
    plt.scatter(data1,data2,c='orchid',edgecolors='none',alpha=0.7)
    if pointNames!=None:
        for i in range(len(pointNames)):
            ax.annotate(pointNames[i],(data1[i],data2[i]),size=6)
    plt.plot(data1,data2Regressed,color='darkorchid',linewidth=2)
    plt.title(plotTitle)
    plt.savefig(outputPath, transparent=False, dpi=300, bbox_inches='tight')
    return plt

def drawLongitudinalCorrelationPlot(data1,data2,r,p,data1Label,data2Label,plotTitle,outputPath,pointNames=None):
    numTimepoints=len(data1[0])
    numSubjects=len(data1)
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_xlabel(data1Label)
    ax.set_ylabel(data2Label)
    data1Flat=data1.flatten()
    data2Flat=data2.flatten()
    coeff=np.polyfit(data1Flat,data2Flat,1)
    data2FlatRegressed=np.polyval(coeff,data1Flat)

    mpl.rcParams['font.sans-serif'] = "Times New Roman"
    mpl.rcParams['font.family'] = "serif"

#    ax.text(0.05,0.9, 'r='+str("{0:.3f}".format(r))+'\np='+("{0:.3f}".format(p)),
#            horizontalalignment='left',verticalalignment='center',
#            transform = ax.transAxes, color='grey',
#            bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    colors=['#2096BA','#FF7E5F','#613A43','#6C4F70','#849974']
    for i in range(numTimepoints):
        plt.scatter(data1[:,i],data2[:,i],c=colors[i],edgecolors='none',alpha=0.7)
    
    for i in range(numSubjects):
        if pointNames!=None:
            for j in range(numTimepoints):
                ax.annotate(pointNames[i]+"_"+str(j),(data1[i,j],data2[i,j]),size=6)
        clr=(np.random.rand(),np.random.rand(),np.random.rand())
        plt.plot(data1[i,:],data2[i,:],color=clr)
#    plt.plot(data1Flat,data2FlatRegressed,color='darkorchid',linewidth=2)
    plt.title(plotTitle)
    plt.savefig(outputPath, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()
    return plt

def drawCorrelationPlotWithConfidenceInterval(data1,data2,r,p,data1Label,data2Label,plotTitle,outputPath,pointNames=None,logScale=False):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_xlabel(data1Label,fontsize=28)
    ax.set_ylabel(data2Label,fontsize=28)
    ax.tick_params(axis='both',labelsize=18)
    
    mpl.rcParams['font.sans-serif'] = "Times New Roman"
    mpl.rcParams['font.family'] = "serif"
    ax.text(0.84,0.87, 'r='+str("{0:.3f}".format(r))+'\np='+("{0:.3f}".format(p)),
            horizontalalignment='center',verticalalignment='center',fontsize=24,
            transform = ax.transAxes, color='grey',
            bbox={'facecolor':'white', 'alpha':0.7, 'pad':10})
    if(logScale):
        ax.set_yscale('log',basey=2)
    
    plt.scatter(data1,data2,c='orchid',alpha=0.7,edgecolors='none')
    if pointNames!=None:
        for i in range(len(pointNames)):
            ax.annotate(pointNames[i],(data1[i],data2[i]),size=6)
    
    ###OLS linear regression
    #http://markthegraph.blogspot.com/2015/05/using-python-statsmodels-for-ols-linear.html
    x=np.array(data1)
    y=np.array(data2)
    xMin=x.min()
    xMax=x.max()
    
    #regression line
    x = sm.add_constant(x)
    model = sm.OLS(y,x)
    fitted = model.fit()
    x_pred=np.linspace(xMin,xMax,50)
    x_pred2=sm.add_constant(x_pred)
    y_pred=fitted.predict(x_pred2)
    
    plt.plot(x_pred,y_pred,color='darkorchid',linewidth=2)
    
    #confidence interval for the regression
    y_hat = fitted.predict(x)
    y_err = y - y_hat
    mean_x = x.T[1].mean()
    n = len(x)
    dof = n - fitted.df_model - 1
    t = stt.t.ppf(1-0.025, df=dof)
    s_err = np.sum(np.power(y_err, 2))
    conf = t * np.sqrt((s_err/(n-2))*
                       (1.0/n + (np.power((x_pred-mean_x),2)/((np.sum(np.power(x_pred,2))) - n*(np.power(mean_x,2))))))
    upper = y_pred + abs(conf)
    lower = y_pred - abs(conf)
    #ax.fill_between(x_pred, lower, upper, color='#888888', alpha=0.4)
    
    #prediction interval
    from statsmodels.sandbox.regression.predstd import wls_prediction_std
    sdev, lower, upper = wls_prediction_std(fitted, exog=x_pred2, alpha=0.05)
    ax.fill_between(x_pred, lower, upper, color='#888888', alpha=0.1)
    
    plt.title(plotTitle)
    
    plt.savefig(outputPath, transparent=False, dpi=300, bbox_inches='tight')

def drawHistogram2Dataset(data1,data2,r,p,scoreLabel1, scoreLabel2,data1Label,data2Label,color1,color2,xLabel,yLabel,plotTitle,outputPath):
    fig=plt.figure()
    plt.hist(data1, bins=32, alpha=0.8,label=data1Label,color=color1,edgecolor='none')
    plt.hist(data2, bins=32, alpha=0.5,label=data2Label,color=color2,edgecolor='none')
#    line1 = plt.Line2D(range(1),range(1),linewidth=8,color=color1,alpha=1.0)
#    line2 = plt.Line2D(range(10),range(10),linewidth=8,color=color2,alpha=1.0)
#    plt.legend((line1,line2),(data1Label,data2Label),loc='upper right', fancybox=True, framealpha=0.5)
    mpl.rcParams['font.sans-serif'] = "Times New Roman"
    mpl.rcParams['font.family'] = "serif"
    plt.legend(loc='upper right', fancybox=True, framealpha=0.5)
    ax=fig.add_subplot(111)
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    ax.text(0.775,0.75, str(scoreLabel1)+' : '+str("{0:.3f}".format(r))+'\n'+ str(scoreLabel2)+' : '+("{0:.3f}".format(p)),
            horizontalalignment='left',verticalalignment='center',
            transform = ax.transAxes, color='black',
            bbox={'facecolor':'white', 'alpha':0.5, 'pad':8})
    plt.title(plotTitle)
    plt.savefig(outputPath, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

def drawHistogramSingleDataset(data,dataLabel,xLabel,yLabel,plotTitle,outputPath):
    mpl.rcParams['font.sans-serif'] = "Times New Roman"
    mpl.rcParams['font.family'] = "serif"
    fig=plt.figure()
    plt.hist(data, bins=32, alpha=0.5,label=dataLabel)
    plt.legend(loc='upper right', fancybox=True, framealpha=0.5)
    ax=fig.add_subplot(111)
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    plt.title(plotTitle)
    plt.savefig(outputPath, transparent=True, dpi=300, bbox_inches='tight')
    plt.close()
    
def drawBarChart(objects,values,xLabel,yLabel,plotTitle,outputPath):
    mpl.rcParams['font.sans-serif'] = "Times New Roman"
    mpl.rcParams['font.family'] = "serif"
    fig=plt.figure()
    y_pos = np.arange(len(objects))
    plt.bar(y_pos,values, align='center',alpha=0.5)
    plt.xticks(y_pos,objects)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.title(plotTitle)
    plt.savefig(outputPath, transparent=True, dpi=300, bbox_inches='tight')
    plt.close()
