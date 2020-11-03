# -*- coding: utf-8 -*-
"""
Created on Sat May  2 17:10:58 2020

@author: zyll1
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import Lasso
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
import statsmodels.api as sm
from scipy.stats import f
import scipy.stats as stat
import math
from sklearn import metrics
from sklearn.metrics import auc 
import scipy as sp
import os

## returns the list minimum index
def indexofMin(arr):
	minindex = 0
	currentindex = 1
	while currentindex < len(arr):
		if arr[currentindex] < arr[minindex]:
			minindex = currentindex
		currentindex += 1
	return minindex

# sorting 
def sor_g(gene_1,gene_2):
    com = []
    for j in range(len(gene_1)):
        tup = (gene_1[j],gene_2[j])
        com.append(tup)
    sort_1 = sorted(com, key=lambda x:x[0])#Sort by gene_1
    sort_2 = sorted(com, key=lambda x:x[1])#Sort by gene_2
    a1 = []
    b1 = []
    for k in range(len(gene_1)):
        a1.append(sort_1[k][0])
        b1.append(sort_1[k][1])   
    return a1,b1
#Lagged regression
def granger(gene_1,gene_2):
	com = []
	for j in range(len(samples1)):
		tup = (gene_1[j],gene_2[j])
		com.append(tup)
	sort_1 = sorted(com, key=lambda x:x[0])
	# granger test based on gene sorting, 
	a1 = []
	b1 = []
	for k in range(len(samples1)):
		a1.append(sort_1[k][0])
		b1.append(sort_1[k][1])
	a1_cha = []
	b1_cha = []
	for k in range(len(samples1)):
		a1_cha.append(300*a1[k]/90001)
		#a1_cha.append(3*a1[k]/90001)
		b1_cha.append(300*b1[k]/90001)
		#b1_cha.append(3*b1[k]/90001)
	# Lag phase 1
	a1_t_1 = np.array(a1_cha[0:len(samples1)-1])
	a1_t1_1= np.array(a1_cha[1:len(samples1)])
	b1_t_1 = np.array(b1_cha[0:len(samples1)-1])
	b1_t1_1 = np.array(b1_cha[1:len(samples1)])   
	b1_t1_t2_1= sm.add_constant(b1_t1_1)
	a1_t1_t2_b1_1 = np.vstack((a1_t1_1,b1_t1_1)).T
	a1_t1_t2_b1_1 = sm.add_constant(a1_t1_t2_b1_1)
	 #F-test
	model = LinearRegression()
	model.fit(a1_t1_t2_b1_1,b1_t_1)
	b1_rssu_1 = np.sum((model.predict(a1_t1_t2_b1_1)-b1_t_1)**2)
	model.fit(b1_t1_t2_1,b1_t_1)
	b1_rssr_1 = np.sum((model.predict(b1_t1_t2_1)-b1_t_1)**2)
	f1_1 = ((b1_rssr_1-b1_rssu_1)/1)/(b1_rssu_1/(len(samples1)-3))
	p1_1 = f.sf(f1_1,1,len(samples1)-3) 
	 #Lag phase 2
	a1_t = np.array(a1_cha[0:len(samples1)-2])
	a1_t1 = np.array(a1_cha[1:len(samples1)-1])
	a1_t2 = np.array(a1_cha[2:len(samples1)])
	b1_t = np.array(b1_cha[0:len(samples1)-2])
	b1_t1 = np.array(b1_cha[1:len(samples1)-1])
	b1_t2 = np.array(b1_cha[2:len(samples1)])		
	a1_t1_t2 =  np.vstack((a1_t1,a1_t2)).T
	a1_t1_t2 = sm.add_constant(a1_t1_t2)		
	b1_t1_t2 = np.vstack((b1_t1,b1_t2)).T
	b1_t1_t2 = sm.add_constant(b1_t1_t2)		
	a1_t1_t2_b1 = np.vstack((a1_t1,a1_t2,b1_t1,b1_t2)).T
	a1_t1_t2_b1 = sm.add_constant(a1_t1_t2_b1)
	#F-test
	model.fit(a1_t1_t2_b1,b1_t)
	b1_rssu_2 = np.sum((model.predict(a1_t1_t2_b1)-b1_t)**2)
	model.fit(b1_t1_t2,b1_t)
	b1_rssr_2 = np.sum((model.predict(b1_t1_t2)-b1_t)**2) 
	f1_2 = ((b1_rssr_2-b1_rssu_2)/2)/(b1_rssu_2/(len(samples1)-2-2-1))
	p1_2= f.sf(f1_2,2,len(samples1)-2-2-1)
	#Lag phase 3
	a1_t_3 = np.array(a1_cha[0:len(samples1)-3])
   
	a1_t1_3= np.array(a1_cha[1:len(samples1)-2])
	a1_t1_t3_3=np.array(a1_cha[2:len(samples1)-1])
	a1_t2_3= np.array(a1_cha[3:len(samples1)])
	
	b1_t_3= np.array(b1_cha[0:len(samples1)-3])
	
	b1_t1_3= np.array(b1_cha[1:len(samples1)-2])
	b1_t1_t3_3=np.array(b1_cha[2:len(samples1)-1])
	b1_t2_3= np.array(b1_cha[3:len(samples1)])  
	  
	a1_t1_t2_3=  np.vstack((a1_t1_3,a1_t2_3,a1_t1_t3_3)).T
	a1_t1_t2_3= sm.add_constant(a1_t1_t2_3) 
	   
	b1_t1_t2_3= np.vstack((b1_t1_3,b1_t2_3,b1_t1_t3_3)).T
	b1_t1_t2_3= sm.add_constant(b1_t1_t2_3)		
	a1_t1_t2_b1_3 = np.vstack((a1_t1_3,a1_t2_3,a1_t1_t3_3,b1_t1_3,b1_t2_3,b1_t1_t3_3)).T
	a1_t1_t2_b1_3 = sm.add_constant(a1_t1_t2_b1_3)
	#F-test
	model = LinearRegression()
	model.fit(a1_t1_t2_b1_3,b1_t_3)
	b1_rssu_3 = np.sum((model.predict(a1_t1_t2_b1_3)-b1_t_3)**2)
	model.fit(b1_t1_t2_3,b1_t_3)
	b1_rssr_3 = np.sum((model.predict(b1_t1_t2_3)-b1_t_3)**2)
	f1_3 = ((b1_rssr_3-b1_rssu_3)/3)/(b1_rssu_3/(len(samples1)-7))
	p1_3 = f.sf(f1_3,3,len(samples1)-7)
	p=[p1_1,p1_2,p1_3]
	p_1=min(p)
	return p_1

def LassoRegression(degree, alpha):
	return Pipeline([
		("poly", PolynomialFeatures(degree=degree)),
		("std_scaler", StandardScaler()),
		("lasso_reg", Lasso(alpha=alpha))])

#Lasso regression:[degree:Highest power exponent; regularization] [parameters:different types of data for different training values]
def lasso_regress(gene_1,gene_2,degree,alpha):
	le_g1=len(gene_1)#长度
	co=max(gene_1)-min(gene_1)#极差
	st=co/le_g1 #步长 
	t=np.linspace(min(gene_1), max(gene_1), num=le_g1)
	#t= np.arange(min(gene_1), max(gene_1), st) #根据极差生成即将映射的y值   
	a1,b1=sor_g(gene_1,gene_2)#排序后的g1和g2
	a1_1=np.array(a1)#np下的g1
	a1_1= a1_1.reshape(-1, 1)
	b1_1=np.array(b1)
	b1_1=b1_1.reshape(-1,1)
	lasso1_reg = LassoRegression(degree,alpha) 
	lasso1_reg.fit(a1_1, t)
	gene2=lasso1_reg.predict(b1_1)
	gene1=lasso1_reg.predict(a1_1)
	return gene1,gene2

#Read data
os.chdir('*/GNIPLR') 

ax_nl=0.04
n_d4='BLCA_cancer.txt'
train_s=pd.read_csv(+str(n_d4),sep='\t')
train_s.T.to_csv('y1_10_d2_t.txt', sep='\t')
fw=open(str(n_d4[:-4])+'_gold.txt')
rels=[]
for p in fw:
	t1=p.split()
	rels.append(t1)
fw.close()  
fla=0
cancer1 = {}
samples1 = []
file1 = open("y1_10_d2_t.txt")
for p in file1:
	fla+=1
	t1 = p.split()
	n1=len(t1)
	samples1=t1[1:]
	cancer1[t1[0]]=[t1[i] for i in range(1,n1)]
if '0' in cancer1:
	del cancer1['0']
if 'wt' in cancer1:
	del cancer1['wt']
if 'strain' in cancer1:
	del cancer1['strain']
if 'Time' in cancer1:
	del cancer1['Time']
file1.close()
le=len(cancer1)
for sub in cancer1:
	cancer1[sub]=[float(x) for x in cancer1[sub]]   
#calculate and return the prediction
fw=open('e1_100_d2_p1.txt','w')
for i in range(len(rels)):
	entry_11=rels[i][0]
	entry_21=rels[i][1]
	gene_1,gene_2=lasso_regress(cancer1[entry_11],cancer1[entry_21],30,ax_nl) 
	p=granger(gene_1,gene_2)
	fw.write(str(entry_11)+'\t'+str(entry_21)+'\t'+str(p)+'\t'+str(rels[i][2])+'\n')
fw.close()   
#############################
relust1=pd.read_csv('e1_100_d2_p1.txt',sep='\t',names = range(4))
bbb=relust1.ix[:,[0,1,2,3]]
bbb.to_csv(str(n_d4[:-4])+'_lasso.txt', sep='\t',index=False,header=False)

