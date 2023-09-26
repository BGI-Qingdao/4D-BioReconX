#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'


data = pd.read_csv('pcg.csv',sep='\t',header=0)
pcgs_info = np.loadtxt(sys.argv[1], delimiter='\t', dtype=str)
classes = pcgs_info[:,2]
A_data = []
P_data = []

time_list2 = ['WT', '0hpa','12hpa','36hpa','3dpa','5dpa','7dpa','10dpa','14dpa']
time_list = ['WT','0hpa1','12hpa2','36hpa2','3dpa2','5dpa1','7dpa2','10dpa1','14dpa1']

for name, smes, c in pcgs_info:
    A_list = []
    P_list = []
    for indv in time_list:
        ddd0 = data[(data['gene']==name)&(data['AP']=='AP_0')&(data['sample']==indv)]
        ddd1 = data[(data['gene']==name)&(data['AP']=='AP_1')&(data['sample']==indv)]
        meanv = np.mean([ddd0['scaled_density'].to_numpy()[0], ddd1['scaled_density'].to_numpy()[0]])
        A_list.append(meanv)
        ddd0 = data[(data['gene']==name)&(data['AP']=='AP_8')&(data['sample']==indv)]
        ddd1 = data[(data['gene']==name)&(data['AP']=='AP_9')&(data['sample']==indv)]
        meanv = np.mean([ddd0['scaled_density'].to_numpy()[0], ddd1['scaled_density'].to_numpy()[0]])
        P_list.append(meanv)
    A_data.append(A_list)
    P_data.append(P_list)


#--------------------
cm = 1/2.54
fntsz = 12
height = 22.5
width = 10
#--------------------
# draw A
#--------------------
A_array = np.array(A_data)
plt.figure(figsize=(10*cm,height*cm))
cg = sns.heatmap(A_array,cmap="coolwarm",vmin=-2,vmax=2,linewidths=.5,cbar=False)
plt.xticks(ticks = np.arange(len(time_list))+0.5, labels=time_list2,rotation=45,fontsize=fntsz)
plt.yticks(ticks = np.arange(pcgs_info.shape[0])+0.5, labels=pcgs_info[:,0],rotation=0,fontsize=fntsz)
plt.tight_layout()
plt.box(False)
plt.savefig('Ascale.pdf', format='pdf')
plt.close()


#--------------------
# draw P
#--------------------
P_array = np.array(P_data)
plt.figure(figsize=(10*cm,height*cm))
cg = sns.heatmap(P_array,cmap="coolwarm",vmin=-2,vmax=2,linewidths=.5,cbar=False)
plt.xticks(ticks = np.arange(len(time_list))+0.5, labels=time_list2,rotation=45,fontsize=fntsz)
plt.yticks(ticks = np.arange(pcgs_info.shape[0])+0.5, labels=pcgs_info[:,0],rotation=0,fontsize=fntsz)
plt.box(False)
plt.tight_layout()
plt.savefig('Pscale.pdf', format='pdf')
plt.close()

'''
AP_label = ['AP_0','AP_0','AP_2','AP_3','AP_4','AP_5','AP_6','AP_7','AP_8','AP_9']

all_in_one = []
for name, smes, c in pcgs_info:
    AP_list = []
    for indv in time_list:
        ddd = data[(data['gene']==name)&(data['sample']==indv)]
        AP_list.append(ddd['scaled_density'].to_numpy())
    AP_array = np.array(AP_list)
    a = AP_array.reshape(-1)
    #plt.figure(figsize=(6,4))
    #sns.heatmap(AP_array,cmap="coolwarm",vmin=-2,vmax=3,linewidths=.5)
    #plt.xticks(ticks = np.arange(len(AP_label))+0.5,labels=AP_label,rotation=45)
    #plt.yticks(ticks = np.arange(len(time_list))+0.5, labels=time_list,rotation=0)
    #plt.title(name)
    #plt.tight_layout()
    #name = name.replace("/","_")
    #plt.savefig(f'{name}.{smes}.png')
    #plt.close()
    #all_in_one.append(AP_array.reshape(-1))
    all_in_one.append(a)


All_array = np.array(all_in_one)
plt.figure(figsize=(30,30))
sns.heatmap(All_array,cmap="coolwarm",vmin=-2,vmax=2,linewidths=.5)
plt.xticks(ticks = np.arange(len(time_list))*10 +5 , labels=time_list2,rotation=0)
#plt.hlines(30,xmin=0,xmax=10)
#plt.yticks(ticks = np.arange(pcgs_info.shape[0])+0.5, fontstyle='italic',labels=yt,rotation=0)
plt.yticks(ticks = np.arange(pcgs_info.shape[0])+0.5, labels=pcgs_info[:,0],rotation=0)
plt.xticks(fontstyle='italic')
plt.box(False)
plt.title('$\it{yt[0]}$')
plt.tight_layout()
plt.savefig('AP_all.pdf', format='pdf')
plt.close()
'''


