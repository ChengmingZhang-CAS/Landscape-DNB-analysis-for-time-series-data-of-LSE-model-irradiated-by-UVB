# 修正p-value和ref_num
import os
import scipy.stats as stat
import math
import numpy as np
import random
import time
import multiprocessing
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import decomposition as skldec  # 用于主成分分析降维的包


def ssn_score(deta,pcc,nn):
    if pcc==1:
        pcc=0.99999999
    if pcc==-1:
        pcc=-0.99999999
    z=deta/((1-pcc*pcc)/(nn-1))
    return z
    

def parallel_procedure(normal, disease, ref, j, title, refnum):
    print("present sample: ", title[j])

    ssn_z = pd.Series(index=ref.keys())
    for p in ref.keys():
        t=p.split()
        r1=ref[p]
        r2=stat.pearsonr(normal[t[0]]+[disease[t[0]][j]],normal[t[1]]+[disease[t[1]][j]])[0]
        r=r2-r1
        z = ssn_score(r, r1, refnum)
        ssn_z[p] = z

    return ssn_z


if __name__ == "__main__":

    refnum=0  
    normal={}
    f=open("Reference_samples.csv")
    flag=0
    for p in f:
        flag+=1
        if flag==1:
            t=p.split(',')
            refnum=len(t)-1
            continue
        t=p.split(',')
        normal[t[0]]=[float(t[i]) for i in range(1,len(t))]
    f.close()

    sd_mean={}
    for key in normal.keys():
        sd_mean[key]=[np.std(normal[key]),np.mean(normal[key])]

    f=open("reference_network.txt")
    network={}
    ref={}
    for p in f:
        t=p.split()
        ref[t[0]+"\t"+t[1]]=float(t[2])

    f.close()

    file = 'CTRL_UVB_data.csv'
    f = open(file)
    disease={}
    title=[]
    flag=0
    for p in f:
        flag+=1
        t=p.split(sep=',')
        if flag==1:
            a = t
            title=[t[k] for k in range(1,len(t))]
            title[-1] = title[-1].strip('\n')
            continue
        disease[t[0]]=[float(t[k]) for k in range(1,len(t))]

    f.close()

    SSN_Z = pd.DataFrame(index=ref.keys(), columns=title)
    for j in range(len(title)):
        ssn_z = parallel_procedure(normal, disease, ref, j, title, refnum)
        SSN_Z.loc[:, title[j]] = ssn_z

    SSN_Z.to_csv('SSN_Z.csv')

    # PCA 可视化Z
    data = SSN_Z.T
    # data = (data - data.mean()) / data.std()  # Z-score标准化
    # 根据两个最大的主成分进行绘图
    pca = skldec.PCA(n_components=0.95, whiten=True)  # 选择方差95%的占比
    pca.fit(data)  # 主成分分析时每一行是一个输入数据
    result = pca.transform(data)  # 计算结果
    pc_df = pd.DataFrame(data=result[:, 0:2], index=data.index, columns=['PC1', 'PC2'])
    print(result.shape)

    plt.figure()
    plt.scatter(pc_df.iloc[0:12, 0], pc_df.iloc[0:12, 1], c='g', s=100, label='CRTL')
    plt.scatter(pc_df.iloc[12:, 0], pc_df.iloc[12:, 1], c='r', s=100, label='UVB')
    TP_idx = ['day-6-UVB-1', 'day-6-UVB-2', 'day-6-UVB-3']
    pc_TP = pc_df.loc[TP_idx, :]
    plt.scatter(pc_TP.iloc[:, 0], pc_TP.iloc[:, 1], c='b', s=100, label='TP-UVB')

    for i in range(result[:, 0].size):
        plt.text(pc_df.iloc[i, 0] + 0, pc_df.iloc[i, 1], pc_df.index[i][4])  # 在每个点边上绘制数据名称

    plt.xlabel('PC1', fontsize=20)  # 绘制x轴标签
    plt.ylabel('PC2', fontsize=20)  # 绘制y轴标签
    plt.title('PCA based on SSN statistic Z', fontsize=30)
    plt.legend()
    plt.show()

