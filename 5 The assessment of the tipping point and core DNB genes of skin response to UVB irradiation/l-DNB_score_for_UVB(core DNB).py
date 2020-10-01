# === l-NDB(Top K DNB modules)===
# 修正p-value和ref_num
import os
import scipy.stats as stat
import math
import numpy as np
import random
import time
import multiprocessing
import pandas as pd


def ssn_score(deta,pcc,nn):
    if pcc==1:
        pcc=0.99999999
    if pcc==-1:
        pcc=-0.99999999
    z=deta/((1-pcc*pcc)/(nn-1))
    return z
    

def parallel_procedure(stage,normal,disease,ref,sd_mean,j,refnum,pvalue,uvb_dnb):
    begin=time.strftime("%Y-%m-%d:%H-%M-%S",time.localtime(time.time()))

    print("Stage: ",stage," Sample: ",j+1)

    network={}
    ssn={}
    full_ssn = pd.DataFrame(index=np.arange(len(ref.keys())), columns=['node1', 'node2', 'PCC1', 'PCC2', 'deltaPCC', 'p_value'])
    i=0
    for p in ref.keys():
        t=p.split()
        r1=ref[p]
        r2=stat.pearsonr(normal[t[0]]+[disease[t[0]][j]],normal[t[1]]+[disease[t[1]][j]])[0]
        r=r2-r1
        z = ssn_score(r,r1,refnum)
        p_value=2*(1-stat.norm.cdf(abs(z)))
        full_ssn.iloc[i, :] = [t[0], t[1], r1, r2, r, p_value]
        i = i + 1
        if p_value <= pvalue:
            r=r if r>0 else -r
            ssn[p]=r
            ssn[t[1]+"\t"+t[0]]=r
            
            if t[0] not in network.keys():
                network[t[0]]=[]
            network[t[0]].append(t[1])

            if t[1] not in network.keys():
                network[t[1]]=[]
            network[t[1]].append(t[0])
    selected_ssn = full_ssn[full_ssn['p_value'] < 0.05]
    selected_ssn.to_csv("SSN for {} in sample {}.csv".format(stage.strip('.csv'), j+1), index=False)

    ci = pd.DataFrame(columns=['score', 'sd', 'pcc_in', 'pcc_out'])
    for p in network.keys():
        if len(network[p])<3:
            continue
        
        sd=abs(disease[p][j]-sd_mean[p][1])/sd_mean[p][0]
        pcc_in=0
        pcc_out=0
        count=0
        for q in network[p]:
            sd+=abs(disease[q][j]-sd_mean[q][1])/sd_mean[q][0]
            pcc_in+=ssn[p+"\t"+q]

            for m in network[q]:
                if (m != p) & (m not in network[p]):
                    pcc_out+=ssn[q+"\t"+m]
                    count+=1
        sd/=len(network[p])+1
        pcc_in/=len(network[p])
        if count==0:
            continue
        pcc_out/=count
        if pcc_out==0:
            continue

        ci.loc[p, :] = [sd*pcc_in/pcc_out, sd, pcc_in, pcc_out]

    ci = ci.loc[uvb_dnb, :]
    ci = ci.sort_values(by='score', ascending=False)
    ci.index.name = 'gene_name'
    print('DNB gene num: ', ci.shape[0])
    ci.to_csv("Max_dnb_score_module in {} for sample {}.csv".format(stage, j + 1))

    # l-DNB mean
    score_list = list(ci['score'])
    score_mean = np.mean(score_list)
    module_num = len(score_list)
    print('module num: ', module_num)
    # print("Begin time is "+begin)
    # print("End time is "+time.strftime("%Y-%m-%d:%H-%M-%S",time.localtime(time.time())))
    print()
    return score_mean


if __name__ == "__main__":
   
    pvalue=1  #p-value threshold is set
    
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

    ESSN_node = pd.read_csv('Expanded_SSN_node_info.csv', header=0, index_col=0)
    uvb_dnb = ESSN_node[(ESSN_node['node_type1'] == 'dnb-marker') & (ESSN_node['node_type2'] != 'plain')].index
    print('dnb-marker strict num:', len(uvb_dnb))

    # uvb_dnb = ESSN_node[ESSN_node['node_type1'] == 'dnb-marker'].index
    # uvb_dnb = ESSN_node[(ESSN_node['node_type1'] == 'in-marker') & (ESSN_node['node_type2'] != 'plain')]

    stages = ["UVB_T1", "UVB_T2", "UVB_T3", "UVB_T4", "UVB_T5", "UVB_T6", "UVB_T7", "UVB_T8", "UVB_T9", "UVB_T10",
              "UVB_T11"]
    pool=multiprocessing.Pool(3)

    l_dnb_score = {}
    for stage in stages:
        file=stage + '.csv'
        f=open(file)
        disease={}
        title=[]
        flag=0
        for p in f:
            flag+=1
            t=p.split(sep=',')
            if flag==1:
                title=[str(k) for k in range(1,len(t))]
                continue
            disease[t[0]]=[float(t[k]) for k in range(1,len(t))]

        f.close()

        for j in range(len(title)):
            # pool.apply_async(parallel_procedure,(stage,normal,disease,title,ref,sd_mean,j,refnum,pvalue,))
            l_dnb_score[stage + '-' + str(j+1)] = \
                parallel_procedure(stage, normal, disease, ref, sd_mean, j, refnum, pvalue, uvb_dnb)
    pool.close()
    pool.join()
    # print(l_dnb_score)
    L_DNB_score = pd.DataFrame(l_dnb_score.values(), index=l_dnb_score.keys(), columns=['mean_score'])
    L_DNB_score.to_excel('L-DNB_score_UVB(dnb-marker strict).xlsx')

    L_DNB_result = pd.DataFrame(index=['UVB_T1', 'UVB_T2', 'UVB_T3', 'UVB_T4', 'UVB_T5', 'UVB_T6', 'UVB_T7', 'UVB_T8',
                                       'UVB_T9', 'UVB_T10', 'UVB_T11'], columns=['stage', 'score'])
    L_DNB_result.loc['UVB_T1', ['stage', 'score']] = ['D1-1h', L_DNB_score.loc[['UVB_T1-1', 'UVB_T1-2', 'UVB_T1-3'], 'mean_score'].mean()]
    L_DNB_result.loc['UVB_T2', ['stage', 'score']] = ['D1-4h', L_DNB_score.loc[['UVB_T2-1', 'UVB_T2-2'], 'mean_score'].mean()]
    L_DNB_result.loc['UVB_T3', ['stage', 'score']] = ['D1-8h', L_DNB_score.loc[['UVB_T3-1', 'UVB_T3-2'], 'mean_score'].mean()]
    L_DNB_result.loc['UVB_T4', ['stage', 'score']] = ['D2-1h', L_DNB_score.loc[['UVB_T4-1', 'UVB_T4-2'], 'mean_score'].mean()]
    L_DNB_result.loc['UVB_T5', ['stage', 'score']] = ['D3-1h', L_DNB_score.loc[['UVB_T5-1', 'UVB_T5-2'], 'mean_score'].mean()]
    L_DNB_result.loc['UVB_T6', ['stage', 'score']] = ['D4-1h', L_DNB_score.loc[['UVB_T6-1', 'UVB_T6-2', 'UVB_T6-3'], 'mean_score'].mean()]
    L_DNB_result.loc['UVB_T7', ['stage', 'score']] = ['D5-1h', L_DNB_score.loc[['UVB_T7-1', 'UVB_T7-2'], 'mean_score'].mean()]
    L_DNB_result.loc['UVB_T8', ['stage', 'score']] = ['D6-1h', L_DNB_score.loc[['UVB_T8-1', 'UVB_T8-2', 'UVB_T8-3'], 'mean_score'].mean()]
    L_DNB_result.loc['UVB_T9', ['stage', 'score']] = ['D6-4h', L_DNB_score.loc[['UVB_T9-1', 'UVB_T9-2'], 'mean_score'].mean()]
    L_DNB_result.loc['UVB_T10', ['stage', 'score']] = ['D6-8h', L_DNB_score.loc[['UVB_T10-1', 'UVB_T10-2'], 'mean_score'].mean()]
    L_DNB_result.loc['UVB_T11', ['stage', 'score']] = ['D7-1h', L_DNB_score.loc[['UVB_T11-1', 'UVB_T11-2'], 'mean_score'].mean()]

    L_DNB_result.to_excel('L_DNB_result_UVB(dnb-marker strict).xlsx')

