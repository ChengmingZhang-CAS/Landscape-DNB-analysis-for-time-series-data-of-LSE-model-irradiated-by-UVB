# ---DEGs for UVB vs CON and S1 vs UVB---
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import warnings

# 载入数据
df = pd.read_excel('RNA-seq_data_FPKM_symbol.xlsx', header=0, index_col=0)
print(df.head())
print(df.shape)


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def calculate_DEGs(group1, group2, group_info):
    print(group_info)
    data = group1.join(group2)
    gene_info = pd.DataFrame(index=data.index, columns=['p_value', 'q_value', 'fold_change'])
    for gene in gene_info.index:
        if (group1.loc[gene, :].sum() != 0) & (group2.loc[gene, :].sum() != 0) & (np.std(data.loc[gene, :]) != 0):
            gene_info.loc[gene, 'fold_change'] = group2.loc[gene, :].mean() / group1.loc[gene, :].mean()
            gene_info.loc[gene, 'p_value'] = stats.ttest_ind(group1.loc[gene, :], group2.loc[gene, :], nan_policy='omit')[1]
        elif (group1.loc[gene, :].sum() != 0) & (group2.loc[gene, :].sum() == 0):
            # print("d-0", gene)
            # print(group1.loc[gene, :], group2.loc[gene, :])
            gene_info.loc[gene, 'fold_change'] = 0
            gene_info.loc[gene, 'p_value'] = stats.ttest_ind(group1.loc[gene, :], group2.loc[gene, :], nan_policy='omit')[1]
        elif (group1.loc[gene, :].sum() == 0) & (group2.loc[gene, :].sum() != 0):
            # print("0-d", gene)
            # print(group1.loc[gene, :], group2.loc[gene, :])
            gene_info.loc[gene, 'fold_change'] = 99999999
            gene_info.loc[gene, 'p_value'] = stats.ttest_ind(group1.loc[gene, :], group2.loc[gene, :], nan_policy='omit')[1]
        else:
            # print("0-0", gene)
            # print(group1.loc[gene, :], group2.loc[gene, :])
            gene_info.loc[gene, 'fold_change'] = 1
            gene_info.loc[gene, 'p_value'] = 0.99999999

    gene_info.loc[:, 'q_value'] = p_adjust_bh(gene_info.loc[:, 'p_value'])
    DEGs = gene_info[(gene_info['p_value'] < 0.05) &
                     (((gene_info['fold_change'] > 0) & (gene_info['fold_change'] < 1/1.2)) |
                      ((gene_info['fold_change'] < 99999999) & (gene_info['fold_change'] > 1.2)))]

    print('DEGs num:', len(DEGs))
    write = pd.ExcelWriter(group_info + '.xlsx')
    DEGs.to_excel(write, sheet_name='DEGs')
    gene_info.to_excel(write, sheet_name='gene_info')
    write.close()

    return len(DEGs)


# Tipping point 1
#  day1-1h day1-4h day1-8h vs day2-1h day3-1h
TP1_group1 = df.loc[:, ['day-1-UVB-1', 'day-1-UVB-2', 'day-1-UVB-3', 'day-1-UVB-4', 'day-1-UVB-5', 'day-1-UVB-6',
                        'day-1-UVB-7']]
TP1_group2 = df.loc[:, ['day-2-UVB-1', 'day-2-UVB-2', 'day-3-UVB-1', 'day-3-UVB-2']]

calculate_DEGs(TP1_group1, TP1_group2, 'TP1-DEGs (day1 vs (day2-1h & day3-1h))')

# Tipping point 2
#  day6-1h vs day4-1h day5-1h
TP2_group1 = df.loc[:, ['day-6-UVB-1', 'day-6-UVB-2', 'day-6-UVB-3']]
TP2_group2 = df.loc[:, ['day-4-UVB-1', 'day-4-UVB-2', 'day-4-UVB-3', 'day-5-UVB-1', 'day-5-UVB-2']]

calculate_DEGs(TP2_group1, TP2_group2, 'TP2-DEGs (day6-1h vs (day4-1h & day5-1h))')

#  day6-1h vs day6-4h day6-8h
TP2_group1 = df.loc[:, ['day-6-UVB-1', 'day-6-UVB-2', 'day-6-UVB-3']]
TP2_group2 = df.loc[:, ['day-6-UVB-4', 'day-6-UVB-5', 'day-6-UVB-6', 'day-6-UVB-7']]

calculate_DEGs(TP2_group1, TP2_group2, 'TP2-DEGs (day6-1h vs (day6-4h & day6-8h))')










def correct_pvalues_for_multiple_testing(pvalues, correction_type="Benjamini-Hochberg"):
    """
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = pvalues.shape[0]
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues

