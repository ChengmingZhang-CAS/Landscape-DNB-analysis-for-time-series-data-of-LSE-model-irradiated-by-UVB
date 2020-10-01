# Public domain UV induced gene change
"""
    计算DNB gene是否显著富集到已知的一些UV相关的marker gene，或者说查看DNB gene和marker gene的
overlap是否显著，同时绘制Venn图
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

df = pd.read_excel('RNA-seq_data_FPKM_symbol.xlsx', header=0)
all_gene = set(df.iloc[:, 0])


# 提前DNB基因
def read_csv(stage_id, sample_num):
    file = "Max_dnb_score_module in UVB_T{} for sample {}.csv".format(stage_id, sample_num)
    df = pd.read_csv(file, header=0, index_col=0)
    return set(df.index[0:600])


# 第一个关键时间点T1 & T2（D1-1h & D1-4h）
T1_dnb1 = read_csv(1, 1)
T1_dnb2 = read_csv(1, 2)
T1_dnb3 = read_csv(1, 3)
T1_dnb = (T1_dnb1 & T1_dnb2) | (T1_dnb1 & T1_dnb3) | (T1_dnb2 & T1_dnb3)

T2_dnb1 = read_csv(2, 1)
T2_dnb2 = read_csv(2, 2)
T2_dnb = T2_dnb1 & T2_dnb2

TP1_dnb = T1_dnb | T2_dnb

print('T1(D1-1h) dnb len: ', len(T1_dnb))
print('T2(D1-4h) dnb len: ', len(T2_dnb))
print('TP1 dnb len: ', len(TP1_dnb))

# 第二个关键时间点T8（D6-1h）
T8_dnb1 = read_csv(8, 1)
T8_dnb2 = read_csv(8, 2)
T8_dnb3 = read_csv(8, 3)

T8_dnb = (T8_dnb1 & T8_dnb2) | (T8_dnb1 & T8_dnb3) | (T8_dnb2 & T8_dnb3)
TP2_dnb = T8_dnb
print('T8(D6-1h) dnb1 len: ', len(T8_dnb1))
print('T8(D6-1h) dnb2 len: ', len(T8_dnb2))
print('T8(D6-1h) dnb3 len: ', len(T8_dnb3))
print('TP2 dnb len: ', len(TP2_dnb))

marker = pd.read_excel('Public domain UV induced gene change.xlsx', header=0)
marker = set(marker.iloc[:, 0])

# fisher exact test
n11 = len(marker & TP2_dnb)
n12 = len(marker) - n11
n21 = len(TP2_dnb) - n11
n22 = len(all_gene | marker) - len(marker | TP2_dnb)

stats.fisher_exact([[n11, n12], [n21, n22]], alternative='greater')

venn2([TP2_dnb, marker], set_labels=('DNB gene', 'Marker gene'))
plt.title('Overlap of DNB and marker gene')
plt.show()

neighbor = []

# load the PPI network
f = open("background900.txt", "r")
line = f.readlines()
for l in line:
    t = l.split()
    if t[0] in marker:
        neighbor.append(t[1])
    if t[1] in marker:
        neighbor.append(t[0])

f.close()

neighbor = set(neighbor)
expand_marker = marker | neighbor
pd.Series(list(expand_marker)).to_csv("skin related signature.csv")
# fisher exact test
n11 = len(expand_marker & TP2_dnb)
n12 = len(expand_marker) - n11
n21 = len(TP2_dnb) - n11
n22 = len(all_gene | expand_marker) - len(expand_marker | TP2_dnb)

plt.figure()
stats.fisher_exact([[n11, n12], [n21, n22]], alternative='greater')

venn2([TP2_dnb, expand_marker], set_labels=('DNB gene', 'Expanded marker'))
plt.title('Overlap of DNB and skin related signature')
plt.show()

#
#
# neighbor = []
# # load the PPI network
# f = open("background900.txt", "r")
# line = f.readlines()
# for l in line:
#     t = l.split()
#     if t[0] in TP2_dnb:
#         neighbor.append(t[1])
#     if t[1] in TP2_dnb:
#         neighbor.append(t[0])
#
# f.close()
#
# neighbor = set(neighbor)
# expand_DNB = TP2_dnb | neighbor
# len(expand_DNB & marker)
