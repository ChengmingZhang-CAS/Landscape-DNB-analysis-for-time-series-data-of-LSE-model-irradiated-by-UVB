# 绘制DNB的landscape
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d


# 提前DNB基因
def read_csv(stage_id, sample_num):
    file = "Max_dnb_score_module in UVB_T{} for sample {}.csv".format(stage_id, sample_num)
    df = pd.read_csv(file, header=0, index_col=0)
    return df.index


def key_gene(stage_id, sample_num):
    """
        计算目标时间点key gene，这里仅计算至少在两个样本中出现的基因
    """
    if sample_num == 3:
        file = "Max_dnb_score_module in UVB_T{} for sample {}.csv".format(stage_id, 1)
        df1 = pd.read_csv(file, header=0, index_col=0)
        idx1 = df1.index

        file = "Max_dnb_score_module in UVB_T{} for sample {}.csv".format(stage_id, 2)
        df2 = pd.read_csv(file, header=0, index_col=0)
        idx2 = df2.index

        file = "Max_dnb_score_module in UVB_T{} for sample {}.csv".format(stage_id, 3)
        df3 = pd.read_csv(file, header=0, index_col=0)
        idx3 = df3.index

        idx = (idx1 & idx2) | (idx1 & idx3) | (idx2 & idx3)

    else:
        file = "Max_dnb_score_module in UVB_T{} for sample {}.csv".format(stage_id, 1)
        df1 = pd.read_csv(file, header=0, index_col=0)
        idx1 = df1.index

        file = "Max_dnb_score_module in UVB_T{} for sample {}.csv".format(stage_id, 2)
        df2 = pd.read_csv(file, header=0, index_col=0)
        idx2 = df2.index

        idx = idx1 & idx2

    return idx


sample_num_lst = [2, 2, 2, 2, 3, 2, 3, 2, 2, 2]  # T2-T11
gene = key_gene(1, 3)
hub = gene
for i in range(len(sample_num_lst)):
    print('当前stage为：', i+1, len(hub))
    stage_gene = key_gene(i+2, sample_num_lst[i])
    hub = hub & stage_gene
print('当前stage为：', i+2, len(hub))


