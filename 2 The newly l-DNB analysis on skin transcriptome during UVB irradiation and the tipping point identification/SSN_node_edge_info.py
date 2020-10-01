# 绘制DNB的landscape
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D   # 需要导入模块Axes3D

l_score = pd.read_excel('L-DNB_score_UVB top600.xlsx')

is_plot = False
if is_plot:
    plt.plot(np.arange(l_score.shape[0]), l_score['score'], 'r', linewidth=5)
    # plt.plot([1, 1], [0, l_score['score'][1]], 'r', linewidth=1)
    plt.xticks(np.arange(l_score.shape[0]), l_score['stage'], fontsize=15)
    plt.ylim(0, 2.5)
    # plt.ylim(0, np.max(l_score['score']+0.5))
    plt.xlabel('time point', fontsize=20)
    plt.ylabel('DNB score', fontsize=20)
    plt.title('DNB score landscape for UVB', fontsize=40)
    plt.show()


# 提前DNB基因
def read_csv(stage_id, sample_num):
    file = "Max_dnb_score_module in UVB_T{} for sample {}.csv".format(stage_id, sample_num)
    df = pd.read_csv(file, header=0, index_col=0)
    return df.index[0:600]


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


# 绘制DNB单样本网络图
def get_SSN_edge(stage_id, sample_num):
    """
        计算目标时间点所有样本的SSN：具体操作为保留至少出现在两个样本的边
    """
    df = pd.DataFrame(columns=['node1', 'node2'])
    if sample_num == 3:
        ssn1 = pd.read_csv('SSN for UVB_T' + str(stage_id) + ' in sample 1.csv')
        ssn1_edge = ssn1.iloc[:, 0] + '_' + ssn1.iloc[:, 1]

        ssn2 = pd.read_csv('SSN for UVB_T' + str(stage_id) + ' in sample 2.csv')
        ssn2_edge = ssn2.iloc[:, 0] + '_' + ssn2.iloc[:, 1]

        ssn3 = pd.read_csv('SSN for UVB_T' + str(stage_id) + ' in sample 3.csv')
        ssn3_edge = ssn3.iloc[:, 0] + '_' + ssn3.iloc[:, 1]

        ssn_edge = (set(ssn1_edge) & set(ssn2_edge)) | (set(ssn1_edge) & set(ssn3_edge)) | (set(ssn2_edge) & set(ssn3_edge))
    else:
        ssn1 = pd.read_csv('SSN for UVB_T' + str(stage_id) + ' in sample 1.csv')
        ssn1_edge = ssn1.iloc[:, 0] + '_' + ssn1.iloc[:, 1]

        ssn2 = pd.read_csv('SSN for UVB_T' + str(stage_id) + ' in sample 2.csv')
        ssn2_edge = ssn2.iloc[:, 0] + '_' + ssn2.iloc[:, 1]

        ssn_edge = set(ssn1_edge) & set(ssn2_edge)

    for edge in ssn_edge:
        df = df.append(pd.DataFrame({'node1': [edge.split('_')[0]], 'node2': [edge.split('_')[1]]}), ignore_index=True)
    df.to_csv('SSN for UVB_T' + str(stage_id) + '.csv', header=True, index=False)
    print('ssn edge num: ', len(ssn_edge))
    return ssn_edge


sample_num_lst = [3, 2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 3]  # 每个stage的样本个数
total_ssn = pd.DataFrame(columns=['node1', 'n1_isTP1_dnb', 'n1_isTP2_dnb', 'node2', 'n2_isTP1_dnb', 'n2_isTP2_dnb',
                                  'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12'])
for i in range(12):
    print('当前stage为：', i+1)
    temp_ssn_edge = get_SSN_edge(i+1, sample_num_lst[i])
    for edge in temp_ssn_edge:
        if edge not in total_ssn.index:
            total_ssn.loc[edge, ['node1', 'node2', 'T' + str(i+1)]] = [edge.split('_')[0], edge.split('_')[1], 1]
        else:
            total_ssn.loc[edge, ['T' + str(i+1)]] = 1

print('ssn background edge num: ', total_ssn.shape[0])

total_ssn.loc[total_ssn.loc[:, 'node1'].isin(TP1_dnb), 'n1_isTP1_dnb'] = 1
total_ssn.loc[total_ssn.loc[:, 'node1'].isin(TP2_dnb), 'n1_isTP2_dnb'] = 1
total_ssn.loc[total_ssn.loc[:, 'node2'].isin(TP1_dnb), 'n2_isTP1_dnb'] = 1
total_ssn.loc[total_ssn.loc[:, 'node2'].isin(TP2_dnb), 'n2_isTP2_dnb'] = 1

total_ssn.to_csv('total_ssn_for_UVB.csv')
ssn_edge_info = total_ssn.drop(columns=['n1_isTP1_dnb', 'n1_isTP2_dnb', 'n2_isTP1_dnb', 'n2_isTP2_dnb'])
ssn_edge_info.fillna(0).to_csv('ssn_edge_info.csv', index=False)

# 设置background network node的属性
ssn_node_info = pd.DataFrame(columns=['node_id', 'node_type'])
ssn_node_info.loc[:, 'node_id'] = list(set(total_ssn.loc[:, 'node1']) | set(total_ssn.loc[:, 'node2']))
ssn_node_info.loc[:, ['node_type']] = 0
ssn_node_info.loc[ssn_node_info.loc[:, 'node_id'].isin(TP1_dnb), 'node_type'] = 1
ssn_node_info.loc[ssn_node_info.loc[:, 'node_id'].isin(TP2_dnb), 'node_type'] = 2
ssn_node_info.loc[ssn_node_info.loc[:, 'node_id'].isin(TP1_dnb & TP2_dnb), 'node_type'] = 12


def stage_dnb_score(stage_id, sample_num):
    """
        计算目标时间点DNB基因的score值，这里仅计算至少在两个样本中出现的DNB基因
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
        df = pd.DataFrame(index=idx, columns=['score1', 'score2', 'score3'])
        df.loc[idx & idx1, 'score1'] = df1.loc[idx & idx1, 'score']
        df.loc[idx & idx2, 'score2'] = df2.loc[idx & idx2, 'score']
        df.loc[idx & idx3, 'score3'] = df3.loc[idx & idx3, 'score']

        df.loc[:, 'mean_score'] = df.mean(axis=1)
    else:
        file = "Max_dnb_score_module in UVB_T{} for sample {}.csv".format(stage_id, 1)
        df1 = pd.read_csv(file, header=0, index_col=0)
        idx1 = df1.index

        file = "Max_dnb_score_module in UVB_T{} for sample {}.csv".format(stage_id, 2)
        df2 = pd.read_csv(file, header=0, index_col=0)
        idx2 = df2.index

        idx = idx1 & idx2
        df = pd.DataFrame(index=idx, columns=['score1', 'score2'])
        df.loc[idx & idx1, 'score1'] = df1.loc[idx & idx1, 'score']
        df.loc[idx & idx2, 'score2'] = df2.loc[idx & idx2, 'score']

        df.loc[:, 'mean_score'] = df.mean(axis=1)

    return df


# 计算每个时间解阶段的dnb分值
sample_num_lst = [3, 2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 3]  # 每个stage的样本个数
for i in range(12):
    print('当前stage为：', i+1)
    temp_stage_score = stage_dnb_score(i+1, sample_num_lst[i])
    in_set = set(ssn_node_info['node_id']) & set(temp_stage_score.index)
    ssn_node_info.loc[ssn_node_info.loc[:, 'node_id'].isin(temp_stage_score.index), 'T'+str(i+1)] = temp_stage_score.loc[in_set, 'mean_score'].values

ssn_node_info = ssn_node_info.fillna(0)
ssn_node_info.to_csv('ssn_node_info.csv', index=False)

