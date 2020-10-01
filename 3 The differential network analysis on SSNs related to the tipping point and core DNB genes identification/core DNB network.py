# 构建core network，提取 core dnb的一阶邻居
import numpy as np
import pandas as pd
# DNB gene
dnb = pd.read_excel('l-DNB_gene_for_UVB.xlsx', sheet_name=3, header=0)
dnb = set(dnb.iloc[:, 0])

# get the core dnb(dnb and skin signature)
key_dnb = pd.read_csv('key_dnb_strict.csv', header=0, index_col=0)
core_dnb = key_dnb[(key_dnb['node_type1'] == 'dnb-marker') & (key_dnb['node_type2'] != 'plain')].index
core_dnb = set(core_dnb)
candidate_dnb = key_dnb[(key_dnb['node_type1'] == 'mDNB1') & (key_dnb['node_type2'] != 'plain')].index
candidate_dnb = set(candidate_dnb)

# SSN Delta PCC
Delta_PCC = pd.read_csv('SSN_delta_PCC.csv', header=0, index_col=0)
# Delta_PCC.index = Delta_PCC['node1'] + '\t' + Delta_PCC['node2']
TP_col = ['day-6-UVB-1', 'day-6-UVB-2', 'day-6-UVB-3']
TP_delta_pcc = Delta_PCC.loc[:, TP_col].mean(axis=1)

# get core net
T8_SSN = pd.read_csv('SSN for UVB_T8.csv', header=0, index_col=None)

# core dnb network node information
node = set(T8_SSN.iloc[:, 0]) | set(T8_SSN.iloc[:, 1])
core_net_node = pd.DataFrame(0, index=node, columns=['label'])
for idx in core_net_node.index:
    if idx in core_dnb:
        core_net_node.loc[idx, 'label'] = 111
        print('111')
    elif idx in candidate_dnb:
        core_net_node.loc[idx, 'label'] = 11
    elif idx in dnb:
        core_net_node.loc[idx, 'label'] = 1
    else:
        core_net_node.loc[idx, 'label'] = 0

# core dnb network edge information
edge = T8_SSN.iloc[:, 0] + '\t' + T8_SSN.iloc[:, 1]
core_net_edge = pd.DataFrame(index=edge, columns=['node1', 'node2', 'delta_pcc'])
core_net_edge.loc[:, 'node1'] = T8_SSN.loc[:, 'node1'].values
core_net_edge.loc[:, 'node2'] = T8_SSN.loc[:, 'node2'].values
core_net_edge.loc[:, 'delta_pcc'] = np.abs(TP_delta_pcc[edge].values)


core_net_node.to_csv('core_net_node.csv', index=True)
core_net_edge.to_csv('core_net_edge.csv', index=False)
