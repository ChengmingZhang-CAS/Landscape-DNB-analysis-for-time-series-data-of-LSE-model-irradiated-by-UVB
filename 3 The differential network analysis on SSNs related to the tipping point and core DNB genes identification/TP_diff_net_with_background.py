# 构建差异网络
import pandas as pd

TP1_DNB = pd.read_excel('l-DNB_gene_for_UVB.xlsx', sheet_name=0, header=0, index_col=None)
TP2_DNB = pd.read_excel('l-DNB_gene_for_UVB.xlsx', sheet_name=3, header=0, index_col=None)
TP1_DNB = set(TP1_DNB.iloc[:, 0])
TP2_DNB = set(TP2_DNB.iloc[:, 0])


# T7_SSN
T7_SSN = pd.read_csv('SSN for UVB_T7.csv', header=0, index_col=None)
t7_ssn = T7_SSN.iloc[:, 0] + '\t' + T7_SSN.iloc[:, 1]
t7_edge = set(t7_ssn)
t7_node = set(T7_SSN.iloc[:, 0]) | set(T7_SSN.iloc[:, 1])

# T8_SSN
T8_SSN = pd.read_csv('SSN for UVB_T8.csv', header=0, index_col=None)
t8_ssn = T8_SSN.iloc[:, 0] + '\t' + T8_SSN.iloc[:, 1]
tp2_edge = set(t8_ssn)
tp2_node = set(T8_SSN.iloc[:, 0]) | set(T8_SSN.iloc[:, 1])


# T9_SSN
T9_SSN = pd.read_csv('SSN for UVB_T9.csv', header=0, index_col=None)
t9_ssn = T9_SSN.iloc[:, 0] + '\t' + T9_SSN.iloc[:, 1]
t9_edge = set(t9_ssn)
t9_node = set(T9_SSN.iloc[:, 0]) | set(T9_SSN.iloc[:, 1])

bg_node = t7_node | tp2_node | t9_node
bg_edge = t7_edge | tp2_edge | t9_edge


diff_bg_edge = pd.DataFrame(index=bg_edge, columns=['node1', 'node2', 'TP2_T7', 'TP2_T9', 'T7', 'T8', 'T9'])
for idx in diff_bg_edge.index:
    diff_bg_edge.loc[idx, ['node1', 'node2']] = idx.split('\t')

diff_bg_edge.loc[t7_edge, 'T7'] = "1"
diff_bg_edge.loc[tp2_edge, 'T8'] = "1"
diff_bg_edge.loc[t9_edge, 'T9'] = "1"


t7_unique78 = t7_edge - tp2_edge
t8_unique78 = tp2_edge - t7_edge
diff_bg_edge.loc[t7_unique78, 'TP2_T7'] = "1-0"
diff_bg_edge.loc[t8_unique78, 'TP2_T7'] = "0-1"

t8_unique89 = tp2_edge - t9_edge
t9_unique89 = t9_edge - tp2_edge
diff_bg_edge.loc[t8_unique89, 'TP2_T9'] = "1-0"
diff_bg_edge.loc[t9_unique89, 'TP2_T9'] = "0-1"
diff_bg_edge = diff_bg_edge.fillna("0")


# 网络的node信息
diff_bg_node = pd.DataFrame(0, index=bg_node, columns=['type'])
diff_bg_node.loc[TP2_DNB, :] = 11


diff_bg_edge.to_csv('diff_bg_edge(T7-T8-T9).csv', index=False, header=True)
diff_bg_node.to_csv('diff_bg_node(T7-T8-T9).csv', index=True, header=True)


