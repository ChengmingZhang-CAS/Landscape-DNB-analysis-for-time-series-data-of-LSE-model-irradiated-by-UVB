# 绘制DNB的landscape
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D   # 需要导入模块Axes3D

l_score = pd.read_excel('L_DNB_result_UVB top600.xlsx')

plt.plot(np.arange(l_score.shape[0]), l_score['score'], 'r', linewidth=5)
# plt.plot([1, 1], [0, l_score['score'][1]], 'r', linewidth=1)
plt.xticks(np.arange(l_score.shape[0]), l_score['stage'], fontsize=15)
plt.ylim(0, 3)
# plt.ylim(0, np.max(l_score['score']+0.5))
plt.xlabel('time point', fontsize=20)
plt.ylabel('DNB score', fontsize=20)
plt.title('l-DNB score for UVB top600', fontsize=40)
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


write = pd.ExcelWriter('l-DNB_gene_for_UVB.xlsx')
pd.DataFrame(TP1_dnb, columns=['gene_name']).to_excel(write, sheet_name='D1-1h_4h_gene', index=None)
pd.DataFrame(T1_dnb, columns=['gene_name']).to_excel(write, sheet_name='D1-1h_DNB_gene', index=None)
pd.DataFrame(T2_dnb, columns=['gene_name']).to_excel(write, sheet_name='D1-4h_DNB_gene', index=None)
pd.DataFrame(TP2_dnb, columns=['gene_name']).to_excel(write, sheet_name='D6-1h_DNB_gene', index=None)

write.close()


# 绘制DNB gene表达landscape
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


# ********************************************************************************************************************
# 计算每个时间解阶段的dnb分值(Tipping point 1)
init_df = pd.DataFrame(index=TP1_dnb, columns=['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12'])
l_dnb_score = init_df.fillna(0)

sample_num_lst = [3, 2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 3]  # 每个stage的样本个数
for i in range(len(sample_num_lst)):
    print('当前stage为：', i+1)
    temp_stage_score = stage_dnb_score(i+1, sample_num_lst[i])
    in_set = TP1_dnb & set(temp_stage_score.index)
    l_dnb_score.loc[in_set, 'T'+str(i+1)] = temp_stage_score.loc[in_set, 'mean_score']


x = np.arange(len(sample_num_lst))
y = np.arange(len(TP1_dnb))
x, y = np.meshgrid(x, y)
z = np.array(l_dnb_score)

X, Y, Z = x.ravel(), y.ravel(), z.ravel()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.bar3d(X, Y, np.zeros_like(Z), dx=0.1, dy=0.4, dz=Z, color='blue')
# ax.bar3d(np.ones(len(dnb)), np.arange(len(dnb)), np.zeros_like(dnb_Z), dx=0.2, dy=0.4, dz=dnb_Z, color='red')
# ax.set_zlim(0, 30)

ax.set_xlabel('Stage')
ax.set_ylabel('DNB gene')
ax.set_zlabel('Score')
plt.xticks(np.arange(12), l_score['stage'], rotation=10)
plt.title('DNB module score bar plot for UVB(TP1)', fontsize=40)
plt.show()


# 彩图
fig = plt.figure()  # 定义图像窗口
ax = Axes3D(fig)   # 在窗口上添加3D坐标轴
# 将colormap ranbow填充颜色，之后将三维图像投影到XY平面做等高线图，其中rstride和cstride表示row和column的宽度
surf = ax.plot_surface(x, y, z, cmap='cool')  # rstride表示图像中分割线的跨图
ax.set_zlim(0, 20)
ax.set_xlabel('Stage')
ax.set_ylabel('DNB gene')
ax.set_zlabel('Score')
plt.xticks(np.arange(12), l_score['stage'], rotation=10)
plt.title('DNB module score surface plot for UVB(TP1)', fontsize=40)
fig.colorbar(surf, shrink=0.6)
plt.show()

# ********************************************************************************************************************
# 计算每个时间解阶段的dnb分值(Tipping point 2)
init_df = pd.DataFrame(index=TP2_dnb, columns=['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12'])
l_dnb_score = init_df.fillna(0)

sample_num_lst = [3, 2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 3]  # 每个stage的样本个数
for i in range(len(sample_num_lst)):
    print('当前stage为：', i+1)
    temp_stage_score = stage_dnb_score(i+1, sample_num_lst[i])
    in_set = TP2_dnb & set(temp_stage_score.index)
    l_dnb_score.loc[in_set, 'T'+str(i+1)] = temp_stage_score.loc[in_set, 'mean_score']


x = np.arange(len(sample_num_lst))
y = np.arange(len(TP2_dnb))
x, y = np.meshgrid(x, y)
z = np.array(l_dnb_score)

X, Y, Z = x.ravel(), y.ravel(), z.ravel()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.bar3d(X, Y, np.zeros_like(Z), dx=0.1, dy=0.4, dz=Z, color='blue')
# ax.bar3d(np.ones(len(dnb)), np.arange(len(dnb)), np.zeros_like(dnb_Z), dx=0.2, dy=0.4, dz=dnb_Z, color='red')
# ax.set_zlim(0, 30)

ax.set_xlabel('Stage')
ax.set_ylabel('DNB gene')
ax.set_zlabel('Score')
plt.xticks(np.arange(12), l_score['stage'], rotation=10)
plt.title('DNB module score bar plot for UVB(TP2)', fontsize=40)
plt.show()


# 彩图
fig = plt.figure()  # 定义图像窗口
ax = Axes3D(fig)   # 在窗口上添加3D坐标轴
# 将colormap ranbow填充颜色，之后将三维图像投影到XY平面做等高线图，其中rstride和cstride表示row和column的宽度
surf = ax.plot_surface(x, y, z, cmap='cool')  # rstride表示图像中分割线的跨图
ax.set_zlim(0, 20)
ax.set_xlabel('Stage')
ax.set_ylabel('DNB gene')
ax.set_zlabel('Score')
plt.xticks(np.arange(12), l_score['stage'], rotation=10)
plt.title('DNB module score surface plot for UVB(TP2)', fontsize=40)
fig.colorbar(surf, shrink=0.6)
plt.show()


