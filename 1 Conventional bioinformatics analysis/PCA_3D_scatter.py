# 提取CTRL和UVB组样本进行可视化
import pandas as pd
from collections import OrderedDict
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import decomposition as skldec  # 用于主成分分析降维的包

df = pd.read_excel('RNA-seq_data_FPKM_symbol.xlsx', header=0, index_col=0)

CTRL_df = df.loc[:, ['day-1-con-1', 'day-1-con-2', 'day-1-con-3',
                     'day-4-con-1', 'day-4-con-2', 'day-4-con-3',
                     'day-6-con-1', 'day-6-con-2', 'day-6-con-3',
                     'day-8-con-1', 'day-8-con-2', 'day-8-con-3']]

UVB_df = df.loc[:, ['day-1-UVB-1', 'day-1-UVB-2', 'day-1-UVB-3', 'day-1-UVB-4', 'day-1-UVB-5', 'day-1-UVB-6',
                    'day-1-UVB-7', 'day-2-UVB-1', 'day-2-UVB-2', 'day-3-UVB-1', 'day-3-UVB-2', 'day-4-UVB-1',
                    'day-4-UVB-2', 'day-4-UVB-3', 'day-5-UVB-1', 'day-5-UVB-2', 'day-6-UVB-1', 'day-6-UVB-2',
                    'day-6-UVB-3', 'day-6-UVB-4', 'day-6-UVB-5', 'day-6-UVB-6', 'day-6-UVB-7', 'day-7-UVB-1',
                    'day-7-UVB-2', 'day-8-UVB-1', 'day-8-UVB-2', 'day-8-UVB-3']]

print('CTRL sample: ', CTRL_df.shape[1])
print('UVB sample: ', UVB_df.shape[1])

data = df.loc[:, CTRL_df.columns | UVB_df.columns].T
data = data.drop(columns=data.columns[data.std() == 0])
data = (data - data.mean()) / data.std()  # Z-score标准化

# 根据3个最大的主成分进行绘图
pca = skldec.PCA(n_components=0.95, whiten=True)    # 选择方差95%的占比
pca.fit(data)   # 主成分分析时每一行是一个输入数据
result = pca.transform(data)  # 计算结果
pc_df = pd.DataFrame(data=result[:, 0:3], index=data.index, columns=['PC1', 'PC2', 'PC3'])
print(result.shape)


fig = plt.figure(dpi=150)
ax = Axes3D(fig)
# -------------------CTRL group 3d scatter-----------------------------
# 1 73 121 169
# ['day-1-con-1', 'day-1-con-2', 'day-1-con-3']
ax.scatter(pc_df.loc[CTRL_df.columns[0:3], 'PC1'],
           pc_df.loc[CTRL_df.columns[0:3], 'PC2'],
           pc_df.loc[CTRL_df.columns[0:3], 'PC3'], c='b', s=50+1/2, label='Control')
# ['day-4-con-1', 'day-4-con-2', 'day-4-con-3']
ax.scatter(pc_df.loc[CTRL_df.columns[3:6], 'PC1'],
           pc_df.loc[CTRL_df.columns[3:6], 'PC2'],
           pc_df.loc[CTRL_df.columns[3:6], 'PC3'], c='b', s=50+73/2, label='Control')
# ['day-6-con-1', 'day-6-con-2', 'day-6-con-3']
ax.scatter(pc_df.loc[CTRL_df.columns[6:9], 'PC1'],
           pc_df.loc[CTRL_df.columns[6:9], 'PC2'],
           pc_df.loc[CTRL_df.columns[6:9], 'PC3'], c='b', s=50+121/2, label='Control')
# ['day-8-con-1', 'day-8-con-2', 'day-8-con-3']
ax.scatter(pc_df.loc[CTRL_df.columns[9:12], 'PC1'],
           pc_df.loc[CTRL_df.columns[9:12], 'PC2'],
           pc_df.loc[CTRL_df.columns[9:12], 'PC3'], c='b', s=50+169/2, label='Control')

# -------------------UVB group 3d scatter-----------------------------
# 1 4 8 25 49 73 97 121 124 128 145 169
# ['day-1-UVB-1', 'day-1-UVB-2', 'day-1-UVB-3']
ax.scatter(pc_df.loc[UVB_df.columns[0:3], 'PC1'],
           pc_df.loc[UVB_df.columns[0:3], 'PC2'],
           pc_df.loc[UVB_df.columns[0:3], 'PC3'], c='r', s=50+1/2, label='UVB')
# ['day-1-UVB-4', 'day-1-UVB-5']
ax.scatter(pc_df.loc[UVB_df.columns[3:5], 'PC1'],
           pc_df.loc[UVB_df.columns[3:5], 'PC2'],
           pc_df.loc[UVB_df.columns[3:5], 'PC3'], c='r', s=50+4/2, label='UVB')
# ['day-1-UVB-6', 'day-1-UVB-7']
ax.scatter(pc_df.loc[UVB_df.columns[5:7], 'PC1'],
           pc_df.loc[UVB_df.columns[5:7], 'PC2'],
           pc_df.loc[UVB_df.columns[5:7], 'PC3'], c='r', s=50+8/2, label='UVB')
# ['day-2-UVB-1', 'day-2-UVB-2']
ax.scatter(pc_df.loc[UVB_df.columns[7:9], 'PC1'],
           pc_df.loc[UVB_df.columns[7:9], 'PC2'],
           pc_df.loc[UVB_df.columns[7:9], 'PC3'], c='r', s=50+25/2, label='UVB')
# ['day-3-UVB-1', 'day-3-UVB-2']
ax.scatter(pc_df.loc[UVB_df.columns[9:11], 'PC1'],
           pc_df.loc[UVB_df.columns[9:11], 'PC2'],
           pc_df.loc[UVB_df.columns[9:11], 'PC3'], c='r', s=50+49/2, label='UVB')
# ['day-4-UVB-1', 'day-4-UVB-2', 'day-4-UVB-3']
ax.scatter(pc_df.loc[UVB_df.columns[11:14], 'PC1'],
           pc_df.loc[UVB_df.columns[11:14], 'PC2'],
           pc_df.loc[UVB_df.columns[11:14], 'PC3'], c='r', s=50+73/2, label='UVB')
# ['day-5-UVB-1', 'day-5-UVB-2']
ax.scatter(pc_df.loc[UVB_df.columns[14:16], 'PC1'],
           pc_df.loc[UVB_df.columns[14:16], 'PC2'],
           pc_df.loc[UVB_df.columns[14:16], 'PC3'], c='r', s=50+97/2, label='UVB')
# ['day-6-UVB-1', 'day-6-UVB-2', 'day-6-UVB-3']
ax.scatter(pc_df.loc[UVB_df.columns[16:19], 'PC1'],
           pc_df.loc[UVB_df.columns[16:19], 'PC2'],
           pc_df.loc[UVB_df.columns[16:19], 'PC3'], c='r', s=50+121/2, label='UVB')
# ['day-6-UVB-4', 'day-6-UVB-5']
ax.scatter(pc_df.loc[UVB_df.columns[19:21], 'PC1'],
           pc_df.loc[UVB_df.columns[19:21], 'PC2'],
           pc_df.loc[UVB_df.columns[19:21], 'PC3'], c='r', s=50+124/2, label='UVB')
# ['day-6-UVB-6', 'day-6-UVB-7']
ax.scatter(pc_df.loc[UVB_df.columns[21:23], 'PC1'],
           pc_df.loc[UVB_df.columns[21:23], 'PC2'],
           pc_df.loc[UVB_df.columns[21:23], 'PC3'], c='r', s=50+128/2, label='UVB')
# ['day-7-UVB-1', 'day-7-UVB-2']
ax.scatter(pc_df.loc[UVB_df.columns[23:25], 'PC1'],
           pc_df.loc[UVB_df.columns[23:25], 'PC2'],
           pc_df.loc[UVB_df.columns[23:25], 'PC3'], c='r', s=50+145/2, label='UVB')
# ['day-8-UVB-1', 'day-8-UVB-2', 'day-8-UVB-3']
ax.scatter(pc_df.loc[UVB_df.columns[25:28], 'PC1'],
           pc_df.loc[UVB_df.columns[25:28], 'PC2'],
           pc_df.loc[UVB_df.columns[25:28], 'PC3'], c='r', s=50+169/2, label='UVB')


for i in range(result[:, 0].size):
    ax.text(pc_df.iloc[i, 0]+0.02, pc_df.iloc[i, 1], pc_df.iloc[i, 2], pc_df.index[i][4], size=15)     # 在每个点边上绘制数据名称

# 设置X、Y、Z轴的名字显示，用刺眼的红色
ax.set_xlabel('PC1', fontdict={'size': 30, 'color': 'black'})
ax.set_ylabel('PC2', fontdict={'size': 30, 'color': 'black'})
ax.set_zlabel('PC3', fontdict={'size': 30, 'color': 'black'})
ax.grid(True)

ax.legend()
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())
# plt.title('PCA visualization of samples')
plt.show()


