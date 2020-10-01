import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_excel("L error bar.xlsx", header=0, index_col=0)
x = df.index
y = df['mean']
sd = df['SD']
# plt.style.use('ggplot')
plt.figure(figsize=(9, 6))
plt.errorbar(x, y, yerr=sd, fmt='o:', ecolor='black', color='red', markersize=10, mfc='black', mec='black', elinewidth=2, capsize=4)
# fmt: 'o' ',' '.' 'x' '+' 'v' '^' '<' '>' 's' 'd' 'p'
# plt.xticks(rotation=45)
plt.xlabel("Time", fontsize=12)
plt.ylabel("L value", fontsize=12)
# plt.tight_layout()
plt.show()
#
# plt.figure(figsize=(9, 6))
# df = pd.read_excel("L error bar.xlsx", header=0, index_col=0)
# x = df['time']
# y = df['mean']
# sd = df['SD']
# plt.errorbar(x, y, yerr=sd, fmt='o:', ecolor='black', color='red', elinewidth=2, capsize=4)
# # fmt: 'o' ',' '.' 'x' '+' 'v' '^' '<' '>' 's' 'd' 'p'
# # plt.xticks(x, df.index, rotation=45)
# plt.xlabel("Time(hour)", fontsize=12)
# plt.ylabel("L value", fontsize=12)
# plt.show()
#
# plt.errorbar(x, y, yerr=sd, fmt='o:', ecolor='hotpink',
#              elinewidth=3, ms=5, mfc='wheat', mec='salmon', capsize=3)
