# Construct Reference Samples & Stage samples
import numpy as np
import pandas as pd

df = pd.read_excel('RNA-seq_data_FPKM_symbol.xlsx', header=0, index_col=0).T
# df = df[(df != 0).astype(int).sum(axis=1) > 57*0.1]
df = (df - df.mean()) / df.std()  # Z-score标准化
df = df.T
print(df.head())

print(df.head())

Reference_samples = df.loc[:, ['day-1-con-1', 'day-1-con-2', 'day-1-con-3',
                               'day-4-con-1', 'day-4-con-2', 'day-4-con-3',
                               'day-6-con-1', 'day-6-con-2', 'day-6-con-3',
                               'day-8-con-1', 'day-8-con-2', 'day-8-con-3']]
Reference_samples.to_csv('Reference_samples.csv')

UVB_T1 = df.loc[:, ['day-1-UVB-1', 'day-1-UVB-2', 'day-1-UVB-3']]
UVB_T1.to_csv('UVB_T1.csv')

UVB_T2 = df.loc[:, ['day-1-UVB-4', 'day-1-UVB-5']]
UVB_T2.to_csv('UVB_T2.csv')

UVB_T3 = df.loc[:, ['day-1-UVB-6', 'day-1-UVB-7']]
UVB_T3.to_csv('UVB_T3.csv')

UVB_T4 = df.loc[:, ['day-2-UVB-1', 'day-2-UVB-2']]
UVB_T4.to_csv('UVB_T4.csv')

UVB_T5 = df.loc[:, ['day-3-UVB-1', 'day-3-UVB-2']]
UVB_T5.to_csv('UVB_T5.csv')

UVB_T6 = df.loc[:, ['day-4-UVB-1', 'day-4-UVB-2', 'day-4-UVB-3']]
UVB_T6.to_csv('UVB_T6.csv')

UVB_T7 = df.loc[:, ['day-5-UVB-1', 'day-5-UVB-2']]
UVB_T7.to_csv('UVB_T7.csv')

UVB_T8 = df.loc[:, ['day-6-UVB-1', 'day-6-UVB-2', 'day-6-UVB-3']]
UVB_T8.to_csv('UVB_T8.csv')

UVB_T9 = df.loc[:, ['day-6-UVB-4', 'day-6-UVB-5']]
UVB_T9.to_csv('UVB_T9.csv')

UVB_T10 = df.loc[:, ['day-6-UVB-6', 'day-6-UVB-7']]
UVB_T10.to_csv('UVB_T10.csv')

UVB_T11 = df.loc[:, ['day-7-UVB-1', 'day-7-UVB-2']]
UVB_T11.to_csv('UVB_T11.csv')

UVB_T12 = df.loc[:, ['day-8-UVB-1', 'day-8-UVB-2', 'day-8-UVB-3']]
UVB_T12.to_csv('UVB_T12.csv')


















