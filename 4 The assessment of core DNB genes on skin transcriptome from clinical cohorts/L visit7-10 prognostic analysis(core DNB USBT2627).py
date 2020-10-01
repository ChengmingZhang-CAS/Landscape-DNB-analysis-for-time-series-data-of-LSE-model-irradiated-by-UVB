# Get the data for survival analysis
import pandas as pd
import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as plt
from collections import OrderedDict

# gene expression matrix
df = pd.read_excel('nomalization data Unilever with mean item.xlsx', index_col=0, header=0)

#  sample group information
group_info = pd.read_excel('groups information from Unilever to Prof Chen.xlsx', index_col=0, header=0)
# plc_idx = group_info[group_info['Group Code'] == 'A'].index
sbt_idx = group_info[group_info['Group Code'] == 'B'].index
col_id = sbt_idx


# 3D model core DNB
ESSN_node = pd.read_csv('Expanded_SSN_node_info.csv', header=0, index_col=0)
uvb_dnb = ESSN_node[(ESSN_node['node_type1'] == 'dnb-marker') & (ESSN_node['node_type2'] != 'plain')].index
# uvb_dnb = ESSN_node[(ESSN_node['node_type1'] == 'dnb-marker')].index

data = df.loc[uvb_dnb, col_id].T

# clinic data phenotype
# CFB(change from baseline)
L_df = pd.read_excel('Clinical phenotype data.xlsx', sheet_name=0, index_col=0, header=0)
# visit 7 SBT group
idx = (L_df.loc[:, 'Test_regimes'] == 'DEG') & (L_df.loc[:, 'Product'] == 'B') & (L_df.loc[:, 'Visit'] == 7)
visit7_L = L_df.loc[idx, 'CFB']

# visit 8 SBT group
idx = (L_df.loc[:, 'Test_regimes'] == 'Phenotype') & (L_df.loc[:, 'Product'] == 'B') & (L_df.loc[:, 'Visit'] == 8)
visit8_L = L_df.loc[idx, 'CFB']

# visit 9 SBT group
idx = (L_df.loc[:, 'Test_regimes'] == 'Phenotype') & (L_df.loc[:, 'Product'] == 'B') & (L_df.loc[:, 'Visit'] == 9)
visit9_L = L_df.loc[idx, 'CFB']

# visit 10 SBT group
idx = (L_df.loc[:, 'Test_regimes'] == 'Phenotype') & (L_df.loc[:, 'Product'] == 'B') & (L_df.loc[:, 'Visit'] == 10)
visit10_L = L_df.loc[idx, 'CFB']

data.insert(loc=0, column='visit10_L*', value=visit10_L.values)
data.insert(loc=0, column='visit09_L*', value=visit9_L.values)
data.insert(loc=0, column='visit08_L*', value=visit8_L.values)
data.insert(loc=0, column='visit07_L*', value=visit7_L.values)


def draw_surv_plot(prob1, prob2, gene, p_value, ratio, is_save_fig=True):

    plt.figure()
    # plt.plot([0, 0.3], [prob1[0], prob1[0]], color='r', linestyle='solid', label='high')  # 添加水平直线
    plt.plot([0.3, 1], [prob1[0], prob1[0]], linewidth=3, color='r',  label='high')  # 添加水平直线
    plt.plot([1, 3], [prob1[1], prob1[1]], linewidth=3, color='r',  label='high')  # 添加水平直线
    plt.plot([3, 10], [prob1[2], prob1[2]], linewidth=3, color='r',  label='high')  # 添加水平直线
    plt.vlines(1, prob1[0], prob1[1], linewidth=3, colors='r',label='high')  # 添加水平直线
    plt.vlines(3, prob1[1], prob1[2], linewidth=3, colors='r',  label='high')  # 添加水平直线
    plt.vlines(10, prob1[2], prob1[3], linewidth=3, colors='r', label='high')  # 添加水平直线

    plt.plot([0.3, 1], [prob2[0], prob2[0]], color='g', linewidth=3, label='low')  # 添加水平直线
    plt.plot([1, 3], [prob2[1], prob2[1]], color='g', linewidth=3, label='low')  # 添加水平直线
    plt.plot([3, 10], [prob2[2], prob2[2]], color='g', linewidth=3, label='low')  # 添加水平直线
    plt.vlines(1, prob2[0], prob2[1], colors='g', linewidth=3, label='low')  # 添加水平直线
    plt.vlines(3, prob2[1], prob2[2], colors='g', linewidth=3, label='low')  # 添加水平直线
    plt.vlines(10, prob2[2], prob2[3], colors='g', linewidth=3, label='low')  # 添加水平直线
    plt.xlim(0, 11)
    plt.ylim(0, 1)
    plt.xlabel('Days after UV irradiation')
    plt.ylabel('L* significant change proportion')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.text(1, 0.8, "p value: ${}$".format(round(p_value, 4)))
    title = gene + "(DNB) prognostic analysis"
    plt.title(title)
    # plt.tight_layout()
    fig_name = 'ratio_{}_{}_{}'.format(str(ratio), gene, str(round(p_value, 3)))
    if is_save_fig:
        plt.savefig(fig_name + '.png')
    plt.ion()
    plt.pause(3)
    plt.ioff()
    plt.close()


def draw_box_plot(data, gene, ratio, is_save_fig=True):

    data = data.sort_values(by=gene, ascending=False)
    data.insert(loc=4, column='DNB_Expression', value=0)
    data.iloc[0:9, 4] = 'high'
    data.iloc[9:, 4] = 'low'
    df = data.iloc[:, 0:5]

    # plt.style.use('ggplot')
    df.boxplot(by='DNB_Expression', patch_artist=True, notch=True, medianprops={'linewidth': 2, 'color': 'black'})
    title = gene + "(DNB) boxplot(grouped by expression)"
    plt.suptitle(title)
    # plt.tight_layout()
    fig_name = 'ratio_{}_{}_boxplot'.format(str(ratio), gene)
    if is_save_fig:
        # gcf: Get Current Figure
        fig = plt.gcf()
        fig.savefig(fig_name + '.png', dpi=100)
    plt.ion()
    plt.pause(3)
    plt.ioff()
    plt.close()


def get_pvalue(gene, data, ratio=0.5, isplot=True):
    lifespan_table = pd.DataFrame(index=['visit7', 'visit8', 'visit9', 'visit10'],
                                  columns=['obs_num1', 'prob1', 'expt_num1', 'obs_num2', 'expt_num2', 'prob2', 'obs_num'])

    gene_median = np.median(data.loc[:, gene])
    group_high = data.loc[data.loc[:, gene] >= gene_median]
    group_low = data.loc[data.loc[:, gene] < gene_median]

    # cutoff = np.quantile(data.iloc[:, 0:4], quantile)
    cutoff = np.min(np.min(data.iloc[:, 0:4])) * ratio
    n1 = 9  # sample num in each group
    n2 = 10  # sample num in each group
    n = n1 + n2
    # ---visit7---
    obs_num1 = np.size(np.where(group_high.loc[:, 'visit07_L*'] < cutoff))
    obs_num2 = np.size(np.where(group_low.loc[:, 'visit07_L*'] < cutoff))
    prob1 = obs_num1 / n1
    prob2 = obs_num2 / n2

    obs_num = obs_num1 + obs_num2
    expt_num1 = obs_num*(n1/n)
    expt_num2 = obs_num*(n2/n)
    lifespan_table.loc['visit7', :] = [obs_num1, prob1, expt_num1, obs_num2, expt_num2, prob2, obs_num]

    # ---visit8---
    obs_num1 = np.size(np.where(group_high.loc[:, 'visit08_L*'] < cutoff))
    obs_num2 = np.size(np.where(group_low.loc[:, 'visit08_L*'] < cutoff))
    prob1 = obs_num1 / n1
    prob2 = obs_num2 / n2
    obs_num = obs_num1 + obs_num2
    expt_num1 = obs_num*(n1/n)
    expt_num2 = obs_num*(n2/n)
    lifespan_table.loc['visit8', :] = [obs_num1, prob1, expt_num1, obs_num2, expt_num2, prob2, obs_num]

    # ---visit9---
    obs_num1 = np.size(np.where(group_high.loc[:, 'visit09_L*'] < cutoff))
    obs_num2 = np.size(np.where(group_low.loc[:, 'visit09_L*'] < cutoff))
    prob1 = obs_num1 / n1
    prob2 = obs_num2 / n2
    obs_num = obs_num1 + obs_num2
    expt_num1 = obs_num*(n1/n)
    expt_num2 = obs_num*(n2/n)
    lifespan_table.loc['visit9', :] = [obs_num1, prob1, expt_num1, obs_num2, expt_num2, prob2, obs_num]

    # ---visit10---
    obs_num1 = np.size(np.where(group_high.loc[:, 'visit10_L*'] < cutoff))
    obs_num2 = np.size(np.where(group_low.loc[:, 'visit10_L*'] < cutoff))
    prob1 = obs_num1 / n1
    prob2 = obs_num2 / n2
    obs_num = obs_num1 + obs_num2
    expt_num1 = obs_num*(n1/n)
    expt_num2 = obs_num*(n2/n)
    lifespan_table.loc['visit10', :] = [obs_num1, prob1, expt_num1, obs_num2, expt_num2, prob2, obs_num]

    O1 = lifespan_table.loc[:, 'obs_num1'].sum()
    O2 = lifespan_table.loc[:, 'obs_num2'].sum()
    E1 = lifespan_table.loc[:, 'expt_num1'].sum()
    E2 = lifespan_table.loc[:, 'expt_num2'].sum()

    X = np.power(O1-E1, 2)/E1 + np.power(O2-E2, 2)/E2
    p_value = 1 - chi2.cdf(X, df=1)
    print(gene, p_value)
    if p_value < 0.1 and isplot:
        prob1 = lifespan_table.loc[:, 'prob1']
        prob2 = lifespan_table.loc[:, 'prob2']
        draw_surv_plot(prob1, prob2, gene, p_value, ratio)
        draw_box_plot(data, gene, ratio)

    return p_value


ratio_list = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
# ratio_list = [0.8,0.9]

p_values = pd.DataFrame(index=uvb_dnb, columns=ratio_list)

for ratio in ratio_list:
    for gene in uvb_dnb:
        p_values.loc[gene, ratio] = get_pvalue(gene, data, ratio=ratio)

print(p_values)
p_values.to_excel('P_value L dnb-marker strict.xlsx')
