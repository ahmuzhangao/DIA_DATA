# import seaborn as sns
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.colors import LinearSegmentedColormap
#
# # 创建数据矩阵
# data = np.array([
#     [-0.186942273, 0.077272783, 0.164449749, -0.331603677],  # Emotional Valence
#     [0.229177721, 0.373034034, -0.167812917, -0.199114027]   # Emotional Arousal
# ])
#
# # 行和列的标签
# row_labels = ['Emotional Valence', 'Emotional Arousal']
# col_labels = ['Number', 'Rate', 'Duration', 'Amplitude']
# colors = ['#D93F49', '#E28187', '#EBBFC2', '#D5E1E3', '#AFC9CF', '#8FB4BE']
# colors.reverse()
# cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
# # 创建 Pandas DataFrame
# # 绘制热力图
# plt.figure(figsize=(10, 6))
# sns.heatmap(data, annot=True, xticklabels=col_labels, yticklabels=row_labels, cmap=cmap, center=0)
#
# # 添加标题
# plt.title('The Correlation Analysis Between Shape-changing Interface Design Elements and Emotion(Valence and Arousal)')
# #
# # 显示图像
# plt.savefig('1.png')
# plt.savefig('1.svg',format='svg')
# plt.show()
#
# import seaborn as sns
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# from matplotlib.colors import LinearSegmentedColormap
#
# # 创建数据框
# data = {
#     'Number': [-0.068760781, -0.062111801, -0.026527905, -0.045750901, -0.270249251, 0.059949065, 0.191837008, 0.062111801, 0.215188342, 0.16936196, -0.042561756, -0.050535836],
#     'Rate': [0.320618551, 0.103302383, 0.189087314, -0.163053111, -0.445902784, -0.054125782, 0.014243627, -0.109205377, 0.124485024, 0.025149773, -0.064045712, -0.138081492],
#     'Duration': [-0.099419482, 0.004954808, 0.049377878, 0.200731264, 0.221572776, -0.081298777, -0.119557025, -0.054502884, -0.223905538, -0.059108046, 0.107516362, 0.151176191],
#     'Amplitude': [-0.234847983, -0.038369457, -0.315145524, -0.076091452, 0.181928336, 0.173772249, 0.105402839, 0.174138303, 0.195619323, -0.015089864, -0.064045712, -0.042024802]
# }
#
# emotions = ['Excitement', 'Happiness', 'Joy', 'Satisfaction', 'Calmness', 'Frustration', 'Sadness', 'Melancholy', 'Anger', 'Fear', 'Others', 'Unclear']
# colors = ['#D93F49', '#E28187', '#EBBFC2', '#D5E1E3', '#AFC9CF', '#8FB4BE']
# colors.reverse()
# cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
# # 创建 Pandas DataFrame
# df = pd.DataFrame(data, index=emotions)
#
# # 绘制热力图
# plt.figure(figsize=(10, 8))
# sns.heatmap(df, annot=True, cmap=cmap, center=0, linewidths=.5)
#
# # 添加标题
# plt.title('The Correlation Analysis Between Shape-changing Interface Design Elements and Emotion(Type)')
#
# plt.savefig('2.png')
# plt.savefig('2.svg',format='svg')
# # 显示图像
# plt.show()
#
# import seaborn as sns
# import matplotlib.pyplot as plt
# import pandas as pd
# from matplotlib.colors import LinearSegmentedColormap
# # 创建 DataFrame
# data = {
#     'Emotional Valence': [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7],
#     'Emotional Arousal': [1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7],
#     'Rating': [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 5, 2, 0, 0, 0, 0, 2, 0, 0]
# }
#
# df = pd.DataFrame(data)
# colors =['#D93F49', '#E28187', '#EBBFC2', '#D5E1E3', '#AFC9CF', '#8FB4BE']
# colors.reverse()
# cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
# # 重新排列为矩阵形式
# heatmap_data = df.pivot('Emotional Valence', 'Emotional Arousal', 'Rating')
#
# # 绘制热力图
# plt.figure(figsize=(10, 6))
# sns.heatmap(heatmap_data, annot=True, cmap=cmap, linewidths=.5)
#
# # 添加标题
# plt.title('Users Emotional Expression Habits When Using Cupid-Echo')
# plt.savefig('3.png')
# plt.savefig('3.svg',format='svg')
# # 显示图像
# plt.show()
#
# import seaborn as sns
# import matplotlib.pyplot as plt
# import pandas as pd
# from matplotlib.colors import LinearSegmentedColormap
#
# # 创建数据
# data = {
#     'Poke': [0.263224017, -0.006179439],
#     'Pinch': [-0.411328191, 0.182925465],
#     'Pat': [0.077455028, 0.076675156],
#     'Hammer': [-0.257676294, 0.253389112],
#     'Stroke': [0.306208373, -0.274170386],
#     'Other': [-0.048041231, -0.168157998]
# }
#
# index = ['Emotional Valence', 'Emotional Arousal']
#
# # 创建 DataFrame
# df = pd.DataFrame(data, index=index)
# colors =['#D93F49','#E28187','#EBBFC2','#D5E1E3','#AFC9CF','#8FB4BE']
# colors.reverse()
# cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
# # 使用新的颜色配色 (YlGnBu)
# plt.figure(figsize=(10, 6))
# sns.heatmap(df, annot=True, cmap=cmap, linewidths=.5, center=0)
#
# # 添加标题
# plt.title('Correlation Analysis of Emotional Stimulus Effects and Behavior')
# plt.savefig('4.png')
# plt.savefig('4.svg',format='svg')
# # 显示图像
# plt.show()
#
# import matplotlib.pyplot as plt
# import pandas as pd
#
# # 创建数据
# data = {
#     'Material': ['Material 1', 'Material 2', 'Material 3', 'Material 4', 'Material 5', 'Material 6', 'Material 7', 'Material 8'],
#     'Valence': [3.478, 4.783, 3.565, 2.609, 4.348, 1.913, 5.304, 4.391],
#     'Arousal': [3.913, 5.652, 5.304, 5.739, 5.087, 5.652, 5.826, 5.565]
# }
#
# # 创建 DataFrame
# df = pd.DataFrame(data)
#
# # 绘制散点图
# plt.figure(figsize=(8, 6))
# plt.scatter(df['Valence'], df['Arousal'], color='#1f77b4', s=10)  # 选择好看的颜色
# for i, txt in enumerate(df['Material']):
#     plt.annotate(txt, (df['Valence'][i], df['Arousal'][i]), fontsize=9, ha='left', va='bottom', color='black')
#
#
# # 设置轴标签和标题
# plt.xlim(2,7)
# plt.xlabel('Valence')
# plt.ylabel('Arousal')
# plt.grid(True, which='both', linestyle='--', linewidth=0.5)  # 增加网格线的密度
# plt.title('Participants perceived emotions for the material are displayed in the emotion axis')
#
# plt.savefig('5.png')
# plt.savefig('5.svg',format='svg')
# plt.show()

import matplotlib.pyplot as plt
plt.rc('font',family='Times New Roman')
import numpy as np

# 创建数据
labels = ['Happy', 'Pleasant', 'Sad', 'Angry']
precisionrecall = [1, 0.5, 1, 1]
recall = [0.333, 1, 1, 1]
f1_score = [0.5, 0.667, 1, 1]

# 设置柱状图的宽度
width = 0.1

# 生成位置数组
x = np.arange(len(labels))

# 设置字体大小
fontsize = 30

# 设置图形大小
plt.figure(figsize=(20, 10))  # 设置图形大小
colors = ['#e38241', '#f2ca6a', '#edbbb9', '#f0dfdd']
# 绘制柱状图
bar1 = plt.bar(x - width, precisionrecall, width, label='PrecisionRecall',color=colors[0])
bar2 = plt.bar(x, recall, width, label='Recall',color=colors[1])
bar3 = plt.bar(x + width, f1_score, width, label='F1-Score',color=colors[2])

# 添加标签和标题
plt.xlabel('Emotion', fontsize=fontsize)
plt.ylabel('Scores', fontsize=fontsize)
plt.title('F1 Score, Recall, and Accuracy for each Emotion', fontsize=fontsize)
plt.xticks(x, labels, fontsize=fontsize)
plt.yticks(fontsize=fontsize)
# 移除上框线和右框线
ax = plt.gca()  # 获取当前的子图
ax.spines['top'].set_visible(False)  # 移除上框线
ax.spines['right'].set_visible(False)  # 移除右框线

# 添加网格线
plt.grid(True, axis='y', linestyle='--', alpha=1, linewidth=1.2)

# 将图例移到右框线外
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=fontsize)


# 设置紧凑布局
plt.tight_layout()
plt.savefig('6.png')
plt.savefig('6.svg',format='svg')
# 显示图形
plt.show()
