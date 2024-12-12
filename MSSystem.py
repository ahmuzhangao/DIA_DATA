# -*- coding: utf-8 -*-


VALUE_MAX_SCAN = 2000000
VALUE_APP_SCAN = 100000
VALUE_ILLEGAL = -7.16  # 搞成整数有实际意义，不行

COUNT_MARK = 200

MARK_LABEL_INFO = ('NONE', 'AA:', 'MOD:')  # 判断的时候忽略大小写

IO_NAME_FILE_EXPORT_Flow1 = ('Precursor_Curve.png', 'Precursor_Curve.svg', 'Precursor_Curve3D.png', 'Precursor_Curve3D.svg')

IO_NAME_FILE_EXPORT_Flow2 = ('Fragment_Curve.png', 'UltraVisual_Fragment_Curve.svg')

IO_NAME_FILE_EXPORT_Flow3 = ('DDA_Precursor_Curve.png', )

IO_NAME_FILE_EXPORT_Flow5 = ('DIA_Check.png', 'DIA_Check.svg')

IO_NAME_FILE_EXPORT_Flow6 = ('DIA_MS1_Signal.txt', 'DIA_MS2_Signal.txt')

IO_NAME_FILE_CONFIG = ('DIA_XLC_cfg_DIANN.txt',)

IO_NAME_FILE_EXPORT_XIC_MATRIX = ('out_XIC.txt')

CFG_TYPE_MS1 = {'MS1': 0}  # 可能有其他的，可能性不大

CFG_TYPE_MS2 = {'MS2': 0, 'MGF': 1}

CFG_TYPE_IDENTIFICATION_RESULT = {'DIA-NN':1,'MaxDIA': 5,'MSFragger':7,'MSFragger-pin':8}  # 在CTaskReadID中用到，判断鉴定结果文件类型，用对应的func进行读取

UNIMOID_TO_STANDARD_MOD = {'UniMod:4': "Carbamidomethyl[C]", 'UniMod:35': "Oxidation[AnyC-termG]",
                           'UniMod:1': "Acetyl[AnyN-term]",
                           'UniMod:385': "Ammonia-loss[AnyN-termC]", 'UniMod:27': "Glu->pyro-Glu[AnyN-termE]",
                           'UniMod:28': "Ammonia-loss[AnyN-termC]",
                           '57.0215': "Carbamidomethyl[C]", '15.9949': "Oxidation[M]",
                           '42.0106': "Acetyl[AnyN-term]"}

CFG_TYPE_VISUAL = {'DIA_XIC': 0, 'Fragment': 1, 'DDA_Precursor': 2, 'Rerank': 3, 'DIA_Check': 4, 'DIA_Phenomenon': 5}  # 在Staff中用到，判断进入不同的Flow

CFG_TYPE_ACCURACY_HALF_WIN_PEAK = {'PPM': 0, 'Da': 1}

CFG_TYPE_PLOT_NUM = {'ALL': 'all', 'CUSTOM': 'custom'}

CFG_TYPE_RERANK_MODEL = {'SVM_SemiSupervise': 0, 'DNN_SemiSupervise': 1, 'DNN_Supervise': 2, 'Resnet_SemiSupervise': 3, 'Merge_SemiSupervise': 4, 'DNN_Offline': -1}


FRAGMENT_LOSS_TYPE_COMP = {'NOLOSS': '', 'H2O': 'H(2)O(1)', 'NH3': 'H(3)N(1)', 'CO': 'C(1)O(1)'}  # 统一所有library结果的by离子loss格式

PLOT_COLOR = ('#00AEFF', '#FF8618', '#CEDB29', '#00A652', '#FF0000', '#8A2BE2', '#EECBAD', '#E0EEE0', '#E0FFFF')

DRAW_COLOR = ('#0000FF', '#5353FF', '#8282FF', '#B3B3FF', '#DADAFF', '#F2F7FF', '#FFFF00', '#B8860B', '#CD5C5C',
              '#F5DEB3', '#B22222', '#FA8072', '#FF8C00', '#FF69B4', '#B03060', '#FF00FF', '#A020F0', '#8B8989', '#FFEFDB',
              '#0000FF', '#5353FF', '#8282FF', '#B3B3FF', '#DADAFF', '#F2F7FF', '#FFFF00', '#B8860B', '#CD5C5C',
              '#F5DEB3', '#B22222', '#FA8072', '#FF8C00', '#FF69B4', '#B03060', '#FF00FF', '#A020F0', '#8B8989', '#FFEFDB',
              '#0000FF', '#5353FF', '#8282FF', '#B3B3FF', '#DADAFF', '#F2F7FF', '#FFFF00', '#B8860B', '#CD5C5C',
              '#F5DEB3', '#B22222', '#FA8072', '#FF8C00', '#FF69B4', '#B03060', '#FF00FF', '#A020F0', '#8B8989', '#FFEFDB',
              )
