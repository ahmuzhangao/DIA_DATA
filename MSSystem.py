import os
import sys

VALUE_ILLEGAL = -7.16  # 搞成整数有实际意义，不行
COUNT_MARK = 500
MASS_PROTON_MONO = 1.00727645224

def defineSoftwareName(inputCMD):
    software_name = ".".join(os.path.basename(inputCMD).split(".")[:-1])

    return software_name

try:
    SOFTWARE_NAME = defineSoftwareName(sys.argv[0])
except:
    SOFTWARE_NAME = "ApuQuant"

IO_NAME_FILE_CONFIG = ('{}_cfg.txt'.format(SOFTWARE_NAME),)

IO_NAME_FILE_EXPORT_Flow1 = ('{}.spectra.list'.format(SOFTWARE_NAME),
					   '{}.protein.list'.format(SOFTWARE_NAME),
                       '{}.modification.list'.format(SOFTWARE_NAME),
                        '{}.precursor.list'.format(SOFTWARE_NAME))

IO_NAME_FILE_EXPORT_Flow2 = ('DIA_Check.png', 'DIA_Check.svg')

IO_NAME_FILE_EXPORT_Flow3 = ('positive_paired_samples.pkl','negative_paired_samples.pkl')

CFG_TYPE_MS1 = {'MS1': 0, 'RAW': 1, 'PKL': 2, 'd': 3}  # 可能有其他的，可能性不大

CFG_TYPE_MS2 = {'MS2': 0, 'RAW': 1, 'PKL': 2, 'd': 3}

CFG_TYPE_RAW = {'RAW': 0}

CFG_TYPE_FLOW = {'DIA_Label_Free': 0, 'DIA_Evidence_Check': 1, 'DIA_Rerank':2}  # 在Staff中用到，判断进入不同的Flow

CFG_TYPE_IDENTIFICATION_RESULT = {'DIA-NN': 0, 'ApuQuant': 1, 'Spectronaut': 2, 'MaxDIA': 3,'DIA-NN_Lib':4,'MSFragger':5,'MSFragger-pin':6}  # 在CTaskReadID中用到，判断鉴定结果文件类型，用对应的func进行读取

CFG_TYPE_PREDICTMS2_MODEL = {'peptdeep': 0}  # 在CTaskReadID中用到，判断鉴定结果文件类型，用对应的func进行读取

UNIMOID_TO_STANDARD_MOD = {'UniMod:4': "Carbamidomethyl[C]", 'UniMod:35': "Oxidation[AnyC-termG]", 'UniMod:1': "Acetyl[AnyN-term]",
                           'UniMod:385': "Ammonia-loss[AnyN-termC]", 'UniMod:27': "Glu->pyro-Glu[AnyN-termE]",
                           'UniMod:28': "Ammonia-loss[AnyN-termC]", '57.0215': "Carbamidomethyl[C]", '15.9949': "Oxidation[M]",
                           '42.0106': "Acetyl[AnyN-term]"}

UNIMOID_TO_STANDARD_MOD2 = {"Carbamidomethyl[C]":'Carbamidomethyl@C' , "Oxidation[AnyC-termG]":'Oxidation@AnyC-termG', "Acetyl[AnyN-term]":"Acetyl@AnyN-term" ,
                           "Ammonia-loss[AnyN-termC]":"Ammonia-loss@AnyN-termC" , "Glu->pyro-Glu[AnyN-termE]": "Glu->pyro-Glu@AnyN-termE",
                            "Oxidation[M]":"Oxidation@M"}

MARK_FLAG_INFER = {'Yes': 0, 'No_WithSameMass': 1}

CFG_TYPE_EXPORT = {'Simple': 0, 'AppendFromID': 1, 'Test': -1}

FRAGMENT_LOSS_TYPE_COMP = {'NOLOSS': '', 'H2O': 'H(2)O(1)', 'NH3': 'H(3)N(1)', 'CO': 'C(1)O(1)'}  # 统一所有library结果的by离子loss格式

NUM_CALIB_CHECK_POINT = 135

TIME_EXPIRATION = {'Year': 2124, 'Month': 12, 'Day': 31}

CFG_DDALF_TYPE_MODEL = {'Model0': 0, 'Model1': 1, 'Model2': 2}

CFG_FLAG_LIBRARY_STRATEGY_FOR_DDALF = {'None': 0, 'Spectra': 1, 'Intensity': 2, "Retention Time": 3}

MARKER_SPLIT_MOD_LIST = "$"

MARK_TYPE_PROTEIN = ('Leading', 'Sameset', 'Subset')

CFG_TYPE_ACCURACY_HALF_WIN_PEAK = {'PPM': 0, 'Da': 1}

IO_NAME_FILE_EXPORT_Flow4 = ('ModelN1.pkl', 'ModelN2.pkl')

DRAW_COLOR = ('#CBF43F', '#1FEA12', '#F15A93', '#F4FF4D', '#14F99A', '#FCA202', '#A3A42A', '#612CFD', '#8EA823', '#4032D2', '#5CB7CC', '#1ED416', '#15D1C0', '#249D3B', '#DFBC76', '#DB8138', '#63AB5C', '#BD09EE', '#6EDE5D', '#BC0B84', '#A40E04', '#2948BF', '#084402', '#86B06D', '#A0011B', '#623633', '#98D5D6', '#9C7935', '#0171C1', '#8020F4', '#F60177', '#B0D704', '#879743', '#614A71', '#BD4654', '#769934', '#E15BAD', '#A80822', '#0423D2', '#80124B', '#765ADF', '#9DC508', '#05CF27', '#6FC007', '#6C4BA8', '#13679A', '#686972', '#12432E', '#818F93', '#59552C', '#2638B5', '#0C561D', '#4C0CD6', '#63F94D', '#DA9780', '#3CC4D2', '#756898', '#1424CA', '#E0A393', '#B221CD', '#657275', '#0B2F5E', '#70254C', '#056284', '#74E3B5', '#EDAF71', '#6E90FC', '#66B41C', '#373800', '#949344', '#A4D096', '#2E5BC2', '#955DED', '#960270', '#FEB024', '#6E4093', '#C6D3E2', '#6E8591', '#49D11E', '#8F06E7', '#6E4EBC', '#8CEDA9', '#8A13E9', '#DF67AE', '#2C8F6F', '#3CDF89', '#2AEA9A', '#AAC77A', '#2259D1', '#3C855F', '#AD6BD3', '#B93932', '#7A5EC3', '#964C32', '#14ABD0', '#1016AE', '#F07E3D', '#DB9FDF', '#970CF8', '#05E4A7', '#ED3E53', '#8BFA20', '#5322CF', '#9038C5', '#2A4520', '#DD221C', '#6BBA92', '#56494A', '#5B519A', '#00CDD6', '#58C094', '#A270CA', '#B3D520', '#305F07', '#76E1A3', '#CBBC9E', '#B864DA', '#86F998', '#866687', '#99201A', '#93C9E6', '#82827A', '#27C109', '#B22787', '#9F4261', '#454F5A', '#A26692', '#DE8FB2', '#AB610E', '#06B3E6', '#1789E8', '#4CEF76', '#BEC97A', '#58052A', '#CE5B6B', '#334786', '#505B61', '#A0DBCE', '#393F61', '#B0F7BB', '#C20999', '#96D6B5', '#EA2A06', '#CC0EC9', '#82EDEC', '#D75B8B', '#18877E', '#22C080', '#371CE6', '#A9B8F3', '#893A22', '#F2B10E', '#569ED4', '#C64816', '#5202E7', '#E8A564', '#E40818', '#AED73B', '#ACF550', '#B2F2B2', '#4B10B5', '#AA614D', '#A7CD40', '#829487', '#2D7A50', '#BB6663', '#691CB6', '#EFF22F', '#1B6687', '#5D4CE7', '#BFC76A', '#61D8DA', '#86BBB8', '#6AE8A2', '#0CFF3E', '#DD9B77', '#FA6984', '#E73EA1', '#BB69D2', '#4CD2BC', '#56A9E5', '#362015', '#22D849', '#EE64FC', '#F0E82F', '#8FC088', '#F1E42A', '#E847F0', '#D5C5B3', '#089A9F', '#17444B', '#B5576D', '#AC6955', '#146F52', '#45F672', '#2B11DF', '#6A24B3', '#08AAAD', '#E9AFEB', '#CA7F93')

TOP_N = 6