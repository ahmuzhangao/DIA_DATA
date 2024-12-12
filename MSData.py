# -*- coding: utf-8 -*-


#  主要放和生物相关的一些文件
class CINI:

    # 这几个东东，只要爱因斯坦不被批斗，不太可能会变
    MASS_ELECTRON = 0.0005485799
    MASS_PROTON_MONO = 1.00727645224  # 1.00782503214-0.0005485799
    MASS_PROTON_ARVG = 1.0025

    DICT0_ELEMENT_MASS = {}
    DICT0_ELEMENT_ABDC = {}

    DICT1_AA_COM = {}
    DICT2_MOD_COM = {}


class Config:

    I0_INI_PATH_ELEMENT = 'element.ini'
    I1_INI_PATH_AA = 'aa.ini'
    I2_INI_PATH_MOD = 'modification.ini'

    A0_FLAG_USE_EXIST_INDEX = 1
    A1_PATH_MS1 = ''  # 里面有多个文件
    A2_TYPE_MS1 = 0  # 这里之所有不搞为string，是因为ms1的这个后缀可能有大小写等等各种变形；
    A3_PATH_MS2 = ''
    A4_TYPE_MS2 = 0

    B1_PATH_IDENTIFICATION_RESULT = ''
    B2_TYPE_IDENTIFICATION_RESULT = 0   # 0:DIA-NN library 1:pFind 2:Spectronaut 3:pQuant
    B3_THRESHOLD_FDR = 0.01
    B4_PATH_LIBRARY_RESULT = r'C:\Users\Lenovo\Desktop\DIA\Olddraw\draw_report-lib.tsv'


    C0_TYPE_FLOW = 2  # 0: DIA色谱曲线，1:DIA色谱曲线并对每个结果单独打分 2:DDA母离子色谱曲线 3:DIA-NN重排序
    C1_DIA_RT_HALF_WIN_IN_MIN = 2  # RT窗口
    C2_DIA_ACCURACY_HALF_WIN_PEAK = 20
    C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK = 0  # 0: ppm ，1: Da

    C4_DIA_PRECURSOR_SEQUENCE = ''
    C5_DIA_PRECURSOR_CHARGE = ''
    C6_DIA_PRECURSOR_RAW = ''
    C7_PLOT_NUM = 'all'  # 'all':画鉴定结果文件所有图，'custom':用户自己输入画图序列到config
    C8_DDA_PRECURSOR_SCAN = ''
    C9_DDA_PRECURSOR_RAW = ''
    C10_DDA_ISOTOPE = ''
    C11_DDA_PRECURSOR_SEQUENCE = ''
    C12_DDA_PRECURSOR_CHARGE = ''
    C13_DDA_RT = ''
    C14_CURVE_SMOOTH = 0
    C15_INT_NUM = -1
    C16_TYPE_RERANK_MODEL = 3
    C17_DRAW_FLAG = 1
    C18_FRAG_NUM = 5 # 0代表提出所有理论碎片离子，1代表提强度最高，2代表提强度前二，依次类推

    D1_PATH_EXPORT = 'D:\\test_score\\'


class CFileMS2:  # 注意：这是按列存储，每个属性都是list。如果按列搞，这个行是一个对象，不太好管理。

    INDEX_SCAN = []   # 这个大小和大家不一样，方便快速索引
    INDEX_RT = []

    LIST_RET_TIME = []  # 2019.08.31，MS1的改了，MS2的应该也要改。为DIA做准备。
    LIST_ION_INJECTION_TIME = []
    LIST_ACTIVATION_CENTER = []
    LIST_PRECURSOR_SCAN = []

    MATRIX_PEAK_MOZ = []  # 每一行是个list
    MATRIX_PEAK_INT = []


class CFileMS1:

    INDEX_SCAN = []   # 这两个大小和大家不一样，方便快速索引
    INDEX_RT = []

    LIST_ION_INJECTION_TIME = []

    MATRIX_PEAK_MOZ = []  # 每一行是个list
    MATRIX_PEAK_INT = []


class CSeedDIACheck:

    MID_SCAN = 0
    MID_RT = 0

    MOZ_FRAGMENT = []
    TYPE_FRAGMENT = []
    MOZ_PRECURSOR = []
    MOZ_LIBRARY = []
    TYPE_LIBRARY = []


class CEvidenceDIACheck:

    MID_IDNEX = -1

    LIST_RET_TIME = []
    LIST_SCAN = []

    MATRIX_MS1_PEAK_MOZ = []
    MATRIX_MS1_PEAK_INT = []
    MATRIX_MS2_PEAK_MOZ = []
    MATRIX_MS2_PEAK_INT = []

    RT_START = -1
    RT_END = -1


class CEvidenceRerank:

    MID_IDNEX = -1
    MATRIX_MS1_PEAK_MOZ = []
    MATRIX_MS1_PEAK_INT = []
    MATRIX_MS2_PEAK_MOZ = []
    MATRIX_MS2_PEAK_INT = []


class CFeature:

    FEATURE_MATRIX = []

    FEATURE_MATRIX_NONORM = []

    FEATURE_MATRIX_PROFILE = []

    LABEL_LIST = []

    ID_INDEX_LIST = []


class CSeed:  # 为得到Evidence准备的

    MID_SCAN = 0  # 根据需求不同，使用scan或rt
    MID_RT = 0

    MOZ_CLC = 0.0


class CSeed_FLOW1:  # 为得到Evidence准备的

    MID_SCAN = 0  # 根据需求不同，使用scan或rt
    MID_RT = 0
    MOZ_PREC = []
    MOZ_FRAG = []
    PREC_TYPE = []
    FRAG_TYPE = []


class CEvidence_FLOW1:

    MATRIX_PROFILE_PREC = []
    MATRIX_PROFILE_FRAG = []

    MATRIX_MASS_DEV_PREC = []
    MATRIX_MASS_DEV_FRAG = []

    LIST_RET_TIME = []   # 曲线上每个点的保留时间
    LIST_SCAN = []  # 曲线上每个点的scan

    FRAG_TYPE = []
    FRAG_MOZ = []
    PREC_TYPE = []
    PREC_MOZ = []

    RT_START = 0  # 色谱曲线的起始点和终止点(存放在rt列表的位置)
    RT_END = 0
    POINT_APEX = 0


class CEvidence:

    # 这是第二批要填的信息
    MATRIX_PROFILE = []  # 一行曲线，强度变化

    # 这是第三批要填的信息
    LIST_RET_TIME = []   # 曲线上每个点的保留时间
    LIST_SCAN = []  # 曲线上每个点的scan

    I_START = 0  # 最后确定的起止。
    I_END = -1


class CSeed2:

    MID_SCAN = 0  # 根据需求不同，使用scan或rt
    MID_RT = 0
    START_RT = 0
    END_RT = 0
    MOZ_PREC = 0
    MOZ_CLC = []
    ION_TYPE = []


class CEvidence2:

    MATRIX_PROFILE = []  # 每一行曲线表示一个子离子的强度变化
    MATRIX_MASS_DEV = []
    LIST_ION_TYPE = []
    MOZ_CLC = []

    LIST_RET_TIME = []   # 曲线上每个点的保留时间
    LIST_SCAN = []  # 曲线上每个点的scan

    RT_START = 0  # 色谱曲线的起始点和终止点(存放在rt列表的位置)
    RT_END = 0
    POINT_APEX = 0


class CFileIDForCheck:

    N_ID = 0

    ID0_RAW_NAME = []
    ID1_SEQ = []
    ID2_SCAN_NO = []
    ID3_RT = []

    ID4_MOD = []
    ID5_PRECURSOR_MOZ_CLC = []
    ID6_CHARGE = []
    ID7_PRECURSOR_ID = []

    ID8_SCORE0 = []  # q value
    ID9_SCORE1 = []  # first DIA-NN score before Rerank
    ID10_TARGET = []
    ID11_RTSTART = []
    ID12_RTEND = []
    ID13_FRAGMENT_TYPE = []  # 每个元素是一个列表，存放library中该母离子的所有子离子类型：如b2++,y5+
    ID14_FRAGMENT_LOSS_TYPE = []  # 每个元素是一个列表，存放该母离子下所有子离子的loss type
    ID15_FRAGMENT_MOZ = []  # 每个元素是一个列表，存放library中该母离子的所有子离子moz


class CFileIDForDIAPhenomenon:

    N_ID = 0

    ID0_RAW_NAME = []
    ID1_SEQ = []
    ID2_SCAN_NO = []
    ID3_RT = []
    ID4_SEQ_WITH_MOD = []
    ID5_MOD = []
    ID6_CHARGE = []
    ID7_RPECURSOR_MOZ_CAL = []

    ID8_SCORE0 = []


class CSeedForDDA:  # 为得到Evidence准备的

    MID_SCAN = 0  # 根据需求不同，使用scan或rt
    MID_RT = 0
    RAW_NAME = 0.0
    PRECURSOR_SQE = 0.0
    PRECURSOR_MOD = 0.0
    PRECURSOR_CHARGE = 0.0
    MOZ_CLC = 0.0


class CFileID:  # 注意：这是按列存储，每个属性都是list。如果按行搞，这个行是一个对象，不太好管理。

    N_ID = 0

    ID0_SEQ = []
    ID1_RAW_NAME = []
    ID2_SCAN_ID = []
    ID3_RT = []

    ID4_SEQ_WITH_MOD = []

    ID5_PRECURSOR_MOZ_EXP = []
    ID6_PRECURSOR_MOZ_CLC = []
    ID7_CHARGE = []
    ID8_PRECURSOR_ID = []  # 存放母离子唯一ID（raw文件+序列+修饰+电荷），方便快速索引

    ID9_SCORE0 = []  # 专门放FDR或q-value

    # 存放by离子有关数据
    ID10_FRAGMENT_TYPE = []  # 每个元素是一个列表，存放library中该母离子的所有子离子类型：如b2++,y5+
    ID11_FRAGMENT_MOZ_EXP = []  # 每个元素是一个列表
    ID12_FRAGMENT_MOZ_CLC = []  # 每个元素是一个列表
    ID13_FRAGMENT_LOSS_TYPE = []  # 每个元素是一个列表，存放该母离子下所有子离子的loss type
    # 起始时间和终止时间
    ID14_RT_BEGIN = []
    ID15_RT_END = []

    ID16_MOD = []

    ID17_TARGET = []


class CFileDDAID:  # 注意：这是按列存储，每个属性都是list。如果按行搞，这个行是一个对象，不太好管理。

    N_PSM = 0

    PSM1_RAW_NAME = []  # 这个可以在鉴定结果中读到，但是，没有路径
    PSM2_SCAN_ID = []
    PSM3_RT = []

    PSM4_SEQ = []
    PSM5_MOD = []  # 每一行都是一个形如6,Carbamidomethyl[C];14,Carbamidomethyl[C];的字符串

    PSM6_PRECURSOR_ID = []  # 用于判断config中的结果是否在ID中，'raw_name + ScanNo'

    PSM7_MH_EXP = []
    PSM8_MH_CLC = []
    PSM9_CHARGE = []

    PSM10_SCORE0 = []  # 专门放FDR
    PSM11_SCORE1 = []  # 各引擎的打分
    PSM12_SCORE2 = []

    PSM13_RT_START = []
    PSM14_RT_END = []


class CResult:
    N_RESULT = 0

    RE0_RAW_NAME = []
    RE1_SEQ = []
    RE2_SCAN_NO = []
    RE3_RT = []

    RE4_MOD = []
    RE5_PRECURSOR_MOZ_CLC = []
    RE6_CHARGE = []
    RE7_PRECURSOR_ID = []

    RE8_SCORE0 = []  # finall score
    RE9_SCORE1 = []  # first DIA-NN score before Rerank
    RE10_SCORE_RERANK = []
    RE11_TARGET = []
    RE12_FEATURE = []


class CDataPack:  # 这个类必须放到最后

    myCFG = Config()
    myINI = CINI()

    #  从config里面搞出来的
    LIST_PATH_MS1 = []
    LIST_PATH_MS2 = []

    LIST_PATH_ID = []

    #  需要全周期维护的
    myID = CFileID()  # 从鉴定结果到seed，seed到evidence，这两个就当临时变量了。如果需要输出，直接写到文件里，不全周期维护了。
    myDDAID = CFileDDAID()
    myIDForCheck = CFileIDForCheck()
    myFeature = CFeature()
    myResult = CResult()
