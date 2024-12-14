class CINI:
    # 这几个东东，只要爱因斯坦不被批斗，不太可能会变
    MASS_ELECTRON = 0.0005485799
    MASS_PROTON_MONO = 1.00727645224  # 1.00782503214-0.0005485799
    MASS_PROTON_ARVG = 1.0025

    DICT0_ELEMENT_MASS = {}#索引是元素名称，值是若干个质量（同位素）
    DICT0_ELEMENT_ABDC = {}#索引是元素名称，值是若干个比例（同位素）

    DICT1_AA_COM = {}#索引是氨基酸名称，值是元素组成和质量
    DICT2_MOD_COM = {}#索引是修饰名称，值是元素组成和质量


class CFileID:  # 注意：这是按列存储，每个属性都是list。如果按行搞，这个行是一个对象，不太好管理。

    N_ID = 0

    ID0_SEQ = []
    ID1_RAW_NAME = []
    ID2_SCAN_ID = []
    ID3_RT = []

    ID4_SEQ_WITH_MOD = [] # 算碎片离子质量的时候比较方便

    ID5_PRECURSOR_MOZ_EXP = []
    ID6_PRECURSOR_MOZ_CLC = []
    ID7_CHARGE = []

    ID8_PRECURSOR_ID = []  # 存放母离子唯一ID（raw文件+序列+修饰+电荷），方便快速索引

    ID9_SCORE0 = []  # 专门放FDR或q-value

    # 蛋白质相关
    ID10_Protein_Groups = []
    ID11_Protein_Ids = []
    ID12_Protein_Names = []
    ID13_Genes = []

    # 起始时间和终止时间
    ID14_RT_BEGIN = []
    ID15_RT_END = []

    ID16_MOD = []

    ID17_SEQ_WITH_MOD_AND_CHARGE = [] # 后续肽段建库需要用到

    ID18_MOD_SITE = []
    ID19_MOD_TYPE = []

    ID20_LIB_MINOR_INFO = []

    MOD1_DICT_PSM = {}
    PRE1_DICT_PSM = {}

    PRO1_DICT_PSM = {}  # 这个用dict是为了方便建倒排
    PRO2_DICT_SEQ = {}  # 这个是为了得到group弄的，里面seq是去冗余的
    PRO3_DICT_N_SEQ = {}  # 这个是为了得到group弄的
    PRO4_LIST_AC = []
    PRO5_LIST_DE = []
    PRO6_LIST_ID = []
    PRO7_LIST_TYPE = []

    PSM_LINE_T = []  # 一般就一行，就是一个字符串。纯是为了实现定量结果在定性的基础上往后加。
    PSM_LINE_C = []

class CFileIDForDIACheck:  # 注意：这是按列存储，每个属性都是list。如果按行搞，这个行是一个对象，不太好管理。

    N_ID = 0
    ID00_RUN = []
    ID0_SEQ = []
    ID1_RAW_NAME = []
    ID2_SCAN_ID = []
    ID3_RT_APEX = []

    ID4_SEQ_WITH_MOD = []  # 算碎片离子质量的时候比较方便

    ID5_PRECURSOR_MOZ_EXP = []
    ID6_PRECURSOR_MOZ_CLC = []
    ID7_CHARGE = []

    ID8_PRECURSOR_ID = []  # 存放母离子唯一ID（raw文件+序列+修饰+电荷），方便快速索引

    ID9_SCORE0 = []  # 专门放FDR或q-value

    # 蛋白质相关
    ID10_Protein_Groups = []
    ID11_Protein_Ids = []
    ID12_Protein_Names = []
    ID13_Genes = []

    # 起始时间和终止时间
    ID14_RT_BEGIN = []
    ID15_RT_END = []

    ID16_MOD = []

    ID17_SEQ_WITH_MOD_AND_CHARGE = []  # 后续肽段建库需要用到

    ID18_FRAGMENT_TYPE = []  # 每个元素是一个列表，存放library中该母离子的所有子离子类型：如b2++,y5+
    ID19_FRAGMENT_MOZ_EXP = []  # 每个元素是一个列表
    ID20_FRAGMENT_MOZ_CLC = []  # 每个元素是一个列表
    ID21_FRAGMENT_LOSS_TYPE = []

    ID22_TARGET = []

    MOD1_DICT_PSM = {}

    PRO1_DICT_PSM = {}  # 这个用dict是为了方便建倒排
    PRO2_DICT_SEQ = {}  # 这个是为了得到group弄的，里面seq是去冗余的
    PRO3_DICT_N_SEQ = {}  # 这个是为了得到group弄的
    PRO4_LIST_AC = []
    PRO5_LIST_DE = []
    PRO6_LIST_ID = []
    PRO7_LIST_TYPE = []

    PSM_LINE_T = []  # 一般就一行，就是一个字符串。纯是为了实现定量结果在定性的基础上往后加。
    PSM_LINE_C = []

class Config:
    I0_INI_PATH_ELEMENT = 'element.ini'
    I1_INI_PATH_AA = 'aa.ini'
    I2_INI_PATH_MOD = 'modification.ini'

    A0_FLAG_USE_EXIST_INDEX = 1
    A1_PATH_MS1 = ''  # 里面有多个文件
    A2_TYPE_MS1 = 0
    A3_PATH_MS2 = ''
    A4_TYPE_MS2 = 0
    A5_PATH_FILE = ''
    A6_TYPE_FILE = 0

    #鉴定结果相关
    B1_PATH_IDENTIFICATION_RESULT = ''
    B2_TYPE_IDENTIFICATION_RESULT = 0
    B3_THRESHOLD_FDR = 0.01
    B4_PREDICTMS2_TYPE_MODEL = 0
    B5_FLAG_DYNAMIC_LIBRARY = 0
    B6_FLAG_DRAW = 0
    B7_FRAG_NUM = 13

    #定量相关
    C0_TYPE_FLOW = 0
    C1_DIA_RT_HALF_WIN_IN_MIN = 0.1  # RT窗口
    C2_DIA_ACCURACY_HALF_WIN_PEAK = 20
    C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK = 0  # 0: ppm ，1: Da
    C4_DIALF_TYPE_MODEL = 0
    C5_DIALF_FLAG_CALIBRATE_RT = 0
    C5_DIALF_PATH_CALIBRATE_RT = '' #如果不需要校正，直接就不填
    C6_FLAG_DECOY = 0
    C7_PROCESS_NUM = 0

    D1_PATH_EXPORT = './MSTestResults/'
    D2_TYPE_EXPORT = 1


class CFileMS1:
    INDEX_SCAN = []  # 这两个大小和大家不一样，方便快速索引
    INDEX_RT = []

    LIST_RET_TIME = []
    LIST_ION_INJECTION_TIME = []

    MATRIX_PEAK_MOZ = []  # 每一行是个list
    MATRIX_PEAK_INT = []
    DICT_PEAK_MOZ = [] # 每一行是个字典


class CFileMS2:  # 注意：这是按列存储，每个属性都是list。如果按列搞，这个行是一个对象，不太好管理。

    INDEX_SCAN = []  # 这个大小和大家不一样，方便快速索引
    INDEX_RT = []

    LIST_RET_TIME = []  # 2019.08.31，MS1的改了，MS2的应该也要改。为DIA做准备。
    LIST_ION_INJECTION_TIME = []

    LIST_ACTIVATION_CENTER = []
    LIST_PRECURSOR_SCAN = []
    LIST_PRECURSOR_MZ = []

    MATRIX_PEAK_MOZ = []  # 每一行是个list
    MATRIX_PEAK_INT = []
    DICT_PEAK_MOZ = []  # 每一行是个字典


class CQuantInfo:
    MATRIX_QUANT_INT_PSM = []
    MATRIX_QUANT_INT_MOD = []
    MATRIX_QUANT_INT_PRE = []
    MATRIX_MATCH_SCORE = []  # 记录match between runs的打分
    LIST_FLAG_INFER_PSM = []
    MATRIX_RT_LENGTH_PSM = []
    MATRIX_RT_START_PSM = []
    MATRIX_RT_MID_PSM = []


class CSeed:  # 为得到Evidence准备的，Labeling只需要这就可以了。LabelFree还需要个Reference。


    MID_SCAN = 0  # 根据需求不同，使用scan或rt
    MID_RT = 0

    DICT_COMPOSITION = {}  # 这个有时也需要

    DIS_ISO_MOZ_CLC = []  # 母离子理论质量
    DIS_ISO_INT_CLC = []  # 母离子理论强度

    DIS_FRA_MOZ_CLC = []
    DIS_FRA_INT_CLC = [] # 用Prosit预测理论强度
    DIS_FRA_TYPE_CLC = [] # 例：b3++（目前考虑一价和二价）

    INDEX_MONO = 0  # 上面几根峰，哪根是mono

    CenterMoz = -1

class CEvidence:
    # 这是第二批要填的信息
    MATRIX_PROFILE_PREC = []# 每一行是个曲线，一样长
    MATRIX_PROFILE_FRAG = []

    MATRIX_MASS_DEV_PREC = []#质量差值
    MATRIX_MASS_DEV_FRAG = []

    # 这是第三批要填的信息
    LIST_RET_TIME = []  # 曲线上每个点的保留时间
    LIST_SCAN = []  # 曲线上每个点的scan

    LIST_I_START_PRE = []  # 每条曲线单独确定起止
    LIST_I_END_PRE = []

    LIST_I_START_FRA = []  # 每条曲线单独确定起止
    LIST_I_END_FRA = []

    POINT_APEX = 0

    I_START = 0  # 最后确定总的起止。在老pQuant里面，这一步很关键。所有profile要取最短的。（目前都是用的mono的起止）
    I_END = -1

    # 这是第四批要填的信息
    DIS_ISO_MOZ_EXP = []  # 实验质量
    DIS_ISO_INT_EXP = []  # 实验强度
    DIS_FRAG_MOZ_EXP = []  # 实验质量
    DIS_FRAG_INT_EXP = []  # 实验强度
    PROFILE_ALL = []  # 一行；是所有同位素峰加和的曲线

    FRAG_TYPE = []
    FRAG_MOZ = []

    SCORE_IS_PEPTIDE = 1.0


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

class CRTBetweenExp:

    LIST_EXPERIMENT = []
    DIC_EXPERIMENT_INDEX = {}

    DIC_MBR_WIDTH = {}
    DIC_RT_SHIFT = {}

    DIC_RT_CORR = {}

class CXLCInfo:
    LIST_PRE = []
    LIST_FRAG = []

class CDataPack:
    myCFG = Config()
    myINI = CINI()
    myRTShift = CRTBetweenExp()
    myQUANT = CQuantInfo()
    myID = CFileID()  # 从鉴定结果到seed，seed到evidence，这两个就当临时变量了。如果需要输出，直接写到文件里，不全周期维护了。
    myIDForDIACHeck = CFileIDForDIACheck()
    myXLCInfo = CXLCInfo()

    LIST_PATH_MS = [] # 单纯存一下raw名
    LIST_PATH_MS1 = []  # ms1文件
    LIST_PATH_MS2 = []  # ms2文件
    LIST_PATH_ID = [] # 鉴定结果文件

    #用于构建正负样本对的
    DIC_PSM_PAIR = {}
    DIC_EVIDENCE = {}
    DIC_EVIDENCE_DECOY = {}

    #用于谱图预测的
    DIC_PRECURSOR_TO_IONS = {}

    NUM_OF_PAIRS = 0

    MS1_START = 380
    MS1_END = 980
    MS2_START = 150
    MS2_END = 2000

    MZ_BIN = 1