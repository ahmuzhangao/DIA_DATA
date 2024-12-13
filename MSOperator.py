from MSLogging import logGetError
from MSTool import toolCountCharInString, toolGetWord, toolGetNameFromPath, toolFindNeighborFromSortedList1
from MSData import CFileMS1, CFileMS2, CINI, CEvidenceRerank, CSeedDIACheck, CEvidenceDIACheck
from MSSystem import VALUE_MAX_SCAN, VALUE_ILLEGAL,UNIMOID_TO_STANDARD_MOD,FRAGMENT_LOSS_TYPE_COMP, \
    CFG_TYPE_ACCURACY_HALF_WIN_PEAK
from MSData import CResult
import os
import numpy as np


def op_INIT_CFILE_MS1(dataMS1: CFileMS1):
    dataMS1.INDEX_SCAN = []
    dataMS1.INDEX_RT = []

    dataMS1.LIST_ION_INJECTION_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN

    dataMS1.MATRIX_PEAK_INT = [[] * 1] * VALUE_MAX_SCAN
    dataMS1.MATRIX_PEAK_MOZ = [[] * 1] * VALUE_MAX_SCAN


def op_INIT_CFILE_MS2(dataMS2: CFileMS2):
    dataMS2.INDEX_SCAN = []
    dataMS2.INDEX_RT = []

    dataMS2.LIST_RET_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS2.LIST_ION_INJECTION_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS2.LIST_PRECURSOR_SCAN = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS2.LIST_ACTIVATION_CENTER = [VALUE_ILLEGAL] * VALUE_MAX_SCAN

    dataMS2.MATRIX_PEAK_MOZ = [[] * 1] * VALUE_MAX_SCAN
    dataMS2.MATRIX_PEAK_INT = [[] * 1] * VALUE_MAX_SCAN


def op_INIT_CSEED_DIACHECK(CSeed: CSeedDIACheck):

    CSeed.MID_RT = VALUE_ILLEGAL
    CSeed.MID_SCAN = VALUE_ILLEGAL
    CSeed.TYPE_FRAGMENT = []
    CSeed.MOZ_FRAGMENT = []
    CSeed.MOZ_PRECURSOR = []
    CSeed.MOZ_LIBRARY = []
    CSeed.TYPE_LIBRARY = []


def op_ININT_CEVIDENCE_DIACHECK(CEvidence: CEvidenceDIACheck):

    CEvidence.MID_IDNEX = VALUE_ILLEGAL
    CEvidence.RT_START = VALUE_ILLEGAL
    CEvidence.RT_END = VALUE_ILLEGAL
    CEvidence.LIST_RET_TIME = []
    CEvidence.LIST_SCAN = []
    CEvidence.MATRIX_MS1_PEAK_MOZ = []
    CEvidence.MATRIX_MS1_PEAK_INT = []
    CEvidence.MATRIX_MS2_PEAK_MOZ = []
    CEvidence.MATRIX_MS2_PEAK_INT = []



def op_INIT_CEVIDENCE_RERANK(CEvidence: CEvidenceRerank):

    CEvidence.MID_IDNEX = VALUE_ILLEGAL
    CEvidence.MATRIX_MS2_PEAK_INT = []
    CEvidence.MATRIX_MS2_PEAK_MOZ = []
    CEvidence.MATRIX_MS1_PEAK_INT = []
    CEvidence.MATRIX_MS1_PEAK_MOZ = []


def op_INIT_RESULT(inputResult: CResult):

    inputResult.N_RESULT = 0

    inputResult.RE0_RAW_NAME = []
    inputResult.RE1_SEQ = []
    inputResult.RE2_SCAN_NO = []
    inputResult.RE3_RT = []

    inputResult.RE4_MOD = []
    inputResult.RE5_PRECURSOR_MOZ_CLC = []
    inputResult.RE6_CHARGE = []
    inputResult.RE7_PRECURSOR_ID = []

    inputResult.RE8_SCORE0 = []  # finall score
    inputResult.RE9_SCORE1 = []  # first DIA-NN score before Rerank
    inputResult.RE10_SCORE_RERANK = []
    inputResult.RE11_TARGET = []
    inputResult.RE12_FEATURE = []


def op_FILL_LIST_PATH_ID(inputPath, inputList):
    separator = '|'

    if len(inputPath) < 1:
        logGetError("MSOperator.py, op_FILL_LIST_PATH_ID, MK_1: Path for identification results is empty!")

    if inputPath[-1] == separator:
        pass
    else:
        inputPath = inputPath + separator

    nFile = toolCountCharInString(inputPath, separator)

    for i in range(nFile):
        inputList.append(toolGetWord(inputPath, i, separator))


def op_FILL_LIST_PATH_MS(inputPath, inputList, inputExt):

    separator = '|'
    listStrPath = []

    if len(inputPath) < 1:
        logGetError("MSOperator.py, op_FILL_LIST_PATH_MS, MK_1: Path for MS is empty!")

    if inputPath[-1] == separator:
        pass
    else:
        inputPath = inputPath + separator

    nFile = toolCountCharInString(inputPath, separator)

    for i in range(nFile):
        listStrPath.append(toolGetWord(inputPath, i, separator))

    for strPath in listStrPath:

        if os.path.isdir(strPath):

            for maindir, subdir, file_name_list in os.walk(strPath):

                for filename in file_name_list:

                    tmpPath = os.path.join(maindir, filename)  # 合并成一个完整路径

                    ext = os.path.splitext(tmpPath)[1]  # 获取文件后缀 [0]获取的是除了文件名以外的内容

                    if ext in inputExt:
                        inputList.append(tmpPath)

        elif os.path.isfile(strPath):

            inputList.append(strPath)

        else:

            logGetError("MSOperator.py, op_FILL_LIST_PATH_MS, MK_2: Path for MS is illegal!")


def opGetStartAndEndForProfile(input_profile, input_seed, input_cutoff, input_n_hole):  # 就这个名字格式特殊

    nHoleLeft = 0
    nHoleRight = 0

    i_middle = input_seed

    if i_middle > 0:
        i_left = i_middle - 1
    else:
        i_left = i_middle

    i_right = i_middle

    result = [i_left, i_right]  # start and end

    int_left = input_profile[i_left]  # int ->当前左边的强度
    int_right = input_profile[i_right]

    int_max = int_left
    int_thr = int_max * input_cutoff

    walkLeft = True
    walkRight = True

    while True:

        if walkLeft or walkRight:
            pass
        else:
            break

        if i_left > 0:

            if walkLeft:
                i_left = i_left - 1

        else:

            walkLeft = False

        if i_right < len(input_profile) - 1:

            if walkRight:
                i_right = i_right + 1

        else:

            walkRight = False

        int_left = input_profile[i_left]
        int_right = input_profile[i_right]

        # max
        if int_max < int_left:
            int_max = int_left

        if int_max < int_right:
            int_max = int_right

        int_thr = int_max * input_cutoff

        # hole
        if int_left < int_thr and walkLeft:
            nHoleLeft = nHoleLeft + 1

        if int_right < int_thr and walkRight:
            nHoleRight = nHoleRight + 1

        if nHoleLeft == input_n_hole:

            walkLeft = False

        if nHoleRight == input_n_hole:

            walkRight = False

    result[0] = i_left
    result[1] = i_right

    return result


def op_STANDARDIZE_MOD_BY_PRECURSOR(input_precursor: str):
    # 将DIA-NN带修饰的母离子转化为标准格式的修饰类型母离子:(UniMod:4)-> Carbamidomethyl
    seperator_start, seperator_end = '(', ')'
    new_precursor = ''
    site_start, site_end = 0, 0
    flag_in_mod = 0
    for i, char in enumerate(input_precursor):
        if char == seperator_start:
            site_start = i
            flag_in_mod = 1
        elif char == seperator_end:
            site_end = i
            mod_type = input_precursor[site_start + 1:site_end]

            try:
                mod_type_standard = UNIMOID_TO_STANDARD_MOD[mod_type]
            except KeyError:
                print(UNIMOID_TO_STANDARD_MOD)
                print(mod_type)
                info = 'DIA-NN mod type has wrong, please check!'
                logGetError(info)

            new_precursor = new_precursor + '(' + mod_type_standard + ')'
            flag_in_mod = 0
        elif flag_in_mod == 0:
            new_precursor = new_precursor + char

    return new_precursor

def op_STANDARDIZE_MOD_MSFraggerPin_BY_PRECURSOR(input_precursor: str):
    # 将DIA-NN带修饰的母离子转化为标准格式的修饰类型母离子:(UniMod:4)-> Carbamidomethyl
    seperator_start, seperator_end = '[', ']'
    new_precursor = ''
    site_start, site_end = 0, 0
    flag_in_mod = 0
    for i, char in enumerate(input_precursor):
        if char == seperator_start:
            site_start = i
            flag_in_mod = 1
        elif char == seperator_end:
            site_end = i
            mod_type = input_precursor[site_start + 1:site_end]

            try:
                mod_type_standard = UNIMOID_TO_STANDARD_MOD[mod_type]
            except KeyError:
                print(UNIMOID_TO_STANDARD_MOD)
                print(mod_type)
                info = 'DIA-NN mod type has wrong, please check!'
                logGetError(info)

            new_precursor = new_precursor + '(' + mod_type_standard + ')'
            flag_in_mod = 0
        elif flag_in_mod == 0:
            new_precursor = new_precursor + char

    return new_precursor


def op_DIVIDE_MOD_FROM_PRECURSOR(input_precursor: str):
    seperator_start, seperator_end = '(', ')'
    precursor_sequence = ''
    mod_dic = {}
    site_start, site_end = 0, 0
    flag_in_mod = 0
    count = 0
    for i, char in enumerate(input_precursor):
        if char == seperator_start:
            site_start = i
            flag_in_mod = 1
        elif char == seperator_end:
            site_end = i
            mod_type = input_precursor[site_start + 1:site_end]
            mod_dic[count] = mod_type
            flag_in_mod = 0
        elif flag_in_mod == 0:
            count = count + 1
            precursor_sequence = precursor_sequence + char

    return precursor_sequence, mod_dic


def op_DIVIDE_ION_TYPE(input_ion_type: str):
    # 拆分by离子类型（如：b2++)到[b/y类型，离子位点，离子电荷数]（如：[b,2,2])
    ion_by = input_ion_type[0]
    site_index = 1
    while input_ion_type[site_index] != '+':
        site_index += 1
    ion_site = int(input_ion_type[1:site_index])
    ion_charge = len(input_ion_type[site_index:])

    return ion_by, ion_site, ion_charge


def op_DIVIDE_ION_TYPE_WITHLOSS(input_ion_type: str):
    # 拆分by离子类型（如：b2++_H2O)到[b/y类型，离子位点，离子电荷数]（如：[b,2,2])
    # 先拆分离子类型和损失
    ion_type_list = input_ion_type.split('_')
    ion_type = ion_type_list[0] if len(ion_type_list) == 1 else ion_type_list[0]
    return op_DIVIDE_ION_TYPE(ion_type)


def op_CAL_MOLECULAR_MASS_FROM_FRAGMENT_LOSS(input_INI_ELEMENT: CINI.DICT0_ELEMENT_MASS, input_molecular:str):
    mass_molecular = 0.0
    molecular_comp = FRAGMENT_LOSS_TYPE_COMP[input_molecular]
    atom_left = 0
    atom_right = 0
    for i, char in enumerate(molecular_comp):
        if char == '(':
            atom_right = i
            atom = molecular_comp[atom_left:atom_right]
            num_left = i + 1
        elif char == ')':
            num_right = i
            atom_left = i + 1
            atom_num = int(molecular_comp[num_left:num_right])
            mass_molecular = mass_molecular + input_INI_ELEMENT[atom][0] * atom_num

    return mass_molecular


def op_CAL_PERCURSOR_MOZ(input_INI: CINI, input_precursor: str, input_charge):
    precursor_sequence, mod_dic = op_DIVIDE_MOD_FROM_PRECURSOR(input_precursor)
    ion_num = len(precursor_sequence)
    mass_aa_sum = 0.0
    for i in range(ion_num):
        if (i + 1) in mod_dic.keys():
            mass_aa_mod = float(input_INI.DICT2_MOD_COM[mod_dic[i + 1]][1])  # 修饰质量
        else:
            mass_aa_mod = 0

        mass_aa = input_INI.DICT1_AA_COM[precursor_sequence[i]][1] + mass_aa_mod
        mass_aa_sum = mass_aa_sum + mass_aa
    mass_precursor = (mass_aa_sum + 18.0105633 + input_INI.MASS_PROTON_MONO * input_charge) / input_charge
    return mass_precursor


def op_CAL_FRAGMENT_MOZ(input_INI: CINI, input_precursor: str, input_ion_list: list, input_ion_loss_list: list):
    precursor_sequence, mod_dic = op_DIVIDE_MOD_FROM_PRECURSOR(input_precursor)
    ion_num = len(precursor_sequence)
    # print('precursor_sequence', precursor_sequence)
    # print('mod_dic', mod_dic)
    # 构建by离子索引
    ion_dic = {}
    for i, ion in enumerate(input_ion_list):
        tmp_ion_by, tmp_ion_site, tmp_ion_charge = op_DIVIDE_ION_TYPE(ion)
        if tmp_ion_site == 0:
            break
        tmp_ion_loss = input_ion_loss_list[i]
        index_key = tmp_ion_by + str(tmp_ion_site)
        if index_key not in ion_dic:
            ion_dic[index_key] = [[tmp_ion_charge, tmp_ion_loss]]
        else:
            ion_dic[index_key].append([tmp_ion_charge, tmp_ion_loss])
    mass_b_aa_sum = 0.0
    mass_y_aa_sum = 0.0
    for i in range(ion_num):
        # b离子
        if (i + 1) in mod_dic.keys():
            mass_aa_mod = float(input_INI.DICT2_MOD_COM[mod_dic[i + 1]][1])  # 修饰质量
        else:
            mass_aa_mod = 0
        mass_aa = input_INI.DICT1_AA_COM[precursor_sequence[i]][1] + mass_aa_mod
        mass_b_aa_sum = mass_b_aa_sum + mass_aa
        # 去by离子索引中计算指定b离子m/z
        index_key_b = 'b' + str(i + 1)
        if index_key_b in ion_dic:
            for j, value in enumerate(ion_dic[index_key_b]):
                charge = value[0]
                loss_type = value[1]
                mass_loss = op_CAL_MOLECULAR_MASS_FROM_FRAGMENT_LOSS(input_INI.DICT0_ELEMENT_MASS, loss_type)
                moz_b_ion = (mass_b_aa_sum + input_INI.MASS_PROTON_MONO * charge - mass_loss) / charge
                ion_dic[index_key_b][j].append(moz_b_ion)

        # y离子
        if (ion_num - i) in mod_dic.keys():
            mass_aa_mod = float(input_INI.DICT2_MOD_COM[mod_dic[ion_num - i]][1])  # 修饰质量
        else:
            mass_aa_mod = 0
        mass_aa = input_INI.DICT1_AA_COM[precursor_sequence[ion_num - i - 1]][1] + mass_aa_mod
        mass_y_aa_sum = mass_y_aa_sum + mass_aa
        # 去by离子索引中计算指定y离子m/z
        index_key_y = 'y' + str(i + 1)
        if index_key_y in ion_dic:
            for j, value in enumerate(ion_dic[index_key_y]):
                charge = value[0]
                loss_type = value[1]
                mass_loss = op_CAL_MOLECULAR_MASS_FROM_FRAGMENT_LOSS(input_INI.DICT0_ELEMENT_MASS, loss_type)
                moz_y_ion = (mass_y_aa_sum + input_INI.MASS_PROTON_MONO * charge + 18.0105633 - mass_loss) / charge
                ion_dic[index_key_y][j].append(moz_y_ion)
    ionTypeList = []
    ionMozList = []
    for key in ion_dic.keys():
        for value in ion_dic[key]:
            if value[1] == 'NOLOSS':
                loss_type = ''
            else:
                loss_type = '_' + value[1] + 'loss'
            ionTypeList.append(key + '+' * (value[0]) + loss_type)
            ionMozList.append(value[2])

    return ionTypeList, ionMozList


def op_FILL_LIST_DRAW_PRECURSOR(input_sequence: str, input_charge: str, input_raw: str):
    separator = '|'
    listStrSequence = []

    # if len(input_sequence) < 1:
    #     logGetError("MSOperator.py, op_FILL_LIST_DRAW_PRECURSOR, MK_1: Sequence for precursor is empty!")
    #
    # if input_sequence[-1] == separator:
    #     pass
    # else:
    #     input_sequence = input_sequence + separator
    #
    # nSequence = toolCountCharInString(input_sequence, separator)
    #
    # for i in range(nSequence):
    #     listStrSequence.append(toolGetWord(input_sequence, i, separator))
    listStrSequence = input_sequence.split('|')
    listStrSequence = [op_STANDARDIZE_MOD_BY_PRECURSOR(i) for i in listStrSequence]

    listStrCharge = []

    # if len(input_charge) < 1:
    #     logGetError("MSOperator.py, op_FILL_LIST_DRAW_PRECURSOR, MK_1: charge for precursor is empty!")
    #
    # if input_charge[-1] == separator:
    #     pass
    # else:
    #     input_charge = input_charge + separator
    # nCharge = toolCountCharInString(input_charge, separator)
    #
    # for i in range(nCharge):
    #     listStrCharge.append(toolGetWord(input_charge, i, separator))
    listStrCharge = input_charge.split('|')

    listStrRaw = []

    # if len(input_raw) < 1:
    #     logGetError("MSOperator.py, op_FILL_LIST_DRAW_PRECURSOR, MK_1: charge for precursor is empty!")
    #
    # if input_raw[-1] == separator:
    #     pass
    # else:
    #     input_raw = input_raw + separator
    #
    # nRaw = toolCountCharInString(input_raw, separator)
    #
    # for i in range(nRaw):
    #     listStrRaw.append(toolGetWord(input_raw, i, separator))
    listStrRaw = input_raw.split('|')
    listStrRaw = [toolGetNameFromPath(i) for i in listStrRaw]
    print(listStrRaw)
    listStrPrecursor = [listStrRaw[i] + '_' + listStrSequence[i] + listStrCharge[i] for i in range(len(listStrCharge))]

    return listStrPrecursor


def op_FILL_LIST_DRAW_DDA_PRECURSOR(input_ScanNo: str, input_raw: str):
    # 去除注释';*'和空格符
    input_ScanNo = input_ScanNo.split(';')[0].split()[0]
    input_raw = input_raw.split(';')[0].split()[0]

    input_ScanNo = input_ScanNo[:-1] if input_ScanNo[-1] == '|' else input_ScanNo
    input_raw = input_raw[:-1] if input_raw[-1] == '|' else input_raw

    input_ScanNo = input_ScanNo.split('|')
    input_raw = input_raw.split('|')
    input_raw = [toolGetNameFromPath(raw) for raw in input_raw]

    precursor_id_list = [input_raw[i]+'_'+input_ScanNo[i] for i in range(len(input_ScanNo))]

    return precursor_id_list


def op_StandardList(l_list):
    # 对数组进行标准化
    mean = np.mean(l_list)
    std = np.std(l_list)
    l_list = (l_list-mean)/(std + 1e-8)
    return l_list


def op_GetCycleWinNum(inputDataMS2: CFileMS2, MidIndex):

    centerMoz = inputDataMS2.LIST_ACTIVATION_CENTER[inputDataMS2.INDEX_SCAN[MidIndex]]
    try:
        num = 1
        while centerMoz != inputDataMS2.LIST_ACTIVATION_CENTER[inputDataMS2.INDEX_SCAN[MidIndex + num]]:
            num = num + 1
    except IndexError:  # 超出激活中心列表的边界，从左边计算
        num = 1
        while centerMoz != inputDataMS2.LIST_ACTIVATION_CENTER[inputDataMS2.INDEX_SCAN[MidIndex - num]]:
            num = num + 1

    return num


def op_GetIndexListFromMS2(inputDataMS1: CFileMS1, inputDataMS2: CFileMS2, Midindex, RT_start, RT_end):

    cycle = op_GetCycleWinNum(inputDataMS2, Midindex)
    MidRT = inputDataMS2.INDEX_RT[Midindex]

    RT_left = MidRT
    tmp_index = Midindex
    count_left = 0
    left_list_ms2 = []
    left_list_ms1 = []
    # while (RT_left >= RT_start * 60) & (tmp_index >= cycle):
    while (count_left < 9) & (tmp_index >= cycle):
        tmp_index -= cycle
        count_left += 1
        RT_left = inputDataMS2.INDEX_RT[tmp_index]
        tmp_index1 = toolFindNeighborFromSortedList1(inputDataMS1.INDEX_RT, RT_left)
        if tmp_index1 > tmp_index:
            tmp_index1 -= 1
        left_list_ms1.append(tmp_index1)
        left_list_ms2.append(tmp_index)
    left_list_ms2.reverse()
    left_list_ms1.reverse()
    RT_right = MidRT
    tmp_index = Midindex
    right_list_ms2 = []
    right_list_ms1 = []
    count_right = 0
    # while (RT_right <= RT_end * 60) & (tmp_index + cycle < len(inputDataMS2.INDEX_RT)):
    while (count_right < 10) & (tmp_index + cycle < len(inputDataMS2.INDEX_RT)):
        tmp_index += cycle
        count_right += 1
        RT_right = inputDataMS2.INDEX_RT[tmp_index]
        tmp_index1 = toolFindNeighborFromSortedList1(inputDataMS1.INDEX_RT, RT_right)
        if tmp_index1 > tmp_index:
            tmp_index1 -= 1
        right_list_ms1.append(tmp_index1)
        right_list_ms2.append(tmp_index)

    Midindex1 = toolFindNeighborFromSortedList1(inputDataMS1.INDEX_RT, inputDataMS2.INDEX_RT[Midindex])
    if Midindex1 > Midindex:
        Midindex1 -= 1

    index_list_ms2 = left_list_ms2 + [Midindex] + right_list_ms2
    index_list_ms1 = left_list_ms1 + [Midindex1] + right_list_ms1
    iMid = len(left_list_ms2)

    return index_list_ms1, index_list_ms2, iMid


def op_GetTopFromProfile(profileMatrix):
    '''
    :param profileMatrix: b/y离子的色谱曲线
    :return: 返回top6强度的b/y离子和色谱曲线矩阵中匹配到的b/y离子数目
    '''
    profileMatrix = np.array(profileMatrix)
    # 还没想好怎么鉴定色谱曲线，先卡匹配到的数目吧
    threshold = 0.7
    topList = np.zeros(shape=(profileMatrix.shape[0], 1))
    for i in range(profileMatrix.shape[0]):
        if sum(profileMatrix[i, :] > 0) / len(profileMatrix[i, :]) >= threshold:
            topList[i] = (sum(profileMatrix[i, :]))
    topList[topList < 0] = 0
    # 去除干扰的b/y离子色谱曲线

    #
    # 按照强度进行排序
    argsortList = np.argsort(-1 * topList, axis=0)
    top6List = argsortList[:6] if len(argsortList) >= 6 else argsortList
    numMatch = sum(topList > 0)

    return top6List, numMatch, topList


def op_GetColor(moz_type):
    '''
    设置离子颜色
    :param moz_type:离子类型
    :return: 颜色标识
    '''

    if moz_type == 'b':
        return 'blue'
    elif moz_type == 'y':
        return 'orange'
    elif moz_type == 'a':
        return 'green'
    elif moz_type == 'x':
        return 'cyan'
    elif moz_type == 'u':
        return 'red'
    elif moz_type == 'v':
        return 'brown'
    elif moz_type == 'p':
        return 'violet'
    elif moz_type == 'M':
        return 'deeppink'
    elif moz_type == 'yb':
        return 'indigo'
    elif moz_type == 'by':
        return 'purple'
    else:
        return 'yellow'


def op_GetFont(moz_type, size):
    # 设置顶部离子标注字典
    color = op_GetColor(moz_type)
    font2 = {'family': 'SimHei',
             'color': color,
             'weight': 'normal',
             'size': size,
             }
    return font2


def op_PlotMozText(clc_tag, plot_moz, moz_start):
    '''
    获取顶部肽段离子标签
    :param clc_tag: 匹配上的离子标签list
    :param plot_moz: 对应质荷比List
    :param moz_start: 需要空的空格字符串，str
    :return:
    '''
    # 顶部离子标注标签
    moz_text = '' + moz_start
    moz_underline_text = '' + moz_start
    moz_type, _, _ = op_DIVIDE_ION_TYPE_WITHLOSS(clc_tag[0])
    if moz_type == 'b':  # b离子
        for pep_index in range(len(plot_moz)):
            if plot_moz[pep_index] != 0:
                if pep_index < 9:
                    moz_text = moz_text + moz_type + "%s" % (pep_index + 1) + '  '
                else:
                    moz_text = moz_text + moz_type + '%s' % (pep_index + 1) + ' '
                moz_underline_text = moz_underline_text + '  __'
            else:
                moz_text = moz_text + '    '
                moz_underline_text = moz_underline_text + '    '
    elif moz_type == 'y':  # y离子
        moz_text = moz_text + '    '
        moz_underline_text = moz_underline_text + '    '
        for pep_index in range(len(plot_moz)):
            if plot_moz[len(plot_moz) - pep_index - 1] != 0:
                if (len(plot_moz) - pep_index - 1) < 9:
                    moz_text = moz_text + moz_type + "%s" % (len(plot_moz) - pep_index) + '  '
                else:
                    moz_text = moz_text + moz_type + '%s' % (len(plot_moz) - pep_index) + ' '
                moz_underline_text = moz_underline_text + '__  '
            else:
                moz_text = moz_text + '    '
                moz_underline_text = moz_underline_text + '    '

    return moz_text, moz_underline_text


def op_PlotMozTopText(clc_tag, plot_moz, ax, moz_start, pep_type):
    '''
    顶部肽段离子标注
    :param clc_tag: 匹配上离子标签字符串,str
    :param plot_moz: 匹配上离子质荷比list
    :param ax: 图2
    :param moz_start: 匹配上离子标签空格字符串，str
    :param pep_type: 肽段类型，alpha,beta，str
    :return:
    '''
    # 顶部离子标注
    moz_type, _, _ = op_DIVIDE_ION_TYPE_WITHLOSS(clc_tag[0])
    moz_text, moz_underline_text = op_PlotMozText(clc_tag, plot_moz, moz_start)
    if moz_type == 'b':
        if pep_type == 'alpha':
            loc1 = 0.75
            loc2 = 0.78
            ax.text(0.100, loc1, moz_text, fontdict=op_GetFont(moz_type, 15))
            ax.text(0.090, loc2, moz_underline_text, fontdict=op_GetFont(moz_type, 15))
        elif pep_type == 'beta':
            loc1 = 0.64
            loc2 = 0.67
            ax.text(0.100, loc1, moz_text, fontdict=op_GetFont(moz_type, 15))
            ax.text(0.090, loc2, moz_underline_text, fontdict=op_GetFont(moz_type, 15))
        elif pep_type == 'single':
            loc1 = 0
            loc2 = 14.5
            ax.text(100, loc1, moz_text, fontdict=op_GetFont(moz_type, 10))
            ax.text(80, loc2, moz_underline_text, fontdict=op_GetFont(moz_type, 10))
        else:
            loc1 = 0
            loc2 = 0

    elif moz_type == 'y':
        if pep_type == 'alpha':
            loc1 = 0.81
            loc2 = 0.81
            ax.text(0.130, loc1, moz_text, fontdict=op_GetFont(moz_type, 15))
            ax.text(0.120, loc2, moz_underline_text, fontdict=op_GetFont(moz_type, 15))
        elif pep_type == 'beta':
            loc1 = 0.70
            loc2 = 0.70
            ax.text(0.130, loc1, moz_text, fontdict=op_GetFont(moz_type, 15))
            ax.text(0.120, loc2, moz_underline_text, fontdict=op_GetFont(moz_type, 15))
        elif pep_type == 'single':
            loc1 = 30.5
            loc2 = 30.5
            ax.text(100, loc1, moz_text, fontdict=op_GetFont(moz_type, 10))
            ax.text(80, loc2, moz_underline_text, fontdict=op_GetFont(moz_type, 10))
        else:
            loc1 = 0
            loc2 = 0
    else:
        pass
    # elif moz_type == 'M':x


def op_FindSame(exe, data, peak, bias_type, bias):
    '''
    查找最近的质荷比
    :param exe: 理论质荷比list
    :param data: 实际质荷比list
    :param peak: 峰强list
    :param bias_type: 偏差类型，ppm或Da
    :param bias: 偏差大小，float
    :return: outputL:满足误差范围内的质荷比list, output_peak：对应强度的List, output_ppm：对应的偏差list
    '''
    # 查找接近的质荷比
    # exe理论质荷比，data实际质荷比，peak实际峰值(已换算成百分比)，bias偏差
    output = list()
    output_ppm = list()
    output_peak = list()

    for index in exe:  # 遍历理论数据
        first = 0
        last = len(data) - 1
        if (index <= data[last]) & (index >= data[first]):
            flag = 0
            while first <= last:  # 二分查找
                mid = (first + last) // 2
                if index > data[mid]:
                    first = mid + 1
                elif index < data[mid]:
                    last = mid - 1
                else:
                    first = mid
                    last = mid - 1
                    break
        else:
            flag = -2
        if flag == 0:
            #Da偏差
            if bias_type == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['Da']:
                if abs(data[first] - index) < bias:
                    loc = first
                    if peak[loc] == 0:
                        output.append(0)
                        output_peak.append(0)
                        output_ppm.append(0)
                    else:
                        output.append(data[loc])
                        output_peak.append(peak[loc])
                        output_ppm.append((- index + data[loc]) / index * 1000000)

                elif abs(data[last] - index) < bias:
                    loc = last
                    if peak[loc] == 0:
                        output.append(0)
                        output_peak.append(0)
                        output_ppm.append(0)
                    else:
                        output.append(data[loc])
                        output_peak.append(peak[loc])
                        output_ppm.append((- index + data[loc]) / index * 1000000)

                else:
                    flag = -2

            elif bias_type == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['PPM']:
                if (abs(data[first] - index) / index * 1000000) < bias:
                    loc = first
                    if peak[loc] == 0:
                        output.append(0)
                        output_peak.append(0)
                        output_ppm.append(0)
                    else:
                        #output.append(data[loc])
                        output.append(index)
                        output_peak.append(peak[loc])
                        output_ppm.append((- index + data[loc]) / index * 1000000)

                elif (abs(data[last] - index) / index * 1000000) < bias:
                    loc = last
                    if peak[loc] == 0:
                        output.append(0)
                        output_peak.append(0)
                        output_ppm.append(0)
                    else:
                        #output.append(data[loc])
                        output.append(index)
                        output_peak.append(peak[loc])
                        output_ppm.append((- index + data[loc]) / index * 1000000)
                else:
                    flag = -2
        if flag == -2:
            output.append(0)
            output_peak.append(0)
            output_ppm.append(0)
    return output, output_peak, output_ppm


def op_Data_fill_plot(clc_moz, clc_tag, exp_moz, exp_intensity, type_bias, bias, ax1, ax2, ax3_half, ax3, ax4, pep_moz_start,
                      type_pep, record_used_moz):
    '''
    :param clc_moz: 理论匹配上的子离子质量
    :param clc_tag: 理论子离子标签
    :param exp_moz: 实际谱图质荷比
    :param exp_intensity: 实际谱图强度
    :param type_bias: 偏差类型
    :param bias: 偏差大小
    :param ax2: 图2
    :param ax3: 图3
    :param ax4: 图4
    :param pep_moz_start: 肽段开始的空格
    :param type_pep: 肽段类型
    :return:
    '''
    # 将输入的离子分类
    clc_moz_type = []
    left_loc = 0

    if len(clc_tag) > 0:
        flag1, _, flag2 = op_DIVIDE_ION_TYPE_WITHLOSS(clc_tag[0])

        for clc_tag_index in range(len(clc_tag)):

            next_flag1, _, next_flag2 = op_DIVIDE_ION_TYPE_WITHLOSS(clc_tag[clc_tag_index])
            if (next_flag1 != flag1) | (flag2 != next_flag2):
                right_loc = clc_tag_index
                clc_moz_type.append([left_loc, right_loc])
                flag1, flag2 = next_flag1, next_flag2
                left_loc = clc_tag_index
        clc_moz_type.append([left_loc, len(clc_tag)])
        moz_type_num = len(clc_moz_type)

        plot_moz_list, plot_peak_list, ppm_list, clc_tag_list = [], [], [], []
        for math_index in range(moz_type_num):
            clc_tag_one = clc_tag[clc_moz_type[math_index][0]:clc_moz_type[math_index][1]]
            clc_moz_one = clc_moz[clc_moz_type[math_index][0]:clc_moz_type[math_index][1]]
            plot_moz, plot_peak, ppm = op_FindSame(clc_moz_one, exp_moz, exp_intensity, type_bias,
                                                   bias)
            plot_moz_list.append(plot_moz)
            plot_peak_list.append(plot_peak)
            ppm_list.append(ppm)
            clc_tag_list.append(clc_tag_one)

        max_match_peak = max(np.max(plot_peak_list), 1)
        exp_intensity_match = np.array(exp_intensity) / max_match_peak * 100.
        exp_intensity_match[exp_intensity_match > 1.2 * 100] = 1.2 * 100.
        info_base_peak = '{:.3g}'.format(max_match_peak)

        ax1.text(0, 5, 'BasePeak:', color='black', fontsize=10)
        ax1.text(8, 5, info_base_peak, color='red', fontsize=10)

        ax3.bar(exp_moz, exp_intensity_match, width=2, facecolor='#BEBEBE', alpha=0.5)

        for plot_moz, plot_peak, ppm, clc_tag_one in zip(plot_moz_list, plot_peak_list, ppm_list, clc_tag_list):
            plot_peak_match = list(np.array(plot_peak) / max_match_peak * 100)
            if type_pep == 'precursor':
                pass
            else:
                print(clc_tag_one)
                print(plot_moz)
                op_PlotMozTopText(clc_tag_one, plot_moz, ax2, pep_moz_start, type_pep)  # 标顶部离子

            # 谱图离子标注
            i = 1
            moz_type, _, valence = op_DIVIDE_ION_TYPE_WITHLOSS(clc_tag_one[0])

            ax3.bar(plot_moz, plot_peak_match, width=2, facecolor=op_GetColor(moz_type))
            ax3_half.bar(plot_moz, plot_peak, width=2, facecolor=op_GetColor(moz_type))
            for xx, yy in zip(plot_moz, plot_peak_match):
                if plot_peak_match[i - 1] != 0:
                    if type_pep == 'alpha':
                        p_str_start = 'α '
                    elif type_pep == 'beta':
                        p_str_start = 'β '
                    elif type_pep == 'single':
                        p_str_start = ''
                    else:
                        p_str_start = ''

                    #p_str = p_str_start + op_GetMozType(clc_tag_one[0]) + str(i) + op_GetMozValenceText(valence) + loss_text
                    if '.' in clc_tag_one[i - 1]:
                        if clc_tag_one[i - 1].split('.')[0][0:2] == 'by' or clc_tag_one[i - 1].split('.')[0][2:] == '0':
                            p_str = p_str_start + 'β ' + clc_tag_one[i - 1].split('.')[1]
                        else:
                            p_str = p_str_start + clc_tag_one[i - 1].split('.')[1]
                    else:
                        p_str = p_str_start + clc_tag_one[i - 1]
                    if yy > 80:
                        if xx in record_used_moz:
                            xx1_add = 10
                            yy1_add = 20
                        else:
                            xx1_add = 0
                            yy1_add = 5
                        ax3.text(xx + 10 + xx1_add, 80 + yy1_add, p_str, color=op_GetColor(moz_type), fontsize=6, ha='center',
                                 rotation='vertical')
                    else:
                        if xx in record_used_moz:
                            xx1_add = 4
                            yy2_add = 20
                        else:
                            xx1_add = 4
                            yy2_add = 10
                        ax3.text(xx + xx1_add, yy + yy2_add, p_str, color=op_GetColor(moz_type), fontsize=6, ha='center',
                                 rotation='vertical')
                        record_used_moz.append(xx)
                i = i + 1

            # 偏差图标注
            i = 1
            for xx, yy in zip(plot_moz, ppm):
                if plot_moz[i - 1] != 0:
                    if plot_peak_match[i - 1] != 0:
                        ax4.scatter(xx, yy, marker=',', color=op_GetColor(moz_type), s=1)
                i = i + 1
        # ax3.set_ylim([0, max_match_peak * 1.2])
        # ytick_ax3 = [max_peak * (i / 5) for i in range(6)]
        # ax3.set_yticks(ytick_ax3, ['0', '20', '40', '60', '80', '100'])
    else:
        pass
