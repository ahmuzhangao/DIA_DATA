import csv
import re
import numpy as np
from MSLogging import logGetError
import os
from MSTool import toolCountCharInString, toolGetWord
from MSSystem import UNIMOID_TO_STANDARD_MOD, CFG_TYPE_EXPORT, FRAGMENT_LOSS_TYPE_COMP, MARKER_SPLIT_MOD_LIST, \
    CFG_TYPE_ACCURACY_HALF_WIN_PEAK
from MSData import CFileMS1, CFileMS2, CINI, CEvidenceDIACheck, CSeed, CDataPack, Config, CFileID, CQuantInfo, CRTBetweenExp

VALUE_ILLEGAL = -7.16
VALUE_MAX_SCAN = 2000000

#这里，为了适应spectronaut的修饰形式:把_删去，把.删去，把[]转化为()，把()转化为[],再把空格去掉。（Oxidation (M) -> Oxidation[M]）
def op_convert_spectronaut(input_string):
    input_string = input_string.replace('_', '').replace('.', '').replace(' ', '')

    # 将方括号替换为临时符号（避免干扰替换）
    input_string = input_string.replace('[', '{').replace(']', '}')

    # 将圆括号替换为方括号
    input_string = input_string.replace('(', '[').replace(')', ']')

    # 将临时符号（原方括号）替换为圆括号
    input_string = input_string.replace('{', '(').replace('}', ')')

    return input_string

#这里，为了适应MaxDIA的修饰形式:把_删去
def op_convert_MaxDIA(input_string):

    input_string = input_string.replace('_', '')

    return input_string

#这里，为了适应peptdeep的修饰形式:首先替换 [ 为 @，然后删除 ](Oxidation[M]->Oxidation@M)
def op_convert_peptdeep(input_string):

    return re.sub(r'\]', '', re.sub(r'\[', '@', input_string))

def op_INIT_CFILE_MS1(dataMS1: CFileMS1):
    #初始化MS1的数据结构
    dataMS1.INDEX_SCAN = []
    dataMS1.INDEX_RT = []

    dataMS1.LIST_RET_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS1.LIST_ION_INJECTION_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN

    dataMS1.MATRIX_PEAK_INT = [[] * 1] * VALUE_MAX_SCAN
    dataMS1.MATRIX_PEAK_MOZ = [[] * 1] * VALUE_MAX_SCAN




def op_INIT_CFILE_MS2(dataMS2: CFileMS2):
    # 初始化MS2的数据结构
    dataMS2.INDEX_SCAN = []
    dataMS2.INDEX_RT = []

    dataMS2.LIST_RET_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS2.LIST_ION_INJECTION_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS2.LIST_PRECURSOR_SCAN = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    dataMS2.LIST_ACTIVATION_CENTER = [VALUE_ILLEGAL] * VALUE_MAX_SCAN

    dataMS2.MATRIX_PEAK_MOZ = [[] * 1] * VALUE_MAX_SCAN
    dataMS2.MATRIX_PEAK_INT = [[] * 1] * VALUE_MAX_SCAN

def op_preallocate_matrix(num_scans, start, end, step):
    """
    预分配一个二维矩阵作为存储空间。
    :param num_scans: 扫描总数
    :param start: mz 范围起始值
    :param end: mz 范围结束值
    :param step: mz 范围步长
    :return: 一个 NumPy 二维数组，大小为 (num_scans, num_bins, 2)
    """
    num_bins = int((end - start) // step)
    # 初始化矩阵，所有值为 -1
    return np.full((num_scans + 1 , num_bins, 2), -1, dtype=int)

def op_update_matrix(matrix, scan_idx, mzs, start, step):
    """
    更新二维矩阵中的范围索引。
    :param matrix: 预分配好的二维矩阵
    :param scan_idx: 当前扫描的索引（对应矩阵的第一维）
    :param mzs: 当前扫描的 mz 值数组
    :param start: mz 范围起始值
    :param step: mz 范围步长
    """
    # 计算所有 m/z 值的 range 索引
    range_indices = ((mzs - start) // step).astype(int)

    # 找到第一次出现的索引
    unique_bins, first_indices = np.unique(range_indices, return_index=True)

    # 找到最后一次出现的索引
    _, last_indices = np.unique(range_indices[::-1], return_index=True)
    last_indices = len(mzs) - 1 - last_indices  # 转换为正向索引

    # 更新矩阵
    matrix[scan_idx, unique_bins, 0] = first_indices
    matrix[scan_idx, unique_bins, 1] = last_indices


def op_FILL_LIST_PATH_MS(inputPath, inputList, inputExt):
    # 将输入的多个文件路径inputPath（以'|'为分割符），转化为多个文件路径的List
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


def op_FILL_LIST_PATH_ID(inputPath, inputList):
    # 将输入的多个文件路径inputPath（以'|'为分割符），转化为多个文件路径的List
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



def op_STANDARDIZE_MOD_BY_PRECURSOR(input_precursor: str, flag='()'):
    # 将DIA-NN带修饰的母离子转化为标准格式的修饰类型母离子:(UniMod:4)-> (Carbamidomethyl)
    if flag == '[]':
        seperator_start, seperator_end = '[', ']'
    else:
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

def op_STANDARDIZE_MOD_BY_PRECURSOR_SPECTORNAUT(input_precursor: str):
    # 将DIA-NN带修饰的母离子转化为标准格式的修饰类型母离子:(UniMod:4)-> (Carbamidomethyl)
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
                mod_type_standard = op_convert(mod_type)
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

def op_PSM2STR(fileID, i_PSM, isTitle, inputStyle):

    result = ''

    if inputStyle == CFG_TYPE_EXPORT['AppendFromID']:

        if isTitle:  # 写title
            # TODO 加新列
            result = result + fileID.PSM_LINE_T[0]
            result = result + '\tSupp_Info'
            result = result + '\tEmpty_Separator\t'

        else:
            result = result + fileID.PSM_LINE_C[i_PSM] + '\t'
            result = result + 'null\t'

    else:

        logGetError("MSOperator, MK227: "+ str(inputStyle))

    return result

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

def op_CAL_FRAGMENT_MOZ(input_INI: CINI, input_precursor: str):
    #顺序是：b1+ b2+ ... b1++ b2++ ... y1+ y2+ ... y1++ y2++ ...
    precursor_sequence, mod_dic = op_DIVIDE_MOD_FROM_PRECURSOR(input_precursor)
    ion_num = len(precursor_sequence)


    ionTypeList = []
    ionMozList = []
    chargMax = 2
    for j in range(chargMax):
        mass_b_aa_sum = 0.0
        for i in range(ion_num - 1):
            # b离子
            if (i + 1) in mod_dic.keys():
                mass_aa_mod = float(input_INI.DICT2_MOD_COM[mod_dic[i + 1]][1])  # 修饰质量
            else:
                mass_aa_mod = 0
            mass_aa = input_INI.DICT1_AA_COM[precursor_sequence[i]][1] + mass_aa_mod
            mass_b_aa_sum = mass_b_aa_sum + mass_aa
            index_key_b = 'b' + str(i + 1)
            charge = j + 1
            moz_b_ion = (mass_b_aa_sum + input_INI.MASS_PROTON_MONO * charge ) / charge
            if i >1:
                ionTypeList.append(index_key_b + '+' * (charge))
                ionMozList.append(moz_b_ion)

    for j in range(chargMax):
        mass_y_aa_sum = 0.0
        for i in range(ion_num - 1):
            # y离子
            if (ion_num - i) in mod_dic.keys():
                mass_aa_mod = float(input_INI.DICT2_MOD_COM[mod_dic[ion_num - i]][1])  # 修饰质量
            else:
                mass_aa_mod = 0
            mass_aa = input_INI.DICT1_AA_COM[precursor_sequence[ion_num - i - 1]][1] + mass_aa_mod
            mass_y_aa_sum = mass_y_aa_sum + mass_aa
            index_key_y = 'y' + str(i + 1)
            charge = j + 1
            moz_y_ion = (mass_y_aa_sum + input_INI.MASS_PROTON_MONO * charge + 18.0105633) / charge
            if i > 1:
                ionTypeList.append(index_key_y + '+' * (charge))
                ionMozList.append(moz_y_ion)

    return ionTypeList, ionMozList

def opGetStartAndEndForProfile(input_profile, input_seed, input_cutoff, input_n_hole):  # 就这个名字格式特殊

    nHoleLeft = 0
    nHoleRight = 0

    i_middle = input_seed

    if i_middle > 0:
        i_left = i_middle
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

def op_MOD2STR(fileID, strMod, isTitle, inputStyle):

    result = ''

    if inputStyle == CFG_TYPE_EXPORT['Simple'] or inputStyle == CFG_TYPE_EXPORT['AppendFromID']:

        if isTitle:

            # result = result + "Protein:Sequence:Modification:" + fileID.PSM6_GLC_T + '\t'
            result = result + "Protein" + '\t'
            result = result + "Sequence" + '\t'
            result = result + "Modification" + '\t'
            result = result + 'Empty_Separator\t'


        else:
            result = result + toolGetWord(strMod, 0, MARKER_SPLIT_MOD_LIST) + '\t'
            result = result + toolGetWord(strMod, 1, MARKER_SPLIT_MOD_LIST) + '\t'
            result = result + toolGetWord(strMod, 2, MARKER_SPLIT_MOD_LIST) + '\t'




    return result

def op_PRE2STR(fileID, strMod, isTitle, inputStyle):

    result = ''

    if inputStyle == CFG_TYPE_EXPORT['Simple'] or inputStyle == CFG_TYPE_EXPORT['AppendFromID']:

        if isTitle:

            # result = result + "Protein:Sequence:Modification:" + fileID.PSM6_GLC_T + '\t'
            result = result + "Protein" + '\t'
            result = result + "Sequence" + '\t'
            result = result + "Modification" + '\t'
            result = result + "Charge" + '\t'
            result = result + 'Empty_Separator\t'


        else:
            result = result + toolGetWord(strMod, 0, MARKER_SPLIT_MOD_LIST) + '\t'
            result = result + toolGetWord(strMod, 1, MARKER_SPLIT_MOD_LIST) + '\t'
            result = result + toolGetWord(strMod, 2, MARKER_SPLIT_MOD_LIST) + '\t'
            result = result + toolGetWord(strMod, 3, MARKER_SPLIT_MOD_LIST) + '\t'

    return result
def op_PRO2STR(fileID, iPRO, isTitle, inputStyle):

    result = ''

    if inputStyle == CFG_TYPE_EXPORT['Simple'] or inputStyle == CFG_TYPE_EXPORT['AppendFromID']:

        if isTitle:

            result = result + "ProteinName" + '\t'
            result = result + "GroupType" + '\t'
            result = result + "GroupID" + '\t'
            result = result + "Number_PSM_ID" + '\t'
            result = result + 'Empty_Separator\t'

        else:

            result = result + fileID.PRO4_LIST_AC[iPRO] + '\t'
            result = result + fileID.PRO7_LIST_TYPE[iPRO] + '\t'
            result = result + fileID.PRO6_LIST_ID[iPRO] + '\t'
            result = result + str(len(fileID.PRO1_DICT_PSM[fileID.PRO4_LIST_AC[iPRO]])) + '\t'

    return result

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

def op_DIVIDE_ION_TYPE_WITHLOSS(input_ion_type: str):
    # 拆分by离子类型（如：b2++_H2O)到[b/y类型，离子位点，离子电荷数]（如：[b,2,2])
    # 先拆分离子类型和损失
    ion_type_list = input_ion_type.split('_')
    ion_type = ion_type_list[0] if len(ion_type_list) == 1 else ion_type_list[0]
    return op_DIVIDE_ION_TYPE(ion_type)

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

def op_GetFont(moz_type, size):
    # 设置顶部离子标注字典
    color = op_GetColor(moz_type)
    font2 = {'family': 'SimHei',
             'color': color,
             'weight': 'normal',
             'size': size,
             }
    return font2

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

def op_INIT_CSEED_DIACHECK(CSeed):

    CSeed.MID_RT = VALUE_ILLEGAL
    CSeed.MID_SCAN = VALUE_ILLEGAL
    CSeed.DIS_ISO_MOZ_CLC = []  # 母离子理论质量
    CSeed.DIS_ISO_INT_CLC = []  # 母离子理论强度

    CSeed.DIS_FRA_MOZ_CLC = []
    CSeed.DIS_FRA_INT_CLC = []  # 留着，预计用Prosit预测理论强度
    CSeed.DIS_FRA_TYPE_CLC = []  # 例：b3++（目前考虑一价和二价）

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

def op_GENERATE_DECOY_FORDIA(inputSeed : CSeed):

    # list_choice = range(5,31,1)
    inputSeed.DIS_FRA_MOZ_CLC = [min(i_moz + 300.0, 1800) for i_moz in inputSeed.DIS_FRA_MOZ_CLC]
    inputSeed.DIS_ISO_MOZ_CLC = [min(i_moz + 300.0, 1800) for i_moz in inputSeed.DIS_ISO_MOZ_CLC]

    return inputSeed

def op_SAVE_MATRIX_TO_CSV(matrix, filename):
    with open(filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for row in matrix:
            csvwriter.writerow(row)


def op_INIT_CDataPack(inputDP: CDataPack):

    inputDP.LIST_PATH_MS1 = []
    inputDP.LIST_PATH_MS = []
    inputDP.LIST_PATH_MS2 = []
    inputDP.LIST_PATH_MS = []
    inputDP.LIST_PATH_ID = []

    inputDP.DIC_PSM_PAIR = {}
    inputDP.DIC_EVIDENCE = {}
    inputDP.DIC_EVIDENCE_DECOY = {}
    inputDP.DIC_PRECURSOR_TO_IONS = {}
    inputDP.NUM_OF_PAIRS = 0

    inputDP.myCFG = Config()
    inputDP.myID = CFileID()
    inputDP.myINI = CINI()
    inputDP.myQUANT = CQuantInfo()
    inputDP.myRTShift = CRTBetweenExp()

def op_ASSIGN_CDataPack(inputDP: CDataPack, outputDP: CDataPack):

    outputDP.LIST_PATH_MS1 = inputDP.LIST_PATH_MS1
    outputDP.LIST_PATH_ID = inputDP.LIST_PATH_ID
    outputDP.LIST_PATH_MS2 = inputDP.LIST_PATH_MS2
    outputDP.LIST_PATH_MS = inputDP.LIST_PATH_MS
    outputDP.DIC_PRECURSOR_TO_IONS = inputDP.DIC_PRECURSOR_TO_IONS

    # assign ID
    outputDP.myID.N_ID = inputDP.myID.N_ID
    outputDP.myID.ID0_SEQ = inputDP.myID.ID0_SEQ
    outputDP.myID.ID1_RAW_NAME = inputDP.myID.ID1_RAW_NAME
    outputDP.myID.ID2_SCAN_ID = inputDP.myID.ID2_SCAN_ID
    outputDP.myID.ID3_RT = inputDP.myID.ID3_RT
    outputDP.myID.ID4_SEQ_WITH_MOD = inputDP.myID.ID4_SEQ_WITH_MOD
    outputDP.myID.ID5_PRECURSOR_MOZ_EXP = inputDP.myID.ID5_PRECURSOR_MOZ_EXP
    outputDP.myID.ID6_PRECURSOR_MOZ_CLC = inputDP.myID.ID6_PRECURSOR_MOZ_CLC
    outputDP.myID.ID7_CHARGE = inputDP.myID.ID7_CHARGE
    outputDP.myID.ID8_PRECURSOR_ID = inputDP.myID.ID8_PRECURSOR_ID
    outputDP.myID.ID9_SCORE0 = inputDP.myID.ID9_SCORE0
    outputDP.myID.ID10_Protein_Groups = inputDP.myID.ID10_Protein_Groups
    outputDP.myID.ID11_Protein_Ids = inputDP.myID.ID11_Protein_Ids
    outputDP.myID.ID12_Protein_Names = inputDP.myID.ID12_Protein_Names
    outputDP.myID.ID13_Genes = inputDP.myID.ID13_Genes
    outputDP.myID.ID14_RT_BEGIN = inputDP.myID.ID14_RT_BEGIN
    outputDP.myID.ID15_RT_END = inputDP.myID.ID15_RT_END
    outputDP.myID.ID16_MOD = inputDP.myID.ID16_MOD
    outputDP.myID.ID17_SEQ_WITH_MOD_AND_CHARGE = inputDP.myID.ID17_SEQ_WITH_MOD_AND_CHARGE
    outputDP.myID.ID18_MOD_SITE = inputDP.myID.ID18_MOD_SITE
    outputDP.myID.ID19_MOD_TYPE = inputDP.myID.ID19_MOD_TYPE

    outputDP.myID.MOD1_DICT_PSM = inputDP.myID.MOD1_DICT_PSM
    outputDP.myID.PRE1_DICT_PSM = inputDP.myID.PRE1_DICT_PSM

    outputDP.myID.PRO1_DICT_PSM = inputDP.myID.PRO1_DICT_PSM
    outputDP.myID.PRO2_DICT_SEQ = inputDP.myID.PRO2_DICT_SEQ
    outputDP.myID.PRO3_DICT_N_SEQ = inputDP.myID.PRO3_DICT_N_SEQ
    outputDP.myID.PRO4_LIST_AC = inputDP.myID.PRO4_LIST_AC
    outputDP.myID.PRO5_LIST_DE = inputDP.myID.PRO5_LIST_DE
    outputDP.myID.PRO6_LIST_ID = inputDP.myID.PRO6_LIST_ID
    outputDP.myID.PRO7_LIST_TYPE = inputDP.myID.PRO7_LIST_TYPE

    # assign cfg
    outputDP.myCFG.I0_INI_PATH_ELEMENT = inputDP.myCFG.I0_INI_PATH_ELEMENT
    outputDP.myCFG.I1_INI_PATH_AA = inputDP.myCFG.I1_INI_PATH_AA
    outputDP.myCFG.I2_INI_PATH_MOD = inputDP.myCFG.I2_INI_PATH_MOD

    outputDP.myCFG.A0_FLAG_USE_EXIST_INDEX = inputDP.myCFG.A0_FLAG_USE_EXIST_INDEX
    outputDP.myCFG.A1_PATH_MS1 = inputDP.myCFG.A1_PATH_MS1
    outputDP.myCFG.A2_TYPE_MS1 = inputDP.myCFG.A2_TYPE_MS1
    outputDP.myCFG.A3_PATH_MS2 = inputDP.myCFG.A3_PATH_MS2
    outputDP.myCFG.A4_TYPE_MS2 = inputDP.myCFG.A4_TYPE_MS2
    outputDP.myCFG.A5_PATH_FILE = inputDP.myCFG.A5_PATH_FILE
    outputDP.myCFG.A6_TYPE_FILE = inputDP.myCFG.A6_TYPE_FILE

    outputDP.myCFG.B1_PATH_IDENTIFICATION_RESULT = inputDP.myCFG.B1_PATH_IDENTIFICATION_RESULT
    outputDP.myCFG.B2_TYPE_IDENTIFICATION_RESULT = inputDP.myCFG.B2_TYPE_IDENTIFICATION_RESULT
    outputDP.myCFG.B3_THRESHOLD_FDR = inputDP.myCFG.B3_THRESHOLD_FDR
    outputDP.myCFG.B4_PREDICTMS2_TYPE_MODEL = inputDP.myCFG.B4_PREDICTMS2_TYPE_MODEL

    outputDP.myCFG.C0_TYPE_FLOW = inputDP.myCFG.C0_TYPE_FLOW
    outputDP.myCFG.C1_DIA_RT_HALF_WIN_IN_MIN = inputDP.myCFG.C1_DIA_RT_HALF_WIN_IN_MIN
    outputDP.myCFG.C2_DIA_ACCURACY_HALF_WIN_PEAK = inputDP.myCFG.C2_DIA_ACCURACY_HALF_WIN_PEAK
    outputDP.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK = inputDP.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK
    outputDP.myCFG.C4_DIALF_TYPE_MODEL = inputDP.myCFG.C4_DIALF_TYPE_MODEL
    outputDP.myCFG.C5_DIALF_PATH_CALIBRATE_RT = inputDP.myCFG.C5_DIALF_PATH_CALIBRATE_RT
    outputDP.myCFG.C6_FLAG_DECOY = inputDP.myCFG.C6_FLAG_DECOY
    outputDP.myCFG.C7_PROCESS_NUM = inputDP.myCFG.C7_PROCESS_NUM

    outputDP.myCFG.D1_PATH_EXPORT = inputDP.myCFG.D1_PATH_EXPORT
    outputDP.myCFG.D2_TYPE_EXPORT = inputDP.myCFG.D2_TYPE_EXPORT

    # assign INI
    outputDP.myINI.MASS_ELECTRON = inputDP.myINI.MASS_ELECTRON
    outputDP.myINI.MASS_PROTON_MONO = inputDP.myINI.MASS_PROTON_MONO
    outputDP.myINI.MASS_PROTON_ARVG = inputDP.myINI.MASS_PROTON_ARVG

    outputDP.myINI.DICT0_ELEMENT_MASS = inputDP.myINI.DICT0_ELEMENT_MASS
    outputDP.myINI.DICT0_ELEMENT_ABDC = inputDP.myINI.DICT0_ELEMENT_ABDC

    outputDP.myINI.DICT1_AA_COM = inputDP.myINI.DICT1_AA_COM
    outputDP.myINI.DICT2_MOD_COM = inputDP.myINI.DICT2_MOD_COM

    # assign quant
    outputDP.myQUANT.MATRIX_QUANT_INT_PSM = inputDP.myQUANT.MATRIX_QUANT_INT_PSM
    outputDP.myQUANT.MATRIX_QUANT_INT_MOD = inputDP.myQUANT.MATRIX_QUANT_INT_MOD
    outputDP.myQUANT.MATRIX_QUANT_INT_PRE = inputDP.myQUANT.MATRIX_QUANT_INT_PRE
    outputDP.myQUANT.MATRIX_MATCH_SCORE = inputDP.myQUANT.MATRIX_MATCH_SCORE

    outputDP.myQUANT.LIST_FLAG_INFER_PSM = inputDP.myQUANT.LIST_FLAG_INFER_PSM
    outputDP.myQUANT.MATRIX_RT_LENGTH_PSM = inputDP.myQUANT.MATRIX_RT_LENGTH_PSM
    outputDP.myQUANT.MATRIX_RT_START_PSM = inputDP.myQUANT.MATRIX_RT_START_PSM
    outputDP.myQUANT.MATRIX_RT_MID_PSM = inputDP.myQUANT.MATRIX_RT_MID_PSM

    # assign rt_info
    outputDP.myRTShift.LIST_EXPERIMENT = inputDP.myRTShift.LIST_EXPERIMENT
    outputDP.myRTShift.DIC_EXPERIMENT_INDEX = inputDP.myRTShift.DIC_EXPERIMENT_INDEX

    outputDP.myRTShift.DIC_MBR_WIDTH = inputDP.myRTShift.DIC_MBR_WIDTH
    outputDP.myRTShift.DIC_RT_SHIFT = inputDP.myRTShift.DIC_RT_SHIFT
    outputDP.myRTShift.DIC_RT_CORR = inputDP.myRTShift.DIC_RT_CORR