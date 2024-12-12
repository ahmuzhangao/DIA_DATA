import os
import pickle
import re

from MSLogging import logGetError
from MSData import Config
from MSTool import toolGetIndexByWord, toolGetWord, toolGetNameFromPath, toolStr2List, toolGetWord1, \
    toolCountCharInString
from MSData import CDataPack, CFileMS1, CFileMS2
from MSOperator import op_INIT_CFILE_MS1, op_INIT_CFILE_MS2, \
    op_STANDARDIZE_MOD_BY_PRECURSOR, op_STANDARDIZE_MOD_MSFraggerPin_BY_PRECURSOR
from MSSystem import VALUE_ILLEGAL, UNIMOID_TO_STANDARD_MOD


class CFunctionConfig:

    def config2file(self, path, config):

        with open(path, 'w') as f:
            f.write("#[Ini files]\n")
            f.write("INI_PATH_ELEMENT=" + config.I0_INI_PATH_ELEMENT + '\n')
            f.write("INI_PATH_AA=" + config.I1_INI_PATH_AA + '\n')
            f.write("INI_PATH_MOD=" + config.I2_INI_PATH_MOD + '\n')

            f.write("\n#[Data]\n")
            f.write("PATH_MS1=" + config.A1_PATH_MS1 + '\n')
            f.write("TYPE_MS1=" + str(config.A2_TYPE_MS1) + '\n')
            f.write("PATH_MS2=" + config.A3_PATH_MS2 + '\n')
            f.write("TYPE_MS2=" + str(config.A4_TYPE_MS2) + '\n')

            f.write("\n#[Identification results]\n")
            f.write("PATH_IDENTIFICATION_RESULT=" + config.B1_PATH_IDENTIFICATION_RESULT + '\n')
            f.write("TYPE_IDENTIFICATION_RESULT=" + str(config.B2_TYPE_IDENTIFICATION_RESULT) + '\n')
            f.write("THRESHOLD_FDR=" + str(config.B3_THRESHOLD_FDR) + '\n')

            f.write('\n# [Curve Visualization]\n')
            f.write('TYPE_FLOW=' + str(config.C0_TYPE_FLOW) + '\n')
            f.write('DIA_RT_HALF_WIN_IN_MIN=' + str(config.C1_DIA_RT_HALF_WIN_IN_MIN) + '\n')
            f.write('DIA_ACCURACY_HALF_WIN_PEAK=' + str(config.C2_DIA_ACCURACY_HALF_WIN_PEAK) + '\n')
            f.write('DIA_TYPE_ACCURACY_HALF_WIN_PEAK=' + str(config.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK) + '\n')
            f.write("CURVE_SMOOTH=" + str(config.C14_CURVE_SMOOTH) + '\n')

            f.write('\n#[Filter]\n')
            f.write('TYPE_PLOT_NUM=' + str(config.C7_PLOT_NUM) + '\n')
            f.write('DIA_PRECURSOR_RAW=' + str(config.C6_DIA_PRECURSOR_RAW) + '\n')
            f.write('DIA_PRECURSOR_SEQUENCE=' + str(config.C4_DIA_PRECURSOR_SEQUENCE) + '\n')
            f.write('DIA_PRECURSOR_CHARGE=' + str(config.C5_DIA_PRECURSOR_CHARGE) + '\n')
            # f.write('DDA_SCAN_NO=' + '\n')
            # f.write('DDA_RAW=' + '\n')
            f.write('DDA_ISOTOPE=' + '\n')
            f.write("\n#[Export]\n")
            f.write("PATH_EXPORT=" + str(config.D1_PATH_EXPORT) + '\n')

    def file2config(self, path, config: Config):
        with open(path, 'r') as f:
            for line in f.readlines():
                if line.startswith("#"):
                    continue
                p_EqualSign = line.find('=')
                if -1 == p_EqualSign:
                    pass
                else:
                    subLine = toolGetWord(line, 0, ';')  # ;后面的是注释
                    self.__soldierParseLine(subLine, config)

    def __soldierParseLine(self, line, cfg):

        str_name = toolGetWord(line, 0, '=')
        str_value = toolGetWord(line, 1, '=').replace("\n", "")

        if "PATH_MS1" == str_name:

            cfg.A1_PATH_MS1 = str_value

        elif "TYPE_MS1" == str_name:

            cfg.A2_TYPE_MS1 = int(str_value)

        elif "PATH_MS2" == str_name:

            cfg.A3_PATH_MS2 = str_value

        elif "TYPE_MS2" == str_name:

            cfg.A4_TYPE_MS2 = int(str_value)

        elif "PATH_IDENTIFICATION_RESULT" == str_name:

            cfg.B1_PATH_IDENTIFICATION_RESULT = str_value

        elif "TYPE_IDENTIFICATION_RESULT" == str_name:

            cfg.B2_TYPE_IDENTIFICATION_RESULT = int(str_value)

        elif "THRESHOLD_FDR" == str_name:

            cfg.B3_THRESHOLD_FDR = float(str_value)

        elif 'PATH_LIBRARY_RESULT' == str_name:

            cfg.B4_PATH_LIBRARY_RESULT = str_value

        elif "TYPE_FLOW" == str_name:

            cfg.C0_TYPE_FLOW = int(str_value)

        elif 'DIA_RT_HALF_WIN_IN_MIN' == str_name:

            cfg.C1_DIA_RT_HALF_WIN_IN_MIN = float(str_value)

        elif 'DIA_ACCURACY_HALF_WIN_PEAK' == str_name:

            cfg.C2_DIA_ACCURACY_HALF_WIN_PEAK = float(str_value)

        elif 'DIA_TYPE_ACCURACY_HALF_WIN_PEAK' == str_name:

            cfg.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK = int(str_value)

        elif 'DIA_PRECURSOR_SEQUENCE' == str_name:

            cfg.C4_DIA_PRECURSOR_SEQUENCE = str_value  # 此处不进行操作，保存为字符串，画图时再将多个序列分开

        elif 'DIA_PRECURSOR_CHARGE' == str_name:

            cfg.C5_DIA_PRECURSOR_CHARGE = str_value  # 此处不进行操作，保存为字符串，画图时再将多个电荷分开

        elif 'DDA_SCAN_NO' == str_name:

            cfg.C8_DDA_PRECURSOR_SCAN = str_value  # 此处不进行操作，保存为字符串，画图时再将多个SCAN号分开

        elif 'DDA_RAW' == str_name:

            cfg.C9_DDA_PRECURSOR_RAW = str_value  # 此处不进行操作，保存为字符串，画图时再分开

        elif 'DDA_ISOTOPE' == str_name:

            cfg.C10_DDA_ISOTOPE = str_value

        elif 'DDA_PRECURSOR_SEQUENCE' == str_name:

            cfg.C11_DDA_PRECURSOR_SEQUENCE = str_value

        elif 'DDA_PRECURSOR_CHARGE' == str_name:

            cfg.C12_DDA_PRECURSOR_CHARGE = str_value

        elif 'DDA_RT' == str_name:

            cfg.C13_DDA_RT = str_value

        elif 'CURVE_SMOOTH' == str_name:

            cfg.C14_CURVE_SMOOTH = int(str_value)

        elif 'TOP_INT_NUM' == str_name:

            cfg.C15_INT_NUM = int(str_value)

        elif 'DRAW_FLAG' == str_name:

            cfg.C17_DRAW_FLAG = int(str_value)  # 是否要画图，画图会慢很多，不画就导出色谱曲线矩阵

        elif 'FRAG_NUM' == str_name:

            cfg.C18_FRAG_NUM = int(str_value)

        elif 'PATH_EXPORT' == str_name:

            cfg.D1_PATH_EXPORT = str_value

        elif "INI_PATH_ELEMENT" == str_name:

            cfg.I0_INI_PATH_ELEMENT = str_value

        elif "INI_PATH_AA" == str_name:

            cfg.I1_INI_PATH_AA = str_value

        elif "INI_PATH_MOD" == str_name:

            cfg.I2_INI_PATH_MOD = str_value

        elif "DIA_PRECURSOR_RAW" == str_name:

            cfg.C6_DIA_PRECURSOR_RAW = str_value

        elif 'TYPE_PLOT_NUM' == str_name:

            cfg.C7_PLOT_NUM = str_value
        else:
            print(str_name)
            info = "MSFunction, MK449, " + str_name + " is all right?"
            logGetError(info)


class CFunctionINI:

    def __init__(self, inputDP):
        self.dp = inputDP

    def file2ini(self):

        self.__captainFile2Element(self.dp.myCFG.I0_INI_PATH_ELEMENT)
        self.__captainFile2AA(self.dp.myCFG.I1_INI_PATH_AA)
        self.__captainFile2Mod(self.dp.myCFG.I2_INI_PATH_MOD)

    def ini2file(self):

        self.__captainElement2File(self.dp.myCFG.I0_INI_PATH_ELEMENT)

    def __captainFile2Mod(self, path):

        with open(path, 'r', encoding='utf8') as f:
            for line in f.readlines():

                if len(line) > 1:

                    if line[-1] == '\n':
                        line = line[0:-1]  # 把最后一个\n干掉

                    if line.startswith("@"):
                        continue

                    if line.startswith("name"):

                        str_name = toolGetWord1(line, '=', ' ')

                    else:

                        nBlank = toolCountCharInString(line, ' ')

                        if nBlank > 5:

                            str_comp = toolGetWord(line, 7, ' ')
                            str_mass = toolGetWord(line, 2, ' ')

                        else:

                            str_comp = toolGetWord(line, 5, ' ')
                            str_mass = toolGetWord(line, 2, ' ')

                        self.dp.myINI.DICT2_MOD_COM[str_name] = (str_comp, str_mass)

    def __captainElement2File(self, path):

        pass

    def __captainFile2AA(self, path):

        with open(path, 'r', encoding='utf8') as f:

            for line in f.readlines():

                if len(line) > 1:
                    str_name = toolGetWord(line, 0, '|')
                    str_comp = toolGetWord(line, 1, '|')
                    str_mass = self.__CalculateAminoacidMass(str_comp, self.dp.myINI.DICT0_ELEMENT_MASS)

                    self.dp.myINI.DICT1_AA_COM[str_name] = (str_comp, str_mass)

    def __captainFile2Element(self, path):

        with open(path, 'r', encoding='utf8') as f:

            for line in f.readlines():

                if len(line) > 1:
                    str_name = toolGetWord(line, 0, '|')
                    str_mass = toolGetWord(line, 1, '|')
                    str_abdc = toolGetWord(line, 2, '|')

                    list_mass = toolStr2List(str_mass, ',')
                    list_abdc = toolStr2List(str_abdc, ',')

                    self.dp.myINI.DICT0_ELEMENT_MASS[str_name] = list_mass
                    self.dp.myINI.DICT0_ELEMENT_ABDC[str_name] = list_abdc

    def __CalculateAminoacidMass(self, inputMolecule, inputAtonMass):  # 计算氨基酸相对分子质量
        index = 0
        molecule_mass = 0
        atom_left = 0
        for i in inputMolecule:
            if (i == '(') | (i == '（'):
                atom_right = index
                num_left = index + 1
                atom = inputMolecule[atom_left: atom_right]
            elif (i == ')') | (i == '）'):
                atom_left = index + 1
                num_right = index
                atom_num = int(inputMolecule[num_left: num_right])
                molecule_mass = molecule_mass + atom_num * inputAtonMass[atom][0]
            index = index + 1
        return molecule_mass


class CFunctionParseIDForDIANN:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soldierParseMOD(self, inputStr):  # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理

        return op_STANDARDIZE_MOD_BY_PRECURSOR(inputStr)

    def read(self, pathID):

        with open(pathID, 'r') as f:
            content = f.readline()
            lines = f.readlines()

        i_RAW_NAME = toolGetIndexByWord(content, 'FileName', '\t')
        i_RT = toolGetIndexByWord(content, 'Tr_recalibrated', '\t')

        i_SEQ_WITH_MOD = toolGetIndexByWord(content, 'ModifiedPeptide', '\t')
        i_SEQ = toolGetIndexByWord(content, 'PeptideSequence', '\t')
        i_PRECURSOR_MOZ_EXP = toolGetIndexByWord(content, 'PrecursorMz', '\t')
        i_CHARGE = toolGetIndexByWord(content, 'PrecursorCharge', '\t')
        i_PRECURSOR = toolGetIndexByWord(content, 'transition_group_id', '\t')

        i_SCORE0 = toolGetIndexByWord(content, 'QValue', '\t')

        i_FRAGMENT_TYPE_1 = toolGetIndexByWord(content, 'FragmentType', '\t')
        i_FRAGMENT_TYPE_2 = toolGetIndexByWord(content, 'FragmentSeriesNumber', '\t')
        i_FRAGMENT_TYPE_3 = toolGetIndexByWord(content, 'FragmentCharge', '\t')

        i_FRAGMENT_MOZ_EXP = toolGetIndexByWord(content, 'ProductMz', '\t')
        i_FRAGMENT_LOSS_TYPE = toolGetIndexByWord(content, 'FragmentLossType', '\t')
        pre_transition_group_id = '-1'  # 用于判断是否读到下一个transition_group
        fragment_moz_list = []  # 暂存当前母离子对应子离子m/z的列表
        fragment_type_list = []  # 暂存当前母离子对应子离子类型的列表
        fragment_loss_list = []  # 暂存当前母离子对应子离子的损失类型的列表

        for line in lines:

            tmpFDR = float(toolGetWord(line, i_SCORE0, '\t'))

            if tmpFDR <= self.dp.myCFG.B3_THRESHOLD_FDR:

                transition_group_id = toolGetWord(line, i_RAW_NAME, '\t') + toolGetWord(line, i_PRECURSOR, '\t')

                if pre_transition_group_id != transition_group_id:

                    pre_transition_group_id = transition_group_id

                    if len(fragment_moz_list) != 0:
                        self.dp.myID.ID10_FRAGMENT_TYPE.append(fragment_type_list)
                        self.dp.myID.ID11_FRAGMENT_MOZ_EXP.append(fragment_moz_list)
                        self.dp.myID.ID13_FRAGMENT_LOSS_TYPE.append(fragment_loss_list)
                        self.dp.myID.ID12_FRAGMENT_MOZ_CLC.append([])
                        self.dp.myID.ID14_RT_BEGIN.append(VALUE_ILLEGAL)
                        self.dp.myID.ID15_RT_END.append(VALUE_ILLEGAL)

                    fragment_type_list = []
                    fragment_moz_list = []
                    fragment_loss_list = []
                    fragment_type_1 = toolGetWord(line, i_FRAGMENT_TYPE_1, '\t')
                    fragment_type_2 = toolGetWord(line, i_FRAGMENT_TYPE_2, '\t')
                    fragment_type_3 = toolGetWord(line, i_FRAGMENT_TYPE_3, '\t')
                    fragment_type_list.append(fragment_type_1 + fragment_type_2 + '+' * int(fragment_type_3))
                    fragment_moz_list.append(float(toolGetWord(line, i_FRAGMENT_MOZ_EXP, '\t')))
                    loss_type = toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t')
                    if loss_type == 'noloss':
                        fragment_loss_list.append('NOLOSS')  # 格式与system中LOSS_TYPE保持一致
                    else:
                        fragment_loss_list.append(toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t'))

                    self.dp.myID.ID1_RAW_NAME.append(toolGetWord(line, i_RAW_NAME, '\t'))
                    self.dp.myID.ID2_SCAN_ID.append(VALUE_ILLEGAL)
                    self.dp.myID.ID3_RT.append(float(toolGetWord(line, i_RT, '\t')))

                    self.dp.myID.ID4_SEQ_WITH_MOD.append(
                        self.__soldierParseMOD(toolGetWord(line, i_SEQ_WITH_MOD, '\t')))
                    self.dp.myID.ID0_SEQ.append(toolGetWord(line, i_SEQ, '\t'))
                    self.dp.myID.ID17_TARGET.append(1)
                    self.dp.myID.ID5_PRECURSOR_MOZ_EXP.append(float(toolGetWord(line, i_PRECURSOR_MOZ_EXP, '\t')))
                    self.dp.myID.ID6_PRECURSOR_MOZ_CLC.append(VALUE_ILLEGAL)
                    self.dp.myID.ID7_CHARGE.append(int(toolGetWord(line, i_CHARGE, '\t')))

                    self.dp.myID.ID8_PRECURSOR_ID.append(toolGetNameFromPath(toolGetWord(line, i_RAW_NAME, '\t')) +
                                                         '_' + self.__soldierParseMOD(
                        toolGetWord(line, i_PRECURSOR, '\t')))

                    self.dp.myID.ID9_SCORE0.append(tmpFDR)
                    self.dp.myID.N_ID = self.dp.myID.N_ID + 1
                else:

                    fragment_type_1 = toolGetWord(line, i_FRAGMENT_TYPE_1, '\t')
                    fragment_type_2 = toolGetWord(line, i_FRAGMENT_TYPE_2, '\t')
                    fragment_type_3 = toolGetWord(line, i_FRAGMENT_TYPE_3, '\t')
                    fragment_type_list.append(fragment_type_1 + fragment_type_2 + '+' * int(fragment_type_3))
                    fragment_moz_list.append(float(toolGetWord(line, i_FRAGMENT_MOZ_EXP, '\t')))
                    loss_type = toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t')
                    if loss_type == 'noloss':
                        fragment_loss_list.append('NOLOSS')
                    else:
                        fragment_loss_list.append(toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t'))

        if len(fragment_moz_list) != 0:
            self.dp.myID.ID10_FRAGMENT_TYPE.append(fragment_type_list)
            self.dp.myID.ID11_FRAGMENT_MOZ_EXP.append(fragment_moz_list)
            self.dp.myID.ID12_FRAGMENT_MOZ_CLC.append([])
            self.dp.myID.ID13_FRAGMENT_LOSS_TYPE.append(fragment_loss_list)


class CFunctionParseIDForMSFragger:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soldierParseMOD(self, inputStr):  # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理

        return op_STANDARDIZE_MOD_BY_PRECURSOR(inputStr)

    def read(self, pathID):

        with open(pathID, 'r') as f:
            content = f.readline()
            lines = f.readlines()

        i_RT = toolGetIndexByWord(content, 'AverageExperimentalRetentionTime', '\t')

        i_SEQ_WITH_MOD = toolGetIndexByWord(content, 'ModifiedPeptideSequence', '\t')
        i_SEQ = toolGetIndexByWord(content, 'PeptideSequence', '\t')
        i_PRECURSOR_MOZ_EXP = toolGetIndexByWord(content, 'PrecursorMz', '\t')
        i_CHARGE = toolGetIndexByWord(content, 'PrecursorCharge', '\t')
        # i_SCORE0 = toolGetIndexByWord(content, 'QValue', '\t')

        i_FRAGMENT_TYPE_1 = toolGetIndexByWord(content, 'FragmentType', '\t')
        i_FRAGMENT_TYPE_3 = toolGetIndexByWord(content, 'FragmentCharge', '\t')
        i_FRAGMENT_TYPE_2 = toolGetIndexByWord(content, 'FragmentSeriesNumber', '\t')
        i_FRAGMENT_MOZ_EXP = toolGetIndexByWord(content, 'ProductMz', '\t')
        i_FRAGMENT_LOSS_TYPE = toolGetIndexByWord(content, 'Annotation', '\t')
        pre_transition_group_id = '-1'  # 用于判断是否读到下一个transition_group
        fragment_moz_list = []  # 暂存当前母离子对应子离子m/z的列表
        fragment_type_list = []  # 暂存当前母离子对应子离子类型的列表
        fragment_loss_list = []  # 暂存当前母离子对应子离子的损失类型的列表

        for line in lines:

            line = line.strip('\n')
            # tmpFDR = float(toolGetWord(line, i_SCORE0, '\t'))
            tmpFDR = 0

            if tmpFDR <= self.dp.myCFG.B3_THRESHOLD_FDR:

                transition_group_id = toolGetWord(line, i_SEQ_WITH_MOD, '\t') + '-' + toolGetWord(line, i_CHARGE, '\t')

                if pre_transition_group_id != transition_group_id:

                    pre_transition_group_id = transition_group_id

                    if len(fragment_moz_list) != 0:
                        self.dp.myID.ID10_FRAGMENT_TYPE.append(fragment_type_list)
                        self.dp.myID.ID11_FRAGMENT_MOZ_EXP.append(fragment_moz_list)
                        self.dp.myID.ID13_FRAGMENT_LOSS_TYPE.append(fragment_loss_list)
                        self.dp.myID.ID12_FRAGMENT_MOZ_CLC.append([])
                        self.dp.myID.ID14_RT_BEGIN.append(VALUE_ILLEGAL)
                        self.dp.myID.ID15_RT_END.append(VALUE_ILLEGAL)

                    fragment_type_list = []
                    fragment_moz_list = []
                    fragment_loss_list = []
                    fragment_type_1 = toolGetWord(line, i_FRAGMENT_TYPE_1, '\t')
                    fragment_type_2 = toolGetWord(line, i_FRAGMENT_TYPE_2, '\t')
                    fragment_type_3 = toolGetWord(line, i_FRAGMENT_TYPE_3, '\t')
                    fragment_type_list.append(fragment_type_1 + fragment_type_2 + '+' * int(fragment_type_3))
                    fragment_moz_list.append(float(toolGetWord(line, i_FRAGMENT_MOZ_EXP, '\t')))
                    loss_type = toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t')
                    if '-' not in loss_type:
                        fragment_loss_list.append('NOLOSS')
                    else:
                        fragment_loss_list.append(toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t'))
                    raw_name = self.dp.myCFG.A1_PATH_MS1.split('\\')[-1].replace('.ms1', '')
                    self.dp.myID.ID1_RAW_NAME.append(raw_name)
                    self.dp.myID.ID2_SCAN_ID.append(VALUE_ILLEGAL)

                    self.dp.myID.ID3_RT.append(float(toolGetWord(line, i_RT, '\t')) / 60)
                    self.dp.myID.ID4_SEQ_WITH_MOD.append(
                        self.__soldierParseMOD(toolGetWord(line, i_SEQ_WITH_MOD, '\t')))
                    self.dp.myID.ID0_SEQ.append(toolGetWord(line, i_SEQ, '\t'))
                    self.dp.myID.ID5_PRECURSOR_MOZ_EXP.append(float(toolGetWord(line, i_PRECURSOR_MOZ_EXP, '\t')))
                    self.dp.myID.ID6_PRECURSOR_MOZ_CLC.append(VALUE_ILLEGAL)
                    self.dp.myID.ID7_CHARGE.append(int(toolGetWord(line, i_CHARGE, '\t')))

                    self.dp.myID.ID8_PRECURSOR_ID.append(raw_name +
                                                         '_' + self.__soldierParseMOD(
                        toolGetWord(line, i_SEQ_WITH_MOD, '\t')))
                    # 修饰需要处理一下

                    self.dp.myID.ID17_TARGET.append(1)

                    self.dp.myID.ID9_SCORE0.append(tmpFDR)
                    self.dp.myID.N_ID = self.dp.myID.N_ID + 1
                else:

                    fragment_type_1 = toolGetWord(line, i_FRAGMENT_TYPE_1, '\t')
                    fragment_type_2 = toolGetWord(line, i_FRAGMENT_TYPE_2, '\t')
                    fragment_type_3 = toolGetWord(line, i_FRAGMENT_TYPE_3, '\t')
                    fragment_type_list.append(fragment_type_1 + fragment_type_2 + '+' * int(fragment_type_3))
                    fragment_moz_list.append(float(toolGetWord(line, i_FRAGMENT_MOZ_EXP, '\t')))
                    loss_type = toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t')
                    if '-' not in loss_type:
                        fragment_loss_list.append('NOLOSS')
                    else:
                        fragment_loss_list.append(toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t'))

        if len(fragment_moz_list) != 0:
            self.dp.myID.ID10_FRAGMENT_TYPE.append(fragment_type_list)
            self.dp.myID.ID11_FRAGMENT_MOZ_EXP.append(fragment_moz_list)
            self.dp.myID.ID12_FRAGMENT_MOZ_CLC.append([])
            self.dp.myID.ID13_FRAGMENT_LOSS_TYPE.append(fragment_loss_list)


class CFunctionParseIDForMSFraggerPin:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soldierParseMOD(self, inputStr):  # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理

        inputStr = re.sub(r'[a-z]', '', inputStr)

        return op_STANDARDIZE_MOD_MSFraggerPin_BY_PRECURSOR(inputStr)

    def read(self, pathID):

        with open(pathID, 'r') as f:
            content = f.readline()
            lines = f.readlines()

        i_RAW_NAME = toolGetIndexByWord(content, 'SpecId', '\t')
        i_SCAN = toolGetIndexByWord(content, 'ScanNr', '\t')
        i_RT = toolGetIndexByWord(content, 'retentiontime', '\t')


        i_SEQ_WITH_MOD = toolGetIndexByWord(content, 'Peptide', '\t')
        i_SEQ = toolGetIndexByWord(content, 'Peptide', '\t')
        # i_PRECURSOR_MOZ_EXP = toolGetIndexByWord(content, 'PrecursorMz', '\t')

        i_CHARGE_1 = toolGetIndexByWord(content, 'charge_1', '\t')
        i_CHARGE_7 = toolGetIndexByWord(content, 'charge_7_or_more', '\t')

        i_DECOY_FLAG = toolGetIndexByWord(content, 'Label', '\t')

        pre_transition_group_id = '-1'  # 用于判断是否读到下一个transition_group

        for line in lines:

            # tmpFDR = float(toolGetWord(line, i_SCORE0, '\t'))

            tmpFDR = 0

            transition_group_id = toolGetWord(line, i_SEQ_WITH_MOD, '\t')
            charge = 0
            for i in range(i_CHARGE_1, i_CHARGE_7):
                charge += 1
                if toolGetWord(line, i, '\t') == '1':
                    break
            if pre_transition_group_id != transition_group_id:

                pre_transition_group_id = transition_group_id

                raw_name = toolGetWord(line, i_RAW_NAME, '\t').split('.')[0] + '.raw'
                self.dp.myID.ID1_RAW_NAME.append(raw_name)
                self.dp.myID.ID2_SCAN_ID.append(int(toolGetWord(line, i_SCAN, '\t')))
                self.dp.myID.ID3_RT.append(float(toolGetWord(line, i_RT, '\t')) / 60)

                self.dp.myID.ID17_TARGET.append(int(toolGetWord(line, i_DECOY_FLAG, '\t')))

                seq_with_mod = toolGetWord(line, i_SEQ_WITH_MOD, '\t')
                mod_dict = self.__captainParseMOD(seq_with_mod[2:-3])
                seq = ''.join(filter(str.isupper, transition_group_id[2:-3]))
                precursor_moz_clc = self.__captainCalMOZ(seq, mod_dict, charge)
                a =  self.__soldierParseMOD(toolGetWord(line, i_SEQ_WITH_MOD, '\t')[2:-3])
                self.dp.myID.ID4_SEQ_WITH_MOD.append(
                    self.__soldierParseMOD(toolGetWord(line, i_SEQ_WITH_MOD, '\t')[2:-3]))
                self.dp.myID.ID0_SEQ.append(seq)
                self.dp.myID.ID16_MOD.append(mod_dict)
                self.dp.myID.ID5_PRECURSOR_MOZ_EXP.append(precursor_moz_clc)
                self.dp.myID.ID6_PRECURSOR_MOZ_CLC.append(VALUE_ILLEGAL)
                self.dp.myID.ID7_CHARGE.append(charge)

                self.dp.myID.ID8_PRECURSOR_ID.append(raw_name +
                                                     '_' + self.__soldierParseMOD(
                    toolGetWord(line, i_SEQ_WITH_MOD, '\t') + str(charge)))
                # 修饰需要处理一下
                self.dp.myID.ID9_SCORE0.append(tmpFDR)
                self.dp.myID.N_ID = self.dp.myID.N_ID + 1


    def __captainParseMOD(self, inputStr):  # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
        # 从DIA-NN格式的带修饰母离子字符串中提取修饰列表，每个元素存放修饰位点和类型
        trantab = str.maketrans('()', '[]')
        seperator_start, seperator_end = '[', ']'
        mod_dic = {}
        mod_type = ''
        site_start, site_end = 0, 0
        flag_in_mod = 0
        count = 0
        if 'n' in inputStr:
            inputStr.replace('n', '')
        for i, char in enumerate(inputStr):
            if char == seperator_start:
                site_start = i
                flag_in_mod = 1
            elif char == seperator_end:
                site_end = i
                try:
                    mod_type = inputStr[site_start + 1:site_end].translate(trantab).replace(' ', '')
                except KeyError:
                    info = 'DIA-NN modify type has wrong, please check!'
                    logGetError(info)
                mod_dic[count-1] = UNIMOID_TO_STANDARD_MOD[mod_type]
                flag_in_mod = 0
            elif flag_in_mod == 0:
                count = count + 1
        return mod_dic

    def __captainCalMOZ(self, inputSeq: str, inputMod: dict, inputCharge: int):

        # cal MOD mass
        mod_mass = 0.0
        for i in list(inputMod.values()):
            mod_mass += float(self.dp.myINI.DICT2_MOD_COM[i][1])
        # cal seq mass
        seq_mass = 0.0
        for i in inputSeq:
            seq_mass += self.dp.myINI.DICT1_AA_COM[i][1]
        precursor_mass = seq_mass + 18.0105633 + mod_mass + inputCharge * self.dp.myINI.MASS_PROTON_MONO
        precursor_moz = precursor_mass / inputCharge

        return precursor_moz



class CFunctionParseIDForMaxDIA:
    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __captainParseMOD(self, inputMod: str, inputLength: int) -> str:
        inputDicFixMod = {'C': ['Carbamidomethyl[C]']}
        dic_FixMOD = inputDicFixMod
        str_ModSequence = inputMod.strip('_')
        str_ModSequence = str_ModSequence.replace(' ', '')
        str_ModSequence = str_ModSequence.replace('(', '[')
        str_ModSequence = str_ModSequence.replace(')', ']')
        MOD_Standard = ''
        i_position = 0
        flag_mod = 0
        stack_mod = 0
        i_mod_name = ''
        modified_sequence = ''

        # 处理N-Term修饰
        if 'N-Term' in dic_FixMOD:
            i_mod_list = dic_FixMOD['N-Term']
            for i_mod_name in i_mod_list:
                modified_sequence += f'{i_mod_name}-'

        for i in range(len(str_ModSequence)):
            if str_ModSequence[i] != '[' and flag_mod == 0:
                i_AA = str_ModSequence[i]
                i_position += 1
                modified_sequence += i_AA  # 先将氨基酸添加到结果序列中
                if i_AA in dic_FixMOD:
                    i_mod_list = dic_FixMOD[i_AA]
                    for i_mod_name in i_mod_list:
                        modified_sequence += f'({i_mod_name})'  # 在氨基酸后面直接加上修饰
                i_mod_name = ''
            elif flag_mod == 0 and str_ModSequence[i] == '[':
                flag_mod = 1
                stack_mod += 1
            else:
                if str_ModSequence[i] == ']':
                    stack_mod -= 1
                elif str_ModSequence[i] == '[':
                    stack_mod += 1
                if stack_mod == 0:
                    flag_mod = 0
                    modified_sequence += f'({i_mod_name})'  # 将修饰加到氨基酸后面
                else:
                    i_mod_name += str_ModSequence[i]

        # 处理C-Term修饰
        if 'C-Term' in dic_FixMOD:
            i_mod_list = dic_FixMOD['C-Term']
            for i_mod_name in i_mod_list:
                modified_sequence += f'-{i_mod_name}'

        return f'{modified_sequence}'

    def __captainJudgeTarget(self, i_Target: str):

        if i_Target == '+':
            return 0
        else:
            return 1


    def read(self, pathID, path_lib):
        with open(pathID, 'r') as f:
            content = f.readline()
            lines = f.readlines()

        content = content.strip().split('\t')
        for i in range(len(content)):
            if content[i] == 'Modified sequence':
                i_SEQ_WITH_MOD = i
            elif content[i] == 'm/z':
                i_PRECURSOR_MOZ_EXP = i
            elif content[i] == 'Charge':
                i_CHARGE = i
            elif content[i] == 'Q-value':
                i_SCORE0 = i
            elif content[i] == 'Sequence':
                i_SEQ = i
            elif content[i] == 'Raw file':
                i_RAW_NAME = i
            elif content[i] == 'Reverse':
                i_DECOY_FLAG = i

        tmp_seq = []
        for line in lines:
            value = line.split('\t')
            # tmpFDR = float(value[i_SCORE0])
            tmpFDR = 0

            if tmpFDR <= self.dp.myCFG.B3_THRESHOLD_FDR:
                mod_seq = self.__captainParseMOD(toolGetWord(line, i_SEQ_WITH_MOD, '\t'),
                                                 toolGetWord(line, i_SEQ, '\t'))
                self.dp.myID.ID4_SEQ_WITH_MOD.append(mod_seq)
                self.dp.myID.ID5_PRECURSOR_MOZ_EXP.append(float(toolGetWord(line, i_PRECURSOR_MOZ_EXP, '\t')))
                self.dp.myID.ID7_CHARGE.append(int(toolGetWord(line, i_CHARGE, '\t')))
                raw_name = toolGetWord(line, i_RAW_NAME, '\t')
                precursor_ID = toolGetWord(line, i_SEQ_WITH_MOD, '\t') + '_'.join(
                    [toolGetNameFromPath(raw_name), toolGetWord(line, i_CHARGE, '\t')])
                self.dp.myID.ID8_PRECURSOR_ID.append(precursor_ID)
                self.dp.myID.ID9_SCORE0.append(tmpFDR)
                tmp_seq.append(toolGetWord(line, i_SEQ, '\t'))
                target = self.__captainJudgeTarget(toolGetWord(line, i_DECOY_FLAG, '\t'))
                self.dp.myID.ID17_TARGET.append(int(target))


        with open(path_lib, 'r') as f:
            content = f.readline()
            lines = f.readlines()

        i_RAW_NAME = toolGetIndexByWord(content, 'Raw file', '\t')
        i_RT = toolGetIndexByWord(content, 'Retention time', '\t')

        i_SEQ = toolGetIndexByWord(content, 'Sequence', '\t')
        i_PRECURSOR = toolGetIndexByWord(content, 'Peptide ID', '\t')
        i_PRECURSOR2 = toolGetIndexByWord(content, 'Mod. peptide ID', '\t')
        i_PRECURSOR3 = toolGetIndexByWord(content, 'Evidence ID', '\t')

        i_FRAGMENT_TYPE_1 = toolGetIndexByWord(content, 'Annotation', '\t')
        i_FRAGMENT_TYPE_2 = toolGetIndexByWord(content, 'Annotation', '\t')
        i_FRAGMENT_TYPE_3 = toolGetIndexByWord(content, 'Charge', '\t')

        i_FRAGMENT_MOZ_EXP = toolGetIndexByWord(content, 'm/z', '\t')
        i_FRAGMENT_LOSS_TYPE = toolGetIndexByWord(content, 'Annotation', '\t')
        pre_transition_group_id = '-1'  # 用于判断是否读到下一个transition_group
        fragment_moz_list = []  # 暂存当前母离子对应子离子m/z的列表
        fragment_type_list = []  # 暂存当前母离子对应子离子类型的列表
        fragment_loss_list = []  # 暂存当前母离子对应子离子的损失类型的列表

        for line in lines:

            tmpFDR = 0

            if tmpFDR <= self.dp.myCFG.B3_THRESHOLD_FDR:

                transition_group_id = toolGetWord(line, i_RAW_NAME, '\t') + toolGetWord(line, i_PRECURSOR,
                                                                                        '\t') + toolGetWord(line,
                                                                                                            i_PRECURSOR2,
                                                                                                            '\t') + toolGetWord(
                    line, i_PRECURSOR3, '\t')

                if pre_transition_group_id != transition_group_id:

                    pre_transition_group_id = transition_group_id

                    if len(fragment_moz_list) != 0:
                        self.dp.myID.ID10_FRAGMENT_TYPE.append(fragment_type_list)
                        self.dp.myID.ID11_FRAGMENT_MOZ_EXP.append(fragment_moz_list)
                        self.dp.myID.ID13_FRAGMENT_LOSS_TYPE.append(fragment_loss_list)
                        self.dp.myID.ID12_FRAGMENT_MOZ_CLC.append([])
                        self.dp.myID.ID14_RT_BEGIN.append(VALUE_ILLEGAL)
                        self.dp.myID.ID15_RT_END.append(VALUE_ILLEGAL)

                    fragment_type_list = []
                    fragment_moz_list = []
                    fragment_loss_list = []
                    fragment_type_1 = toolGetWord(line, i_FRAGMENT_TYPE_1, '\t')[0]
                    fragment_type_2 = toolGetWord(line, i_FRAGMENT_TYPE_2, '\t').split('-')[0].split('(')[0][1:]
                    if '(2+)' in toolGetWord(line, i_FRAGMENT_TYPE_2, '\t'):
                        fragment_type_3 = 2
                    else:
                        fragment_type_3 = 1
                    if fragment_type_1 == 'b' or fragment_type_1 == 'y':
                        fragment_type_list.append(fragment_type_1 + fragment_type_2 + '+' * int(fragment_type_3))
                        fragment_moz_list.append(float(toolGetWord(line, i_FRAGMENT_MOZ_EXP, '\t')))
                        loss_type = toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t')
                        if '-' not in loss_type:
                            fragment_loss_list.append('NOLOSS')  # 格式与system中LOSS_TYPE保持一致
                        else:
                            fragment_loss_list.append(loss_type.split('-')[-1])

                    self.dp.myID.ID0_SEQ.append(toolGetWord(line, i_SEQ, '\t'))
                    self.dp.myID.ID1_RAW_NAME.append(toolGetWord(line, i_RAW_NAME, '\t'))
                    self.dp.myID.ID2_SCAN_ID.append(VALUE_ILLEGAL)
                    self.dp.myID.ID3_RT.append(float(toolGetWord(line, i_RT, '\t')))

                    self.dp.myID.ID6_PRECURSOR_MOZ_CLC.append(VALUE_ILLEGAL)

                    self.dp.myID.ID9_SCORE0.append(tmpFDR)
                    self.dp.myID.N_ID = self.dp.myID.N_ID + 1


                else:

                    fragment_type_1 = toolGetWord(line, i_FRAGMENT_TYPE_1, '\t')[0]
                    fragment_type_2 = toolGetWord(line, i_FRAGMENT_TYPE_2, '\t').split('-')[0].split('(')[0][1:]
                    if '(2+)' in toolGetWord(line, i_FRAGMENT_TYPE_2, '\t'):
                        fragment_type_3 = 2
                    else:
                        fragment_type_3 = 1
                    if fragment_type_1 == 'b' or fragment_type_1 == 'y':
                        fragment_type_list.append(fragment_type_1 + fragment_type_2 + '+' * int(fragment_type_3))
                        fragment_moz_list.append(float(toolGetWord(line, i_FRAGMENT_MOZ_EXP, '\t')))
                        loss_type = toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t')
                        if '-' not in loss_type:
                            fragment_loss_list.append('NOLOSS')  # 格式与system中LOSS_TYPE保持一致
                        else:
                            fragment_loss_list.append(loss_type.split('-')[-1])

        if len(fragment_moz_list) != 0:
            self.dp.myID.ID10_FRAGMENT_TYPE.append(fragment_type_list)
            self.dp.myID.ID11_FRAGMENT_MOZ_EXP.append(fragment_moz_list)
            self.dp.myID.ID12_FRAGMENT_MOZ_CLC.append([])
            self.dp.myID.ID13_FRAGMENT_LOSS_TYPE.append(fragment_loss_list)
        tmp_point = 0
        seq_point = 0
        while len(tmp_seq) != len(self.dp.myID.ID0_SEQ):
            if tmp_seq[tmp_point] == self.dp.myID.ID0_SEQ[seq_point]:
                tmp_point += 1
                seq_point += 1
            else:
                self.dp.myID.ID0_SEQ.pop(seq_point)
                self.dp.myID.ID1_RAW_NAME.pop(seq_point)
                self.dp.myID.ID2_SCAN_ID.pop(seq_point)
                self.dp.myID.ID3_RT.pop(seq_point)
                self.dp.myID.ID6_PRECURSOR_MOZ_CLC.pop(seq_point)
                self.dp.myID.ID9_SCORE0.pop(seq_point)
                self.dp.myID.ID10_FRAGMENT_TYPE.pop(seq_point)
                self.dp.myID.ID11_FRAGMENT_MOZ_EXP.pop(seq_point)
                self.dp.myID.ID13_FRAGMENT_LOSS_TYPE.pop(seq_point)
                self.dp.myID.ID12_FRAGMENT_MOZ_CLC.pop(seq_point)
                self.dp.myID.ID14_RT_BEGIN.pop(seq_point)
                self.dp.myID.ID15_RT_END.pop(seq_point)
                self.dp.myID.N_ID = self.dp.myID.N_ID - 1
        print('over')


class CFunctionParseMS1:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def loadPKL(self, pathMS1):

        pathPKL = pathMS1 + '.pkl'
        dataMS1 = CFileMS1()

        f_pkl = open(pathPKL, 'rb')
        dataMS1 = pickle.load(f_pkl)
        f_pkl.close()

        return dataMS1

    def ms1Topkl(self, pathMS1):

        path_pkl = pathMS1 + '.pkl'

        if os.access(path_pkl, os.F_OK) and self.dp.myCFG.A0_FLAG_USE_EXIST_INDEX == 1:

            pass  # pkl文件已存在

        else:

            dataMS1 = CFileMS1()
            op_INIT_CFILE_MS1(dataMS1)

            with open(pathMS1, 'r') as f:
                lines = f.readlines()

            for line in lines:

                if len(line) > 1:

                    if line.startswith('S'):

                        tmpScan = int(toolGetWord(line, 1, '\t'))
                        dataMS1.INDEX_SCAN.append(tmpScan)

                        dataMS1.MATRIX_PEAK_MOZ[tmpScan] = []
                        dataMS1.MATRIX_PEAK_INT[tmpScan] = []

                    elif line.startswith('H'):
                        pass
                    elif line.startswith('I'):

                        if line.startswith('I	RetTime'):
                            dataMS1.INDEX_RT.append(float(toolGetWord(line, 2, '\t')))

                        elif line.startswith('I	IonInjectionTime'):
                            dataMS1.LIST_ION_INJECTION_TIME[tmpScan] = float(toolGetWord(line, 2, '\t'))
                    else:

                        dataMS1.MATRIX_PEAK_MOZ[tmpScan].append(float(toolGetWord(line, 0, ' ')))
                        dataMS1.MATRIX_PEAK_INT[tmpScan].append(float(toolGetWord(line, 1, ' ')))

            f_pkl = open(path_pkl, 'wb')
            pickle.dump(dataMS1, f_pkl)
            f_pkl.close()


class CFunctionParseMS2:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def loadPKL(self, pathMS2):

        pathPKL = pathMS2 + '.pkl'
        dataMS2 = CFileMS2()

        f_pkl = open(pathPKL, 'rb')
        dataMS2 = pickle.load(f_pkl)

        return dataMS2

    def ms2topkl(self, pathMS2):

        path_pkl = pathMS2 + '.pkl'

        if os.access(path_pkl, os.F_OK) and self.dp.myCFG.A0_FLAG_USE_EXIST_INDEX == 1:
            pass
        else:

            dataMS2 = CFileMS2()
            op_INIT_CFILE_MS2(dataMS2)

            with open(pathMS2, 'r') as f:
                lines = f.readlines()

            flag_repeat = 0
            for line in lines:

                if len(line) > 1:

                    if line.startswith('S'):

                        tmpScan = int(toolGetWord(line, 1, '\t'))
                        if tmpScan in dataMS2.INDEX_SCAN:

                            flag_repeat = 1
                        else:

                            dataMS2.INDEX_SCAN.append(tmpScan)

                            dataMS2.MATRIX_PEAK_INT[tmpScan] = []
                            dataMS2.MATRIX_PEAK_MOZ[tmpScan] = []

                            flag_repeat = 0

                    elif line.startswith('H'):
                        pass

                    elif line.startswith('Z'):
                        pass

                    elif flag_repeat == 0:  # 重复的scan号只读一遍

                        if line.startswith('I'):

                            if line.startswith('I	RetTime'):
                                tmpRT = float(toolGetWord(line, 2, '\t'))
                                dataMS2.LIST_RET_TIME[tmpScan] = tmpRT
                                dataMS2.INDEX_RT.append(tmpRT)

                            if line.startswith('I	PrecursorScan'):
                                dataMS2.LIST_PRECURSOR_SCAN[tmpScan] = int(toolGetWord(line, 2, '\t'))

                            if line.startswith('I	ActivationCenter'):
                                dataMS2.LIST_ACTIVATION_CENTER[tmpScan] = float(toolGetWord(line, 2, '\t'))

                            if line.startswith('I	IonInjectionTime'):
                                dataMS2.LIST_ION_INJECTION_TIME[tmpScan] = float(toolGetWord(line, 2, '\t'))

                        else:

                            dataMS2.MATRIX_PEAK_MOZ[tmpScan].append(float(toolGetWord(line, 0, ' ')))
                            dataMS2.MATRIX_PEAK_INT[tmpScan].append(float(toolGetWord(line, 1, ' ')))

            f_pkl = open(path_pkl, 'wb')
            pickle.dump(dataMS2, f_pkl)
            f_pkl.close()
