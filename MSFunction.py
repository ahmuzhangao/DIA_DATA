import re
from tqdm import tqdm
from MSTool import toolGetWord, toolGetWord1, toolCountCharInString, toolStr2List, toolGetIndexByWord, \
    toolGetNameFromPath
from MSData import CFileMS1, CDataPack, CFileMS2
from MSOperator import op_INIT_CFILE_MS1, op_STANDARDIZE_MOD_BY_PRECURSOR, op_INIT_CFILE_MS2, op_preallocate_matrix, op_convert_peptdeep, op_STANDARDIZE_MOD_BY_PRECURSOR_SPECTORNAUT, \
    op_convert_spectronaut, op_convert_MaxDIA
from MSLogging import logGetError
from MSSystem import VALUE_ILLEGAL, UNIMOID_TO_STANDARD_MOD2
from MSSystem import UNIMOID_TO_STANDARD_MOD
import pickle
import os

class CFunctionParseMS1:

    def __init__(self, inputDP):
        self.dp = inputDP

    def loadPKL(self, pathMS1):

        pathPKL = pathMS1 + '.pkl'
        dataMS1 = CFileMS1()

        f_pkl = open(pathPKL, 'rb')
        dataMS1 = pickle.load(f_pkl)
        f_pkl.close()

        return dataMS1

    def ms1TOpkl(self, pathMS1):

        # check
        path_pkl = pathMS1 + ".pkl"
        if os.access(path_pkl, os.F_OK) and 1 == self.dp.myCFG.A0_FLAG_USE_EXIST_INDEX:
            pass
        else:
            # init
            dataMS1 = CFileMS1()
            op_INIT_CFILE_MS1(dataMS1)

            # open
            with open(pathMS1, 'r') as f:

                i_MS1 = -1

                for line in tqdm(f.readlines(),desc="Read MS1"):

                    len_line = len(line)
                    if len_line > 1:

                        if line.startswith("S	"):
                            i_MS1 = i_MS1 + 1
                            tmpScan = int(toolGetWord(line, 1, '	'))
                            dataMS1.INDEX_SCAN.append(tmpScan)

                            dataMS1.MATRIX_PEAK_MOZ[tmpScan] = []  # 这两行不能少
                            dataMS1.MATRIX_PEAK_INT[tmpScan] = []

                        elif line.startswith("H	"):
                            pass
                        elif line.startswith("I	"):
                            if line.startswith("I	IonInjectionTime"):
                                t = toolGetWord(line, 2, '	')
                                dataMS1.LIST_ION_INJECTION_TIME[tmpScan] = float(t)
                            elif line.startswith("I	RetTime	"):
                                t = toolGetWord(line, 2, '	')
                                dataMS1.INDEX_RT.append(float(t))
                                dataMS1.LIST_RET_TIME[tmpScan] = float(t)
                        else:
                            dataMS1.MATRIX_PEAK_MOZ[tmpScan].append(float(toolGetWord(line, 0, ' ')))
                            dataMS1.MATRIX_PEAK_INT[tmpScan].append(float(toolGetWord(line, 1, ' ')))

            # write pkl
            fid_pkl = open(path_pkl, 'wb')
            pickle.dump(dataMS1, fid_pkl)
            fid_pkl.close()

class CFunctionParseMS2:

    def __init__(self, inputDP:CDataPack):

        self.dp = inputDP

    def loadPKL(self, pathMS2):

        pathPKL = pathMS2 + '.pkl'
        dataMS2 = CFileMS2()

        f_pkl = open(pathPKL, 'rb')
        dataMS2 = pickle.load(f_pkl)
        f_pkl.close()

        return dataMS2

    def ms2topkl(self, pathMS2):

        path_pkl = pathMS2 + '.pkl'

        if os.access(path_pkl, os.F_OK) and self.dp.myCFG.A0_FLAG_USE_EXIST_INDEX == 1:
            pass
        else:

            dataMS2 = CFileMS2()
            op_INIT_CFILE_MS2(dataMS2)

            with open(pathMS2, 'r')as f:
                lines = f.readlines()

            prev_scan = None
            flag_repeat = 0
            for line in tqdm(lines, desc="Read MS2"):
            # for line in lines:
                if len(line) > 1:

                    if line.startswith('S'):

                        tmpScan = int(toolGetWord(line, 1, '\t'))
                        if tmpScan == prev_scan:

                            flag_repeat = 1
                        else:
                            prev_scan = tmpScan

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

                        else:#

                            dataMS2.MATRIX_PEAK_MOZ[tmpScan].append(float(toolGetWord(line, 0, ' ')))
                            dataMS2.MATRIX_PEAK_INT[tmpScan].append(float(toolGetWord(line, 1, ' ')))

            f_pkl = open(path_pkl, 'wb')
            pickle.dump(dataMS2, f_pkl)
            f_pkl.close()

class CFunctionConfig:

    def config2file(self, path, config):

        with open(path, 'w') as f:
            f.write("#[Ini files]\n")
            f.write("INI_PATH_ELEMENT=" + config.I0_INI_PATH_ELEMENT + '\n')
            f.write("INI_PATH_AA=" + config.I1_INI_PATH_AA + '\n')
            f.write("INI_PATH_MOD=" + config.I2_INI_PATH_MOD + '\n')

            f.write("\n#[Data]\n")
            # f.write("PATH_MS1=" + config.A1_PATH_MS1 + '\n')
            # f.write("TYPE_MS1=" + str(config.A2_TYPE_MS1) + '\n')
            # f.write("PATH_MS2=" + config.A3_PATH_MS2 + '\n')
            # f.write("TYPE_MS2=" + str(config.A4_TYPE_MS2) + '\n')
            f.write("PATH_FILE=" + config.A5_PATH_FILE + '\n')
            f.write("TYPE_FILE=" + str(config.A6_TYPE_FILE) + '\n')

            f.write("\n#[Type of task]\n")
            f.write("TYPE_FLOW=" + str(config.C0_TYPE_FLOW) + '\n')

            f.write("\n#[Identification results]\n")
            f.write("PATH_IDENTIFICATION_RESULT=" + config.B1_PATH_IDENTIFICATION_RESULT + '\n')
            f.write("TYPE_IDENTIFICATION_RESULT=" + str(config.B2_TYPE_IDENTIFICATION_RESULT) + '\n')
            f.write("THRESHOLD_FDR=" + str(config.B3_THRESHOLD_FDR) + '\n')
            f.write("FLAG_DYNAMIC_LIBRARY=" + str(config.B5_FLAG_DYNAMIC_LIBRARY) + '\n')

            f.write('\n# [Quantitation]\n')
            f.write('DIA_RT_HALF_WIN_IN_MIN=' + str(config.C1_DIA_RT_HALF_WIN_IN_MIN) + '\n')
            f.write('DIA_ACCURACY_HALF_WIN_PEAK=' + str(config.C2_DIA_ACCURACY_HALF_WIN_PEAK) + '\n')
            f.write('DIA_TYPE_ACCURACY_HALF_WIN_PEAK=' + str(config.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK) + '\n')
            f.write("DIALF_TYPE_MODEL=" + str(config.C4_DIALF_TYPE_MODEL) + '\n')
            f.write("FLAG_CALIBRATE_RT=" + str(config.C5_DIALF_FLAG_CALIBRATE_RT) + '\n')
            f.write("DIALF_PATH_CALIBRATE_RT=" + str(config.C5_DIALF_PATH_CALIBRATE_RT) + '\n')
            f.write("DIALF_FLAG_DECOY=" + str(config.C6_FLAG_DECOY) + '\n')
            f.write("DIALF_PROCESS_NUM=" + str(config.C7_PROCESS_NUM) + '\n')

            f.write("\n#[Export]\n")
            f.write("PATH_EXPORT=" + config.D1_PATH_EXPORT + '\n')
            f.write("TYPE_EXPORT=" + str(config.D2_TYPE_EXPORT) + '\n')

    def file2config(self, path, config):

        with open(path, 'r', encoding='utf-8') as f:

            print("\n#################### Your Parameters ####################")

            for line in f.readlines():

                if line.startswith("#") or line.startswith("["):
                    continue

                p_EqualSign = line.find('=')

                if -1 == p_EqualSign:
                    pass
                else:

                    subLine = toolGetWord(line, 0, ';')  # ;后面的是注释
                    self.__soldierParseLine(subLine, config)

            print("#################### Your Parameters ####################\n")

    def __soldierParseLine(self, line, cfg):

        str_name = toolGetWord(line, 0, '=')
        str_value = toolGetWord(line, 1, '=').replace("\n", "")
        print(str_name + ": " + str_value)

        if "PATH_MS1" == str_name:

            cfg.A1_PATH_MS1 = str_value

        elif "TYPE_MS1" == str_name:

            cfg.A2_TYPE_MS1 = int(str_value)

        elif "PATH_MS2" == str_name:

            cfg.A3_PATH_MS2 = str_value

        elif "TYPE_MS2" == str_name:

            cfg.A4_TYPE_MS2 = int(str_value)

        elif "PATH_FILE" == str_name:

            cfg.A5_PATH_FILE = str_value

        elif "TYPE_FILE" == str_name:

            cfg.A6_TYPE_FILE = int(str_value)

        elif "PATH_IDENTIFICATION_RESULT" == str_name:

            cfg.B1_PATH_IDENTIFICATION_RESULT = str_value

        elif "TYPE_IDENTIFICATION_RESULT" == str_name:

            cfg.B2_TYPE_IDENTIFICATION_RESULT = int(str_value)

        elif "THRESHOLD_FDR" == str_name:

            cfg.B3_THRESHOLD_FDR = float(str_value)

        elif 'DIA_RT_HALF_WIN_IN_MIN' == str_name:

            cfg.C1_DIA_RT_HALF_WIN_IN_MIN = float(str_value)

        elif 'DIA_ACCURACY_HALF_WIN_PEAK' == str_name:

            cfg.C2_DIA_ACCURACY_HALF_WIN_PEAK = float(str_value)

        elif 'DIA_TYPE_ACCURACY_HALF_WIN_PEAK' == str_name:

            cfg.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK = int(str_value)

        elif 'DIALF_TYPE_MODEL' == str_name:

            cfg.C4_DIALF_TYPE_MODEL = int(str_value)

        elif 'DIALF_PATH_CALIBRATE_RT' == str_name:

            cfg.C5_DIALF_PATH_CALIBRATE_RT = str_value

        elif 'DIALF_FLAG_DECOY' == str_name:

            cfg.C6_FLAG_DECOY = int(str_value)

        elif 'DIALF_PROCESS_NUM' == str_name:

            cfg.C7_PROCESS_NUM = int(str_value)

        elif 'TYPE_FLOW' == str_name:

            cfg.C0_TYPE_FLOW = int(str_value)

        elif 'PATH_EXPORT' == str_name:

            cfg.D1_PATH_EXPORT = str_value

        elif 'TYPE_EXPORT' == str_name:

            cfg.D2_TYPE_EXPORT = int(str_value)

        elif "INI_PATH_ELEMENT" == str_name:

            cfg.I0_INI_PATH_ELEMENT = str_value

        elif "INI_PATH_AA" == str_name:

            cfg.I1_INI_PATH_AA = str_value

        elif "INI_PATH_MOD" == str_name:

            cfg.I2_INI_PATH_MOD = str_value

        elif "FLAG_DYNAMIC_LIBRARY" == str_name:

            cfg.B5_FLAG_DYNAMIC_LIBRARY = int(str_value)

        elif "FLAG_CALIBRATE_RT" == str_name:

            cfg.C5_DIALF_FLAG_CALIBRATE_RT = int(str_value)

        elif "FLAG_DRAW" == str_name:

            cfg.B6_FLAG_DRAW = int(str_value)

        elif "FRAG_NUM" == str_name:

            cfg.B7_FRAG_NUM = int(str_value)

        else:
            print(str_name)
            info = "MSFunction, MK299, " + str_name + " is all right?"
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

    def __CalculateAminoacidMass(self, inputMolecule, inputAtonMass):
        # 计算氨基酸相对分子质量
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

    def  __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soldierParseMOD(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理

        return op_STANDARDIZE_MOD_BY_PRECURSOR(inputStr)

    def __soldierParseMOD2(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
        # 从DIA-NN格式的带修饰母离子字符串中提取修饰列表，每个元素存放修饰位点和类型
        seperator_start, seperator_end = '(', ')'
        mod_dic = {}
        mod_type = ''
        site_start, site_end = 0, 0
        flag_in_mod = 0
        count = 0
        for i, char in enumerate(inputStr):
            if char == seperator_start:
                site_start = i
                flag_in_mod = 1
            elif char == seperator_end:
                site_end = i
                try:
                    mod_type = UNIMOID_TO_STANDARD_MOD[inputStr[site_start + 1:site_end]]
                except KeyError:
                    print(inputStr[site_start + 1:site_end])
                    print('wrong')
                    info = 'DIA-NN modify type has wrong, please check!'
                    logGetError(info)
                mod_dic[count] = mod_type
                flag_in_mod = 0
            elif flag_in_mod == 0:
                count = count + 1

        return mod_dic

    def __soldierBinary_search(self, lines_Lib, target_precursor, i_PRECURSOR_Lib):

        start, end = -1, -1
        found = False

        low, high = 0, len(lines_Lib) - 1

        while low <= high:
            mid = (low + high) // 2
            current_precursor = toolGetWord(lines_Lib[mid], i_PRECURSOR_Lib, '\t')

            if current_precursor == target_precursor:
                start = mid
                found = True
                break
            elif current_precursor < target_precursor:
                low = mid + 1
            else:
                high = mid - 1

        if found:
            # 找到起始索引后，寻找结束索引
            for i in range(start, len(lines_Lib)):
                if toolGetWord(lines_Lib[i], i_PRECURSOR_Lib, '\t') != target_precursor:
                    end = i - 1
                    break
            else:
                end = len(lines_Lib) - 1

        return start, end

    def read(self, path):

        with open(path, 'r')as f:
            content = f.readline()
            lines = f.readlines()

        self.dp.myID.PSM_LINE_T.append(content.strip('\n'))

        i_RAW_NAME = toolGetIndexByWord(content, 'File.Name', '\t')
        i_RT = toolGetIndexByWord(content, 'RT', '\t')

        i_SEQ_WITH_MOD = toolGetIndexByWord(content, 'Modified.Sequence', '\t')
        i_SEQ = toolGetIndexByWord(content, 'Stripped.Sequence', '\t')
        i_CHARGE = toolGetIndexByWord(content, 'Precursor.Charge', '\t')
        i_PRECURSOR = toolGetIndexByWord(content, 'Precursor.Id', '\t')

        i_SCORE0 = toolGetIndexByWord(content, 'Q.Value', '\t')

        i_Protein_Group = toolGetIndexByWord(content, 'Protein.Group', '\t')
        i_Protein_ID = toolGetIndexByWord(content, 'Protein.Ids', '\t')
        i_Protein_Name = toolGetIndexByWord(content, 'Protein.Names', '\t')
        i_Gene = toolGetIndexByWord(content, 'Genes', '\t')

        i_RT_BEGIN = toolGetIndexByWord(content, 'RT.Start', '\t')
        i_RT_STOP = toolGetIndexByWord(content, 'RT.Stop', '\t')

        for line in lines:

            tmpFDR = float(toolGetWord(line, i_SCORE0, '\t'))

            self.dp.myID.PSM_LINE_C.append(line.strip('\n'))
            self.dp.myID.ID0_SEQ.append(toolGetWord(line, i_SEQ, '\t'))

            self.dp.myID.ID1_RAW_NAME.append(toolGetWord(line, i_RAW_NAME, '\t'))

            #单位是分钟
            self.dp.myID.ID3_RT.append(float(toolGetWord(line, i_RT, '\t')))

            #seq+mod：mod在seq的对应位置中间，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK
            self.dp.myID.ID4_SEQ_WITH_MOD.append(self.__soldierParseMOD(toolGetWord(line, i_SEQ_WITH_MOD, '\t')))

            self.dp.myID.ID5_PRECURSOR_MOZ_EXP.append(VALUE_ILLEGAL)
            self.dp.myID.ID6_PRECURSOR_MOZ_CLC.append(VALUE_ILLEGAL)

            self.dp.myID.ID7_CHARGE.append(int(toolGetWord(line, i_CHARGE, '\t')))

            # raw+seq+mod+charge：raw在最前面，mod在seq的对应位置中间，charge在seq后面，例如：CNCP2023_cpnl_480_hela_HRDIA_2_AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myID.ID8_PRECURSOR_ID.append(toolGetNameFromPath(toolGetWord(line, i_RAW_NAME, '\t')) + '_' + self.__soldierParseMOD(toolGetWord(line, i_PRECURSOR, '\t')))

            self.dp.myID.ID9_SCORE0.append(tmpFDR)

            self.dp.myID.ID10_Protein_Groups.append(toolGetWord(line, i_Protein_Group, '\t'))
            self.dp.myID.ID11_Protein_Ids.append(toolGetWord(line,i_Protein_ID,'\t'))
            self.dp.myID.ID12_Protein_Names.append(toolGetWord(line, i_Protein_Name, '\t'))
            self.dp.myID.ID13_Genes.append(toolGetWord(line, i_Gene, '\t'))

            self.dp.myID.ID14_RT_BEGIN.append(float(toolGetWord(line, i_RT_BEGIN, '\t')))
            self.dp.myID.ID15_RT_END.append(float(toolGetWord(line, i_RT_STOP, '\t')))

            # seq+mod+charge：mod在seq的对应位置中间，charge在seq后面，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myID.ID17_SEQ_WITH_MOD_AND_CHARGE.append(self.__soldierParseMOD(toolGetWord(line, i_PRECURSOR, '\t')))

            seq_with_mod = toolGetWord(line, i_SEQ_WITH_MOD, '\t')
            mod_dict = self.__soldierParseMOD2(seq_with_mod)
            MOD = ''
            mod_type = ''
            mod_site = ''
            for i_key in mod_dict.keys():
                i_value = mod_dict[i_key]
                MOD += str(i_key) + ',' + str(i_value) + ';'
                mod_type += UNIMOID_TO_STANDARD_MOD2[str(i_value)] + ';'
                mod_site += str(i_key) + ';'

            #格式是：modsite,modtype;
            self.dp.myID.ID16_MOD.append(MOD)

            # 格式是：modtype;(这里为了适配pepdeep的修饰形式，需要再转一下)
            self.dp.myID.ID19_MOD_TYPE.append(mod_type)

            # 格式是：modsite;
            self.dp.myID.ID18_MOD_SITE.append(mod_site)

            self.dp.myID.N_ID = self.dp.myID.N_ID + 1

class CFunctionParseIDForSpectronaut:

    def  __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soldierParseMOD(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理

        return op_STANDARDIZE_MOD_BY_PRECURSOR_SPECTORNAUT(inputStr)

    def __soldierParseMOD2(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
        # 从DIA-NN格式的带修饰母离子字符串中提取修饰列表，每个元素存放修饰位点和类型
        seperator_start, seperator_end = '(', ')'
        mod_dic = {}
        mod_type = ''
        site_start, site_end = 0, 0
        flag_in_mod = 0
        count = 0
        for i, char in enumerate(inputStr):
            if char == seperator_start:
                site_start = i
                flag_in_mod = 1
            elif char == seperator_end:
                site_end = i
                mod_type = inputStr[site_start + 1:site_end]
                mod_dic[count] = mod_type
                flag_in_mod = 0
            elif flag_in_mod == 0:
                count = count + 1

        return mod_dic

    def read(self, path):

        with open(path, 'r')as f:
            content = f.readline()
            lines = f.readlines()

        self.dp.myID.PSM_LINE_T.append(content.strip('\n'))

        i_RAW_NAME = toolGetIndexByWord(content, 'R.FileName', '\t')
        i_RT = toolGetIndexByWord(content, 'EG.ApexRT', '\t')

        i_SEQ_WITH_MOD = toolGetIndexByWord(content, 'EG.ModifiedSequence', '\t') # 得去掉下划线
        i_SEQ = toolGetIndexByWord(content, 'PEP.StrippedSequence', '\t')
        i_CHARGE = toolGetIndexByWord(content, 'FG.Charge', '\t')
        i_PRECURSOR = toolGetIndexByWord(content, 'EG.PrecursorId', '\t') # 得去掉下划线

        i_SCORE0 = toolGetIndexByWord(content, 'EG.Qvalue', '\t')

        i_Protein_Group = toolGetIndexByWord(content, 'PG.ProteinGroups', '\t')
        i_Protein_ID = toolGetIndexByWord(content, 'PG.UniProtIds', '\t')
        i_Protein_Name = toolGetIndexByWord(content, 'PG.ProteinNames', '\t')

        i_RT_BEGIN = toolGetIndexByWord(content, 'EG.StartQuantRT', '\t')
        i_RT_STOP = toolGetIndexByWord(content, 'EG.EndQuantRT', '\t')

        for line in lines:

            tmpFDR = float(toolGetWord(line, i_SCORE0, '\t'))

            self.dp.myID.PSM_LINE_C.append(line.strip('\n'))
            self.dp.myID.ID0_SEQ.append(toolGetWord(line, i_SEQ, '\t'))

            self.dp.myID.ID1_RAW_NAME.append(toolGetWord(line, i_RAW_NAME, '\t'))

            #单位是分钟
            self.dp.myID.ID3_RT.append(float(toolGetWord(line, i_RT, '\t')))

            #seq+mod：mod在seq的对应位置中间，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK
            self.dp.myID.ID4_SEQ_WITH_MOD.append(op_convert_spectronaut(toolGetWord(line, i_SEQ_WITH_MOD, '\t')))

            self.dp.myID.ID5_PRECURSOR_MOZ_EXP.append(VALUE_ILLEGAL)
            self.dp.myID.ID6_PRECURSOR_MOZ_CLC.append(VALUE_ILLEGAL)

            self.dp.myID.ID7_CHARGE.append(int(toolGetWord(line, i_CHARGE, '\t')))

            # raw+seq+mod+charge：raw在最前面，mod在seq的对应位置中间，charge在seq后面，例如：CNCP2023_cpnl_480_hela_HRDIA_2_AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myID.ID8_PRECURSOR_ID.append(toolGetNameFromPath(toolGetWord(line, i_RAW_NAME, '\t')) + '_' + op_convert_spectronaut(toolGetWord(line, i_PRECURSOR, '\t')))

            self.dp.myID.ID9_SCORE0.append(tmpFDR)

            self.dp.myID.ID10_Protein_Groups.append(toolGetWord(line, i_Protein_Group, '\t'))
            self.dp.myID.ID11_Protein_Ids.append(toolGetWord(line,i_Protein_ID,'\t'))
            self.dp.myID.ID12_Protein_Names.append(toolGetWord(line, i_Protein_Name, '\t'))
            # self.dp.myID.ID13_Genes.append(toolGetWord(line, i_Gene, '\t'))

            self.dp.myID.ID14_RT_BEGIN.append(float(toolGetWord(line, i_RT_BEGIN, '\t')))
            self.dp.myID.ID15_RT_END.append(float(toolGetWord(line, i_RT_STOP, '\t')))

            # seq+mod+charge：mod在seq的对应位置中间，charge在seq后面，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myID.ID17_SEQ_WITH_MOD_AND_CHARGE.append(op_convert_spectronaut(toolGetWord(line, i_PRECURSOR, '\t')))

            seq_with_mod = op_convert_spectronaut(toolGetWord(line, i_SEQ_WITH_MOD, '\t'))
            mod_dict = self.__soldierParseMOD2(seq_with_mod)
            MOD = ''
            mod_type = ''
            mod_site = ''
            for i_key in mod_dict.keys():
                i_value = mod_dict[i_key]
                MOD += str(i_key) + ',' + str(i_value) + ';'
                if op_convert_peptdeep(str(i_value)) == 'Acetyl@ProteinN-term':
                    mod_type += 'Acetyl@Protein_N-term' + ';'
                else:
                    mod_type += op_convert_peptdeep(str(i_value)) + ';'

                mod_site += str(i_key) + ';'

            #格式是：modsite,modtype;
            self.dp.myID.ID16_MOD.append(MOD)

            # 格式是：modtype;(这里为了适配pepdeep的修饰形式，需要再转一下)
            self.dp.myID.ID19_MOD_TYPE.append(mod_type)

            # 格式是：modsite;
            self.dp.myID.ID18_MOD_SITE.append(mod_site)

            self.dp.myID.N_ID = self.dp.myID.N_ID + 1


class CFunctionParseIDForMaxDIA:

    def  __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soldierParseMOD(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理

        return op_STANDARDIZE_MOD_BY_PRECURSOR_SPECTORNAUT(inputStr)

    def __soldierParseMOD2(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
        # 从DIA-NN格式的带修饰母离子字符串中提取修饰列表，每个元素存放修饰位点和类型
        seperator_start, seperator_end = '(', ')'
        mod_dic = {}
        mod_type = ''
        site_start, site_end = 0, 0
        flag_in_mod = 0
        count = 0
        for i, char in enumerate(inputStr):
            if char == seperator_start:
                site_start = i
                flag_in_mod = 1
            elif char == seperator_end:
                site_end = i
                mod_type = inputStr[site_start + 1:site_end]
                mod_dic[count] = mod_type
                flag_in_mod = 0
            elif flag_in_mod == 0:
                count = count + 1

        return mod_dic

    def read(self, path):

        with open(path, 'r')as f:
            content = f.readline()
            lines = f.readlines()

        self.dp.myID.PSM_LINE_T.append(content.strip('\n'))

        i_RAW_NAME = toolGetIndexByWord(content, 'Raw file', '\t')
        i_RT = toolGetIndexByWord(content, 'EG.ApexRT', '\t')

        i_SEQ_WITH_MOD = toolGetIndexByWord(content, 'Modified sequence', '\t') # 得去掉下划线
        i_SEQ = toolGetIndexByWord(content, 'Sequence', '\t')
        i_CHARGE = toolGetIndexByWord(content, 'Charge', '\t')

        i_SCORE0 = toolGetIndexByWord(content, 'Q-value', '\t')

        i_Protein_Group = toolGetIndexByWord(content, 'Proteins', '\t')
        i_Protein_ID = toolGetIndexByWord(content, 'Leading razor protein', '\t')
        i_Protein_Name = toolGetIndexByWord(content, 'Protein names', '\t')


        for line in lines:

            tmpFDR = float(toolGetWord(line, i_SCORE0, '\t'))

            self.dp.myID.PSM_LINE_C.append(line.strip('\n'))
            self.dp.myID.ID0_SEQ.append(toolGetWord(line, i_SEQ, '\t'))

            self.dp.myID.ID1_RAW_NAME.append(toolGetWord(line, i_RAW_NAME, '\t'))

            #单位是分钟
            self.dp.myID.ID3_RT.append(float(toolGetWord(line, i_RT, '\t')))

            #seq+mod：mod在seq的对应位置中间，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK
            self.dp.myID.ID4_SEQ_WITH_MOD.append(op_convert_MaxDIA(toolGetWord(line, i_SEQ_WITH_MOD, '\t')))

            self.dp.myID.ID5_PRECURSOR_MOZ_EXP.append(VALUE_ILLEGAL)
            self.dp.myID.ID6_PRECURSOR_MOZ_CLC.append(VALUE_ILLEGAL)

            self.dp.myID.ID7_CHARGE.append(int(toolGetWord(line, i_CHARGE, '\t')))

            # raw+seq+mod+charge：raw在最前面，mod在seq的对应位置中间，charge在seq后面，例如：CNCP2023_cpnl_480_hela_HRDIA_2_AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myID.ID8_PRECURSOR_ID.append(toolGetNameFromPath(toolGetWord(line, i_RAW_NAME, '\t')) + '_' + op_convert_MaxDIA(toolGetWord(line, i_SEQ_WITH_MOD, '\t')) + '_' + int(toolGetWord(line, i_CHARGE, '\t')))

            self.dp.myID.ID9_SCORE0.append(tmpFDR)

            self.dp.myID.ID10_Protein_Groups.append(toolGetWord(line, i_Protein_Group, '\t'))
            self.dp.myID.ID11_Protein_Ids.append(toolGetWord(line,i_Protein_ID,'\t'))
            self.dp.myID.ID12_Protein_Names.append(toolGetWord(line, i_Protein_Name, '\t'))
            # self.dp.myID.ID13_Genes.append(toolGetWord(line, i_Gene, '\t'))

            # self.dp.myID.ID14_RT_BEGIN.append(float(toolGetWord(line, i_RT_BEGIN, '\t')))
            # self.dp.myID.ID15_RT_END.append(float(toolGetWord(line, i_RT_STOP, '\t')))

            # seq+mod+charge：mod在seq的对应位置中间，charge在seq后面，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myID.ID17_SEQ_WITH_MOD_AND_CHARGE.append(op_convert_MaxDIA(toolGetWord(line, i_SEQ_WITH_MOD, '\t')) + '_' + int(toolGetWord(line, i_CHARGE, '\t')))

            seq_with_mod = op_convert_spectronaut(toolGetWord(line, i_SEQ_WITH_MOD, '\t'))
            mod_dict = self.__soldierParseMOD2(seq_with_mod)
            MOD = ''
            mod_type = ''
            mod_site = ''
            for i_key in mod_dict.keys():
                i_value = mod_dict[i_key]
                MOD += str(i_key) + ',' + str(i_value) + ';'
                if op_convert_peptdeep(str(i_value)) == 'Acetyl@ProteinN-term':
                    mod_type += 'Acetyl@Protein_N-term' + ';'
                else:
                    mod_type += op_convert_peptdeep(str(i_value)) + ';'

                mod_site += str(i_key) + ';'

            #格式是：modsite,modtype;
            self.dp.myID.ID16_MOD.append(MOD)

            # 格式是：modtype;(这里为了适配pepdeep的修饰形式，需要再转一下)
            self.dp.myID.ID19_MOD_TYPE.append(mod_type)

            # 格式是：modsite;
            self.dp.myID.ID18_MOD_SITE.append(mod_site)

            self.dp.myID.N_ID = self.dp.myID.N_ID + 1

class CFunctionReadApuQuantIDForCheck:

    def __init__(self, inputDP: CDataPack):
        self.dp = inputDP

    def __captainGetRawName(self, Intensity_Column: str):

        separator_start, separator_end = '(', ')'
        start = Intensity_Column.index(separator_start) + 1
        end = len(Intensity_Column) - list(reversed(Intensity_Column)).index(separator_end) - 1
        raw_name = Intensity_Column[start:end]

        return raw_name

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

    def __soldierParseMOD(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理

        return op_STANDARDIZE_MOD_BY_PRECURSOR(inputStr)

    def __soldierParseMOD2(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
        # 从DIA-NN格式的带修饰母离子字符串中提取修饰列表，每个元素存放修饰位点和类型
        seperator_start, seperator_end = '(', ')'
        mod_dic = {}
        mod_type = ''
        site_start, site_end = 0, 0
        flag_in_mod = 0
        count = 0
        for i, char in enumerate(inputStr):
            if char == seperator_start:
                site_start = i
                flag_in_mod = 1
            elif char == seperator_end:
                site_end = i
                try:
                    mod_type = UNIMOID_TO_STANDARD_MOD[inputStr[site_start + 1:site_end]]
                except KeyError:
                    print(inputStr[site_start + 1:site_end])
                    print('wrong')
                    info = 'DIA-NN modify type has wrong, please check!'
                    logGetError(info)
                mod_dic[count] = mod_type
                flag_in_mod = 0
            elif flag_in_mod == 0:
                count = count + 1

        return mod_dic

    def __captainParseMOD(self, inputStr):  # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
        # 从DIA-NN格式的带修饰母离子字符串中提取修饰列表，每个元素存放修饰位点和类型
        seperator_start, seperator_end = '(', ')'
        mod_dic = {}
        mod_type = ''
        site_start, site_end = 0, 0
        flag_in_mod = 0
        count = 0
        for i, char in enumerate(inputStr):
            if char == seperator_start:
                site_start = i
                flag_in_mod = 1
            elif char == seperator_end:
                site_end = i
                try:
                    mod_type = UNIMOID_TO_STANDARD_MOD[inputStr[site_start + 1:site_end]]
                except KeyError:
                    print('wrong')
                    info = 'DIA-NN modify type has wrong, please check!'
                    logGetError(info)
                mod_dic[count] = mod_type
                flag_in_mod = 0
            elif flag_in_mod == 0:
                count = count + 1

        return mod_dic

    def read(self, path_Report):

        with open(path_Report, 'rb')as f:
            lines = f.read().decode(encoding='utf-8').split('\r\n')
        lines = lines[: -1] if lines[-1] == '' else lines
        content = lines[0]
        lines = lines[1:]

        i_Empty_Separator = toolGetIndexByWord(content, 'Empty_Separator', '\t')
        i_SEQ_WITH_MOD = toolGetIndexByWord(content, 'Modified.Sequence', '\t')
        i_SEQ = toolGetIndexByWord(content, 'Stripped.Sequence', '\t')
        i_RUN = toolGetIndexByWord(content, 'Run', '\t')

        n_raw_name = 0  # 统计raw文件的数目
        raw_name_list = []  # 存放raw文件名字的列表
        i_RT_Start_list = []  # 存放n个(n=raw文件数目)结果的保留时间起点
        i_RT_Length_list = []  # 存放n个(n=raw文件数目)结果的保留时间长度
        i_RT_Mid_list = []
        i_MH_list = []  # 存放n个(n=raw文件数目)结果的母离子质量mass
        i_Charge_list = []  # 存放n个(n=raw文件数目)结果的母离子电荷
        i_Score_list = []  # 存放n个(n=raw文件数目)结果的MBR打分
        content = content.split('\t')

        # 统计raw文件
        for i in content[(i_Empty_Separator + 1):]:
            if i.startswith('Intensity'):
                raw_name_list.append(self.__captainGetRawName(i))
                n_raw_name += 1
            else:
                break
        # 计算每个raw文件对应的保留时间起点和保留时间长度的位置
        position_RT = i_Empty_Separator + n_raw_name + 1
        for i, rt_start in enumerate(content[position_RT: (position_RT + 3 * n_raw_name): 3]):
            if rt_start.startswith('RT_Start_Sec'):
                i_RT_Start_list.append(i_Empty_Separator + 1 + n_raw_name + 3 * i)
        for i, i_rt_start in enumerate(i_RT_Start_list):
            i_rt_length = i_rt_start + 1
            if content[i_rt_length].startswith('RT_Length_Sec'):
                i_RT_Length_list.append(i_rt_length)
        for i, i_rt_start in enumerate(i_RT_Start_list):
            i_rt_length = i_rt_start + 2
            if content[i_rt_length].startswith('RT_Mid_Sec'):
                i_RT_Mid_list.append(i_rt_length)

        # 计算每个raw文件对应的质荷比和电荷数的位置
        position_MH = i_Empty_Separator + 4 * n_raw_name + 1
        for i, MH in enumerate(content[position_MH: (position_MH + 2 * n_raw_name): 2]):
            if MH.startswith('MZ'):
                i_MH_list.append(position_MH + 2 * i)
        for i, i_MH in enumerate(i_MH_list):
            i_charge = i_MH + 1
            if content[i_charge].startswith('Charge'):
                i_Charge_list.append(i_charge)

        position_Score = i_Empty_Separator + 6 * n_raw_name + 1
        for i, Score in enumerate(content[position_Score: (position_Score + n_raw_name): 1]):
            if Score.startswith('Match Score'):
                i_Score_list.append(position_Score + i)

        for number, line in enumerate(lines):
            line = line.split('\t')
            seq_with_mod = line[i_SEQ_WITH_MOD]
            mod_dict = self.__soldierParseMOD2(seq_with_mod)
            MOD = ''
            for i_key in mod_dict.keys():
                i_value = mod_dict[i_key]
                MOD += str(i_key) + ',' + str(i_value) + ';'

            for i in range(n_raw_name):
                raw_name = raw_name_list[i]
                rt_start = float(line[i_RT_Start_list[i]]) * 60.0
                rt_length = float(line[i_RT_Length_list[i]]) * 60.0
                rt_mid = float(line[i_RT_Mid_list[i]]) * 60.0
                Moz = float(line[i_MH_list[i]])
                Charge = float(line[i_Charge_list[i]])
                Score = float(line[i_Charge_list[i]])

                self.dp.myIDForDIACHeck.ID00_RUN.append(line[i_RUN])
                self.dp.myIDForDIACHeck.ID1_RAW_NAME.append(raw_name)
                self.dp.myIDForDIACHeck.ID3_RT_APEX.append(rt_mid)
                self.dp.myIDForDIACHeck.ID4_SEQ_WITH_MOD.append(self.__soldierParseMOD(seq_with_mod))
                self.dp.myIDForDIACHeck.ID0_SEQ.append(line[i_SEQ])
                self.dp.myIDForDIACHeck.ID16_MOD.append(mod_dict)
                self.dp.myIDForDIACHeck.ID8_PRECURSOR_ID.append(raw_name + '_' + self.__soldierParseMOD(seq_with_mod) + '_' + str(Charge))
                self.dp.myIDForDIACHeck.ID5_PRECURSOR_MOZ_EXP.append(Moz)
                self.dp.myIDForDIACHeck.ID7_CHARGE.append(int(Charge))
                self.dp.myIDForDIACHeck.ID14_RT_BEGIN.append(rt_start)
                self.dp.myIDForDIACHeck.ID15_RT_END.append(rt_start + rt_length)

                self.dp.myIDForDIACHeck.ID9_SCORE0.append(Score)
                self.dp.myIDForDIACHeck.N_ID = self.dp.myIDForDIACHeck.N_ID + 1

class CFunctionReadDIANNIDForCheck:

    def  __init__(self, inputDP: CDataPack):

        self.dp = inputDP

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

    def __soldierParseMOD(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理

        return op_STANDARDIZE_MOD_BY_PRECURSOR(inputStr)

    def __soldierParseMOD2(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
        # 从DIA-NN格式的带修饰母离子字符串中提取修饰列表，每个元素存放修饰位点和类型
        seperator_start, seperator_end = '(', ')'
        mod_dic = {}
        mod_type = ''
        site_start, site_end = 0, 0
        flag_in_mod = 0
        count = 0
        for i, char in enumerate(inputStr):
            if char == seperator_start:
                site_start = i
                flag_in_mod = 1
            elif char == seperator_end:
                site_end = i
                try:
                    mod_type = UNIMOID_TO_STANDARD_MOD[inputStr[site_start + 1:site_end]]
                except KeyError:
                    print(inputStr[site_start + 1:site_end])
                    print('wrong')
                    info = 'DIA-NN modify type has wrong, please check!'
                    logGetError(info)
                mod_dic[count] = mod_type
                flag_in_mod = 0
            elif flag_in_mod == 0:
                count = count + 1

        return mod_dic

    def __soldierBinary_search(self, lines_Lib, target_precursor, i_PRECURSOR_Lib):

        start, end = -1, -1
        found = False

        low, high = 0, len(lines_Lib) - 1

        while low <= high:
            mid = (low + high) // 2
            current_precursor = toolGetWord(lines_Lib[mid], i_PRECURSOR_Lib, '\t')

            if current_precursor == target_precursor:
                start = mid
                found = True
                break
            elif current_precursor < target_precursor:
                low = mid + 1
            else:
                high = mid - 1

        if found:
            # 找到起始索引后，寻找结束索引
            for i in range(start, len(lines_Lib)):
                if toolGetWord(lines_Lib[i], i_PRECURSOR_Lib, '\t') != target_precursor:
                    end = i - 1
                    break
            else:
                end = len(lines_Lib) - 1

        return start, end

    def read(self, path):

        with open(path, 'r')as f:
            content = f.readline()
            lines = f.readlines()

        self.dp.myIDForDIACHeck.PSM_LINE_T.append(content.strip('\n'))

        i_RAW_NAME = toolGetIndexByWord(content, 'File.Name', '\t')
        i_RT = toolGetIndexByWord(content, 'RT', '\t')

        i_SEQ_WITH_MOD = toolGetIndexByWord(content, 'Modified.Sequence', '\t')
        i_SEQ = toolGetIndexByWord(content, 'Stripped.Sequence', '\t')
        i_CHARGE = toolGetIndexByWord(content, 'Precursor.Charge', '\t')
        i_PRECURSOR = toolGetIndexByWord(content, 'Precursor.Id', '\t')

        i_SCORE0 = toolGetIndexByWord(content, 'Q.Value', '\t')

        i_Protein_Group = toolGetIndexByWord(content, 'Protein.Group', '\t')
        i_Protein_ID = toolGetIndexByWord(content, 'Protein.Ids', '\t')
        i_Protein_Name = toolGetIndexByWord(content, 'Protein.Names', '\t')
        i_Gene = toolGetIndexByWord(content, 'Genes', '\t')

        i_RT_BEGIN = toolGetIndexByWord(content, 'RT.Start', '\t')
        i_RT_STOP = toolGetIndexByWord(content, 'RT.Stop', '\t')

        for line in lines:

            tmpFDR = float(toolGetWord(line, i_SCORE0, '\t'))

            seq_with_mod = toolGetWord(line, i_SEQ_WITH_MOD, '\t')

            seq = toolGetWord(line, i_SEQ, '\t')
            charge = int(toolGetWord(line, i_CHARGE, '\t'))
            mod_dict = self.__soldierParseMOD2(seq_with_mod)

            self.dp.myIDForDIACHeck.PSM_LINE_C.append(line.strip('\n'))

            self.dp.myIDForDIACHeck.ID0_SEQ.append(toolGetWord(line, i_SEQ, '\t'))

            self.dp.myIDForDIACHeck.ID00_RUN.append(toolGetNameFromPath(toolGetWord(line, i_RAW_NAME, '\t')))

            self.dp.myIDForDIACHeck.ID1_RAW_NAME.append(toolGetWord(line, i_RAW_NAME, '\t'))

            #单位是分钟
            self.dp.myIDForDIACHeck.ID3_RT_APEX.append(float(toolGetWord(line, i_RT, '\t'))*60)

            #seq+mod：mod在seq的对应位置中间，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK
            self.dp.myIDForDIACHeck.ID4_SEQ_WITH_MOD.append(self.__soldierParseMOD(toolGetWord(line, i_SEQ_WITH_MOD, '\t')))

            self.dp.myIDForDIACHeck.ID5_PRECURSOR_MOZ_EXP.append(self.__captainCalMOZ(seq,mod_dict,charge))
            self.dp.myIDForDIACHeck.ID6_PRECURSOR_MOZ_CLC.append(VALUE_ILLEGAL)

            self.dp.myIDForDIACHeck.ID7_CHARGE.append(int(toolGetWord(line, i_CHARGE, '\t')))

            # raw+seq+mod+charge：raw在最前面，mod在seq的对应位置中间，charge在seq后面，例如：CNCP2023_cpnl_480_hela_HRDIA_2_AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myIDForDIACHeck.ID8_PRECURSOR_ID.append(toolGetNameFromPath(toolGetWord(line, i_RAW_NAME, '\t')) + '_' + self.__soldierParseMOD(toolGetWord(line, i_PRECURSOR, '\t')))

            self.dp.myIDForDIACHeck.ID9_SCORE0.append(tmpFDR)

            self.dp.myIDForDIACHeck.ID10_Protein_Groups.append(toolGetWord(line, i_Protein_Group, '\t'))
            self.dp.myIDForDIACHeck.ID11_Protein_Ids.append(toolGetWord(line,i_Protein_ID,'\t'))
            self.dp.myIDForDIACHeck.ID12_Protein_Names.append(toolGetWord(line, i_Protein_Name, '\t'))
            self.dp.myIDForDIACHeck.ID13_Genes.append(toolGetWord(line, i_Gene, '\t'))

            self.dp.myIDForDIACHeck.ID14_RT_BEGIN.append(float(toolGetWord(line, i_RT_BEGIN, '\t')) * 60)
            self.dp.myIDForDIACHeck.ID15_RT_END.append(float(toolGetWord(line, i_RT_STOP, '\t')) * 60)

            # seq+mod+charge：mod在seq的对应位置中间，charge在seq后面，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myIDForDIACHeck.ID17_SEQ_WITH_MOD_AND_CHARGE.append(self.__soldierParseMOD(toolGetWord(line, i_PRECURSOR, '\t')))



            MOD = ''
            mod_type = ''
            mod_site = ''
            for i_key in mod_dict.keys():
                i_value = mod_dict[i_key]
                MOD += str(i_key) + ',' + str(i_value) + ';'
                mod_type += UNIMOID_TO_STANDARD_MOD2[str(i_value)] + ';'
                mod_site += str(i_key) + ';'

            #格式是：modsite,modtype;
            self.dp.myIDForDIACHeck.ID16_MOD.append(mod_dict)

            self.dp.myIDForDIACHeck.N_ID = self.dp.myIDForDIACHeck.N_ID + 1

class CFunctionReadDIANNLibIDForCheck:

    def  __init__(self, inputDP: CDataPack):

        self.dp = inputDP

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

    def __soldierParseMOD(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理

        return op_STANDARDIZE_MOD_BY_PRECURSOR(inputStr)

    def __soldierParseMOD2(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
        # 从DIA-NN格式的带修饰母离子字符串中提取修饰列表，每个元素存放修饰位点和类型
        seperator_start, seperator_end = '(', ')'
        mod_dic = {}
        mod_type = ''
        site_start, site_end = 0, 0
        flag_in_mod = 0
        count = 0
        for i, char in enumerate(inputStr):
            if char == seperator_start:
                site_start = i
                flag_in_mod = 1
            elif char == seperator_end:
                site_end = i
                try:
                    mod_type = UNIMOID_TO_STANDARD_MOD[inputStr[site_start + 1:site_end]]
                except KeyError:
                    print(inputStr[site_start + 1:site_end])
                    print('wrong')
                    info = 'DIA-NN modify type has wrong, please check!'
                    logGetError(info)
                mod_dic[count] = mod_type
                flag_in_mod = 0
            elif flag_in_mod == 0:
                count = count + 1

        return mod_dic

    def __soldierBinary_search(self, lines_Lib, target_precursor, i_PRECURSOR_Lib):

        start, end = -1, -1
        found = False

        low, high = 0, len(lines_Lib) - 1

        while low <= high:
            mid = (low + high) // 2
            current_precursor = toolGetWord(lines_Lib[mid], i_PRECURSOR_Lib, '\t')

            if current_precursor == target_precursor:
                start = mid
                found = True
                break
            elif current_precursor < target_precursor:
                low = mid + 1
            else:
                high = mid - 1

        if found:
            # 找到起始索引后，寻找结束索引
            for i in range(start, len(lines_Lib)):
                if toolGetWord(lines_Lib[i], i_PRECURSOR_Lib, '\t') != target_precursor:
                    end = i - 1
                    break
            else:
                end = len(lines_Lib) - 1

        return start, end

    def read(self, path):

        with open(path, 'r')as f:
            content = f.readline()
            lines = f.readlines()

        self.dp.myIDForDIACHeck.PSM_LINE_T.append(content.strip('\n'))

        i_RAW_NAME = toolGetIndexByWord(content, 'File.Name', '\t')
        i_RT = toolGetIndexByWord(content, 'RT', '\t')

        i_SEQ_WITH_MOD = toolGetIndexByWord(content, 'Modified.Sequence', '\t')
        i_SEQ = toolGetIndexByWord(content, 'Stripped.Sequence', '\t')
        i_CHARGE = toolGetIndexByWord(content, 'Precursor.Charge', '\t')
        i_PRECURSOR = toolGetIndexByWord(content, 'Precursor.Id', '\t')

        i_SCORE0 = toolGetIndexByWord(content, 'Q.Value', '\t')

        i_Protein_Group = toolGetIndexByWord(content, 'Protein.Group', '\t')
        i_Protein_ID = toolGetIndexByWord(content, 'Protein.Ids', '\t')
        i_Protein_Name = toolGetIndexByWord(content, 'Protein.Names', '\t')
        i_Gene = toolGetIndexByWord(content, 'Genes', '\t')

        i_RT_BEGIN = toolGetIndexByWord(content, 'RT.Start', '\t')
        i_RT_STOP = toolGetIndexByWord(content, 'RT.Stop', '\t')

        for line in lines:

            tmpFDR = float(toolGetWord(line, i_SCORE0, '\t'))

            seq_with_mod = toolGetWord(line, i_SEQ_WITH_MOD, '\t')

            seq = toolGetWord(line, i_SEQ, '\t')
            charge = int(toolGetWord(line, i_CHARGE, '\t'))
            mod_dict = self.__soldierParseMOD2(seq_with_mod)

            self.dp.myIDForDIACHeck.PSM_LINE_C.append(line.strip('\n'))

            self.dp.myIDForDIACHeck.ID0_SEQ.append(toolGetWord(line, i_SEQ, '\t'))

            self.dp.myIDForDIACHeck.ID00_RUN.append(toolGetNameFromPath(toolGetWord(line, i_RAW_NAME, '\t')))

            self.dp.myIDForDIACHeck.ID1_RAW_NAME.append(toolGetWord(line, i_RAW_NAME, '\t'))

            #单位是分钟
            self.dp.myIDForDIACHeck.ID3_RT_APEX.append(float(toolGetWord(line, i_RT, '\t'))*60)

            #seq+mod：mod在seq的对应位置中间，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK
            self.dp.myIDForDIACHeck.ID4_SEQ_WITH_MOD.append(self.__soldierParseMOD(toolGetWord(line, i_SEQ_WITH_MOD, '\t')))

            self.dp.myIDForDIACHeck.ID5_PRECURSOR_MOZ_EXP.append(self.__captainCalMOZ(seq,mod_dict,charge))
            self.dp.myIDForDIACHeck.ID6_PRECURSOR_MOZ_CLC.append(VALUE_ILLEGAL)

            self.dp.myIDForDIACHeck.ID7_CHARGE.append(int(toolGetWord(line, i_CHARGE, '\t')))

            # raw+seq+mod+charge：raw在最前面，mod在seq的对应位置中间，charge在seq后面，例如：CNCP2023_cpnl_480_hela_HRDIA_2_AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myIDForDIACHeck.ID8_PRECURSOR_ID.append(toolGetNameFromPath(toolGetWord(line, i_RAW_NAME, '\t')) + '_' + self.__soldierParseMOD(toolGetWord(line, i_PRECURSOR, '\t')))

            self.dp.myIDForDIACHeck.ID9_SCORE0.append(tmpFDR)

            self.dp.myIDForDIACHeck.ID10_Protein_Groups.append(toolGetWord(line, i_Protein_Group, '\t'))
            self.dp.myIDForDIACHeck.ID11_Protein_Ids.append(toolGetWord(line,i_Protein_ID,'\t'))
            self.dp.myIDForDIACHeck.ID12_Protein_Names.append(toolGetWord(line, i_Protein_Name, '\t'))
            self.dp.myIDForDIACHeck.ID13_Genes.append(toolGetWord(line, i_Gene, '\t'))

            self.dp.myIDForDIACHeck.ID14_RT_BEGIN.append(float(toolGetWord(line, i_RT_BEGIN, '\t')) * 60)
            self.dp.myIDForDIACHeck.ID15_RT_END.append(float(toolGetWord(line, i_RT_STOP, '\t')) * 60)

            # seq+mod+charge：mod在seq的对应位置中间，charge在seq后面，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myIDForDIACHeck.ID17_SEQ_WITH_MOD_AND_CHARGE.append(self.__soldierParseMOD(toolGetWord(line, i_PRECURSOR, '\t')))



            MOD = ''
            mod_type = ''
            mod_site = ''
            for i_key in mod_dict.keys():
                i_value = mod_dict[i_key]
                MOD += str(i_key) + ',' + str(i_value) + ';'
                mod_type += UNIMOID_TO_STANDARD_MOD2[str(i_value)] + ';'
                mod_site += str(i_key) + ';'

            #格式是：modsite,modtype;
            self.dp.myIDForDIACHeck.ID16_MOD.append(mod_dict)

            self.dp.myIDForDIACHeck.N_ID = self.dp.myIDForDIACHeck.N_ID + 1

class CFunctionReadMSFraggerIDForCheck:

    def  __init__(self, inputDP: CDataPack):

        self.dp = inputDP

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

    def __soldierParseMOD(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理

        return op_STANDARDIZE_MOD_BY_PRECURSOR(inputStr)


    def __soldierParseMOD2(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
        # 从DIA-NN格式的带修饰母离子字符串中提取修饰列表，每个元素存放修饰位点和类型
        seperator_start, seperator_end = '(', ')'
        mod_dic = {}
        mod_type = ''
        site_start, site_end = 0, 0
        flag_in_mod = 0
        count = 0
        for i, char in enumerate(inputStr):
            if char == seperator_start:
                site_start = i
                flag_in_mod = 1
            elif char == seperator_end:
                site_end = i
                try:
                    mod_type = UNIMOID_TO_STANDARD_MOD[inputStr[site_start + 1:site_end]]
                except KeyError:
                    print(inputStr[site_start + 1:site_end])
                    print('wrong')
                    info = 'DIA-NN modify type has wrong, please check!'
                    logGetError(info)
                mod_dic[count] = mod_type
                flag_in_mod = 0
            elif flag_in_mod == 0:
                count = count + 1

        return mod_dic

    def __soldierBinary_search(self, lines_Lib, target_precursor, i_PRECURSOR_Lib):

        start, end = -1, -1
        found = False

        low, high = 0, len(lines_Lib) - 1

        while low <= high:
            mid = (low + high) // 2
            current_precursor = toolGetWord(lines_Lib[mid], i_PRECURSOR_Lib, '\t')

            if current_precursor == target_precursor:
                start = mid
                found = True
                break
            elif current_precursor < target_precursor:
                low = mid + 1
            else:
                high = mid - 1

        if found:
            # 找到起始索引后，寻找结束索引
            for i in range(start, len(lines_Lib)):
                if toolGetWord(lines_Lib[i], i_PRECURSOR_Lib, '\t') != target_precursor:
                    end = i - 1
                    break
            else:
                end = len(lines_Lib) - 1

        return start, end

    def read(self, path):

        with open(path, 'r')as f:
            content = f.readline()
            lines = f.readlines()

        self.dp.myIDForDIACHeck.PSM_LINE_T.append(content.strip('\n'))

        i_RT = toolGetIndexByWord(content, 'AverageExperimentalRetentionTime', '\t')

        i_SEQ_WITH_MOD = toolGetIndexByWord(content, 'ModifiedPeptideSequence', '\t')
        i_SEQ = toolGetIndexByWord(content, 'PeptideSequence', '\t')
        i_PRECURSOR_MOZ_EXP = toolGetIndexByWord(content, 'PrecursorMz', '\t')
        i_CHARGE = toolGetIndexByWord(content, 'PrecursorCharge', '\t')
        # i_SCORE0 = toolGetIndexByWord(content, 'QValue', '\t')

        i_Gene = toolGetIndexByWord(content, 'GeneName', '\t')
        i_Protein_Name = toolGetIndexByWord(content, 'ProteinId', '\t')

        i_FRAGMENT_TYPE_1 = toolGetIndexByWord(content, 'FragmentType', '\t')
        i_FRAGMENT_TYPE_3 = toolGetIndexByWord(content, 'FragmentCharge', '\t')
        i_FRAGMENT_TYPE_2 = toolGetIndexByWord(content, 'FragmentSeriesNumber', '\t')
        i_FRAGMENT_MOZ_EXP = toolGetIndexByWord(content, 'ProductMz', '\t')
        i_FRAGMENT_LOSS_TYPE = toolGetIndexByWord(content, 'Annotation', '\t')
        pre_transition_group_id = '-1'  # 用于判断是否读到下一个transition_group
        fragment_moz_list = []  # 暂存当前母离子对应子离子m/z的列表
        fragment_type_list = []  # 暂存当前母离子对应子离子类型的列表
        fragment_loss_list = []  # 暂存当前母离子对应子离子的损失类型的列表

        for line in tqdm(lines, desc="Read MSFragger Library"):
        # for line in lines:

            tmpFDR = 0
            transition_group_id = toolGetWord(line, i_SEQ_WITH_MOD, '\t') + '-' + toolGetWord(line, i_CHARGE, '\t')

            if pre_transition_group_id != transition_group_id:

                pre_transition_group_id = transition_group_id

                if len(fragment_moz_list) != 0:
                    self.dp.myIDForDIACHeck.ID18_FRAGMENT_TYPE.append(fragment_type_list)
                    self.dp.myIDForDIACHeck.ID19_FRAGMENT_MOZ_EXP.append(fragment_moz_list)
                    self.dp.myIDForDIACHeck.ID21_FRAGMENT_LOSS_TYPE.append(fragment_loss_list)
                fragment_type_list = []
                fragment_moz_list = []
                fragment_loss_list = []
                fragment_type_1 = toolGetWord(line, i_FRAGMENT_TYPE_1, '\t')
                fragment_type_2 = toolGetWord(line, i_FRAGMENT_TYPE_2, '\t')
                fragment_type_3 = toolGetWord(line, i_FRAGMENT_TYPE_3, '\t')
                if int(fragment_type_3) >= 3:
                    pass
                else:

                    fragment_moz_list.append(float(toolGetWord(line, i_FRAGMENT_MOZ_EXP, '\t')))
                    loss_type = toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t')
                    if '-' not in loss_type:
                        fragment_loss_list.append('NOLOSS')
                    else:
                        fragment_loss_list.append(toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t'))
                raw_name = self.dp.myCFG.A1_PATH_MS1.split('\\')[-1].replace('.ms1', '')

                seq_with_mod = toolGetWord(line, i_SEQ_WITH_MOD, '\t')
                seq = toolGetWord(line, i_SEQ, '\t')
                charge = int(toolGetWord(line, i_CHARGE, '\t'))

                mod_dict = self.__soldierParseMOD2(seq_with_mod)

                self.dp.myIDForDIACHeck.PSM_LINE_C.append(line.strip('\n'))

                self.dp.myIDForDIACHeck.ID0_SEQ.append(toolGetWord(line, i_SEQ, '\t'))

                self.dp.myIDForDIACHeck.ID00_RUN.append(raw_name)

                self.dp.myIDForDIACHeck.ID1_RAW_NAME.append(raw_name)

                #单位是分钟
                self.dp.myIDForDIACHeck.ID3_RT_APEX.append(float(toolGetWord(line, i_RT, '\t')))

                #seq+mod：mod在seq的对应位置中间，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK
                self.dp.myIDForDIACHeck.ID4_SEQ_WITH_MOD.append(self.__soldierParseMOD(toolGetWord(line, i_SEQ_WITH_MOD, '\t')))

                self.dp.myIDForDIACHeck.ID5_PRECURSOR_MOZ_EXP.append(float(toolGetWord(line, i_PRECURSOR_MOZ_EXP, '\t')))
                self.dp.myIDForDIACHeck.ID6_PRECURSOR_MOZ_CLC.append(VALUE_ILLEGAL)

                self.dp.myIDForDIACHeck.ID7_CHARGE.append(int(toolGetWord(line, i_CHARGE, '\t')))

                # raw+seq+mod+charge：raw在最前面，mod在seq的对应位置中间，charge在seq后面，例如：CNCP2023_cpnl_480_hela_HRDIA_2_AAAGDLGGDHLAFSC(UniMod:4)DVAK3
                self.dp.myIDForDIACHeck.ID8_PRECURSOR_ID.append(raw_name + '_' + self.__soldierParseMOD(seq_with_mod) + str(charge))

                self.dp.myIDForDIACHeck.ID9_SCORE0.append(tmpFDR)

                self.dp.myIDForDIACHeck.ID12_Protein_Names.append(toolGetWord(line, i_Protein_Name, '\t'))
                self.dp.myIDForDIACHeck.ID13_Genes.append(toolGetWord(line, i_Gene, '\t'))

                # seq+mod+charge：mod在seq的对应位置中间，charge在seq后面，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK3
                self.dp.myIDForDIACHeck.ID17_SEQ_WITH_MOD_AND_CHARGE.append(self.__soldierParseMOD(seq_with_mod) + str(charge))

                MOD = ''
                mod_type = ''
                mod_site = ''
                for i_key in mod_dict.keys():
                    i_value = mod_dict[i_key]
                    MOD += str(i_key) + ',' + str(i_value) + ';'
                    mod_type += UNIMOID_TO_STANDARD_MOD2[str(i_value)] + ';'
                    mod_site += str(i_key) + ';'

                #格式是：modsite,modtype;
                self.dp.myIDForDIACHeck.ID16_MOD.append(mod_dict)

                self.dp.myIDForDIACHeck.N_ID = self.dp.myIDForDIACHeck.N_ID + 1

            else:
                fragment_type_1 = toolGetWord(line, i_FRAGMENT_TYPE_1, '\t')
                fragment_type_2 = toolGetWord(line, i_FRAGMENT_TYPE_2, '\t')
                fragment_type_3 = toolGetWord(line, i_FRAGMENT_TYPE_3, '\t')
                if int(fragment_type_3) >= 3:
                    pass
                else:
                    fragment_type_list.append(fragment_type_1 + fragment_type_2 + '+' * int(fragment_type_3))
                    fragment_moz_list.append(float(toolGetWord(line, i_FRAGMENT_MOZ_EXP, '\t')))
                    loss_type = toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t')
                    if '-' not in loss_type:
                        fragment_loss_list.append('NOLOSS')
                    else:
                        fragment_loss_list.append(toolGetWord(line, i_FRAGMENT_LOSS_TYPE, '\t'))

        if len(fragment_moz_list) != 0:
            self.dp.myIDForDIACHeck.ID18_FRAGMENT_TYPE.append(fragment_type_list)
            self.dp.myIDForDIACHeck.ID19_FRAGMENT_MOZ_EXP.append(fragment_moz_list)
            self.dp.myIDForDIACHeck.ID20_FRAGMENT_MOZ_CLC.append([])
            self.dp.myIDForDIACHeck.ID21_FRAGMENT_LOSS_TYPE.append(fragment_loss_list)

class CFunctionReadMSFraggerPinIDForCheck:

    def  __init__(self, inputDP: CDataPack):

        self.dp = inputDP

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

    def __soldierParseMOD(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
        inputStr = re.sub(r'[a-z]', '', inputStr)
        return op_STANDARDIZE_MOD_BY_PRECURSOR(inputStr,flag='[]')

    def __soldierParseMOD2(self, inputStr):
        # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
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

    def __soldierBinary_search(self, lines_Lib, target_precursor, i_PRECURSOR_Lib):

        start, end = -1, -1
        found = False

        low, high = 0, len(lines_Lib) - 1

        while low <= high:
            mid = (low + high) // 2
            current_precursor = toolGetWord(lines_Lib[mid], i_PRECURSOR_Lib, '\t')

            if current_precursor == target_precursor:
                start = mid
                found = True
                break
            elif current_precursor < target_precursor:
                low = mid + 1
            else:
                high = mid - 1

        if found:
            # 找到起始索引后，寻找结束索引
            for i in range(start, len(lines_Lib)):
                if toolGetWord(lines_Lib[i], i_PRECURSOR_Lib, '\t') != target_precursor:
                    end = i - 1
                    break
            else:
                end = len(lines_Lib) - 1

        return start, end

    def read(self, path):

        with open(path, 'r')as f:
            content = f.readline()
            lines = f.readlines()

        self.dp.myIDForDIACHeck.PSM_LINE_T.append(content.strip('\n'))

        i_RAW_NAME = toolGetIndexByWord(content, 'SpecId', '\t')
        i_RT = toolGetIndexByWord(content, 'retentiontime', '\t')

        i_SEQ_WITH_MOD = toolGetIndexByWord(content, 'Peptide', '\t')

        i_Protein_Name = toolGetIndexByWord(content, 'Proteins', '\t')

        i_DECOY_FLAG = toolGetIndexByWord(content, 'Label', '\t')

        for line in tqdm(lines, desc="Read MSFragger Pin"):
        # for line in lines:

            tmpFDR = 0

            seq_with_mod = toolGetWord(line, i_SEQ_WITH_MOD, '\t')
            raw_name = toolGetWord(line, i_RAW_NAME, '\t').split('.')[0] + '.raw'
            seq = ''.join(filter(str.isupper, seq_with_mod[2:-3]))
            charge = int(seq_with_mod[-3])
            mod_dict = self.__soldierParseMOD2(seq_with_mod)

            self.dp.myIDForDIACHeck.PSM_LINE_C.append(line.strip('\n'))

            self.dp.myIDForDIACHeck.ID0_SEQ.append(seq)

            self.dp.myIDForDIACHeck.ID00_RUN.append(raw_name)

            self.dp.myIDForDIACHeck.ID1_RAW_NAME.append(raw_name)

            #单位是分钟
            self.dp.myIDForDIACHeck.ID3_RT_APEX.append(float(toolGetWord(line, i_RT, '\t'))*60)

            #seq+mod：mod在seq的对应位置中间，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK
            self.dp.myIDForDIACHeck.ID4_SEQ_WITH_MOD.append(self.__soldierParseMOD(seq_with_mod[2:-3]))

            self.dp.myIDForDIACHeck.ID5_PRECURSOR_MOZ_EXP.append(self.__captainCalMOZ(seq,mod_dict,charge))
            self.dp.myIDForDIACHeck.ID6_PRECURSOR_MOZ_CLC.append(VALUE_ILLEGAL)

            self.dp.myIDForDIACHeck.ID7_CHARGE.append(charge)

            # raw+seq+mod+charge：raw在最前面，mod在seq的对应位置中间，charge在seq后面，例如：CNCP2023_cpnl_480_hela_HRDIA_2_AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myIDForDIACHeck.ID8_PRECURSOR_ID.append(raw_name + '_' + self.__soldierParseMOD(seq_with_mod[2:-3] + str(charge)))

            self.dp.myIDForDIACHeck.ID9_SCORE0.append(tmpFDR)

            # self.dp.myIDForDIACHeck.ID10_Protein_Groups.append(toolGetWord(line, i_Protein_Group, '\t'))
            # self.dp.myIDForDIACHeck.ID11_Protein_Ids.append(toolGetWord(line,i_Protein_ID,'\t'))
            self.dp.myIDForDIACHeck.ID12_Protein_Names.append(toolGetWord(line, i_Protein_Name, '\t'))
            # self.dp.myIDForDIACHeck.ID13_Genes.append(toolGetWord(line, i_Gene, '\t'))
            self.dp.myIDForDIACHeck.ID22_TARGET.append(int(toolGetWord(line, i_DECOY_FLAG, '\t')))
            # self.dp.myIDForDIACHeck.ID14_RT_BEGIN.append(float(toolGetWord(line, i_RT_BEGIN, '\t')) * 60)
            # self.dp.myIDForDIACHeck.ID15_RT_END.append(float(toolGetWord(line, i_RT_STOP, '\t')) * 60)

            # seq+mod+charge：mod在seq的对应位置中间，charge在seq后面，例如：AAAGDLGGDHLAFSC(UniMod:4)DVAK3
            self.dp.myIDForDIACHeck.ID17_SEQ_WITH_MOD_AND_CHARGE.append(self.__soldierParseMOD(toolGetWord(line, i_SEQ_WITH_MOD, '\t') + str(charge)))

            MOD = ''
            mod_type = ''
            mod_site = ''
            for i_key in mod_dict.keys():
                i_value = mod_dict[i_key]
                MOD += str(i_key) + ',' + str(i_value) + ';'
                mod_type += UNIMOID_TO_STANDARD_MOD2[str(i_value)] + ';'
                mod_site += str(i_key) + ';'

            #格式是：modsite,modtype;
            self.dp.myIDForDIACHeck.ID16_MOD.append(mod_dict)

            self.dp.myIDForDIACHeck.N_ID = self.dp.myIDForDIACHeck.N_ID + 1