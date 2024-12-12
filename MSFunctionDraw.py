import pickle
import time
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.pyplot import subplots
import matplotlib.patches as pc
from scipy.interpolate import make_interp_spline, BSpline, CubicSpline
from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import MultipleLocator
from mpl_toolkits.mplot3d import Axes3D
from MSData import CDataPack, CSeed, CFileMS1, CEvidence, CFileMS2, CSeed2,\
    CEvidence2, CSeed_FLOW1, CEvidence_FLOW1, CSeedForDDA, CINI,\
    CSeedDIACheck, CEvidenceDIACheck
from MSTool import toolFindNeighborFromSortedList1, toolGetNameFromPath, toolPearsonForNumpy, \
    tooFindNeighborListFromSortedList, toolCosineForMatrixAndList, toolMaxContinueBY
from MSLogging import logGetError, logGetWarning, logToUser
from MSOperator import op_FILL_LIST_DRAW_PRECURSOR, op_CAL_FRAGMENT_MOZ,\
    op_CAL_PERCURSOR_MOZ, op_FILL_LIST_DRAW_DDA_PRECURSOR, op_DIVIDE_MOD_FROM_PRECURSOR,\
    op_INIT_CSEED_DIACHECK, op_CAL_MOLECULAR_MASS_FROM_FRAGMENT_LOSS, op_DIVIDE_ION_TYPE,\
    op_Data_fill_plot
from MSFunctionEvidence import CFunctionEvidence, CFunctionEvidence2, CFunctionEvidenceDDA,\
    CFunctionEvidenceForDIACheck
from MSFunction import CFunctionParseMS1, CFunctionParseMS2
from MSFunctionComposition import CFunctionComposition
from MSEmass import CEmass
from MSSystem import VALUE_ILLEGAL, IO_NAME_FILE_EXPORT_Flow2, CFG_TYPE_PLOT_NUM,\
    IO_NAME_FILE_EXPORT_XIC_MATRIX, PLOT_COLOR, IO_NAME_FILE_EXPORT_Flow1, DRAW_COLOR,\
    FRAGMENT_LOSS_TYPE_COMP, IO_NAME_FILE_EXPORT_Flow5, CFG_TYPE_ACCURACY_HALF_WIN_PEAK, \
    IO_NAME_FILE_EXPORT_Flow6


class CFunctionDrawDIAXIC:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soliderCalMOZ(self, iID):

        ionPrecursor = self.dp.myID.ID4_SEQ_WITH_MOD[iID]
        ionCharge = self.dp.myID.ID7_CHARGE[iID]

        ionMoz = op_CAL_PERCURSOR_MOZ(self.dp.myINI, ionPrecursor, ionCharge)

        return ionMoz

    def __soliderCalFragMoz(self, iID):
        ionPrecursor = self.dp.myID.ID4_SEQ_WITH_MOD[iID]

        if self.dp.myID.ID10_FRAGMENT_TYPE:

            ionTypeList = self.dp.myID.ID10_FRAGMENT_TYPE[iID]
            ionLossList = self.dp.myID.ID13_FRAGMENT_LOSS_TYPE[iID]

            ionTypeList, ionMozList = op_CAL_FRAGMENT_MOZ(self.dp.myINI, ionPrecursor, ionTypeList, ionLossList)

        else:

            length = len(self.dp.myID.ID0_SEQ[iID])
            ionLossList = ['NOLOSS'] * (length-1) * 4
            ionTypeList = []
            for i in range(length-1):
                ionTypeList.append('b' + str(i+1) + '+')
                ionTypeList.append('y' + str(i+1) + '+')
                ionTypeList.append('b' + str(i+1) + '++')
                ionTypeList.append('y' + str(i+1) + '++')

            ionTypeList, ionMozList = op_CAL_FRAGMENT_MOZ(self.dp.myINI, ionPrecursor, ionTypeList, ionLossList)

        return ionTypeList, ionMozList

    def __captainFillSeed(self, iID, dataMS2: CFileMS2):

        # 初始化seed
        myseed = CSeed_FLOW1()

        # 根据rt得到ms1中对应的索引
        if self.dp.myID.ID2_SCAN_ID[iID] != VALUE_ILLEGAL:
            indexMS2 = toolFindNeighborFromSortedList1(dataMS2.INDEX_SCAN, self.dp.myID.ID2_SCAN_ID[iID])
        elif self.dp.myID.ID3_RT[iID] != VALUE_ILLEGAL:
            indexMS2 = toolFindNeighborFromSortedList1(dataMS2.INDEX_RT, self.dp.myID.ID3_RT[iID] * 60)
        # elif self.dp.myID.ID2_SCAN_ID[iID] != VALUE_ILLEGAL:
        #     indexMS2 = toolFindNeighborFromSortedList1(dataMS2.INDEX_SCAN, self.dp.myID.ID2_SCAN_ID[iID])
        else:
            indexMS2 = VALUE_ILLEGAL
            info = 'can not index rt or scan'
            logGetError(info)

        myseed.MID_RT = dataMS2.INDEX_RT[indexMS2]
        myseed.MID_SCAN = dataMS2.INDEX_SCAN[indexMS2]

        charge_prec = self.dp.myID.ID7_CHARGE[iID]
        Mono_moz = self.__soliderCalMOZ(iID)
        isotope1_moz = (Mono_moz * charge_prec + 1.0) / charge_prec
        isotope2_moz = (Mono_moz * charge_prec + 2.0) / charge_prec
        isotope3_moz = (Mono_moz * charge_prec + 3.0) / charge_prec
        # isotope4_moz = (Mono_moz * charge_prec + 4.0) / charge_prec

        myseed.MOZ_PREC = [Mono_moz, isotope1_moz, isotope2_moz,isotope3_moz]
        myseed.PREC_TYPE = ['precursor', 'precursor[M+1]', 'precursor[M+2]','precursor[M+3]']

        myseed.FRAG_TYPE, myseed.MOZ_FRAG = self.__soliderCalFragMoz(iID)

        return myseed

    def __captainGetEvidence(self, inputListEvidence, inputListSeed, inputListIndex):
        '''
        :param inputListEvidence: 输入空的evidence列表用于填充，得到的evidence用于画图
        :param inputListSeed: 输入空的seed列表用于填充，得到的seed用于对应evidence的获取
        :param inputListIndex: 输入鉴定结果列表的索引列表，去鉴定结果数据结构中寻址得到rt，母离子等信息
        :return:
        '''
        nMS1 = len(self.dp.LIST_PATH_MS1)
        nMS2 = len(self.dp.LIST_PATH_MS2)

        for iMS1 in range(nMS1):

            # init ms1
            pathMS1 = self.dp.LIST_PATH_MS1[iMS1]
            nameRawMS1 = toolGetNameFromPath(pathMS1)

            # load
            dataMS1 = CFileMS1()
            functionParseMS1 = CFunctionParseMS1(self.dp)
            dataMS1 = functionParseMS1.loadPKL(pathMS1)

            dataMS2 = ''

            for iMS2 in range(nMS2):

                # init ms2
                pathMS2 = self.dp.LIST_PATH_MS2[iMS2]
                nameRawMS2 = toolGetNameFromPath(pathMS2)

                # load
                if nameRawMS2 == nameRawMS1:
                    dataMS2 = CFileMS2()
                    functionParseMS2 = CFunctionParseMS2(self.dp)
                    dataMS2 = functionParseMS2.loadPKL(pathMS2)

                    for i in range(len(inputListIndex)):

                        iID = inputListIndex[i]
                        nameRawID = toolGetNameFromPath(self.dp.myID.ID1_RAW_NAME[iID])

                        if nameRawID == nameRawMS1 or self.dp.myID.ID1_RAW_NAME[iID] == nameRawMS1:

                            tmpSeed = self.__captainFillSeed(iID, dataMS2)

                            inputListSeed[i] = tmpSeed

                            functionEvidence = CFunctionEvidence(self.dp)

                            tmpEvidence = functionEvidence.fillEvidence1(dataMS1, dataMS2, tmpSeed)
                            inputListEvidence[i] = tmpEvidence


    def __captaindrawPrecursor(self, inputEvidence_list: list, inputIndex_list: list):

        fig = plt.figure(figsize=(20, 18))
        gs = GridSpec(16, 12,hspace=0.8)
        text_base_title = 'Base Peak:'
        text_precursor_title = 'Precursor:'
        text_raw_title = 'Raw File:'
        text_qvalue_title = 'q value:'
        text_charge_title = 'Charge:'
        text_base_prec_title = 'Base Peak Prec:'

        for tmp_i in range(len(inputEvidence_list)):

            if tmp_i % 2000 == 0:
                print(tmp_i,'/',len(inputEvidence_list))

            inputEvidence = CEvidence_FLOW1()
            inputEvidence = inputEvidence_list[tmp_i]
            inputIndex = inputIndex_list[tmp_i]
            rt_apex = self.dp.myID.ID3_RT[inputIndex]
            precursor_id = self.dp.myID.ID8_PRECURSOR_ID[inputIndex]
            n_fragment = len(inputEvidence.FRAG_TYPE)
            n_iostope = len(inputEvidence.PREC_TYPE)

            if len(inputEvidence.MATRIX_PROFILE_FRAG) != 0:
                base_intensity = np.max(inputEvidence.MATRIX_PROFILE_FRAG)
                if base_intensity > 0.01:
                    matrix_profile_frag = inputEvidence.MATRIX_PROFILE_FRAG / base_intensity
                else:
                    matrix_profile_frag = inputEvidence.MATRIX_PROFILE_FRAG
                ax_frag = plt.subplot(gs[9:, 0:9])
                frag_type_text = []

                for i in range(n_fragment):
                    x_list = inputEvidence.LIST_RET_TIME
                    x_list = [i / 60.0 for i in x_list]
                    y_list = matrix_profile_frag[i]
                    plt.plot(x_list, y_list, color=DRAW_COLOR[i],)
                    frag_type_text.append(
                        inputEvidence.FRAG_TYPE[i] + ':' + format(inputEvidence.FRAG_MOZ[i], '.2f') + ' m/z')
                    # plt.scatter(x_list, y_list, color=DRAW_COLOR[i])

                xlim = [np.min(x_list), np.max(x_list)]
                plt.xlim(xlim)
                plt.ylim([0, 1])
                # plot rt apex line, rt start line, rt end line
                if inputEvidence.RT_START != VALUE_ILLEGAL:
                    rt_start = inputEvidence.LIST_RET_TIME[inputEvidence.RT_START] / 60.0
                    rt_end = inputEvidence.LIST_RET_TIME[inputEvidence.RT_END] / 60.0
                    # plt.plot([rt_apex, rt_apex], [0, 1], color='b')
                    # plt.plot([rt_start, rt_start], [0, 1], color='r', linestyle=':')
                    # plt.plot([rt_end, rt_end], [0, 1], color='r', linestyle=':')
                fontsize = 18
                plt.xticks(fontsize=fontsize)
                plt.yticks(fontsize=fontsize)
                plt.xlabel('RT (min)',fontsize=fontsize)
                plt.ylabel('Intensity (%)',fontsize=fontsize)

                ax_frag_type = fig.add_subplot(gs[9:, 9:])
                plt.xlim(0, 1)
                plt.ylim(0, 1)
                for i, text in enumerate(frag_type_text):
                    ax_frag_type.add_patch(
                        pc.Rectangle(  # 长方形
                            (0.01, 0.95 - i * 0.05),  # （x,y）
                            0.06,  # 长
                            0.01,  # 宽
                            color=DRAW_COLOR[i]
                        )
                    )
                    plt.text(0.1, 0.95 - i * 0.05, text, fontsize=fontsize)
                plt.axis('off')

            else:
                base_intensity = 0

            base_intensity_prec = np.max(inputEvidence.MATRIX_PROFILE_PREC)
            base_intensity_prec = base_intensity_prec if base_intensity_prec > 0 else 1
            matrix_profile_prec = inputEvidence.MATRIX_PROFILE_PREC/base_intensity_prec

            ax = plt.subplot(gs[0:2, 0:])
            text_base_value = format(base_intensity, '.2e')
            text_base_prec_value = format(base_intensity_prec, '.2e')
            text_precursor_value = self.dp.myID.ID4_SEQ_WITH_MOD[inputIndex]
            text_raw_value = toolGetNameFromPath(self.dp.myID.ID1_RAW_NAME[inputIndex])
            text_charge_value = str(self.dp.myID.ID7_CHARGE[inputIndex])
            text_qvalue_value = format(self.dp.myID.ID9_SCORE0[inputIndex], '.5f')

            plt.xlim((0, 200))
            plt.ylim((0, 20))
            plt.xticks([])
            plt.yticks([])

            plt.text(2, 15, text_base_title, fontdict={'size': fontsize, 'color': 'black'})
            plt.text(23, 15, text_base_value, fontdict={'size': fontsize, 'color': 'red'})
            plt.text(2, 0, text_raw_title, fontdict={'size': fontsize, 'color': 'black'})
            plt.text(23, 0, text_raw_value, fontdict={'size': fontsize, 'color': 'red'})
            plt.text(2, 8, text_precursor_title, fontdict={'size': fontsize, 'color': 'black'})
            plt.text(23, 8, text_precursor_value, fontdict={'size': fontsize, 'color': 'red'})
            plt.text(80, 0, text_charge_title, fontdict={'size': fontsize, 'color': 'black'})
            plt.text(94, 0, text_charge_value, fontdict={'size': fontsize, 'color': 'red'})
            plt.text(102, 0, text_qvalue_title, fontdict={'size': fontsize, 'color': 'black'})
            plt.text(116, 0, text_qvalue_value, fontdict={'size': fontsize, 'color': 'red'})
            plt.text(52, 15, text_base_prec_title, fontdict={'size': fontsize, 'color': 'black'})
            plt.text(80, 15, text_base_prec_value, fontdict={'size': fontsize, 'color': 'red'})
            plt.axis('off')
            ax_prec = plt.subplot(gs[2:9, 0:9])
            prec_type_text = []
            for i in range(n_iostope):
                x_list = inputEvidence.LIST_RET_TIME
                x_list = [i / 60.0 for i in x_list]
                y_list = matrix_profile_prec[i]
                plt.plot(x_list, y_list, color=DRAW_COLOR[i])
                prec_type_text.append(
                    inputEvidence.PREC_TYPE[i] + ':' + format(inputEvidence.PREC_MOZ[i], '.2f') + ' m/z')
                # plt.scatter(x_list, y_list, color=DRAW_COLOR[i])
            plt.xticks(fontsize=fontsize)
            plt.yticks(fontsize=fontsize)
            xlim = [np.min(x_list), np.max(x_list)]
            plt.xlim(xlim)
            # plot rt apex line, rt start line, rt end line
            if inputEvidence.RT_START != VALUE_ILLEGAL:
                rt_start = inputEvidence.LIST_RET_TIME[inputEvidence.RT_START] / 60.0
                rt_end = inputEvidence.LIST_RET_TIME[inputEvidence.RT_END] / 60.0
                # plt.plot([rt_apex, rt_apex], [0, 1], color='b')
                # plt.plot([rt_start, rt_start], [0, 1], color='r', linestyle=':')
                # plt.plot([rt_end, rt_end], [0, 1], color='r', linestyle=':')
            # plt.xlabel('RT (min)')
            plt.ylabel('Intensity (%)',fontsize=fontsize)

            ax_prec_type = fig.add_subplot(gs[2:9, 9:])

            plt.xlim(0, 1)
            plt.ylim(0, 1)
            for i, text in enumerate(prec_type_text):
                ax_prec_type.add_patch(
                    pc.Rectangle(  # 长方形
                        (0.01, 0.95 - i * 0.1),  # （x,y）
                        0.06,  # 长
                        0.01,  # 宽
                        color=DRAW_COLOR[i]
                    )
                )
                plt.text(0.1, 0.95 - i * 0.1, text, fontsize=fontsize)
            plt.axis('off')
            plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + precursor_id + '_' +
                        IO_NAME_FILE_EXPORT_Flow1[0])

            # plt.show()
            plt.clf()
        plt.close()

    def __captainParseMOD(self, inputStr):  # 各个鉴定结果导出的修饰格式不一样，搞成一样的格式，方便后续处理
        # 从DIA-NN格式的带修饰母离子字符串中提取修饰列表，每个元素存放修饰位点和类型
        trantab = str.maketrans('()', '[]')
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
                    mod_type = inputStr[site_start + 1:site_end].translate(trantab).replace(' ', '')
                except KeyError:
                    info = 'DIA-NN modify type has wrong, please check!'
                    logGetError(info)
                if count == 0:
                    mod_dic[count] = mod_type
                else:
                    mod_dic[count-1] = mod_type
                flag_in_mod = 0
            elif flag_in_mod == 0:
                count = count + 1
        return mod_dic

    def __captainOutputXLC(self, listEvidence):
        pre = []
        frg = []

        for j,evi in enumerate(listEvidence):

            if j%2000 == 0:
                print(j,'/',len(listEvidence))

            tmp_frg = []
            tmp_pre = []

            tmp_frg.append(evi.MATRIX_PROFILE_FRAG)
            tmp_frg.append(evi.FRAG_MOZ)
            tmp_frg.append(evi.FRAG_TYPE)
            rt = [i/60 for i in evi.LIST_RET_TIME]
            tmp_frg.append(rt)
            tmp_frg.append(self.dp.myID.ID0_SEQ[j])
            tmp_frg.append(self.dp.myID.ID7_CHARGE[j])
            mod = self.__captainParseMOD(self.dp.myID.ID4_SEQ_WITH_MOD[j])
            tmp_frg.append(mod)
            tmp_frg.append(self.dp.myID.ID17_TARGET[j])

            tmp_id = np.array(evi.PREC_TYPE).argsort()
            tmp_pre.append(evi.MATRIX_PROFILE_PREC[tmp_id])
            tmp_pre.append(np.array(evi.PREC_MOZ)[tmp_id])
            tmp_pre.append(rt)
            tmp_pre.append(self.dp.myID.ID0_SEQ[j])
            tmp_pre.append(self.dp.myID.ID7_CHARGE[j])
            tmp_pre.append(mod)
            tmp_pre.append(self.dp.myID.ID17_TARGET[j])

            frg.append(tmp_frg)
            pre.append(tmp_pre)
        return frg,pre

    def draw(self):
        # 得到要画图的结果列表
        listPrecursor_id = []

        if self.dp.myCFG.C7_PLOT_NUM == CFG_TYPE_PLOT_NUM['ALL']:  # 画identification文件中所有鉴定结果的色谱曲线

            listPrecursor_id = self.dp.myID.ID8_PRECURSOR_ID

        elif self.dp.myCFG.C7_PLOT_NUM == CFG_TYPE_PLOT_NUM['CUSTOM']:  # 根据自定义输入的鉴定结果画色谱曲线
            listPrecursor_id = op_FILL_LIST_DRAW_PRECURSOR(self.dp.myCFG.C4_DIA_PRECURSOR_SEQUENCE,
                                                           self.dp.myCFG.C5_DIA_PRECURSOR_CHARGE,
                                                           self.dp.myCFG.C6_DIA_PRECURSOR_RAW)
        else:
            info = "MSFunctionDraw, MK449, " + self.dp.myCFG.C7_PLOT_NUM + " is all right?"
            logGetError(info)

        listIndex = []
        for precursor_id in listPrecursor_id:

            if precursor_id in self.dp.myID.ID8_PRECURSOR_ID:
                index = self.dp.myID.ID8_PRECURSOR_ID.index(precursor_id)
                listIndex.append(index)
            else:
                pass

        listEvidence = [[] * 1] * len(listIndex)  # 在函数中是全局的
        listSeed = [[] * 1] * len(listIndex)

        self.__captainGetEvidence(listEvidence, listSeed, listIndex)

        frag,pre = self.__captainOutputXLC(listEvidence)

        with open(self.dp.myCFG.D1_PATH_EXPORT+'pre.pkl','wb') as f:
            pickle.dump(pre,f)

        with open(self.dp.myCFG.D1_PATH_EXPORT+'frag.pkl','wb') as f:
            pickle.dump(frag,f)

        if self.dp.myCFG.C17_DRAW_FLAG:
            self.__captaindrawPrecursor(listEvidence, listIndex)




class CFunctionExtractAllScanForMS1:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soldierfillListIsotope(self, inputIsotope: str):

        if inputIsotope[-1] == '|':
            pass
        else:
            inputIsotope = inputIsotope + '|'

        isotopeList = inputIsotope.split('|')[:-1]
        isotopeList = [int(i) for i in isotopeList]

        return isotopeList

    def __captainCalSignal(self, inputTmpListIndex: list, inputdataMS1: CFileMS1):

        functionComposition = CFunctionComposition(self.dp)
        functionEMass = CEmass(self.dp)
        num_ID = len(inputTmpListIndex)
        List_Iostope = self.__soldierfillListIsotope(self.dp.myCFG.C10_DDA_ISOTOPE)
        num_Isotope = len(List_Iostope)
        MATRIX_MS1_Signal = np.zeros(shape=(num_ID * num_Isotope, 1))
        # 计算当前结果列表中全部的moz，并和索引建立标记
        List_MOZ_Order = []
        for tmp_index in inputTmpListIndex:
            tmp_moz = self.dp.myIDForPhenomenon.ID7_RPECURSOR_MOZ_CAL[tmp_index]
            tmp_charge = self.dp.myIDForPhenomenon.ID6_CHARGE[tmp_index]
            # tmp_mass = tmp_moz * tmp_charge - tmp_charge * self.dp.myINI.MASS_PROTON_MONO
            moz_isotope = [[tmp_moz + i / tmp_charge, tmp_index] for i in List_Iostope]
            List_MOZ_Order.extend(moz_isotope)
        List_MOZ = [i[0] for i in List_MOZ_Order]
        List_argsort_Moz = list(np.argsort(List_MOZ))
        # 建立索引词典
        Index_Dic = {}
        for i in range(len(List_MOZ_Order)):
            if List_MOZ_Order[i][1] not in Index_Dic:
                Index_Dic[List_MOZ_Order[i][1]] = [List_argsort_Moz.index(i)]
            else:
                Index_Dic[List_MOZ_Order[i][1]].append(List_argsort_Moz.index(i))

        List_sorted_moz = np.array(List_MOZ)[List_argsort_Moz]
        # 遍历所有ms1谱图，对每张谱图返回母离子moz列表对应的最近索引
        for i, tmp_scan in enumerate(inputdataMS1.INDEX_SCAN):
            List_ms1_moz = np.array(inputdataMS1.MATRIX_PEAK_MOZ[tmp_scan])
            List_ms1_int = np.array(inputdataMS1.MATRIX_PEAK_INT[tmp_scan])
            List_Index = tooFindNeighborListFromSortedList(List_ms1_moz, List_sorted_moz)
            List_Int_Index = List_ms1_int[List_Index]  # 不限制质量偏差下返回的母离子强度列表
            if i == 8:
                print(List_Index)

            if self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['PPM']:
                List_Mass_Dev = np.abs(List_ms1_moz[List_Index] - List_sorted_moz) / List_sorted_moz * 1e6
                List_Int_Index[List_Mass_Dev > self.dp.myCFG.C2_DIA_ACCURACY_HALF_WIN_PEAK] = 0.
            elif self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['Da']:
                List_Mass_Dev = np.abs(List_ms1_moz[List_Index] - List_sorted_moz)
                List_Int_Index[List_Mass_Dev > self.dp.myCFG.C2_DIA_ACCURACY_HALF_WIN_PEAK] = 0.
            if i == 0:
                MATRIX_MS1_Signal = List_Int_Index.reshape(-1, 1)
            else:
                MATRIX_MS1_Signal = np.concatenate((MATRIX_MS1_Signal, List_Int_Index.reshape(-1, 1)), axis=1)

        # 按照原来鉴定结果索引的顺序
        List_Recover_Index = []
        for tmp_index in inputTmpListIndex:

            tmp_index_isotope = Index_Dic[tmp_index]
            List_Recover_Index.extend(tmp_index_isotope)

        MATRIX_MS1_Signal = MATRIX_MS1_Signal[List_Recover_Index]  # ‘修复’顺序后的ms1谱图信号矩阵
        # 对每一个鉴定结果，遍历每一张ms1谱图，计算在该ms1谱图上该鉴定结果是否有存在的证据
        List_Sim = []
        for i in range(len(inputTmpListIndex)):
            # 对每一个结果先调用emass算法计算多个同位素的强度分布情况
            i_sequence_pep = self.dp.myIDForPhenomenon.ID1_SEQ[inputTmpListIndex[i]]
            i_mod_dict = self.dp.myIDForPhenomenon.ID5_MOD[inputTmpListIndex[i]]  # 这个格式和emass需要的格式不太一样
            i_charge = self.dp.myIDForPhenomenon.ID6_CHARGE[inputTmpListIndex[i]]
            # 对齐emass需要的修饰格式，标准形式为6,Carbamidomethyl[C];14,Carbamidomethyl[C];（我的是字典形式）
            i_mod = ''

            for i_key in i_mod_dict.keys():
                i_value = i_mod_dict[i_key]
                i_mod += str(i_value) + ',' + str(i_value) + ';'
            # 调用EMASS算法计算得到鉴定结果的理论同位素峰强度
            i_composition = functionComposition.getStrComposition(i_sequence_pep, i_mod, '', '', self.dp.myINI)
            i_dic_composition = functionComposition.getDictComposition(i_composition)
            i_isotope_peaks_all = functionEMass.getCalculatedIsotopicPeaks(i_dic_composition, i_charge)
            # 得到在config中记录的理论同位素峰强度
            i_isotope_peaks_theory = [i_isotope_peaks_all[0][i] for i in List_Iostope]
            # 对当前一个鉴定结果，计算RT范围上所有ms1谱图对应的实际同位素谱峰强度和理论强度的相似度
            start_ID, end_ID = i * num_Isotope, (i+1) * num_Isotope
            i_MATRIX_MS1_Signal = MATRIX_MS1_Signal[start_ID: end_ID]
            i_List_isotope_peaks_theory = np.array(i_isotope_peaks_theory).reshape(-1, 1)
            print(i_List_isotope_peaks_theory)
            i_List_Sim = toolCosineForMatrixAndList(i_MATRIX_MS1_Signal, i_List_isotope_peaks_theory)
            List_Sim.append(list(i_List_Sim.reshape(-1,)))

        return MATRIX_MS1_Signal, List_Sim

    def __captainWrite(self, inputMatrixSignal, inputListSim: list, inputListIndex: list, inputDataMS1: CFileMS1):

        file_folder = self.dp.myCFG.D1_PATH_EXPORT + '\\' if self.dp.myCFG.D1_PATH_EXPORT[-1] != '\\' else self.dp.myCFG.D1_PATH_EXPORT
        writePath = file_folder + IO_NAME_FILE_EXPORT_Flow6[0]

        table_info = ['No', 'File.Name', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'Isotope', 'Isotope_MOZ']
        table_empty = ['empty_seperator']
        table_intensity = ['Intensity_Scan' + str(i) for i in range(len(inputDataMS1.INDEX_RT))]
        table_scan = ['Scan'] + [''] * 6
        table_rt = ['RT'] + [''] * 6
        info_scan = [str(i) for i in inputDataMS1.INDEX_SCAN]
        info_rt = [str(i) for i in inputDataMS1.INDEX_RT]
        List_Iostope = self.__soldierfillListIsotope(self.dp.myCFG.C10_DDA_ISOTOPE)
        num_Isotope = len(List_Iostope)

        with open(writePath, 'w')as f:
            # 写主标题
            table_main = table_info + table_empty + table_intensity
            f.write('\t'.join(table_main) + '\n')
            # 写副标题
            table_sub1 = table_scan + [''] + info_scan
            f.write('\t'.join(table_sub1) + '\n')
            table_sub2 = table_rt + [''] + info_rt
            f.write('\t'.join(table_sub2) + '\n')
            # 写每一个鉴定结果的信息及母离子同位素峰信号强度
            for i_, i_index in enumerate(inputListIndex):
                file_name = self.dp.myIDForPhenomenon.ID0_RAW_NAME[i_index]
                sequence = self.dp.myIDForPhenomenon.ID1_SEQ[i_index]
                sequence_with_mod = self.dp.myIDForPhenomenon.ID4_SEQ_WITH_MOD[i_index]
                charge = self.dp.myIDForPhenomenon.ID6_CHARGE[i_index]
                precursor_moz = self.dp.myIDForPhenomenon.ID7_RPECURSOR_MOZ_CAL[i_index]
                precursor_moz_isotope_list = [str(precursor_moz + j / charge) for j in List_Iostope]
                info_ID = [file_name, sequence, sequence_with_mod, str(charge)]
                start_position, end_position = i_ * num_Isotope, (i_ + 1) * num_Isotope
                for i, i_isotope in enumerate(List_Iostope):
                    line_ID = [str(i_)] + info_ID + [str(i), precursor_moz_isotope_list[i]]
                    list_intensity = inputMatrixSignal[start_position + i]
                    list_intensity = [str(intensity) for intensity in list_intensity]
                    line_one_isotope = line_ID + [''] + list_intensity
                    f.write('\t'.join(line_one_isotope) + '\n')
                line_score = [str(sim) for sim in inputListSim[i_]]
                line_score = [str(i_)] + [''] * 7 + line_score
                f.write('\t'.join(line_score) + '\n')

    def __captainDraw(self, inputMatrixSignal, inputListSim: list, inputListIndex: list, inputDataMS1: CFileMS1):

        List_scan = np.array(inputDataMS1.INDEX_SCAN)  # scan列表用于画图x轴标签
        List_Iostope = self.__soldierfillListIsotope(self.dp.myCFG.C10_DDA_ISOTOPE)
        num_Isotope = len(List_Iostope)
        # 对每一个鉴定结果画母离子信号的3维图像
        fig = plt.figure(figsize=(15, 12))
        ax = Axes3D(fig)
        # 设置x,y,z轴的背景颜色为空
        # ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        # ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        # ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        # 取消网格线
        ax.grid(False)
        ax.set_xlabel('RT(min)', fontsize=13)
        ax.set_ylabel('M / Z', fontsize=13)
        ax.set_zlabel('Intensity', fontsize=13)
        # 调整输出3D图像的角度（elev是沿着y轴偏斜，azim是沿着z轴旋转角度
        ax.view_init(elev=30, azim=235)
        # ax.set_xticks(scan_list)
        a = plt.gca()
        # 强度显示为指数
        a.zaxis.get_major_formatter().set_powerlimits((0, 1))
        for i_, i_index in enumerate(inputListIndex):
            file_name = self.dp.myIDForPhenomenon.ID0_RAW_NAME[i_index]
            sequence = self.dp.myIDForPhenomenon.ID1_SEQ[i_index]
            sequence_with_mod = self.dp.myIDForPhenomenon.ID4_SEQ_WITH_MOD[i_index]
            charge = self.dp.myIDForPhenomenon.ID6_CHARGE[i_index]
            precursor_moz = self.dp.myIDForPhenomenon.ID7_RPECURSOR_MOZ_CAL[i_index]
            precursor_moz_isotope_list = [precursor_moz + j / charge for j in List_Iostope]
            info_ID = [file_name, sequence, sequence_with_mod, str(charge)]
            start_position, end_position = i_ * num_Isotope, (i_ + 1) * num_Isotope
            matrix_isotope = inputMatrixSignal[start_position: end_position]
            List_index_green = np.array(inputListSim[i_]) >= 0.95
            List_index_black = np.array(inputListSim[i_]) < 0.95
            List_scan_green = List_scan[List_index_green]
            matrix_isotope_green = matrix_isotope[:, List_index_green]
            List_scan_black = List_scan[List_index_black]
            matrix_isotope_black = matrix_isotope[:, List_index_black]

            for i in range(matrix_isotope_green.shape[0]):

                xline = np.array(List_scan_green)
                yline = np.array(matrix_isotope_green[i])
                zline = precursor_moz_isotope_list[i]
                ax.bar(xline, yline, zs=zline, zdir='y', color='#00CD00', linewidth=250, width=250)

            for i in range(matrix_isotope_black.shape[0]):

                xline = np.array(List_scan_black)
                yline = np.array(matrix_isotope_black[i])
                zline = precursor_moz_isotope_list[i]
                ax.bar(xline, yline, zs=zline, zdir='y', color='#000000', linewidth=250, width=250)
            ax.set_zlim(0, 1.6 * np.max(matrix_isotope))
            ax.set_ylim(min(precursor_moz_isotope_list) - 0.6 * (precursor_moz_isotope_list[0] - precursor_moz_isotope_list[1]),
                        max(precursor_moz_isotope_list) + 0.6 * (precursor_moz_isotope_list[0] - precursor_moz_isotope_list[1]))
            # ax.set_zticks([(i / 5) * np.max(matrix_isotope) for i in range(6)])
            ax.set_yticks(list(reversed(precursor_moz_isotope_list)))
            ax.set_yticklabels(['{:.2f}'.format(i) for i in reversed(precursor_moz_isotope_list)])
            plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + 'precursor' + file_name.split('\\')[-1] + sequence_with_mod + '.png', format='png')
            plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + 'precursor' + file_name.split('\\')[-1] + sequence_with_mod + '.svg', format='svg')

    def getSignal(self, inputListIndex):

        nMS1 = len(self.dp.LIST_PATH_MS1)

        for iMS1 in range(nMS1):

            # init ms1
            pathMS1 = self.dp.LIST_PATH_MS1[iMS1]
            nameRawMS1 = toolGetNameFromPath(pathMS1)

            # load
            dataMS1 = CFileMS1()
            functionParseMS1 = CFunctionParseMS1(self.dp)
            dataMS1 = functionParseMS1.loadPKL(pathMS1)

            tmpListIndex = []
            for i in range(len(inputListIndex)):
                iID = inputListIndex[i]
                nameRawID = toolGetNameFromPath(self.dp.myIDForPhenomenon.ID0_RAW_NAME[iID])
                if nameRawID == nameRawMS1:
                    tmpListIndex.append(iID)

            MatrixSignal, ListSim = self.__captainCalSignal(tmpListIndex, dataMS1)

            self.__captainWrite(MatrixSignal, ListSim, tmpListIndex, dataMS1)

            self.__captainDraw(MatrixSignal, ListSim, tmpListIndex, dataMS1)

            del dataMS1
            del MatrixSignal
            del ListSim


class CFunctionExtractAllScanForMS2:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __captainCalWinMS2(self, inputDataMS2: CFileMS2):

        MS2_MOZ_List = inputDataMS2.LIST_ACTIVATION_CENTER
        Scan_Index = inputDataMS2.INDEX_SCAN
        Dic_WinMS2 = {}
        print(Scan_Index)
        for i, i_Index in enumerate(Scan_Index):
            i_moz = MS2_MOZ_List[i_Index]
            if i_moz not in Dic_WinMS2:
                Dic_WinMS2[i_moz] = [i_Index]
            else:
                Dic_WinMS2[i_moz].append(i_Index)

        return Dic_WinMS2

    def __captainCalMozForID(self, inputListIndex: list, inputListMozMS2: list):

        len_MS2 = len(inputListMozMS2)
        ListIndexMoz = []
        for i_Index in inputListIndex:
            i_scan = int(self.dp.myIDForPhenomenon.ID2_SCAN_NO[i_Index])
            order_scan = i_scan % (len_MS2 + 1) - 2
            ListIndexMoz.append(inputListMozMS2[order_scan])

        return ListIndexMoz

    def __captainWrite(self, inputMatrixSignal, inputListSim: list, inputIndex,inputMozFragment: list,
                       inputTypeFragment: list, inputMozList: list, inputDataMS2: CFileMS2):
        print(len(inputMozList))
        file_name = self.dp.myIDForPhenomenon.ID0_RAW_NAME[inputIndex]
        sequence = self.dp.myIDForPhenomenon.ID1_SEQ[inputIndex]
        sequence_with_mod = self.dp.myIDForPhenomenon.ID4_SEQ_WITH_MOD[inputIndex]
        charge = self.dp.myIDForPhenomenon.ID6_CHARGE[inputIndex]

        file_folder = self.dp.myCFG.D1_PATH_EXPORT + '\\' if self.dp.myCFG.D1_PATH_EXPORT[-1] != '\\' else self.dp.myCFG.D1_PATH_EXPORT
        writePath = file_folder + sequence_with_mod + IO_NAME_FILE_EXPORT_Flow6[1]

        table_info = ['File.Name', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'Fragment', 'Fragment_MOZ']
        table_empty = ['empty_seperator']
        table_intensity = ['Intensity_Scan' + str(i) for i in range(len(inputDataMS2.INDEX_RT))]
        table_scan = ['Scan'] + [''] * 5
        table_rt = ['RT'] + [''] * 5
        info_scan = [str(inputDataMS2.LIST_PRECURSOR_SCAN[i]) for i in inputMozList]
        info_rt = [str(inputDataMS2.LIST_RET_TIME[i]) for i in inputMozList]

        with open(writePath, 'w')as f:
            # 写主标题
            table_main = table_info + table_empty + table_intensity
            f.write('\t'.join(table_main) + '\n')
            # 写副标题
            table_sub1 = table_scan + [''] + info_scan
            f.write('\t'.join(table_sub1) + '\n')
            table_sub2 = table_rt + [''] + info_rt
            f.write('\t'.join(table_sub2) + '\n')
            # 鉴定结果的信息及碎片离子信号强度
            fragment_moz_list = inputMozFragment
            fragment_type_list = inputTypeFragment
            info_ID = [file_name, sequence, sequence_with_mod, str(charge)]
            start_position, end_position = 0, len(fragment_type_list)
            for i, i_moz in enumerate(fragment_moz_list):
                line_ID = info_ID + [fragment_type_list[i], str(i_moz)]
                list_intensity = inputMatrixSignal[start_position + i]
                list_intensity = [str(intensity) for intensity in list_intensity]
                line_one_isotope = line_ID + [''] + list_intensity
                f.write('\t'.join(line_one_isotope) + '\n')
            line_score = [str(sim) for sim in inputListSim]
            line_score = [''] * 7 + line_score
            f.write('\t'.join(line_score) + '\n')

    def __captainDraw(self, inputMatrixSignal, inputListSim: list, inputIndex,inputMozFragment: list,
                       inputTypeFragment: list, inputMozList: list, inputDataMS2: CFileMS2):

        List_scan = np.array([inputDataMS2.LIST_PRECURSOR_SCAN[i] for i in inputMozList])  # scan列表用于画图x轴标签
        # 对每一个鉴定结果画碎片离子信号的3维图像
        fig = plt.figure(figsize=(13, 12))
        ax = Axes3D(fig)
        # 设置x,y,z轴的背景颜色为空
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        # 取消网格线
        ax.grid(False)
        ax.set_xlabel('Scan', fontsize=13)
        ax.set_ylabel('M / Z', fontsize=13)
        ax.set_zlabel('Intensity', fontsize=13)
        # 调整输出3D图像的角度（elev是沿着y轴偏斜，azim是沿着z轴旋转角度
        ax.view_init(elev=17, azim=225)
        # ax.set_xticks(scan_list)
        a = plt.gca()
        # 强度显示为指数
        a.zaxis.get_major_formatter().set_powerlimits((0, 1))
        i_index = inputIndex
        file_name = self.dp.myIDForPhenomenon.ID0_RAW_NAME[i_index]
        sequence = self.dp.myIDForPhenomenon.ID1_SEQ[i_index]
        sequence_with_mod = self.dp.myIDForPhenomenon.ID4_SEQ_WITH_MOD[i_index]
        charge = self.dp.myIDForPhenomenon.ID6_CHARGE[i_index]
        inputMozFragment = [float(i) for i in inputMozFragment]
        fragment_moz_list = np.array(inputMozFragment)
        print(inputMozFragment)
        argsort_fragment = np.argsort(fragment_moz_list)
        print(argsort_fragment)
        fragment_moz_list = fragment_moz_list[argsort_fragment]
        info_ID = [file_name, sequence, sequence_with_mod, str(charge)]
        matrix_isotope = inputMatrixSignal[argsort_fragment]
        matrix_isotope = matrix_isotope.reshape(matrix_isotope.shape[0], -1)
        List_index_green = np.array(inputListSim) >= 2.
        List_index_black = np.array(inputListSim) < 2.
        List_scan_green = List_scan[List_index_green]
        matrix_fragment_green = matrix_isotope[:, List_index_green]
        List_scan_black = List_scan[List_index_black]
        matrix_fragment_black = matrix_isotope[:, List_index_black]
        print(fragment_moz_list)
        print(List_scan_green)
        print(List_scan_black)
        for i in range(matrix_fragment_black.shape[0]):

            xline = np.array(List_scan_black)
            yline = np.array(matrix_fragment_black[i])
            zline = fragment_moz_list[i]
            # ax.bar(xline, yline, zs=zline, zdir='y', color='#000000', linewidth=250, width=250)
            ax.scatter(xline, yline, zs=zline, zdir='y', color='#000000', linewidth=0.2, s=1.2, alpha=0.4)

        for i in range(matrix_fragment_green.shape[0]):

            xline = np.array(List_scan_green)
            yline = np.array(matrix_fragment_green[i])
            zline = fragment_moz_list[i]
            ax.bar(xline, yline, zs=zline, zdir='y', color='#00CD00', linewidth=1000, width=1000)

        ax.set_zlim(0., 1.2 * np.max(matrix_isotope))
        max_moz = max(fragment_moz_list)
        min_moz = min(fragment_moz_list)
        ax.set_ylim(min_moz, max_moz)
        # ax.set_zticks([(i / 5) * np.max(matrix_isotope) for i in range(6)])
        ax.set_yticks(list(reversed(fragment_moz_list)))
        ax.set_yticklabels(['{:.2f}'.format(i) for i in reversed(fragment_moz_list)])
        plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + 'fragment' + file_name.split('\\')[-1] + sequence_with_mod + '.png', format='png')
        plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + 'fragment' + file_name.split('\\')[-1] + sequence_with_mod + '.svg', format='svg')
        plt.close()

    def __captainCalSignal(self, inputIndex, inputListMozMS2, inputdataMS2: CFileMS2):
        # 1.计算当前鉴定结果的by离子moz
        Sequence = self.dp.myIDForPhenomenon.ID1_SEQ[inputIndex]
        Charge = self.dp.myIDForPhenomenon.ID6_CHARGE[inputIndex]
        Mod_dict = self.dp.myIDForPhenomenon.ID5_MOD[inputIndex]
        Precursor_moz = self.dp.myIDForPhenomenon.ID7_RPECURSOR_MOZ_CAL[inputIndex]
        MOZ_FRAGMENT = []
        TYPE_FRAGMENT = []
        # 计算b/y离子质荷比用于PSM打分
        tmp_b_mass = 0.0
        b_mass_list = np.zeros(shape=(len(Sequence) - 1, 1))
        Precursor_mass = Precursor_moz * Charge - Charge * self.dp.myINI.MASS_PROTON_MONO
        for i, aa in enumerate(Sequence[:-1]):
            tmp_b_mass += self.dp.myINI.DICT1_AA_COM[aa][1]
            if i in Mod_dict:
                tmp_b_mass += float(self.dp.myINI.DICT2_MOD_COM[Mod_dict[i]][1])
            b_mass_list[i] = tmp_b_mass
        y_mass_list = Precursor_mass - b_mass_list
        y_mass_list = np.sort(y_mass_list.reshape([1, -1])).reshape([-1, 1])
        for tmp_charge in range(1, min(Charge + 1, 3)):
            # b fragment
            tmp_b_moz_list = (b_mass_list + tmp_charge * self.dp.myINI.MASS_PROTON_MONO) / tmp_charge
            MOZ_FRAGMENT.extend(list(tmp_b_moz_list))
            TYPE_FRAGMENT.extend([''.join(['b', str(i + 1), '+' * tmp_charge])
                                         for i in range(len(Sequence) - 1)])
            # y fragment
            tmp_y_moz_list = (y_mass_list + tmp_charge * self.dp.myINI.MASS_PROTON_MONO) / tmp_charge
            MOZ_FRAGMENT.extend(list(tmp_y_moz_list))
            TYPE_FRAGMENT.extend([''.join(['y', str(i + 1), '+' * tmp_charge])
                                         for i in range(len(Sequence) - 1)])
        MOZ_FRAGMENT = np.array(MOZ_FRAGMENT)
        # 2.遍历所有ms2谱图，计算by离子信号矩阵和打分
        MATRIX_MS2_Signal = np.zeros(shape=(len(MOZ_FRAGMENT), 1))
        for i, tmp_scan in enumerate(inputListMozMS2):
            List_ms2_moz = np.array(inputdataMS2.MATRIX_PEAK_MOZ[tmp_scan])
            List_ms2_int = np.array(inputdataMS2.MATRIX_PEAK_INT[tmp_scan])
            List_Index = tooFindNeighborListFromSortedList(List_ms2_moz, MOZ_FRAGMENT)
            List_Int_Index = List_ms2_int[List_Index]  # 不限制质量偏差下返回的母离子强度列表
            # 对二分最近查找返回的结果进行质量精度的控制
            if self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['PPM']:
                List_Mass_Dev = np.abs(List_ms2_moz[List_Index] - MOZ_FRAGMENT) / MOZ_FRAGMENT * 1e6
                List_Int_Index[List_Mass_Dev > 20.] = 0.
            elif self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['Da']:
                List_Mass_Dev = np.abs(List_ms2_moz[List_Index] - MOZ_FRAGMENT)
                List_Int_Index[List_Mass_Dev > 0.02] = 0.
            if i == 0:
                MATRIX_MS2_Signal = List_Int_Index.reshape(-1, 1)
            else:
                MATRIX_MS2_Signal = np.concatenate((MATRIX_MS2_Signal, List_Int_Index.reshape(-1, 1)), axis=1)
        # 对该鉴定结果，遍历每一张ms2谱图，计算在该ms2谱图上该鉴定结果是否有存在的证据（计算打分）
        print(MATRIX_MS2_Signal.shape)
        List_Sim = []
        type_by_list = np.array(['b' if i[0] == 'b' else 'y'for i in TYPE_FRAGMENT])
        type_b_list = (type_by_list == 'b')
        type_y_list = (type_by_list == 'y')
        # 计算分数1.匹配到by离子数目（归一化）
        score0_match_num_List = list(np.sum(MATRIX_MS2_Signal > 100., axis=0) / len(Sequence))
        # 计算分数2，匹配到最长连续b离子数目
        MATRIX_MS2_Signal_B = MATRIX_MS2_Signal[type_b_list] > 100.
        MATRIX_MS2_Signal_B_Simple = []
        for i_charge in range(1, min(Charge + 1, 3)):
            if i_charge == 1:
                MATRIX_MS2_Signal_B_Simple = MATRIX_MS2_Signal_B[(i_charge-1)*(len(Sequence)-1): (i_charge)*(len(Sequence)-1)]
            else:
                MATRIX_MS2_Signal_B_Simple = MATRIX_MS2_Signal_B_Simple + MATRIX_MS2_Signal_B[
                                             (i_charge - 1) * (len(Sequence) - 1): (i_charge) * (len(Sequence) - 1)]
        print(MATRIX_MS2_Signal_B_Simple.shape)
        score1_max_continue_b_list = [toolMaxContinueBY(list(MATRIX_MS2_Signal_B_Simple[:, i])) / len(Sequence)
                                      for i in range(MATRIX_MS2_Signal_B_Simple.shape[1])]
        # 计算分数3，匹配到最长连续y离子数目
        MATRIX_MS2_Signal_Y = MATRIX_MS2_Signal[type_y_list] > 100.
        MATRIX_MS2_Signal_Y_Simple = []
        for i_charge in range(1, min(Charge + 1, 3)):
            if i_charge == 1:
                MATRIX_MS2_Signal_Y_Simple = MATRIX_MS2_Signal_Y[
                                             (i_charge - 1) * (len(Sequence) - 1): (i_charge) * (len(Sequence) - 1)]
            else:
                MATRIX_MS2_Signal_Y_Simple = MATRIX_MS2_Signal_Y_Simple + MATRIX_MS2_Signal_Y[
                                                                   (i_charge - 1) * (len(Sequence) - 1): (i_charge) * (
                                                                               len(Sequence) - 1)]
        print(MATRIX_MS2_Signal_Y_Simple.shape)
        score2_max_continue_y_list = [toolMaxContinueBY(list(MATRIX_MS2_Signal_Y_Simple[:, i])) / len(Sequence)
                                      for i in range(MATRIX_MS2_Signal_Y_Simple.shape[1])]

        List_Sim = list(np.array(score0_match_num_List) + np.array(score1_max_continue_b_list) +
                        np.array(score2_max_continue_y_list))
        TYPE_FRAGMENT = list(TYPE_FRAGMENT)
        MOZ_FRAGMENT = list(MOZ_FRAGMENT)
        return MATRIX_MS2_Signal, List_Sim, MOZ_FRAGMENT, TYPE_FRAGMENT

    def getSignal(self, inputListIndex):

        nMS2 = len(self.dp.LIST_PATH_MS2)

        for iMS2 in range(nMS2):

            # init ms2
            pathMS2 = self.dp.LIST_PATH_MS2[iMS2]
            nameRawMS2 = toolGetNameFromPath(pathMS2)

            # load
            dataMS2 = CFileMS2()
            functionParseMS2 = CFunctionParseMS2(self.dp)
            dataMS2 = functionParseMS2.loadPKL(pathMS2)

            tmpListIndex = []
            for i in range(len(inputListIndex)):
                iID = inputListIndex[i]
                nameRawID = toolGetNameFromPath(self.dp.myIDForPhenomenon.ID0_RAW_NAME[iID])
                if nameRawID == nameRawMS2:
                    tmpListIndex.append(iID)
            # 将MS2文件的所有隔离窗口下的ms2谱图存在一个列表
            Dic_WinMS2 = self.__captainCalWinMS2(dataMS2)
            List_WinMS2 = list(Dic_WinMS2.keys())
            List_WinMS2.sort()
            print(Dic_WinMS2)
            print(List_WinMS2)
            # 对鉴定结果列表中的每个结果计算它的隔离窗口moz
            tmpListIndexMoz = self.__captainCalMozForID(tmpListIndex, List_WinMS2)
            print(tmpListIndexMoz)

            # 对每一个鉴定结果，计算by离子信号并画图
            for i, i_Index in enumerate(tmpListIndex):
                i_IndexMoz = tmpListIndexMoz[i]
                i_ListMozMS2 = Dic_WinMS2[i_IndexMoz]
                MatrixSignal, ListSim, Moz_Fragment, Type_Fragment = self.__captainCalSignal(i_Index,
                                                                                             i_ListMozMS2, dataMS2)

                self.__captainWrite(MatrixSignal, ListSim, i_Index, Moz_Fragment, Type_Fragment,
                                    i_ListMozMS2, dataMS2)

                self.__captainDraw(MatrixSignal, ListSim, i_Index, Moz_Fragment, Type_Fragment,
                                    i_ListMozMS2, dataMS2)

                del MatrixSignal
                del ListSim



class CFunctionDrawAllScanForMS1:  # 暂时不用这个了（直接在内存中画图）

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def draw(self, inputListIndex, inputMatrixSignal, inputListSim, inputDataMS1: CFileMS1):

        pass



class CFunctionDrawDIACheck:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soliderFillSeed(self, inputIndex, inputdataMS2: CFileMS2, LossDic: dict):

        # 初始化seed
        myseed = CSeedDIACheck()
        op_INIT_CSEED_DIACHECK(myseed)

        if self.dp.myIDForCheck.ID3_RT[inputIndex] != VALUE_ILLEGAL:
            indexMS2 = toolFindNeighborFromSortedList1(inputdataMS2.INDEX_RT,
                                                            self.dp.myIDForCheck.ID3_RT[inputIndex] * 60.)
        elif self.dp.myIDForCheck.ID2_SCAN_NO[inputIndex] != VALUE_ILLEGAL:
            indexMS2 = toolFindNeighborFromSortedList1(inputdataMS2.INDEX_SCAN,
                                                              self.dp.myIDForCheck.ID2_SCAN_NO[inputIndex])
        else:
            info = "MSFunctionDIACheck CFunctionDrawDIACheck in Flow5, MK56, RT:{:.2f} " \
                   "or Scan{:d} is right ?".format(self.dp.myIDForCheck.ID3_RT[inputIndex],
                                                   self.dp.myIDForCheck.ID2_SCAN_NO[inputIndex])
            logGetError(info)

        myseed.MID_RT = inputdataMS2.INDEX_RT[indexMS2]
        myseed.MID_SCAN = inputdataMS2.INDEX_SCAN[indexMS2]

        Sequence = self.dp.myIDForCheck.ID1_SEQ[inputIndex]
        Mod_dict = self.dp.myIDForCheck.ID4_MOD[inputIndex]
        Charge = self.dp.myIDForCheck.ID6_CHARGE[inputIndex]
        Precursor_moz = self.dp.myIDForCheck.ID5_PRECURSOR_MOZ_CLC[inputIndex]
        # 计算母离子同位素峰质荷比
        for i in [0, 1, 2]:
            isotope_moz = Precursor_moz + i / Charge
            myseed.MOZ_PRECURSOR.append(isotope_moz)
        # 计算b/y离子质荷比用于PSM打分
        tmp_b_mass = 0.0
        b_mass_list = np.zeros(shape=(len(Sequence) - 1, 1))
        Precursor_mass = Precursor_moz * Charge - Charge * self.dp.myINI.MASS_PROTON_MONO
        for i, aa in enumerate(Sequence[:-1]):
            tmp_b_mass += self.dp.myINI.DICT1_AA_COM[aa][1]
            if i in Mod_dict:
                tmp_b_mass += float(self.dp.myINI.DICT2_MOD_COM[Mod_dict[i]][1])
            b_mass_list[i] = tmp_b_mass
        y_mass_list = Precursor_mass - b_mass_list
        y_mass_list = np.sort(y_mass_list.reshape([1, -1])).reshape([-1, 1])
        for tmp_charge in range(1, min(Charge + 1, 3)):
            # b fragment
            tmp_b_moz_list = (b_mass_list + tmp_charge * self.dp.myINI.MASS_PROTON_MONO) / tmp_charge
            myseed.MOZ_FRAGMENT.extend(list(tmp_b_moz_list))
            myseed.TYPE_FRAGMENT.extend([''.join(['b', str(i + 1), '+' * tmp_charge])
                                         for i in range(len(Sequence) - 1)])
            # y fragment
            tmp_y_moz_list = (y_mass_list + tmp_charge * self.dp.myINI.MASS_PROTON_MONO) / tmp_charge
            myseed.MOZ_FRAGMENT.extend(list(tmp_y_moz_list))
            myseed.TYPE_FRAGMENT.extend([''.join(['y', str(i + 1), '+' * tmp_charge])
                                         for i in range(len(Sequence) - 1)])
            # b fragment loss
            for Loss in LossDic.keys():
                loss_mass = LossDic[Loss]
                tmp_b_moz_list = (b_mass_list - loss_mass + tmp_charge * self.dp.myINI.MASS_PROTON_MONO) / tmp_charge
                myseed.MOZ_FRAGMENT.extend(list(tmp_b_moz_list))
                myseed.TYPE_FRAGMENT.extend([''.join(['b', str(i + 1), '+' * tmp_charge, '_' + Loss])
                                             for i in range(len(Sequence) - 1)])
                tmp_y_moz_list = (y_mass_list - loss_mass + tmp_charge * self.dp.myINI.MASS_PROTON_MONO) / tmp_charge
                myseed.MOZ_FRAGMENT.extend(list(tmp_y_moz_list))
                myseed.TYPE_FRAGMENT.extend([''.join(['y', str(i + 1), '+' * tmp_charge, '_' + Loss])
                                             for i in range(len(Sequence) - 1)])
        # Library中的b/y离子moz
        fragmentTypeList = self.dp.myIDForCheck.ID13_FRAGMENT_TYPE[inputIndex]
        fragmentLossList = self.dp.myIDForCheck.ID14_FRAGMENT_LOSS_TYPE[inputIndex]
        fragmentTypeForCalList = list(map(op_DIVIDE_ION_TYPE, fragmentTypeList))
        fragmentLossForCalLost = ['' if i == 'NOLOSS' else i for i in fragmentLossList]
        fragmentLossList = ['' if i == 'NOLOSS' else '_' + i for i in fragmentLossList]
        fragmentTypeList = [i[0] + i[1] for i in zip(fragmentTypeList, fragmentLossList)]
        fragmentMozList = []
        for i in range(len(fragmentLossForCalLost)):
            ion_by, ion_site, ion_charge = fragmentTypeForCalList[i]
            fragmentLoss = fragmentLossForCalLost[i]
            fragmentLossMass = 0 if fragmentLoss == '' else LossDic[fragmentLoss]
            if ion_by == 'b':
                fragmentMoz = (b_mass_list[ion_site - 1] + ion_charge * self.dp.myINI.MASS_PROTON_MONO
                               - fragmentLossMass) / ion_charge
            else:
                fragmentMoz = (y_mass_list[ion_site - 1] + ion_charge * self.dp.myINI.MASS_PROTON_MONO
                               - fragmentLossMass) / ion_charge
            fragmentMozList.append(fragmentMoz)

        myseed.MOZ_LIBRARY, myseed.TYPE_LIBRARY = fragmentMozList, fragmentTypeList

        return myseed

    def __captainGetEvidence(self, inputListEvidence, inputListSeed, inputListIndex):
        '''
        :param inputListEvidence: 输入空的evidence列表用于填充，得到的evidence用于画图
        :param inputListSeed: 输入空的seed列表用于填充，得到的seed用于对应evidence的获取
        :param inputListIndex: 输入鉴定结果列表的索引列表，去鉴定结果数据结构中寻址得到rt，母离子等信息
        :return:
        '''
        nMS1 = len(self.dp.LIST_PATH_MS1)
        nMS2 = len(self.dp.LIST_PATH_MS2)

        # 计算中性损失的质量
        LossDic = {}
        for loss in FRAGMENT_LOSS_TYPE_COMP.keys():
            if loss == 'NOLOSS':
                pass
            else:
                LossDic[loss] = op_CAL_MOLECULAR_MASS_FROM_FRAGMENT_LOSS(self.dp.myINI.DICT0_ELEMENT_MASS,
                                                                         loss)
        for iMS1 in range(nMS1):
            # init ms1
            pathMS1 = self.dp.LIST_PATH_MS1[iMS1]
            nameRawMS1 = toolGetNameFromPath(pathMS1)

            for iMS2 in range(nMS2):
                # init ms2
                pathMS2 = self.dp.LIST_PATH_MS2[iMS2]
                nameRawMS2 = toolGetNameFromPath(pathMS2)
                # load
                if nameRawMS2 == nameRawMS1:

                    # load ms1
                    dataMS1 = CFileMS1()
                    functionParseMS1 = CFunctionParseMS1(self.dp)
                    dataMS1 = functionParseMS1.loadPKL(pathMS1)

                    dataMS2 = CFileMS2()
                    functionParseMS2 = CFunctionParseMS2(self.dp)
                    dataMS2 = functionParseMS2.loadPKL(pathMS2)

                    for i in range(len(inputListIndex)):

                        iID = inputListIndex[i]
                        nameRawID = toolGetNameFromPath(self.dp.myIDForCheck.ID0_RAW_NAME[iID])
                        if nameRawID == nameRawMS1:
                            tmpSeed = self.__soliderFillSeed(iID, dataMS2, LossDic)
                            inputListSeed[i] = tmpSeed
                            functionEvidence = CFunctionEvidenceForDIACheck(self.dp)
                            tmpEvidence = functionEvidence.fillEvidence(dataMS1, dataMS2, tmpSeed, iID)
                            inputListEvidence[i] = tmpEvidence
                    del dataMS2
            del dataMS1

            # if dataMS2 == '':
            #     info = "MSFunctionDrawPrecursor in Flow1, ms1 file and ms2 file not match. MK449, " + self.dp.myCFG.A3_PATH_MS2 + 'and' + \
            #            self.dp.myCFG.A1_PATH_MS1 + " is all right?"
            #     logGetWarning(info)



    def __GetPlotPepText(self, clc_pep, mod_dic: dict):
        # 顶部肽段标签
        pep_text = ''
        pep_mod_text = ''
        for pep_index in range(len(clc_pep)):
            if pep_index < (len(clc_pep) - 1):
                if pep_index + 1 in mod_dic.keys():
                    pep_text = pep_text + ' ' + ' | '
                    pep_mod_text = pep_mod_text + clc_pep[pep_index] + '   '
                else:
                    pep_text = pep_text + clc_pep[pep_index] + ' | '
                    pep_mod_text = pep_mod_text + ' ' + '   '
            else:
                pep_text = pep_text + clc_pep[pep_index]
                pep_mod_text = pep_mod_text + ' '
        return pep_text, pep_mod_text

    def __soliderDrawDIACurve(self, inputSeed: CSeedDIACheck, inputEvidence: CEvidenceDIACheck, inputIndex: int, gs, fig):

        # 根据seed和evidence计算每一个鉴定结果Lib中的色谱曲线
        moz_frag = inputSeed.MOZ_LIBRARY
        type_frag = inputSeed.TYPE_LIBRARY
        moz_prec = inputSeed.MOZ_PRECURSOR
        n_fragment = len(type_frag)
        n_iostope = len(moz_prec)
        n_rt = len(inputEvidence.LIST_RET_TIME)
        accuracy = self.dp.myCFG.C2_DIA_ACCURACY_HALF_WIN_PEAK

        MATRIX_PROFILE_FRAG = np.zeros(shape=(n_fragment, n_rt))
        MATRIX_PROFILE_PREC = np.zeros(shape=(n_iostope, n_rt))

        for i in range(n_rt):

            tmp_MOZ_list = inputEvidence.MATRIX_MS2_PEAK_MOZ[i]
            tmp_INT_list = inputEvidence.MATRIX_MS2_PEAK_INT[i]
            tmp_MOZ_list_ms1 = inputEvidence.MATRIX_MS1_PEAK_MOZ[i]
            tmp_INT_list_ms1 = inputEvidence.MATRIX_MS1_PEAK_INT[i]

            if self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['PPM']:
                # print(tmp_MOZ_list)
                # print(precursor_MOZ)
                for j, tmp_isotope_moz in enumerate(moz_prec):
                    indexINT_ms1 = toolFindNeighborFromSortedList1(tmp_MOZ_list_ms1, tmp_isotope_moz)
                    mass_dev = abs(tmp_MOZ_list_ms1[indexINT_ms1] - tmp_isotope_moz) / tmp_isotope_moz * 1e6
                    if mass_dev < accuracy:
                        MATRIX_PROFILE_PREC[j, i] = tmp_INT_list_ms1[indexINT_ms1]
                    else:
                        MATRIX_PROFILE_PREC[j, i] = 0.0

                for j, fragment_moz in enumerate(moz_frag):

                    indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, fragment_moz)

                    mass_dev = abs(tmp_MOZ_list[indexINT] - fragment_moz) / fragment_moz * 1e6
                    if mass_dev < accuracy:
                        MATRIX_PROFILE_FRAG[j, i] = tmp_INT_list[indexINT]
                    else:
                        MATRIX_PROFILE_FRAG[j, i] = 0.0

            elif self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['Da']:

                for j, tmp_isotope_moz in enumerate(moz_prec):
                    indexINT_ms1 = toolFindNeighborFromSortedList1(tmp_MOZ_list_ms1, tmp_isotope_moz)
                    mass_dev = abs(tmp_MOZ_list_ms1[indexINT_ms1] - tmp_isotope_moz)
                    if mass_dev < accuracy:
                        MATRIX_PROFILE_PREC[j, i] = tmp_INT_list_ms1[indexINT_ms1]
                    else:
                        MATRIX_PROFILE_PREC[j, i] = 0.0

                for j, fragment_moz in enumerate(moz_frag):

                    indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, fragment_moz)

                    mass_dev = abs(tmp_MOZ_list[indexINT] - fragment_moz)
                    if mass_dev < accuracy:
                        MATRIX_PROFILE_FRAG[j, i] = tmp_INT_list[indexINT]
                    else:
                        MATRIX_PROFILE_FRAG[j, i] = 0.0

            else:
                logGetError('MSFunctionEvidence,TYPE_ACCURACY:' + str(
                    self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK) + 'is all right?')
        #  根据peak rt的强度对b/y离子和isotope 母离子进行排序

        frag_rank_list = np.argsort(-1 * MATRIX_PROFILE_FRAG[:, inputEvidence.MID_IDNEX])
        MATRIX_PROFILE_FRAG = MATRIX_PROFILE_FRAG[frag_rank_list]
        FRAG_MOZ = [moz_frag[i] for i in frag_rank_list]
        FRAG_TYPE = [type_frag[i] for i in frag_rank_list]
        PREC_MOZ = moz_prec
        PREC_TYPE = ['precursor', 'precursor[M+1]', 'precursor[M+2]']

        text_base_title = 'Base Peak:'
        text_precursor_title = 'Precursor:'
        text_raw_title = 'Raw File:'
        text_qvalue_title = 'q value:'
        text_charge_title = 'Charge:'
        text_base_prec_title = 'Base Peak Prec:'

        if len(MATRIX_PROFILE_FRAG) != 0:
            base_intensity = np.max(MATRIX_PROFILE_FRAG)
            if base_intensity > 0.01:
                matrix_profile_frag = MATRIX_PROFILE_FRAG / base_intensity
            else:
                matrix_profile_frag = MATRIX_PROFILE_FRAG
            ax_frag = plt.subplot(gs[9:16, 0:9])
            frag_type_text = []

            x_list = inputEvidence.LIST_RET_TIME
            x_list = [i / 60.0 for i in x_list]

            for i in range(n_fragment):
                y_list = matrix_profile_frag[i]
                plt.plot(x_list, y_list, color=DRAW_COLOR[i], )
                frag_type_text.append(
                    FRAG_TYPE[i] + ':' + format(float(FRAG_MOZ[i]), '.2f') + ' m/z')
                # plt.scatter(x_list, y_list, color=DRAW_COLOR[i])

            xlim = [np.min(x_list), np.max(x_list)]
            plt.xlim(xlim)
            plt.ylim([0, 1])
            # plot rt apex line, rt start line, rt end line
            if inputEvidence.RT_START != VALUE_ILLEGAL:
                rt_start = inputEvidence.LIST_RET_TIME[inputEvidence.RT_START] / 60.0
                rt_end = inputEvidence.LIST_RET_TIME[inputEvidence.RT_END] / 60.0
                plt.plot([rt_start, rt_start], [0, 1], color='r', linestyle=':')
                plt.plot([rt_end, rt_end], [0, 1], color='r', linestyle=':')
            plt.xlabel('RT (min)')
            plt.ylabel('Intensity (%)')

            ax_frag_type = fig.add_subplot(gs[9:16, 9:])
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            for i, text in enumerate(frag_type_text):
                ax_frag_type.add_patch(
                    pc.Rectangle(  # 长方形
                        (0.01, 0.95 - i * 0.05),  # （x,y）
                        0.06,  # 长
                        0.01,  # 宽
                        color=DRAW_COLOR[i]
                    )
                )
                plt.text(0.1, 0.95 - i * 0.05, text, fontsize=8)
            plt.axis('off')

        else:
            base_intensity = 0

        base_intensity_prec = np.max(MATRIX_PROFILE_PREC)
        matrix_profile_prec = MATRIX_PROFILE_PREC / base_intensity_prec

        ax = plt.subplot(gs[0:2, 0:])
        text_base_value = format(base_intensity, '.2e')
        text_base_prec_value = format(base_intensity_prec, '.2e')
        text_precursor_value = self.dp.myIDForCheck.ID7_PRECURSOR_ID[inputIndex]
        text_raw_value = toolGetNameFromPath(self.dp.myIDForCheck.ID0_RAW_NAME[inputIndex])
        text_charge_value = str(self.dp.myIDForCheck.ID6_CHARGE[inputIndex])
        text_qvalue_value = format(self.dp.myIDForCheck.ID8_SCORE0[inputIndex], '.5f')

        plt.xlim((0, 200))
        plt.ylim((0, 20))
        plt.xticks([])
        plt.yticks([])
        plt.text(2, 15, text_base_title, fontdict={'size': 10, 'color': 'black'})
        plt.text(23, 15, text_base_value, fontdict={'size': 10, 'color': 'red'})
        plt.text(2, 0, text_raw_title, fontdict={'size': 10, 'color': 'black'})
        plt.text(23, 0, text_raw_value, fontdict={'size': 10, 'color': 'red'})
        plt.text(2, 8, text_precursor_title, fontdict={'size': 10, 'color': 'black'})
        plt.text(23, 8, text_precursor_value, fontdict={'size': 10, 'color': 'red'})
        plt.text(52, 0, text_charge_title, fontdict={'size': 10, 'color': 'black'})
        plt.text(65, 0, text_charge_value, fontdict={'size': 10, 'color': 'red'})
        plt.text(102, 0, text_qvalue_title, fontdict={'size': 10, 'color': 'black'})
        plt.text(116, 0, text_qvalue_value, fontdict={'size': 10, 'color': 'red'})
        plt.text(52, 15, text_base_prec_title, fontdict={'size': 10, 'color': 'black'})
        plt.text(80, 15, text_base_prec_value, fontdict={'size': 10, 'color': 'red'})
        plt.axis('off')
        ax_prec = plt.subplot(gs[2:9, 0:9])
        prec_type_text = []
        x_list = inputEvidence.LIST_RET_TIME
        x_list = [i / 60.0 for i in x_list]
        for i in range(n_iostope):
            y_list = matrix_profile_prec[i]
            plt.plot(x_list, y_list, color=DRAW_COLOR[i])
            prec_type_text.append(
                PREC_TYPE[i] + ':' + format(float(PREC_MOZ[i]), '.2f') + ' m/z')
            # plt.scatter(x_list, y_list, color=DRAW_COLOR[i])

        xlim = [np.min(x_list), np.max(x_list)]
        plt.xlim(xlim)
        # plot rt apex line, rt start line, rt end line
        # if inputEvidence.RT_START != VALUE_ILLEGAL:
        #     rt_start = inputEvidence.LIST_RET_TIME[inputEvidence.RT_START] / 60.0
        #     rt_end = inputEvidence.LIST_RET_TIME[inputEvidence.RT_END] / 60.0
        #     plt.plot([rt_start, rt_start], [0, 1], color='r')
        #     plt.plot([rt_end, rt_end], [0, 1], color='r')
        plt.xlabel('RT (min)')
        plt.ylabel('Intensity (%)')

        ax_prec_type = fig.add_subplot(gs[2:9, 9:])

        plt.xlim(0, 1)
        plt.ylim(0, 1)
        for i, text in enumerate(prec_type_text):
            ax_prec_type.add_patch(
                pc.Rectangle(  # 长方形
                    (0.01, 0.95 - i * 0.1),  # （x,y）
                    0.06,  # 长
                    0.01,  # 宽
                    color=DRAW_COLOR[i]
                )
            )
            plt.text(0.1, 0.95 - i * 0.1, text, fontsize=8)
        plt.axis('off')

    def __soliderDrawPSMLabel(self, i_label: int, inputSeed: CSeedDIACheck, inputEvidence: CEvidenceDIACheck, inputIndex: int, gs, inputMaxPeak):

        clc_moz = inputSeed.MOZ_FRAGMENT
        clc_tag = inputSeed.TYPE_FRAGMENT
        exp_moz = np.array(inputEvidence.MATRIX_MS2_PEAK_MOZ[i_label + inputEvidence.RT_START])
        exp_int = np.array(inputEvidence.MATRIX_MS2_PEAK_INT[i_label + inputEvidence.RT_START])
        max_int = max(exp_int)
        # 对ms2谱图所有谱峰，获取config中设置的topN强度的峰
        if self.dp.myCFG.C15_INT_NUM == -1:
            exp_int = list(exp_int)
            exp_moz = list(exp_moz)
        elif len(exp_moz) > self.dp.myCFG.C15_INT_NUM:
            top_num_flag = np.sort(exp_int)[len(exp_moz) - self.dp.myCFG.C15_INT_NUM - 1]
            flag_list = exp_int >= top_num_flag
            exp_int = list(exp_int[flag_list])
            exp_moz = list(exp_moz[flag_list])
        else:
            exp_int = list(exp_int)
            exp_moz = list(exp_moz)

        Y1_NAME = 'Relative Intensity(%)'  # 谱图y轴标注
        Y2_NAME = 'ppm'  # 偏差图y标注
        X_NAME = 'm/z'  # 横轴坐标标注

        # Label图框架
        # i_label当前第几张需要标记的ms2谱图，gs是整个check图的网格布局
        gs_start = 17 + i_label * 17

        # 图1.PSM的基本信息
        ax1 = plt.subplot(gs[gs_start+0, :])
        plt.xlim((0, 110))
        plt.ylim((0, 20))
        plt.xticks([])
        plt.yticks([])

        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        info_scan = inputEvidence.LIST_SCAN[i_label + inputEvidence.RT_START]
        info_rt = '{:.3f}'.format(inputEvidence.LIST_RET_TIME[i_label + inputEvidence.RT_START] / 60.)
        info_raw_name = self.dp.myIDForCheck.ID0_RAW_NAME[inputIndex]
        info_base_peak = '{:.3g}'.format(max_int)
        info_moz = '{:.3f}'.format(self.dp.myIDForCheck.ID5_PRECURSOR_MOZ_CLC[inputIndex])
        mod_dic = self.dp.myIDForCheck.ID4_MOD[inputIndex]

        # plt.text(0, 15, 'Title:', color='black', fontsize=10)
        # plt.text(4, 15, info_raw_name, color='red', fontsize=10)
        plt.text(50, 5, 'MaxPeakMS2:', color='black', fontsize=10)
        plt.text(61, 5, info_base_peak, color='red', fontsize=10)
        plt.text(17, 5, 'Moz:', color='black', fontsize=10)
        plt.text(21, 5, info_moz, color='red', fontsize=10)
        plt.text(30, 5, 'RT', color='black', fontsize=10)
        plt.text(33, 5, info_rt, color='red', fontsize=10)
        plt.text(40, 5, 'Scan', color='black', fontsize=10)
        plt.text(44, 5, info_scan, color='red', fontsize=10)

        # 图2.匹配序列
        ax2 = plt.subplot(gs[gs_start+1:gs_start+2, :])
        plt.xlim((0, 2500))
        plt.ylim((0, 30))
        plt.xticks([])
        plt.yticks([])

        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['left'].set_visible(False)

        info_pep = self.dp.myIDForCheck.ID1_SEQ[inputIndex]
        info_charge = self.dp.myIDForCheck.ID6_CHARGE[inputIndex]

        # 顶部标注
        font1 = {'family': 'SimHei',
                 'color': 'black',
                 'weight': 'normal',
                 'size': 10,
                 }
        font2 = {'family': 'SimHei',
                 'color': 'red',
                 'weight': 'normal',
                 'size': 10,
                 }

        pep_text, pep_mod_text = self.__GetPlotPepText(info_pep, mod_dic)

        plt.text(100, 15, pep_text, fontdict=font1)  # 黑色序列
        plt.text(95, 15, pep_mod_text, fontdict=font2)  # 标红修饰的
        plt.text(10, 25, str(info_charge) + '+', color='red', fontsize=10)  # 在序列左上角标注电荷

        # 图3.ms2谱图谱峰标注
        ax_3_half = plt.subplot(gs[gs_start+2:gs_start+6, 0:])
        ax_3_half.spines['top'].set_visible(False)
        ax_3_half.spines['right'].set_visible(False)
        ax_3_half.bar(exp_moz, exp_int, width=2, facecolor='black')
        plt.yticks([])
        x_axis_max = (math.ceil(inputMaxPeak / 100) + 1) * 100  # 横坐标最大值
        if x_axis_max > 2500:
            x_axis_max = 2500
        plt.xlim((0, x_axis_max))
        plt.xticks([])
        # x_axis_max = (math.ceil(max(exp_moz) / 100) + 1) * 100  # 横坐标最大值

        ax3 = plt.subplot(gs[gs_start+6:gs_start+15, 0:])
        # ax3.bar(exp_moz, exp_int, width=2, facecolor='#BEBEBE', alpha=0.5)

        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)

        # x_axis_max = (math.ceil(max(exp_moz) / 100) + 1) * 100  # 横坐标最大值
        # x_axis_max = (math.ceil(inputMaxPeak / 100) + 1) * 100  # 横坐标最大值
        #
        # if x_axis_max > 2500:
        #     x_axis_max = 2500
        plt.xlim((0, x_axis_max))
        plt.xticks([])
        plt.ylim((0, 120))
        plt.yticks([0, 20, 40, 60, 80, 100])
        plt.ylabel(Y1_NAME)
        # plt.title(TITLE, FontProperties=font)

        # 图4. ms2谱图谱峰匹配的质量偏差值
        ax4 = plt.subplot(gs[gs_start+15, 0:])
        plt.xlim((0, x_axis_max))
        # 设置刻度（主刻度为100，副刻度为10）
        ax = plt.gca()
        x_major_locator = MultipleLocator(100)
        x_minor_locator = MultipleLocator(10)
        ax.xaxis.set_minor_locator(x_minor_locator)
        ax.xaxis.set_major_locator(x_major_locator)
        info_type_bias = self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK
        info_bias = float(self.dp.myCFG.C2_DIA_ACCURACY_HALF_WIN_PEAK)

        if info_type_bias == 0:  # 0: ppm 1: Da
            ax4_y = info_bias
        else:
            ax4_y = 20

        plt.ylim((-ax4_y, ax4_y))
        plt.grid(axis="y", linestyle='--')
        plt.xlabel(X_NAME)
        plt.ylabel(Y2_NAME)

        pep_moz_start = ''
        type_pep = 'single'

        op_Data_fill_plot(clc_moz, clc_tag, exp_moz, exp_int, info_type_bias, info_bias, ax1,
                          ax2, ax_3_half, ax3, ax4, pep_moz_start, type_pep, [])

    def __captainDrawDIACheck(self, inputListSeed, inputListEvidence, inputListIndex):

        # 遍历结果列表中所有鉴定结果，对每一个结果分别画结果的library色谱曲线图和在结果rt_start到rt_end保留时间范围内的PSM标图
        for tmp_i in range(len(inputListIndex)):
            if inputListSeed[tmp_i]:  # 有的鉴定结果可能是无效的（如不存在文件名）
                tmp_index = inputListIndex[tmp_i]
                tmp_seed = inputListSeed[tmp_i]
                tmp_evidence = inputListEvidence[tmp_i]
                precursor_id = self.dp.myIDForCheck.ID7_PRECURSOR_ID[tmp_index]

                num_label = tmp_evidence.RT_END - tmp_evidence.RT_START + 1  # 色谱曲线保留时间范围内的ms2谱图个数
                fig_Check = plt.figure(figsize=(15, 17 + num_label * 5))  # 新建画布（色谱曲线 + ms2谱图PSM标图）
                gs = GridSpec(17 + 17 * num_label, 12)
                self.__soliderDrawDIACurve(tmp_seed, tmp_evidence, tmp_index, gs, fig_Check)  # 画母离子和Lib中by离子的色谱曲线
                # 计算rt范围内ms2谱图最大moz
                max_peak = np.max([np.max(i_ms2) for i_ms2 in tmp_evidence.MATRIX_MS2_PEAK_MOZ[tmp_evidence.RT_START:tmp_evidence.RT_END]])
                max_peak = min(max_peak, np.max(tmp_seed.MOZ_FRAGMENT))
                # 下面这行代码是错误的，得到的是最大的一个ms2谱图moz，因为对于不同长度的列表转化成numpy数组，得到的不是二维数组
                # 而是以列表为元素的一维数组，对齐进行max操作返回的是值最大的那个列表而不是所有ms2谱图moz列表中最大的数值
                # max_peak = np.max(np.array(tmp_evidence.MATRIX_MS2_PEAK_MOZ[tmp_evidence.RT_START:tmp_evidence.RT_END]).reshape([-1,]))
                for i in range(num_label):
                    i_label = i
                    self.__soliderDrawPSMLabel(i_label, tmp_seed, tmp_evidence, tmp_index, gs, max_peak)  # 画RT范围上ms2谱图的标图
                print(precursor_id)
                plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + precursor_id + '_' +
                            IO_NAME_FILE_EXPORT_Flow5[0])
                plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + precursor_id + '_' +
                            IO_NAME_FILE_EXPORT_Flow5[1])
                plt.close(fig_Check)
            else:
                pass

    def draw(self):

        # 获取要画图的结果在ID数据结构中的位置
        listIndex = list(range(self.dp.myIDForCheck.N_ID))

        listEvidence = [[] * 1] * len(listIndex)
        listSeed = [[] * 1] * len(listIndex)

        self.__captainGetEvidence(listEvidence, listSeed, listIndex)
        self.__captainDrawDIACheck(listSeed, listEvidence, listIndex)


class CFunctionDrawFragment:

    def __init__(self, inputDP: CDataPack):
        self.dp = inputDP

    def __captainCalMoz(self, iID):
        ionPrecursor = self.dp.myID.ID4_SEQ_WITH_MOD[iID]
        ionTypeList = self.dp.myID.ID10_FRAGMENT_TYPE[iID]
        ionLossList = self.dp.myID.ID13_FRAGMENT_LOSS_TYPE[iID]

        ionTypeList, ionMozList = op_CAL_FRAGMENT_MOZ(self.dp.myINI, ionPrecursor, ionTypeList, ionLossList)

        return ionTypeList, ionMozList

    def __captainCalPRECMOZ(self, iID):

        ionPrecursor = self.dp.myID.ID4_SEQ_WITH_MOD[iID]
        ionCharge = self.dp.myID.ID7_CHARGE[iID]

        ionMoz = op_CAL_PERCURSOR_MOZ(self.dp.myINI, ionPrecursor, ionCharge)
        #return self.dp.myID.ID5_PRECURSOR_MOZ_EXP[iID]  # 这里要根据序列电荷计算m/z，还没做，先标记一下
        return ionMoz

    def __captainFillSeed(self, iID, inputListMS2Scan, inputListMS2RT):
        # init seed
        myseed = CSeed2()
        # 根据rt或者scan号得到在ms2文件中的索引位置
        if self.dp.myID.ID2_SCAN_ID[iID] != VALUE_ILLEGAL:
            indexMS2 = toolFindNeighborFromSortedList1(inputListMS2Scan, self.dp.myID.ID2_SCAN_ID[iID])
        else:
            indexMS2 = toolFindNeighborFromSortedList1(inputListMS2RT, self.dp.myID.ID3_RT[iID] * 60)

        myseed.MID_SCAN = inputListMS2Scan[indexMS2]
        myseed.MID_RT = inputListMS2RT[indexMS2]
        # 根据鉴定结果的起始时间和结束时间得到在ms2文件的索引位置，并得到在ms2中的起始时间和结束时间
        if self.dp.myID.ID14_RT_BEGIN[iID] != VALUE_ILLEGAL:
            indexMS2 = toolFindNeighborFromSortedList1(inputListMS2Scan, self.dp.myID.ID14_RT_BEGIN[iID] * 60)
            myseed.START_RT = inputListMS2RT[indexMS2]
        else:
            myseed.START_RT = VALUE_ILLEGAL
            # myseed.START_RT = max(myseed.MID_RT - 0.07 * 60, 0)
        if self.dp.myID.ID15_RT_END[iID] != VALUE_ILLEGAL:
            indexMS2 = toolFindNeighborFromSortedList1(inputListMS2Scan, self.dp.myID.ID15_RT_END[iID] * 60)
            myseed.END_RT = inputListMS2RT[indexMS2]
        else:
            myseed.END_RT = VALUE_ILLEGAL
            # myseed.END_RT = min(myseed.MID_RT + 0.07 * 60, inputListMS2RT[len(inputListMS2RT) - 1])
        myseed.MOZ_PREC = self.__captainCalPRECMOZ(iID)
        myseed.ION_TYPE, myseed.MOZ_CLC = self.__captainCalMoz(iID)

        return myseed

    def __captainGetEvidence(self, inputListEvidence2, inputListSeed2, inputListIndex):

        nMS2 = len(self.dp.LIST_PATH_MS2)
        for tmp_i in range(nMS2):
            # init ms2
            pathMS2 = self.dp.LIST_PATH_MS2[tmp_i]
            nameRawMS2 = toolGetNameFromPath(pathMS2)
            dataMS2 = CFileMS2()
            functionParseMS2 = CFunctionParseMS2(self.dp)
            dataMS2 = functionParseMS2.loadPKL(pathMS2)
            # init ms1
            pathMS1 = self.dp.LIST_PATH_MS1[tmp_i]
            nameRawMS1 = toolGetNameFromPath(pathMS1)
            if nameRawMS1 != nameRawMS2:
                info = "MSFunctionDrawFragment in Flow2, MK449, " + self.dp.myCFG.A3_PATH_MS2 + 'and' + \
                       self.dp.myCFG.A1_PATH_MS1 + " is all right?"
                logGetError(info)
            else:
                pass
            dataMS1 = CFileMS1()
            functionParseMS1 = CFunctionParseMS1(self.dp)
            dataMS1 = functionParseMS1.loadPKL(pathMS1)

            for i in range(len(inputListIndex)):
                i_ID = inputListIndex[i]
                nameRawID = toolGetNameFromPath(self.dp.myID.ID1_RAW_NAME[i_ID])
                # 是同一raw文件
                if nameRawMS2 == nameRawID:
                    tmpSeed2 = self.__captainFillSeed(i_ID, dataMS2.INDEX_SCAN, dataMS2.INDEX_RT)
                    inputListSeed2[i] = tmpSeed2
                    # print(tmpSeed2.MID_SCAN,tmpSeed2.MID_RT,tmpSeed2.MOZ_CLC)
                    functionEvidence = CFunctionEvidence2(self.dp)
                    tmpEvidence = functionEvidence.fillEvidence2(dataMS1, dataMS2, tmpSeed2)
                    inputListEvidence2[i] = tmpEvidence
                else:
                    continue

    def __captainDrawFragment(self, inputEvidence_list: list, inputIndex_list: list, score_list: list):

        plt.figure(figsize=(18, 12))
        gs = GridSpec(15, 10)
        text_base_title = 'Base Peak:'
        text_precursor_title = 'Precursor:'
        text_raw_title = 'Raw File:'
        text_qvalue_title = 'q value:'
        text_charge_title = 'Charge:'
        text_score_title = 'Score:'
        for tmp_i in range(len(inputEvidence_list)):
            inputEvidence = inputEvidence_list[tmp_i]
            inputIndex = inputIndex_list[tmp_i]
            rt_apex = self.dp.myID.ID3_RT[inputIndex]
            precursor_id = self.dp.myID.ID8_PRECURSOR_ID[inputIndex]
            n_fragment = len(inputEvidence.LIST_ION_TYPE)

            if len(inputEvidence.MATRIX_PROFILE) == 0:
                break
            base_intensity = np.max(inputEvidence.MATRIX_PROFILE[:-1])
            matrix_profile = inputEvidence.MATRIX_PROFILE / base_intensity

            ax = plt.subplot(gs[0:2, 0:])
            text_base_value = format(base_intensity, '.2e')
            text_precursor_value = self.dp.myID.ID4_SEQ_WITH_MOD[inputIndex]
            text_raw_value = toolGetNameFromPath(self.dp.myID.ID1_RAW_NAME[inputIndex])
            text_charge_value = str(self.dp.myID.ID7_CHARGE[inputIndex])
            text_qvalue_value = format(self.dp.myID.ID9_SCORE0[inputIndex], '.5f')
            text_score_value = score_list[tmp_i] if score_list[tmp_i] != VALUE_ILLEGAL else 'None'
            plt.xlim((0, 200))
            plt.ylim((0, 20))
            plt.xticks([])
            plt.yticks([])
            plt.text(2,  15, text_base_title,      fontdict={'size': 10, 'color': 'black'})
            plt.text(16, 15, text_base_value,      fontdict={'size': 10, 'color': 'red'})
            plt.text(2,  5,  text_raw_title,       fontdict={'size': 10, 'color': 'black'})
            plt.text(14, 5,  text_raw_value,       fontdict={'size': 10, 'color': 'red'})
            plt.text(52, 15, text_precursor_title, fontdict={'size': 10, 'color': 'black'})
            plt.text(65, 15, text_precursor_value, fontdict={'size': 10, 'color': 'red'})
            plt.text(52, 5,  text_charge_title,    fontdict={'size': 10, 'color': 'black'})
            plt.text(62, 5,  text_charge_value,    fontdict={'size': 10, 'color': 'red'})
            plt.text(102, 5, text_qvalue_title,    fontdict={'size': 10, 'color': 'black'})
            plt.text(113, 5, text_qvalue_value,    fontdict={'size': 10, 'color': 'red'})
            plt.text(102, 15, text_score_title,    fontdict={'size': 10, 'color': 'black'})
            plt.text(113, 15, text_score_value,    fontdict={'size': 10, 'color': 'red'})
            # print(inputEvidence.LIST_I_START)

            ax2 = plt.subplot(gs[2:, 0:])
            for i in range(n_fragment):
                x_list = inputEvidence.LIST_RET_TIME
                x_list = [i/60.0 for i in x_list]
                y_list = matrix_profile[i]
                plt.plot(x_list, y_list, label=inputEvidence.LIST_ION_TYPE[i] + ':' +
                         format(inputEvidence.MOZ_CLC[i], '.2f') + ' m/z')
                plt.scatter(x_list, y_list)
            plt.legend(loc='upper right', frameon=False, title='Fragment ion type')

            # plot rt apex line, rt start line, rt end line
            if inputEvidence.RT_START != VALUE_ILLEGAL:
                rt_start = inputEvidence.LIST_RET_TIME[inputEvidence.RT_START] / 60.0
                rt_end = inputEvidence.LIST_RET_TIME[inputEvidence.RT_END] / 60.0
                # plt.plot([rt_apex, rt_apex], [0, 1], color='b')
                plt.plot([rt_start, rt_start], [0, 1], color='r')
                plt.plot([rt_end, rt_end], [0, 1], color='r')
            plt.xlabel('RT (min)')
            plt.ylabel('Intensity (%)')
            plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + precursor_id + '_' +
                        IO_NAME_FILE_EXPORT_Flow2[0])
            # plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + text_qvalue_value + '_' + precursor_id + '_' +
            #             IO_NAME_FILE_EXPORT_Flow2[0], type='png')
            # plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + text_qvalue_value + '_' + precursor_id + '_' +
            #             IO_NAME_FILE_EXPORT_Flow2[0], type='svg')
            # plt.close()
            plt.clf()
        plt.close()

    def __captainCalScore(self, inputEvidence_list: list, inputIndex_list: list, score_list: list):

        output_list = []

        for tmp_i in range(len(inputEvidence_list)):

            inputEvidence = inputEvidence_list[tmp_i]

            if inputEvidence.RT_START == VALUE_ILLEGAL:
                score_list[tmp_i] = VALUE_ILLEGAL
            else:
                # 得到母离子在鉴定结果中索引
                inputIndex = inputIndex_list[tmp_i]
                out_raw_name = self.dp.myID.ID1_RAW_NAME[inputIndex]
                out_precursor_id = self.dp.myID.ID4_SEQ_WITH_MOD[inputIndex] + str(self.dp.myID.ID7_CHARGE[inputIndex])
                precursor_seq = self.dp.myID.ID0_SEQ[inputIndex]
                # 根据evidence中带有起始时间和结束时间的色谱曲线计算母离子各个分数
                # 目前实现分数有（色谱曲线相似度，匹配by离子数目，质量偏差）
                rt_start_position = inputEvidence.RT_START
                rt_end_position = inputEvidence.RT_END
                base_intensity = np.max(inputEvidence.MATRIX_PROFILE)
                matrix_profile = inputEvidence.MATRIX_PROFILE / base_intensity
                fragment_curve = inputEvidence.MATRIX_PROFILE[:-1, rt_start_position:rt_end_position+1]
                precursor_curve = inputEvidence.MATRIX_PROFILE[-1, rt_start_position:rt_end_position+1]
                sum_intensity = np.sum(fragment_curve, axis=1)
                index_sort = np.argsort(sum_intensity)[::-1]
                fragment_curve_top = fragment_curve[index_sort[:6]] if index_sort.shape[0] >= 6 else fragment_curve
                match_num_score = fragment_curve.shape[0]
                # 得到best fragment和pearson相似度
                max_num = -1
                max_no = 0
                for i in range(fragment_curve_top.shape[0]):
                    sum = 0
                    for j in range(fragment_curve_top.shape[0]):
                        sum += toolPearsonForNumpy(fragment_curve_top[i], fragment_curve_top[j])
                    if sum > max_num:
                        max_no = i
                        max_num = sum
                max_num = max_num
                rest_corr = 0
                if index_sort.shape[0] > 6:
                    fragment_curve_rest = fragment_curve[index_sort[6:]]
                    for j in range(fragment_curve_rest.shape[0]):
                        rest_corr += toolPearsonForNumpy(fragment_curve[max_no], fragment_curve_rest[j])
                rest_corr = rest_corr / fragment_curve_rest.shape[0] * match_num_score / len(precursor_seq)
                # 计算best fragment和母离子色谱曲线的相似度
                prec_corr = toolPearsonForNumpy(fragment_curve[max_no], precursor_curve)

                score_similarity = max_num + 0.5 * prec_corr + rest_corr
                score_similarity_average = score_similarity / index_sort.shape[0]
                # 计算by离子色谱曲线上的质量偏差
                matrix_mass_dev = inputEvidence.MATRIX_MASS_DEV[:, int((3 * rt_start_position + rt_end_position) / 4)
                                                                :int((3 * rt_end_position + rt_start_position) / 4) + 1]
                matrix_mass_dev_top = matrix_mass_dev[index_sort[:6]] if index_sort.shape[0] >= 6 else matrix_mass_dev
                mass_dev_score_top = np.cos(np.mean(matrix_mass_dev_top) / 20.0 * np.pi)
                mass_dev_score_top = mass_dev_score_top if index_sort.shape[0] >= 6 else mass_dev_score_top * index_sort.shape[0] / 6
                if index_sort.shape[0] > 6:
                    matrix_mass_dev_rest = matrix_mass_dev[index_sort[6:]]
                    mass_dev_score_rest = np.cos(np.mean(matrix_mass_dev_rest) / 20.0 * np.pi)
                else:
                    mass_dev_score_rest = 0.0
                matrix_mass_dev_prec = matrix_mass_dev[-1]
                mass_dev_score_prec = np.cos(np.mean(matrix_mass_dev_prec) / 20.0 * np.pi)
                score_mass_dev = 2 * mass_dev_score_top + 0.5 * mass_dev_score_rest + 0.5 * mass_dev_score_prec
                
                score_weight = score_similarity + score_mass_dev
                score_list[tmp_i] = score_weight

                out_ion = ';'.join([str(i) for i in inputEvidence.LIST_ION_TYPE])
                n_fragment = len(inputEvidence.LIST_ION_TYPE)
                score_1 = str(score_similarity)
                score_2 = str(score_similarity_average)
                score_3 = str(mass_dev_score_top)
                score_4 = str(score_weight)
                score_5 = str(rest_corr)
                out_XIC = ';'.join([','.join([str(j) for j in list(matrix_profile[i][:])]) for i in range(n_fragment)])
                output_list.append([out_raw_name, out_precursor_id, out_ion, out_XIC, score_1, score_2, score_3, score_4, score_5])

        with open(self.dp.myCFG.D1_PATH_EXPORT + IO_NAME_FILE_EXPORT_XIC_MATRIX, 'w')as f:

            f.write('\t'.join(['raw name', 'precursor_id', 'fragment ion', 'fragment XIC', 'similarity score',
                               'average similarity score', 'mass dev score', 'weight score']) + '\n')
            for i in output_list:

                f.write('\t'.join(i) + '\n')

    def draw(self):

        # 得到要画图的结果列表
        if self.dp.myCFG.C7_PLOT_NUM == CFG_TYPE_PLOT_NUM['ALL']:

            listPrecursor_id = self.dp.myID.ID8_PRECURSOR_ID

        elif self.dp.myCFG.C7_PLOT_NUM == CFG_TYPE_PLOT_NUM['CUSTOM']:
            listPrecursor_id = op_FILL_LIST_DRAW_PRECURSOR(self.dp.myCFG.C4_DIA_PRECURSOR_SEQUENCE,
                                                           self.dp.myCFG.C5_DIA_PRECURSOR_CHARGE,
                                                           self.dp.myCFG.C6_DIA_PRECURSOR_RAW)
        else:
            info = "MSFunctionDraw, MK449, " + self.dp.myCFG.C7_PLOT_NUM + " is all right?"
            logGetError(info)
        # 得到结果列表在ID中的索引列表
        listIndex = []
        for precursor_id in listPrecursor_id:
            if precursor_id in self.dp.myID.ID8_PRECURSOR_ID:

                index = self.dp.myID.ID8_PRECURSOR_ID.index(precursor_id)
                listIndex.append(index)
            else:
                print(precursor_id + '  not in ID')
                # logGetError('input precursor was wrong , please check:' + precursor_id)
        # 得到每个结果的evidence
        listEvidence2 = [[] * 1] * len(listIndex)
        listSeed2 = [[] * 1] * len(listIndex)
        self.__captainGetEvidence(listEvidence2, listSeed2, listIndex)
        # 对每个鉴定结果计算一个自定义的分数
        score_custom = [[] * 1] * len(listIndex)
        # self.__captainCalScore(listEvidence2, listIndex, score_custom)
        # 对结果绘制by离子色谱曲线
        self.__captainDrawFragment(listEvidence2, listIndex, score_custom)


class CFunctionDrawDDAPrecursor:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soldierfillListIsotope(self, inputIsotope: str):

        if inputIsotope[-1] == '|':
            pass
        else:
            inputIsotope = inputIsotope + '|'

        isotopeList = inputIsotope.split('|')[:-1]
        isotopeList = [int(i) for i in isotopeList]

        return isotopeList

    def __soldier3Dplot(self, inputEvidence: CEvidence, inputSeed: CSeed, precursor_id):

        matrix_profile = inputEvidence.MATRIX_PROFILE
        rt_list = [i / 60. for i in inputEvidence.LIST_RET_TIME]
        moz_list = list(reversed(inputSeed.MOZ_CLC))

        fig = plt.figure(figsize=(10, 8))
        ax = Axes3D(fig)
        # 设置x,y,z轴的背景颜色为空
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

        for i in range(len(moz_list)):
            xline = np.array(rt_list)
            yline = np.array(matrix_profile[i])
            zline = moz_list[i]
            ax.plot(xline, yline, zs=zline, zdir='y', color=PLOT_COLOR[i], fillstyle='full', linewidth='2.5', alpha=0.8)

        # 取消网格线
        ax.grid(False)
        ax.set_xlabel('RT(min)', fontsize=13)
        ax.set_ylabel('M / Z', fontsize=13)
        ax.set_zlabel('Intensity', fontsize=13)
        ax.set_zlim(0, 2 * np.max(matrix_profile))
        ax.set_ylim(min(moz_list) - 0.6 * (moz_list[0] - moz_list[1]), max(moz_list) + 0.6 * (moz_list[0] - moz_list[1]))
        ax.set_yticks(moz_list)
        ax.set_yticklabels(['{:.2f}'.format(i) for i in reversed(moz_list)])
        # 调整输出3D图像的角度（elev是沿着y轴偏斜，azim是沿着z轴旋转角度
        ax.view_init(elev=30, azim=235)
        # ax.set_xticks(scan_list)
        a = plt.gca()
        # 强度显示为指数
        a.zaxis.get_major_formatter().set_powerlimits((0, 1))
        plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + '3D_' + precursor_id + '_' +
                    IO_NAME_FILE_EXPORT_Flow1[2])
        plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + '3D_' + precursor_id + '_' +
                    IO_NAME_FILE_EXPORT_Flow1[3])
        # plt.show()
        plt.close()

    def __captainCalMOZ(self, iID):

        ionPrecursor = self.dp.myID.ID4_SEQ_WITH_MOD[iID]
        ionCharge = self.dp.myID.ID7_CHARGE[iID]

        ionMoz = op_CAL_PERCURSOR_MOZ(self.dp.myINI, ionPrecursor, ionCharge)
        print('result', ionPrecursor, ionMoz)
        #return self.dp.myID.ID5_PRECURSOR_MOZ_EXP[iID]  # 这里要根据序列电荷计算m/z，还没做，先标记一下
        return ionMoz

    def __soliderCAL_PERCURSOR_MOZ(self, input_INI: CINI, input_precursor: str, input_charge):
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

        return mass_precursor, mod_dic.keys()

    def __soliderCAL_PERCURSOR_MOZ2(self, input_INI: CINI, precursor_sequence: str, input_charge, mod_dic):
        # precursor_sequence, mod_dic = op_DIVIDE_MOD_FROM_PRECURSOR(input_precursor)
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
        return mass_precursor, mod_dic.keys()

    def __captainFillSeed(self, iID, dataMS1: CFileMS1, inputInfo):

        # 初始化seed
        myseed = CSeedForDDA()

        listMS1Scan = dataMS1.INDEX_SCAN
        listMS1RT = dataMS1.INDEX_RT

        #计算母离子m/z
        if inputInfo == VALUE_ILLEGAL:
            # 根据rt得到ms1中对应的索引
            if self.dp.myDDAID.PSM3_RT[iID] != VALUE_ILLEGAL:
                indexMS1 = toolFindNeighborFromSortedList1(listMS1RT, self.dp.myDDAID.PSM3_RT[iID])
                if listMS1RT[indexMS1] > self.dp.myDDAID.PSM3_RT[iID]:
                    indexMS1 -= 1
            else:
                indexMS1 = toolFindNeighborFromSortedList1(listMS1Scan, self.dp.myDDAID.PSM2_SCAN_ID[iID])
                if listMS1Scan[indexMS1] > self.dp.myDDAID.PSM2_SCAN_ID[iID]:
                    indexMS1 -= 1
            myseed.MID_RT = listMS1RT[indexMS1]
            myseed.MID_SCAN = listMS1Scan[indexMS1]

            if self.dp.myDDAID.PSM8_MH_CLC[iID] == VALUE_ILLEGAL:  # 对pFind结果的母离子计算moz
                myseed.PRECURSOR_SQE = self.dp.myDDAID.PSM4_SEQ[iID]
                myseed.PRECURSOR_CHARGE = self.dp.myDDAID.PSM9_CHARGE[iID]
                myseed.PRECURSOR_MOD = self.dp.myDDAID.PSM5_MOD[iID]
                mod = {}
                for i in myseed.PRECURSOR_MOD.split(';')[:-1]:
                    num, mod_str = i.split(',')
                    mod[float(num)] = mod_str

                precursor_moz, mod = self.__soliderCAL_PERCURSOR_MOZ2(self.dp.myINI, myseed.PRECURSOR_SQE, myseed.PRECURSOR_CHARGE, mod)
                # precursor_moz = (self.dp.myDDAID.PSM8_MH_CLC[iID] + self.dp.myDDAID.PSM9_CHARGE[iID] - 1) / self.dp.myDDAID.PSM9_CHARGE[iID]
            else:  # 对于pQuant结果没有读取母离子序列等信息，不计算母离子moz，直接使用结果文件给定moz
                myseed.PRECURSOR_SQE = 'None' if self.dp.myDDAID.PSM4_SEQ[iID] == VALUE_ILLEGAL \
                    else self.dp.myDDAID.PSM4_SEQ[iID]
                myseed.PRECURSOR_CHARGE = self.dp.myDDAID.PSM9_CHARGE[iID]
                myseed.PRECURSOR_MOD = '' if self.dp.myDDAID.PSM5_MOD[iID] == VALUE_ILLEGAL \
                    else self.dp.myDDAID.PSM5_MOD[iID]
                precursor_moz = (self.dp.myDDAID.PSM8_MH_CLC[iID] - 1.00727645224 + myseed.PRECURSOR_CHARGE * 1.00727645224) / self.dp.myDDAID.PSM9_CHARGE[iID]
                
        else:
            # 根据rt得到ms1中对应的索引
            indexMS1 = toolFindNeighborFromSortedList1(listMS1RT, inputInfo[0])
            if listMS1Scan[indexMS1] > self.dp.myDDAID.PSM2_SCAN_ID[iID]:
                indexMS1 -= 1
            myseed.MID_RT = listMS1RT[indexMS1]
            myseed.MID_SCAN = listMS1Scan[indexMS1]
            myseed.PRECURSOR_SQE = inputInfo[1]
            myseed.PRECURSOR_CHARGE = inputInfo[2]
            precursor_moz, mod = self.__soliderCAL_PERCURSOR_MOZ(self.dp.myINI, inputInfo[1], inputInfo[2])
            myseed.PRECURSOR_MOD = ','.join([str(i) for i in list(mod)])
        cfg_isotope = self.dp.myCFG.C10_DDA_ISOTOPE
        isotope_list = self.__soldierfillListIsotope(cfg_isotope)
        precursor_moz_list = [precursor_moz + i * 1.003/self.dp.myDDAID.PSM9_CHARGE[iID] for i in isotope_list]

        myseed.MOZ_CLC = precursor_moz_list
        # myseed.MOZ_CLC = self.__captainCalMOZ(iID)

        return myseed

    def __captainGetEvidence(self, inputListEvidence, inputListSeed, inputListIndex, inputListInfo):

        nMS1 = len(self.dp.LIST_PATH_MS1)

        for iMS1 in range(nMS1):

            # init ms1
            pathMS1 = self.dp.LIST_PATH_MS1[iMS1]
            nameRawMS1 = toolGetNameFromPath(pathMS1)

            # load
            dataMS1 = CFileMS1()
            functionParseMS1 = CFunctionParseMS1(self.dp)
            dataMS1 = functionParseMS1.loadPKL(pathMS1)

            for i in range(len(inputListIndex)):

                iID = inputListIndex[i]
                nameRawID = self.dp.myDDAID.PSM1_RAW_NAME[iID]
                if len(inputListInfo) != 0:
                    tmp_info = inputListInfo[i]
                else:
                    tmp_info = VALUE_ILLEGAL

                if nameRawID == nameRawMS1:

                    tmpSeed = self.__captainFillSeed(iID, dataMS1, tmp_info)
                    tmpSeed.RAW_NAME = nameRawID
                    # if nameRawID != nameRawMS1:
                    #     # print('来自于{}的鉴定结果去{} raw中画色谱曲线'.format(nameRawID, nameRawMS1))
                    #     tmpSeed.RAW_NAME = nameRawMS1
                    # else:
                    #     tmpSeed.RAW_NAME = nameRawID
                    inputListSeed[i] = tmpSeed

                    functionEvidence = CFunctionEvidenceDDA(self.dp)

                    tmpEvidence = functionEvidence.fillEvidence(dataMS1, tmpSeed)
                    inputListEvidence[i] = tmpEvidence

    def __captainDraw(self, inputListEvidence, inputListSeed, inputListIndex):

        ax = plt.figure(figsize=(18, 5))
        gs = GridSpec(15, 14)
        text_base_title = 'Base Peak:'
        text_sequence_title = 'Sequence:'
        text_raw_title = 'Raw File:'
        text_qvalue_title = 'q value:'
        text_charge_title = 'Charge:'
        text_modify_title = 'Modification:'
        text_midScan_titile = 'Mid Scan:'

        for tmp_i in range(len(inputListEvidence)):

            inputSeed = inputListSeed[tmp_i]
            inputEvidence = inputListEvidence[tmp_i]
            inputIndex = inputListIndex[tmp_i]
            # precursor_id = self.dp.myDDAID.PSM1_RAW_NAME[inputIndex] + self.dp.myDDAID.PSM4_SEQ[inputIndex] + \
            #                str(self.dp.myDDAID.PSM9_CHARGE[inputIndex])
            precursor_id = self.dp.myDDAID.PSM6_PRECURSOR_ID[inputIndex]
            rt_start = self.dp.myDDAID.PSM13_RT_START[inputIndex]
            rt_end = self.dp.myDDAID.PSM14_RT_END[inputIndex]
            n_fragment = inputEvidence.MATRIX_PROFILE.shape[0]

            # 输出3D色谱曲线的一个截面图
            # self.__soldier3Dplot(inputEvidence, inputSeed, precursor_id)

            if self.dp.myCFG.C14_CURVE_SMOOTH:
                max_intensity = 0.
                matrix_profile = np.zeros(shape=(inputEvidence.MATRIX_PROFILE.shape[0], 300))
                x_list = [i / 60 for i in inputEvidence.LIST_RET_TIME]
                for i in range(inputEvidence.MATRIX_PROFILE.shape[0]):
                    y_list = inputEvidence.MATRIX_PROFILE[i]
                    # 300 represents number of points to make between T.min and T.max
                    xnew = np.linspace(min(x_list), max(x_list), 300)
                    spl = CubicSpline(x_list, y_list)
                    # spl = make_interp_spline(x_list, y_list, k=3)  # type: BSpline
                    power_smooth = spl(xnew)
                    matrix_profile[i] = power_smooth
                    tmp_max = np.max(power_smooth)
                    if tmp_max > max_intensity:
                        max_intensity = tmp_max
                base_intensity = max_intensity
                matrix_profile = matrix_profile / base_intensity
                matrix_profile[matrix_profile < 0] = 0
            else:
                # base_intensity = np.max(inputEvidence.MATRIX_PROFILE)
                # # print(base_intensity)
                # matrix_profile = inputEvidence.MATRIX_PROFILE / (base_intensity)
                # 根据母离子色谱曲线匹配结果的最大值画图
                x_list = [i / 60 for i in inputEvidence.LIST_RET_TIME]
                start_position = toolFindNeighborFromSortedList1(x_list, rt_start / 60)
                end_position = toolFindNeighborFromSortedList1(x_list, rt_end / 60)
                # base_intensity = np.max(inputEvidence.MATRIX_PROFILE[:, start_position: end_position])
                base_intensity = np.max(inputEvidence.MATRIX_PROFILE)
                #matrix_profile = inputEvidence.MATRIX_PROFILE
                matrix_profile = inputEvidence.MATRIX_PROFILE / max(base_intensity, 1)

            ax = plt.subplot(gs[3:6, 0:])
            text_raw_value = inputSeed.RAW_NAME
            length_raw_name = len(text_raw_value)
            text_midScan_value = inputSeed.MID_SCAN

            plt.xlim((0, 200))
            plt.ylim((0, 10))
            plt.xticks([])
            plt.yticks([])
            plt.text(2, 5, text_raw_title, fontdict={'size': 10, 'color': 'black'})
            plt.text(14, 5, text_raw_value, fontdict={'size': 10, 'color': 'red'})
            plt.text(42 - 5 + length_raw_name, 5, text_midScan_titile, fontdict={'size': 10, 'color': 'black'})
            plt.text(55 - 5 + length_raw_name, 5, text_midScan_value, fontdict={'size': 10, 'color': 'red'})

            plt.axis('off')

            # print(inputEvidence.LIST_I_START)

            ax2 = plt.subplot(gs[6:, 0:])
            for i in range(n_fragment):
                plt.plot([], [], label=format(inputSeed.MOZ_CLC[i], '.4f')+'m/z', color=PLOT_COLOR[i],)
                # plt.scatter(x_list, y_list, s=0.5, color=DRAW_COLOR[i])
            plt.legend(loc=3, frameon=False, bbox_to_anchor=(0.85, 1.1), borderaxespad=0)
            max_intensity = 0.
            for i in range(n_fragment - 1, -1, -1):
                x_list = [i / 60 for i in inputEvidence.LIST_RET_TIME]
                if self.dp.myCFG.C14_CURVE_SMOOTH:
                    x_list = np.linspace(min(x_list), max(x_list), 300)
                y_list = matrix_profile[i]
                # 类似于skyline的fill进行半透明填充，效果不太好
                # plt.fill(x_list, y_list, color=PLOT_COLOR[i], alpha=0.75, linewidth=1.5)
                plt.plot(x_list, y_list, color=PLOT_COLOR[i], alpha=0.75, linewidth=1.5)

            # 如果读取的结果文件有记录起始时间和终止时间（如pQuant的定量结果文件），在色谱曲线中标注出来
            if rt_start != VALUE_ILLEGAL and rt_end != VALUE_ILLEGAL:
                rt_start = rt_start / 60.0
                rt_end = rt_end / 60.0
                plt.plot([rt_start, rt_start], [0.0, 1.0], color='#FF0000', linestyle='--')
                plt.plot([rt_end, rt_end], [0.0, 1.0], color='#FF0000', linestyle='--')

            plt.yticks([0, 0.5, 1])
            plt.ylim(0., 1.)
            # 此处保证pQuant画图时一个结果对应多个raw文件母离子色谱曲线的x轴（rt）严格对齐
            x_limit = [self.dp.myDDAID.PSM3_RT[inputIndex] / 60. - self.dp.myCFG.C1_DIA_RT_HALF_WIN_IN_MIN,
                       self.dp.myDDAID.PSM3_RT[inputIndex] / 60. + self.dp.myCFG.C1_DIA_RT_HALF_WIN_IN_MIN]
            # plt.xlim(x_limit)
            plt.xlabel('RT(minute)')
            plt.ylabel('Intensity (%)')
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)

            ax = plt.subplot(gs[0:3, 0:])
            text_base_value = format(base_intensity, '.2e')
            text_sequence_value = inputSeed.PRECURSOR_SQE
            text_charge_value = str(int(inputSeed.PRECURSOR_CHARGE))
            text_modify_value = inputSeed.PRECURSOR_MOD
            if self.dp.myDDAID.PSM10_SCORE0[inputIndex] == VALUE_ILLEGAL:
                text_qvalue_value = 'None'
            else:
                text_qvalue_value = format(self.dp.myDDAID.PSM10_SCORE0[inputIndex], '.3e')
            if len(text_modify_value) == 0:
                text_modify_value = 'None'
            plt.xlim((0, 200))
            plt.ylim((0, 10))
            plt.axis('off')
            plt.xticks([])
            plt.yticks([])
            plt.text(2, 5, text_base_title, fontdict={'size': 10, 'color': 'black'})
            plt.text(16, 5, text_base_value, fontdict={'size': 10, 'color': 'red'})
            plt.text(2, 0, text_sequence_title, fontdict={'size': 10, 'color': 'black'})
            plt.text(14, 0, text_sequence_value, fontdict={'size': 10, 'color': 'red'})
            plt.text(38, 5, text_charge_title, fontdict={'size': 10, 'color': 'black'})
            plt.text(49, 5, text_charge_value, fontdict={'size': 10, 'color': 'red'})
            plt.text(60, 5, text_qvalue_title, fontdict={'size': 10, 'color': 'black'})
            plt.text(72, 5, text_qvalue_value, fontdict={'size': 10, 'color': 'red'})
            plt.text(42 - 5 + length_raw_name, 0, text_modify_title, fontdict={'size': 10, 'color': 'black'})
            plt.text(57 - 5 + length_raw_name, 0, text_modify_value, fontdict={'size': 10, 'color': 'red'})
            print(self.dp.myCFG.D1_PATH_EXPORT + precursor_id + '_' +
                        IO_NAME_FILE_EXPORT_Flow1[0])
            plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + precursor_id + '_' +
                        IO_NAME_FILE_EXPORT_Flow1[0])
            plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + precursor_id + '_' +
                        IO_NAME_FILE_EXPORT_Flow1[1])
            # # 输出母离子色谱曲线到txt用于反查
            # with open(self.dp.myCFG.D1_PATH_EXPORT + precursor_id + '.txt', 'w')as f:
            #     start_position = toolFindNeighborFromSortedList1(x_list, rt_start)
            #     end_position = toolFindNeighborFromSortedList1(x_list, rt_end)
            #     f.write(str(start_position) + '\t' + str(end_position) + '\n')
            #     for i in range(inputEvidence.MATRIX_PROFILE.shape[0]):
            #         out_list = inputEvidence.MATRIX_PROFILE[i]
            #         out_list = [str(i) for i in out_list]
            #         f.write('\t'.join(out_list) + '\n')

            logToUser('Visualization completed:\t{}'.format(precursor_id))
            # plt.show()
            plt.clf()

    def __captainDrawCompare(self, inputListEvidence, inputListSeed, inputListIndex):  # 暂时弃用

        figsize = (18, len(inputListEvidence) * 3 + 1)
        plt.figure(figsize=figsize)
        gs = GridSpec(figsize[1], 18)
        text_base_title = 'Base Peak:'
        text_sequence_title = 'Sequence:'
        text_raw_title = 'Raw File:'
        text_qvalue_title = 'q value:'
        text_charge_title = 'Charge:'
        text_modify_title = 'Modification:'
        text_midScan_titile = 'Mid Scan'
        inputSeed = inputListSeed[0]
        base_intensity = 0.
        rt_start, rt_end = -1, -1
        for i in range(len(inputListEvidence)):
            tmp_intensity = np.max(inputListEvidence[i].MATRIX_PROFILE)
            if tmp_intensity > base_intensity:
                base_intensity = tmp_intensity
            tmp_rt_start = np.min(inputListEvidence[i].LIST_RET_TIME) / 60
            tmp_rt_end = np.max(inputListEvidence[i].LIST_RET_TIME) / 60
            if tmp_rt_start < rt_start or rt_start == -1:
                rt_start = tmp_rt_start
            if tmp_rt_end > rt_end or rt_end == -1:
                rt_end = tmp_rt_end
        # print(rt_start, rt_end)
        ax = plt.subplot(gs[0:1, 0:])
        text_base_value = format(base_intensity, '.2e')
        text_sequence_value = inputSeed.PRECURSOR_SQE
        text_charge_value = str(inputSeed.PRECURSOR_CHARGE)
        text_modify_value = inputSeed.PRECURSOR_MOD[:]
        text_qvalue_value = format(self.dp.myDDAID.PSM10_SCORE0[inputListIndex[0]], '.3e')
        if len(text_modify_value) == 0:
            text_modify_value = 'None'
        length_mod = len(text_modify_value)
        plt.xlim((0, 200))
        plt.ylim((0, 10))
        plt.axis('off')
        plt.xticks([])
        plt.yticks([])
        plt.text(2, 3, text_base_title, fontdict={'size': 10, 'color': 'black'})
        plt.text(16, 3, text_base_value, fontdict={'size': 10, 'color': 'red'})
        plt.text(92 + length_mod - 3, 3, text_sequence_title, fontdict={'size': 10, 'color': 'black'})
        plt.text(108 + length_mod - 3, 3, text_sequence_value, fontdict={'size': 10, 'color': 'red'})
        plt.text(28, 3, text_charge_title, fontdict={'size': 10, 'color': 'black'})
        plt.text(36, 3, text_charge_value, fontdict={'size': 10, 'color': 'red'})
        plt.text(45, 3, text_qvalue_title, fontdict={'size': 10, 'color': 'black'})
        plt.text(55, 3, text_qvalue_value, fontdict={'size': 10, 'color': 'red'})
        plt.text(68, 3, text_modify_title, fontdict={'size': 10, 'color': 'black'})
        plt.text(83, 3, text_modify_value, fontdict={'size': 10, 'color': 'red'})
        for tmp_i in range(len(inputListEvidence)):
            # print(tmp_i * 3)
            inputSeed = inputListSeed[tmp_i]
            inputEvidence = inputListEvidence[tmp_i]
            inputIndex = inputListIndex[tmp_i]

            n_fragment = inputEvidence.MATRIX_PROFILE.shape[0]

            text_raw_value = inputSeed.RAW_NAME
            length_raw_name = len(text_raw_value)

            matrix_profile = inputEvidence.MATRIX_PROFILE / base_intensity
            ax = plt.subplot(gs[tmp_i * 3 + 1: tmp_i * 3 + 2, 0:])

            text_midScan_value = inputSeed.MID_SCAN

            plt.xlim((0, 200))
            plt.ylim((0, 10))
            plt.xticks([])
            plt.yticks([])
            plt.text(2, 1, text_raw_title, fontdict={'size': 10, 'color': 'black'})
            plt.text(14, 1, text_raw_value, fontdict={'size': 10, 'color': 'red'})
            plt.text(42 - 5 + length_raw_name, 1, text_midScan_titile, fontdict={'size': 10, 'color': 'black'})
            plt.text(57 - 5 + length_raw_name, 1, text_midScan_value, fontdict={'size': 10, 'color': 'red'})


            plt.axis('off')
            # print(inputEvidence.LIST_I_START)

            ax2 = plt.subplot(gs[tmp_i * 3 + 2: tmp_i * 3 + 4, 0:])
            for i in range(n_fragment):

                x_list = [i / 60 for i in inputEvidence.LIST_RET_TIME]
                y_list = matrix_profile[i]
                plt.plot(x_list, y_list, label=format(inputSeed.MOZ_CLC[i], '.2f')+'m/z', color=DRAW_COLOR[i], alpha=0.7)
                # plt.scatter(x_list, y_list, s=0.5, color=DRAW_COLOR[i])
            plt.legend(loc='upper right', frameon=False,)
            plt.yticks([0.0, 0.5, 1.0])
            plt.ylim(0.0, 1.0)
            plt.xlim(rt_start, rt_end)
            plt.xlabel('RT(minute)')
            plt.ylabel('Intensity (%)')
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            # plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + text_qvalue_value + '_' + precursor_id + '_' +
            #             IO_NAME_FILE_EXPORT_Flow2[0], type='png')
            # plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + text_qvalue_value + '_' + precursor_id + '_' +
            #             IO_NAME_FILE_EXPORT_Flow2[0], type='svg')
            # plt.close()
        plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + text_sequence_value + '_' + text_charge_value + '_' +
                    IO_NAME_FILE_EXPORT_Flow1[0], dpi=300)
        plt.show()

        plt.close()

    def draw(self):

        # 得到要画图的结果列表
        if self.dp.myCFG.C7_PLOT_NUM == CFG_TYPE_PLOT_NUM['ALL']:
            # 对config中输入鉴定结果文件的所有结果进行画图
            listIndex = [i for i in range(self.dp.myDDAID.N_PSM)]
            listPrecursor_info = []

        elif self.dp.myCFG.C7_PLOT_NUM == CFG_TYPE_PLOT_NUM['CUSTOM']:
            # 这里是先前版本直接在config文件输入要画图的结果（scan,raw），现已不再使用
            listPrecursor_id = op_FILL_LIST_DRAW_DDA_PRECURSOR(self.dp.myCFG.C8_DDA_PRECURSOR_SCAN,
                                                               self.dp.myCFG.C9_DDA_PRECURSOR_RAW)
            listSEQ = self.dp.myCFG.C11_DDA_PRECURSOR_SEQUENCE.split('|')
            if listSEQ[-1] == '':
                listSEQ = listSEQ[:-1]
            listRT = self.dp.myCFG.C13_DDA_RT.split('|')
            if listRT[-1] == '':
                listRT = listRT[:-1]
            listChar = self.dp.myCFG.C12_DDA_PRECURSOR_CHARGE.split('|')
            if listChar[-1] == '':
                listChar = listChar[:-1]
            listPrecursor_info = [[float(listRT[i]), listSEQ[i], int(listChar[i])] for i in range(len(listChar))]
            # 得到结果列表在ID中的索引列表
            listIndex = []

            for precursor_id in listPrecursor_id:
                position = -1
                num = self.dp.myDDAID.PSM6_PRECURSOR_ID.count(precursor_id)
                # print(self.dp.myDDAID.PSM6_PRECURSOR_ID)
                if num > 0:
                    for i in range(num):
                        index = self.dp.myDDAID.PSM6_PRECURSOR_ID[position+1:].index(precursor_id)
                        listIndex.append(index)
                        position = index
                        print(precursor_id + ' in ID ( position: {:d} )'.format(index))
                else:
                    print(precursor_id + '  not in ID')
        else:
            info = "MSFunctionDraw, MK449, " + self.dp.myCFG.C7_PLOT_NUM + " is all right?"
            logGetError(info)

        # 得到每个结果的evidence
        listEvidence = [[] * 1] * len(listIndex)
        listSeed = [[] * 1] * len(listIndex)
        self.__captainGetEvidence(listEvidence, listSeed, listIndex, listPrecursor_info)

        # 进行母离子色谱曲线的画图
        self.__captainDraw(listEvidence, listSeed, listIndex)

