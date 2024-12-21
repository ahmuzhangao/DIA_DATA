import math
import pickle
from  tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pc
from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import MultipleLocator
from MSData import CDataPack, CSeed, CFileMS2, CEvidenceDIACheck
from MSTool import toolFindNeighborFromSortedList1, toolGetNameFromPath
from MSLogging import logGetError
from MSOperator import op_CAL_FRAGMENT_MOZ,op_Data_fill_plot, op_INIT_CSEED_DIACHECK
from MSFunctionEvidence import CFunctionEvidenceForDIACheck
from MSFunction import CFunctionParseMS1, CFunctionParseMS2
from MSSystem import VALUE_ILLEGAL, CFG_TYPE_ACCURACY_HALF_WIN_PEAK, DRAW_COLOR, \
    IO_NAME_FILE_EXPORT_Flow2


class CFunctionDrawDIACheck:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soliderCalFragMoz(self, iID):
        ionPrecursor = self.dp.myIDForDIACHeck.ID4_SEQ_WITH_MOD[iID]

        ionTypeList, ionMozList = op_CAL_FRAGMENT_MOZ(self.dp.myINI, ionPrecursor)

        return ionTypeList, ionMozList

    def __soliderFillSeed(self, inputIndex, inputdataMS2: CFileMS2):

        # 初始化seed
        myseed = CSeed()
        op_INIT_CSEED_DIACHECK(myseed)

        if self.dp.myIDForDIACHeck.ID3_RT_APEX[inputIndex] != VALUE_ILLEGAL:
            indexMS2 = toolFindNeighborFromSortedList1(inputdataMS2.INDEX_RT,self.dp.myIDForDIACHeck.ID3_RT_APEX[inputIndex])
        else:
            info = "MSFunctionDIACheck CFunctionDrawDIACheck in Flow5, MK56, RT:{:.2f} " "or Scan{:d} is right ?".format(self.dp.myIDForDIACHeck.ID3_RT_APEX[inputIndex],self.dp.myIDForDIACHeck.ID2_SCAN_ID[inputIndex])
            logGetError(info)

        Charge = self.dp.myIDForDIACHeck.ID7_CHARGE[inputIndex]
        Precursor_moz = self.dp.myIDForDIACHeck.ID5_PRECURSOR_MOZ_EXP[inputIndex]
        # 计算母离子同位素峰质荷比
        for i in [0, 1, 2, 3]:
            isotope_moz = Precursor_moz + i * self.dp.myINI.MASS_PROTON_MONO / Charge
            myseed.DIS_ISO_MOZ_CLC.append(isotope_moz)

        myseed.DIS_FRA_TYPE_CLC, myseed.DIS_FRA_MOZ_CLC = self.__soliderCalFragMoz(inputIndex)

        myseed.MID_RT = inputdataMS2.INDEX_RT[indexMS2]
        myseed.MID_SCAN = inputdataMS2.INDEX_SCAN[indexMS2]

        return myseed

    def __soliderGetCycleWinNum(self, inputDataMS2: CFileMS2, MidIndex):

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

            for iMS2 in range(nMS2):
                # init ms2
                pathMS2 = self.dp.LIST_PATH_MS2[iMS2]
                nameRawMS2 = toolGetNameFromPath(pathMS2)

                # load
                if nameRawMS2 == nameRawMS1:
                    functionParseMS1 = CFunctionParseMS1(self.dp)
                    dataMS1 = functionParseMS1.loadPKL(pathMS1)
                    functionParseMS2 = CFunctionParseMS2(self.dp)
                    dataMS2 = functionParseMS2.loadPKL(pathMS2)
                    for i in range(len(inputListIndex)):
                        iID = inputListIndex[i]
                        nameRawID = toolGetNameFromPath(self.dp.myIDForDIACHeck.ID1_RAW_NAME[iID])

                        if nameRawID == nameRawMS1:
                            tmpSeed = self.__soliderFillSeed(iID, dataMS2)
                            inputListSeed[i] = tmpSeed
                            functionEvidence = CFunctionEvidenceForDIACheck(self.dp)
                            tmpEvidence = functionEvidence.fillEvidence(dataMS1, dataMS2, tmpSeed, iID)
                            inputListEvidence[i] = tmpEvidence
                    del dataMS2
            del dataMS1

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

    def __soliderOutputXLC(self,preInfo,fragInfo,inputIndex,rt):
        #[MATRIX_PROFILE_PREC,PREC_MOZ], [MATRIX_PROFILE_FRAG,FRAG_MOZ,FRAG_TYPE]
        tmp_frg = []
        tmp_pre = []

        tmp_frg.append(fragInfo[0])
        tmp_frg.append(fragInfo[1])
        tmp_frg.append(fragInfo[2])
        tmp_frg.append(rt)
        tmp_frg.append(self.dp.myIDForDIACHeck.ID0_SEQ[inputIndex])
        tmp_frg.append(self.dp.myIDForDIACHeck.ID7_CHARGE[inputIndex])
        tmp_frg.append(self.dp.myIDForDIACHeck.ID16_MOD[inputIndex])
        if self.dp.myIDForDIACHeck.ID22_TARGET:
            tmp_frg.append(self.dp.myIDForDIACHeck.ID22_TARGET[inputIndex])
        else:
            tmp_frg.append(1)
        self.dp.myXLCInfo.LIST_FRAG.append(tmp_frg)

        tmp_pre.append(preInfo[0])
        tmp_pre.append(preInfo[1])
        tmp_pre.append(rt)
        tmp_pre.append(self.dp.myIDForDIACHeck.ID0_SEQ[inputIndex])
        tmp_pre.append(self.dp.myIDForDIACHeck.ID7_CHARGE[inputIndex])
        tmp_pre.append(self.dp.myIDForDIACHeck.ID16_MOD[inputIndex])
        if self.dp.myIDForDIACHeck.ID22_TARGET:
            tmp_pre.append(self.dp.myIDForDIACHeck.ID22_TARGET[inputIndex])
        else:
            tmp_pre.append(1)
        self.dp.myXLCInfo.LIST_PRE.append(tmp_pre)

    def __soliderGetDIACurve(self, inputSeed: CSeed, inputEvidence: CEvidenceDIACheck, inputIndex: int):

        # 根据seed和evidence计算每一个鉴定结果Lib中的色谱曲线
        moz_frag = inputSeed.DIS_FRA_MOZ_CLC
        type_frag = inputSeed.DIS_FRA_TYPE_CLC
        moz_prec = inputSeed.DIS_ISO_MOZ_CLC
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

        # frag_rank_list = np.argsort(-1 * MATRIX_PROFILE_FRAG[:, inputEvidence.MID_IDNEX])
        # MATRIX_PROFILE_FRAG = MATRIX_PROFILE_FRAG[frag_rank_list]
        # FRAG_MOZ = [moz_frag[i] for i in frag_rank_list]
        # FRAG_TYPE = [type_frag[i] for i in frag_rank_list]

        frag_rank_list = np.argsort(-1 * MATRIX_PROFILE_FRAG[:, inputEvidence.MID_IDNEX])
        if self.dp.myCFG.B7_FRAG_NUM !=0:
            MATRIX_PROFILE_FRAG = MATRIX_PROFILE_FRAG[frag_rank_list][:self.dp.myCFG.B7_FRAG_NUM]  # 取前12个元素
            FRAG_MOZ = [moz_frag[i] for i in frag_rank_list][:self.dp.myCFG.B7_FRAG_NUM]  # 取前12个元素
            FRAG_TYPE = [type_frag[i] for i in frag_rank_list][:self.dp.myCFG.B7_FRAG_NUM]
            if len(FRAG_TYPE) < self.dp.myCFG.B7_FRAG_NUM:
                TMP_FRAG = (self.dp.myCFG.B7_FRAG_NUM - len(FRAG_TYPE)) * ['*']
                TMP_MOZ = (self.dp.myCFG.B7_FRAG_NUM - len(FRAG_TYPE)) * [0]
                TMP_PROFILE =  np.zeros((self.dp.myCFG.B7_FRAG_NUM - len(FRAG_TYPE), MATRIX_PROFILE_FRAG.shape[1]))

                FRAG_MOZ += TMP_MOZ
                FRAG_TYPE += TMP_FRAG
                MATRIX_PROFILE_FRAG = np.vstack((MATRIX_PROFILE_FRAG, TMP_PROFILE))
        else:
            MATRIX_PROFILE_FRAG = MATRIX_PROFILE_FRAG[frag_rank_list]  # 取前12个元素
            FRAG_MOZ = [moz_frag[i] for i in frag_rank_list] # 取前12个元素
            FRAG_TYPE = [type_frag[i] for i in frag_rank_list]

        PREC_MOZ = moz_prec
        PREC_TYPE = ['precursor', 'precursor[M+1]', 'precursor[M+2]','precursor[M+3]']

        self.__soliderOutputXLC([MATRIX_PROFILE_PREC,PREC_MOZ], [MATRIX_PROFILE_FRAG,FRAG_MOZ,FRAG_TYPE],inputIndex,inputEvidence.LIST_RET_TIME)

    def __soliderDrawDIACurve(self, inputSeed: CSeed, inputEvidence: CEvidenceDIACheck, inputIndex: int, gs, fig):

        # 根据seed和evidence计算每一个鉴定结果Lib中的色谱曲线
        moz_frag = inputSeed.DIS_FRA_MOZ_CLC
        type_frag = inputSeed.DIS_FRA_TYPE_CLC
        moz_prec = inputSeed.DIS_ISO_MOZ_CLC
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

        # frag_rank_list = np.argsort(-1 * MATRIX_PROFILE_FRAG[:, inputEvidence.MID_IDNEX])
        # MATRIX_PROFILE_FRAG = MATRIX_PROFILE_FRAG[frag_rank_list]
        # FRAG_MOZ = [moz_frag[i] for i in frag_rank_list]
        # FRAG_TYPE = [type_frag[i] for i in frag_rank_list]

        frag_rank_list = np.argsort(-1 * MATRIX_PROFILE_FRAG[:, inputEvidence.MID_IDNEX])

        if self.dp.myCFG.B7_FRAG_NUM !=0:
            MATRIX_PROFILE_FRAG = MATRIX_PROFILE_FRAG[frag_rank_list][:self.dp.myCFG.B7_FRAG_NUM]  # 取前12个元素
            FRAG_MOZ = [moz_frag[i] for i in frag_rank_list][:self.dp.myCFG.B7_FRAG_NUM]  # 取前12个元素
            FRAG_TYPE = [type_frag[i] for i in frag_rank_list][:self.dp.myCFG.B7_FRAG_NUM]
        else:
            MATRIX_PROFILE_FRAG = MATRIX_PROFILE_FRAG[frag_rank_list]  # 取前12个元素
            FRAG_MOZ = [moz_frag[i] for i in frag_rank_list] # 取前12个元素
            FRAG_TYPE = [type_frag[i] for i in frag_rank_list]


        PREC_MOZ = moz_prec
        PREC_TYPE = ['precursor', 'precursor[M+1]', 'precursor[M+2]','precursor[M+3]']

        self.__soliderOutputXLC([MATRIX_PROFILE_PREC,PREC_MOZ], [MATRIX_PROFILE_FRAG,FRAG_MOZ,FRAG_TYPE],inputIndex,inputEvidence.LIST_RET_TIME)

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

            for i in range(len(FRAG_TYPE)):
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
        text_precursor_value = self.dp.myIDForDIACHeck.ID8_PRECURSOR_ID[inputIndex]
        text_raw_value = toolGetNameFromPath(self.dp.myIDForDIACHeck.ID1_RAW_NAME[inputIndex])
        text_charge_value = str(self.dp.myIDForDIACHeck.ID7_CHARGE[inputIndex])
        text_qvalue_value = format(self.dp.myIDForDIACHeck.ID9_SCORE0[inputIndex], '.5f')

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
            prec_type_text.append(PREC_TYPE[i] + ':' + format(float(PREC_MOZ[i]), '.2f') + ' m/z')

        xlim = [np.min(x_list), np.max(x_list)]
        plt.xlim(xlim)
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

    def __soliderDrawPSMLabel(self, i_label: int, inputSeed: CSeed, inputEvidence: CEvidenceDIACheck, inputIndex: int, gs, inputMaxPeak):

        clc_moz = inputSeed.DIS_FRA_MOZ_CLC
        clc_tag = inputSeed.DIS_FRA_TYPE_CLC
        exp_moz = np.array(inputEvidence.MATRIX_MS2_PEAK_MOZ[i_label + inputEvidence.RT_START])
        exp_int = np.array(inputEvidence.MATRIX_MS2_PEAK_INT[i_label + inputEvidence.RT_START])
        max_int = max(exp_int)
        # 对ms2谱图所有谱峰，获取config中设置的topN强度的峰

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
        info_raw_name = self.dp.myIDForDIACHeck.ID1_RAW_NAME[inputIndex]
        info_base_peak = '{:.3g}'.format(max_int)
        info_moz = '{:.3f}'.format(self.dp.myIDForDIACHeck.ID5_PRECURSOR_MOZ_EXP[inputIndex])
        mod_dic = self.dp.myIDForDIACHeck.ID16_MOD[inputIndex]

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

        info_pep = self.dp.myIDForDIACHeck.ID0_SEQ[inputIndex]
        info_charge = self.dp.myIDForDIACHeck.ID7_CHARGE[inputIndex]

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
        for tmp_i in tqdm(range(len(inputListIndex)),desc="Get DIA Curve Information"):
        # for tmp_i in range(len(inputListIndex)):
            if inputListSeed[tmp_i]:  # 有的鉴定结果可能是无效的（如不存在文件名）
                tmp_index = inputListIndex[tmp_i]
                tmp_seed = inputListSeed[tmp_i]
                tmp_evidence = inputListEvidence[tmp_i]
                precursor_id = self.dp.myIDForDIACHeck.ID8_PRECURSOR_ID[tmp_index]

                num_label = tmp_evidence.RT_END - tmp_evidence.RT_START + 1  # 色谱曲线保留时间范围内的ms2谱图个数

                if self.dp.myCFG.B6_FLAG_DRAW:
                    fig_Check = plt.figure(figsize=(15, 17 + num_label * 5))  # 新建画布（色谱曲线 + ms2谱图PSM标图）
                    gs = GridSpec(17 + 17 * num_label, 12)
                    self.__soliderDrawDIACurve(tmp_seed, tmp_evidence, tmp_index, gs, fig_Check)  # 画母离子和Lib中by离子的色谱曲线
                    # 计算rt范围内ms2谱图最大moz
                    max_peak = np.max([np.max(i_ms2) for i_ms2 in tmp_evidence.MATRIX_MS2_PEAK_MOZ[tmp_evidence.RT_START:tmp_evidence.RT_END]])
                    max_peak = min(max_peak, np.max(tmp_seed.DIS_FRA_MOZ_CLC))
                    # 下面这行代码是错误的，得到的是最大的一个ms2谱图moz，因为对于不同长度的列表转化成numpy数组，得到的不是二维数组
                    # 而是以列表为元素的一维数组，对齐进行max操作返回的是值最大的那个列表而不是所有ms2谱图moz列表中最大的数值
                    # max_peak = np.max(np.array(tmp_evidence.MATRIX_MS2_PEAK_MOZ[tmp_evidence.RT_START:tmp_evidence.RT_END]).reshape([-1,]))
                    for i in tqdm(range(num_label), desc="Draw DIA Curve"):
                        i_label = i
                        self.__soliderDrawPSMLabel(i_label, tmp_seed, tmp_evidence, tmp_index, gs, max_peak)  # 画RT范围上ms2谱图的标图
                    run_precursor_id = self.dp.myIDForDIACHeck.ID00_RUN[tmp_index] + '_' + self.dp.myIDForDIACHeck.ID8_PRECURSOR_ID[tmp_index]
                    print(run_precursor_id)
                    plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + run_precursor_id + '_' +
                                IO_NAME_FILE_EXPORT_Flow2[0], format='png')
                    plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + run_precursor_id + '_' +
                                IO_NAME_FILE_EXPORT_Flow2[1], format='svg')
                    plt.close(fig_Check)
                else:
                    self.__soliderGetDIACurve(tmp_seed, tmp_evidence, tmp_index)  # 画母离子和Lib中by离子的色谱曲线
            else:
                pass



    def draw(self):

        # 获取要画图的结果在ID数据结构中的位置
        listIndex = list(range(self.dp.myIDForDIACHeck.N_ID))

        listEvidence = [[] * 1] * len(listIndex)
        listSeed = [[] * 1] * len(listIndex)

        self.__captainGetEvidence(listEvidence, listSeed, listIndex)

        self.__captainDrawDIACheck(listSeed, listEvidence, listIndex)

        with open(self.dp.myCFG.D1_PATH_EXPORT + 'pre.pkl', 'wb') as f:
            pickle.dump(self.dp.myXLCInfo.LIST_PRE, f)

        with open(self.dp.myCFG.D1_PATH_EXPORT + 'frag.pkl', 'wb') as f:
            pickle.dump(self.dp.myXLCInfo.LIST_FRAG, f)