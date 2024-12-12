import numpy as np
from MSData import CSeed, CSeedDIACheck, CSeed_FLOW1, CSeed2, \
    CFileMS1, CDataPack,CFileMS2,\
    CEvidence, CEvidence_FLOW1, CEvidenceDIACheck
from MSSystem import CFG_TYPE_ACCURACY_HALF_WIN_PEAK,VALUE_ILLEGAL
from MSTool import toolFindNeighborFromSortedList1,toolFindIndexFromSortedList1
from MSLogging import logGetError, logGetWarning
from MSOperator import opGetStartAndEndForProfile, op_ININT_CEVIDENCE_DIACHECK


class CFunctionEvidence2:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __captainFillEvidence2(self, cycle_win_num, inputwinRT, iMid, inputdataMS1: CFileMS1, inputdataMS2: CFileMS2,
                               inputSeed: CSeed2, flagGetStartAndEnd):

        # rtMid = inputdataMS2.INDEX_RT[iMid]
        rtMid = inputSeed.MID_RT
        rtStart = inputSeed.START_RT
        rtEnd = inputSeed.END_RT
        rtStartIndex = 0
        rtEndIndex = 0
        listScanLeftIndex = []
        listScanRightIndex = []

        iLeft = iMid
        borderRTLeft = max(0, rtMid - inputwinRT * 60)

        while iLeft > cycle_win_num:
            if rtStart != VALUE_ILLEGAL and inputdataMS2.INDEX_RT[iLeft] >= rtStart:
                rtStartIndex = iLeft
            if inputdataMS2.INDEX_RT[iLeft] < borderRTLeft:
                break
            else:
                listScanLeftIndex.append(iLeft)
                iLeft = iLeft - cycle_win_num
        if rtStart != VALUE_ILLEGAL:
            rtStartPosition = list(reversed(listScanLeftIndex)).index(rtStartIndex)
        else:
            rtStartPosition = VALUE_ILLEGAL
        iRight = iMid + cycle_win_num
        borderRTRight = min(rtMid + inputwinRT * 60, inputdataMS2.INDEX_RT[-1])

        while iRight < inputdataMS2.INDEX_SCAN[-1] - cycle_win_num:
            if rtEndIndex != VALUE_ILLEGAL and inputdataMS2.INDEX_RT[iRight] <= rtEnd:
                rtEndIndex = iRight
            if inputdataMS2.INDEX_RT[iRight] > borderRTRight:
                break
            else:
                listScanRightIndex.append(iRight)
                iRight = iRight + cycle_win_num
        if rtEnd != VALUE_ILLEGAL:
            rtEndPosition = listScanRightIndex.index(rtEndIndex) + len(listScanLeftIndex)
        else:
            rtEndPosition = VALUE_ILLEGAL
        # 合并两个scan列表
        listScan = list(reversed(listScanLeftIndex)) + listScanRightIndex
        nScan = len(listScan)
        # print(iMid)
        # print('list scan', listScan)
        # print('nScan', nScan)

        listFragment_MOZ = inputSeed.MOZ_CLC
        accuracy = self.dp.myCFG.C2_DIA_ACCURACY_HALF_WIN_PEAK

        # init
        outputEvidence = CEvidence2()

        outputEvidence.MOZ_CLC = inputSeed.MOZ_CLC
        outputEvidence.LIST_ION_TYPE = inputSeed.ION_TYPE
        outputEvidence.MATRIX_PROFILE = np.zeros(shape=[len(listFragment_MOZ) + 1, nScan])
        outputEvidence.MATRIX_MASS_DEV = np.zeros(shape=[len(listFragment_MOZ) + 1, nScan])
        outputEvidence.LIST_RET_TIME = [VALUE_ILLEGAL] * nScan
        outputEvidence.LIST_SCAN = [VALUE_ILLEGAL] * nScan

        for i, iScan in enumerate(listScan):

            # print(i, iScan)
            tmp_MOZ_list = inputdataMS2.MATRIX_PEAK_MOZ[inputdataMS2.INDEX_SCAN[iScan]]
            tmp_INT_list = inputdataMS2.MATRIX_PEAK_INT[inputdataMS2.INDEX_SCAN[iScan]]
            # 根据ms2的scan索引去ms1中找到对应当前ms2的ms1 scan号
            iScan_ms1 = toolFindNeighborFromSortedList1(inputdataMS1.INDEX_SCAN, inputdataMS2.INDEX_SCAN[iScan])
            iScan_ms1 = iScan_ms1 if inputdataMS1.INDEX_SCAN[iScan_ms1] < inputdataMS2.INDEX_SCAN[
                iScan] else iScan_ms1 - 1
            tmp_MOZ_list_ms1 = inputdataMS1.MATRIX_PEAK_MOZ[inputdataMS1.INDEX_SCAN[iScan_ms1]]
            tmp_INT_list_ms1 = inputdataMS1.MATRIX_PEAK_INT[inputdataMS1.INDEX_SCAN[iScan_ms1]]

            outputEvidence.LIST_SCAN[i] = inputdataMS2.INDEX_SCAN[iScan]
            outputEvidence.LIST_RET_TIME[i] = inputdataMS2.INDEX_RT[iScan]

            if self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['PPM']:
                # print(tmp_MOZ_list)
                # print(precursor_MOZ)
                indexINT_ms1 = toolFindNeighborFromSortedList1(tmp_MOZ_list_ms1, inputSeed.MOZ_PREC)
                mass_dev = abs(tmp_MOZ_list_ms1[indexINT_ms1] - inputSeed.MOZ_PREC) / inputSeed.MOZ_PREC * 1e6
                if mass_dev < accuracy:
                    outputEvidence.MATRIX_PROFILE[-1, i] = tmp_INT_list_ms1[indexINT_ms1]
                    outputEvidence.MATRIX_MASS_DEV[-1, i] = mass_dev
                else:
                    outputEvidence.MATRIX_PROFILE[-1, i] = 0.0
                    outputEvidence.MATRIX_MASS_DEV[-1, i] = accuracy
                for j, fragment_moz in enumerate(listFragment_MOZ):

                    indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, fragment_moz)

                    mass_dev = abs(tmp_MOZ_list[indexINT] - fragment_moz) / fragment_moz * 1e6
                    if mass_dev < accuracy:
                        outputEvidence.MATRIX_PROFILE[j, i] = tmp_INT_list[indexINT]
                        outputEvidence.MATRIX_MASS_DEV[j, i] = mass_dev
                    else:
                        outputEvidence.MATRIX_PROFILE[j, i] = 0.0
                        outputEvidence.MATRIX_MASS_DEV[j, i] = accuracy

            elif self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['Da']:

                indexINT_ms1 = toolFindNeighborFromSortedList1(tmp_MOZ_list_ms1, inputSeed.MOZ_PREC)
                mass_dev = abs(tmp_MOZ_list_ms1[indexINT_ms1] - inputSeed.MOZ_PREC) / inputSeed.MOZ_PREC * 1e6
                if mass_dev < accuracy:
                    outputEvidence.MATRIX_PROFILE[-1, i] = tmp_INT_list_ms1[indexINT_ms1]
                    outputEvidence.MATRIX_MASS_DEV[-1, i] = mass_dev
                else:
                    outputEvidence.MATRIX_PROFILE[-1, i] = 0.0
                    outputEvidence.MATRIX_MASS_DEV[-1, i] = accuracy

                for j, fragment_moz in enumerate(listFragment_MOZ):

                    indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, fragment_moz)

                    mass_dev = abs(tmp_MOZ_list[indexINT] - fragment_moz)
                    if abs(tmp_MOZ_list[indexINT] - fragment_moz) < accuracy:
                        outputEvidence.MATRIX_PROFILE[j, i] = tmp_INT_list[indexINT]
                        outputEvidence.MATRIX_MASS_DEV[j, i] = mass_dev
                    else:
                        outputEvidence.MATRIX_PROFILE[j, i] = 0.0
                        outputEvidence.MATRIX_MASS_DEV[j, i] = accuracy

            else:

                logGetError('MSFunctionEvidence,TYPE_ACCURACY:' + str(
                    self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK) + 'is all right?')
        # 得到色谱曲线的起始点和终止点
        rt_start_list = []
        rt_end_list = []
        for i in range(outputEvidence.MATRIX_PROFILE.shape[0]):
            i_left, i_right = opGetStartAndEndForProfile(outputEvidence.MATRIX_PROFILE[i], len(listScanLeftIndex), 0.1,
                                                         2)
            rt_start_list.append(i_left)
            rt_end_list.append(i_right)
        rtStartPosition = np.argmax(np.bincount(rt_start_list))
        rtEndPosition = np.argmax(np.bincount(rt_end_list))
        outputEvidence.RT_START = rtStartPosition
        outputEvidence.RT_END = rtEndPosition

        # if flagGetStartAndEnd:
        #     outputEvidence.LIST_I_START = [0] * len(listFragment_MOZ)
        #     outputEvidence.LIST_I_END = [len(outputEvidence.LIST_RET_TIME) - 1] * len(listFragment_MOZ)  # 先默认输出全部曲线，这里可能要改，先标记一下-------
        # else:
        #     outputEvidence.LIST_I_START = [0] * len(listFragment_MOZ)
        #     outputEvidence.LIST_I_END = [len(outputEvidence.LIST_RET_TIME) - 1] * len(listFragment_MOZ)
        #
        # outputEvidence.POINT_APEX = len(listScanLeftIndex)

        return outputEvidence

    def __captainGetCycleWinNum(self, inputDataMS2: CFileMS2, MidIndex):
        
        centerMoz = inputDataMS2.LIST_ACTIVATION_CENTER[inputDataMS2.INDEX_SCAN[MidIndex]]
        try:
            num = 1
            while centerMoz != inputDataMS2.LIST_ACTIVATION_CENTER[inputDataMS2.INDEX_SCAN[MidIndex+num]]:
                num = num + 1
        except IndexError:  # 超出激活中心列表的边界，从左边计算
            num = 1
            while centerMoz != inputDataMS2.LIST_ACTIVATION_CENTER[inputDataMS2.INDEX_SCAN[MidIndex-num]]:
                num = num + 1

        return num

    def fillEvidence2(self, inputdataMS1: CFileMS1, inputdataMS2: CFileMS2, inputSeed: CSeed2):

        winRT = self.dp.myCFG.C1_DIA_RT_HALF_WIN_IN_MIN

        iMid = toolFindIndexFromSortedList1(inputdataMS2.INDEX_SCAN, inputSeed.MID_SCAN)

        if inputdataMS2.INDEX_SCAN[iMid] != inputSeed.MID_SCAN:
            logGetError("MSFunctionEvidence, MK142: Can not find scan " + str(
                inputSeed.MID_SCAN) + " in inputDataMS1.INDEX_SCAN. inputDataMS1.INDEX_SCAN[iMid] is " + str(
                inputdataMS2.INDEX_SCAN[iMid]))

        cycle_win_num = self.__captainGetCycleWinNum(inputdataMS2, iMid)
        # print('cycle_win_num', cycle_win_num)

        return self.__captainFillEvidence2(cycle_win_num, winRT, iMid, inputdataMS1, inputdataMS2, inputSeed, True)


class CFunctionEvidence:

    def __init__(self,inputDP:CDataPack):

        self.dp = inputDP

    def __captainGetCycleWinNum(self, inputDataMS2: CFileMS2, MidIndex):

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

    def  __captainGetIndexListFromMS2(self, inputDataMS1: CFileMS1, inputDataMS2: CFileMS2, Midindex, inputWinRT):

        cycle = self.__captainGetCycleWinNum(inputDataMS2, Midindex)
        MidRT = inputDataMS2.INDEX_RT[Midindex]

        RT_left = MidRT
        tmp_index = Midindex
        left_list_ms2 = []
        left_list_ms1 = []
        while (RT_left >= MidRT - inputWinRT * 60) & (tmp_index >= cycle):
            tmp_index -= cycle
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
        while (RT_right <= MidRT + inputWinRT * 60) & (tmp_index + cycle <len(inputDataMS2.INDEX_RT)):
            tmp_index += cycle
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

    def __captainfillEvidence1(self, index_list_ms1, index_list_ms2, mid, inputdataMS1: CFileMS1, inputdataMS2: CFileMS2,
                               inputSeed: CSeed_FLOW1, flagGetStartAndEnd: True):

        listScan = index_list_ms2
        listScan_ms1 = index_list_ms1
        nScan = len(listScan)


        listFragment_MOZ = inputSeed.MOZ_FRAG
        listPrec_MOZ = inputSeed.MOZ_PREC
        accuracy = self.dp.myCFG.C2_DIA_ACCURACY_HALF_WIN_PEAK

        # init

        outputEvidence = CEvidence_FLOW1()
        outputEvidence.PREC_MOZ = inputSeed.MOZ_PREC
        outputEvidence.PREC_TYPE = inputSeed.PREC_TYPE
        outputEvidence.FRAG_MOZ = inputSeed.MOZ_FRAG
        outputEvidence.FRAG_TYPE = inputSeed.FRAG_TYPE
        outputEvidence.MATRIX_PROFILE_FRAG = np.zeros(shape=[len(listFragment_MOZ), nScan])
        outputEvidence.MATRIX_MASS_DEV_FRAG = np.zeros(shape=[len(listFragment_MOZ), nScan])
        outputEvidence.MATRIX_PROFILE_PREC = np.zeros(shape=[len(listPrec_MOZ), nScan])
        outputEvidence.MATRIX_MASS_DEV_PREC = np.zeros(shape=[len(listPrec_MOZ), nScan])
        outputEvidence.LIST_SCAN = [[]] * nScan
        outputEvidence.LIST_RET_TIME = [[]] * nScan
        outputEvidence.POINT_APEX = mid

        for i, iScan in enumerate(listScan):

            # print(i, iScan)
            tmp_MOZ_list = inputdataMS2.MATRIX_PEAK_MOZ[inputdataMS2.INDEX_SCAN[iScan]]
            tmp_INT_list = inputdataMS2.MATRIX_PEAK_INT[inputdataMS2.INDEX_SCAN[iScan]]
            # 根据ms2的scan索引去ms1中找到对应当前ms2的ms1 scan号
            tmp_MOZ_list_ms1 = inputdataMS1.MATRIX_PEAK_MOZ[inputdataMS1.INDEX_SCAN[listScan_ms1[i]]]
            tmp_INT_list_ms1 = inputdataMS1.MATRIX_PEAK_INT[inputdataMS1.INDEX_SCAN[listScan_ms1[i]]]


            outputEvidence.LIST_SCAN[i] = inputdataMS2.INDEX_SCAN[iScan]
            outputEvidence.LIST_RET_TIME[i] = inputdataMS2.INDEX_RT[iScan]

            if self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['PPM']:
                # print(tmp_MOZ_list)
                # print(precursor_MOZ)
                for j, tmp_isotope_moz in enumerate(inputSeed.MOZ_PREC):
                    indexINT_ms1 = toolFindNeighborFromSortedList1(tmp_MOZ_list_ms1, tmp_isotope_moz)
                    mass_dev = abs(tmp_MOZ_list_ms1[indexINT_ms1] - tmp_isotope_moz) / tmp_isotope_moz * 1e6
                    if mass_dev < accuracy:
                        outputEvidence.MATRIX_PROFILE_PREC[j, i] = tmp_INT_list_ms1[indexINT_ms1]
                        outputEvidence.MATRIX_MASS_DEV_PREC[j, i] = mass_dev
                    else:
                        outputEvidence.MATRIX_PROFILE_PREC[j, i] = 0.0
                        outputEvidence.MATRIX_MASS_DEV_PREC[j, i] = accuracy

                for j, fragment_moz in enumerate(listFragment_MOZ):

                    indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, fragment_moz)

                    mass_dev = abs(tmp_MOZ_list[indexINT] - fragment_moz) / fragment_moz * 1e6
                    if mass_dev < accuracy:
                        outputEvidence.MATRIX_PROFILE_FRAG[j, i] = tmp_INT_list[indexINT]
                        outputEvidence.MATRIX_MASS_DEV_FRAG[j, i] = mass_dev
                    else:
                        outputEvidence.MATRIX_PROFILE_FRAG[j, i] = 0.0
                        outputEvidence.MATRIX_MASS_DEV_FRAG[j, i] = accuracy

            elif self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['Da']:

                for j, tmp_isotope_moz in enumerate(inputSeed.MOZ_PREC):
                    indexINT_ms1 = toolFindNeighborFromSortedList1(tmp_MOZ_list_ms1, tmp_isotope_moz)
                    mass_dev = abs(tmp_MOZ_list_ms1[indexINT_ms1] - tmp_isotope_moz)
                    if mass_dev < accuracy:
                        outputEvidence.MATRIX_PROFILE_PREC[j, i] = tmp_INT_list_ms1[indexINT_ms1]
                        outputEvidence.MATRIX_MASS_DEV_PREC[j, i] = mass_dev
                    else:
                        outputEvidence.MATRIX_PROFILE_PREC[j, i] = 0.0
                        outputEvidence.MATRIX_MASS_DEV_PREC[j, i] = accuracy

                for j, fragment_moz in enumerate(listFragment_MOZ):

                    indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, fragment_moz)

                    mass_dev = abs(tmp_MOZ_list[indexINT] - fragment_moz)
                    if mass_dev < accuracy:
                        outputEvidence.MATRIX_PROFILE_FRAG[j, i] = tmp_INT_list[indexINT]
                        outputEvidence.MATRIX_MASS_DEV_FRAG[j, i] = mass_dev
                    else:
                        outputEvidence.MATRIX_PROFILE_FRAG[j, i] = 0.0
                        outputEvidence.MATRIX_MASS_DEV_FRAG[j, i] = accuracy

            else:
                logGetError('MSFunctionEvidence,TYPE_ACCURACY:' + str(
                    self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK) + 'is all right?')
        #  根据peak rt的强度对b/y离子和isotope 母离子进行排序

        frag_rank_list = np.argsort(-1 * outputEvidence.MATRIX_PROFILE_FRAG[:, mid])
        prec_rank_list = np.argsort(-1 * outputEvidence.MATRIX_PROFILE_PREC[:, mid])
        outputEvidence.MATRIX_PROFILE_FRAG = outputEvidence.MATRIX_PROFILE_FRAG[frag_rank_list]
        outputEvidence.MATRIX_MASS_DEV_FRAG = outputEvidence.MATRIX_MASS_DEV_FRAG[frag_rank_list]
        outputEvidence.MATRIX_PROFILE_PREC = outputEvidence.MATRIX_PROFILE_PREC[prec_rank_list]
        outputEvidence.MATRIX_MASS_DEV_PREC = outputEvidence.MATRIX_MASS_DEV_PREC[prec_rank_list]
        FRAG_TYPE = []
        FRAG_MOZ = []
        PREC_TYPE = []
        PREC_MOZ = []
        for i in frag_rank_list:
            FRAG_TYPE.append(outputEvidence.FRAG_TYPE[i])
            FRAG_MOZ.append(outputEvidence.FRAG_MOZ[i])
        for i in prec_rank_list:
            PREC_TYPE.append(outputEvidence.PREC_TYPE[i])
            PREC_MOZ.append(outputEvidence.PREC_MOZ[i])

        if self.dp.myCFG.C18_FRAG_NUM == 0:
            outputEvidence.FRAG_TYPE = FRAG_TYPE
            outputEvidence.FRAG_MOZ = FRAG_MOZ
            outputEvidence.MATRIX_PROFILE_FRAG = outputEvidence.MATRIX_PROFILE_FRAG[frag_rank_list]
            outputEvidence.MATRIX_MASS_DEV_FRAG = outputEvidence.MATRIX_MASS_DEV_FRAG[frag_rank_list]
        else:
            outputEvidence.FRAG_TYPE = FRAG_TYPE[:self.dp.myCFG.C18_FRAG_NUM]
            outputEvidence.FRAG_MOZ = FRAG_MOZ[:self.dp.myCFG.C18_FRAG_NUM]
            outputEvidence.MATRIX_PROFILE_FRAG = outputEvidence.MATRIX_PROFILE_FRAG[frag_rank_list[:self.dp.myCFG.C18_FRAG_NUM]]
            outputEvidence.MATRIX_MASS_DEV_FRAG = outputEvidence.MATRIX_MASS_DEV_FRAG[frag_rank_list[:self.dp.myCFG.C18_FRAG_NUM]]

        outputEvidence.PREC_TYPE = PREC_TYPE
        outputEvidence.PREC_MOZ = PREC_MOZ

        # 得到色谱曲线的起始点和终止点
        rt_start_list = []
        rt_end_list = []
        if len(listFragment_MOZ) > 0:
            for i in range(outputEvidence.MATRIX_PROFILE_FRAG.shape[0]):
                i_left, i_right = opGetStartAndEndForProfile(outputEvidence.MATRIX_PROFILE_FRAG[i], mid, 0.1, 1)
                rt_start_list.append(i_left)
                rt_end_list.append(i_right)
        elif len(listPrec_MOZ) > 0:
            for i in range(outputEvidence.MATRIX_PROFILE_PREC.shape[0]):
                i_left, i_right = opGetStartAndEndForProfile(outputEvidence.MATRIX_PROFILE_PREC[i], mid, 0.1, 1)
                rt_start_list.append(i_left)
                rt_end_list.append(i_right)
        rtStartPosition = np.argmax(np.bincount(rt_start_list))
        rtEndPosition = np.argmax(np.bincount(rt_end_list))
        outputEvidence.RT_START = rtStartPosition
        outputEvidence.RT_END = rtEndPosition

        # if flagGetStartAndEnd:
        #     outputEvidence.LIST_I_START = [0] * len(listFragment_MOZ)
        #     outputEvidence.LIST_I_END = [len(outputEvidence.LIST_RET_TIME) - 1] * len(listFragment_MOZ)  # 先默认输出全部曲线，这里可能要改，先标记一下-------
        # else:
        #     outputEvidence.LIST_I_START = [0] * len(listFragment_MOZ)
        #     outputEvidence.LIST_I_END = [len(outputEvidence.LIST_RET_TIME) - 1] * len(listFragment_MOZ)
        #
        # outputEvidence.POINT_APEX = len(listScanLeftIndex)

        return outputEvidence

    def fillEvidence1(self, inputdataMS1: CFileMS1, inputdataMS2: CFileMS2, inputSeed: CSeed_FLOW1):

        winRT = self.dp.myCFG.C1_DIA_RT_HALF_WIN_IN_MIN

        iMid = toolFindIndexFromSortedList1(inputdataMS2.INDEX_SCAN, inputSeed.MID_SCAN)

        if inputdataMS2.INDEX_SCAN[iMid] != inputSeed.MID_SCAN:
            info = "MSFunctionEvidence, MK142: Can not find scan " + str(
                inputSeed.MID_SCAN) + " in inputDataMS2.INDEX_SCAN. inputDataMS2.INDEX_SCAN[iMid] is " + str(
                inputdataMS2.INDEX_SCAN[iMid])
            logGetError(info)

        # 得到保留时间窗口范围内的scan号索引列表
        index_list_ms1, index_list_ms2, mid = self.__captainGetIndexListFromMS2(inputdataMS1, inputdataMS2, iMid, winRT)

        return self.__captainfillEvidence1(index_list_ms1, index_list_ms2, mid, inputdataMS1, inputdataMS2, inputSeed, True)


class CFunctionEvidenceForDIACheck:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

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

    def __captainGetIndexListFromMS2(self, inputDataMS1: CFileMS1, inputDataMS2: CFileMS2, Midindex: int, inputWinRT: int):

        cycle = self.__soliderGetCycleWinNum(inputDataMS2, Midindex)
        MidRT = inputDataMS2.INDEX_RT[Midindex]

        RT_left = MidRT
        tmp_index = Midindex
        left_list_ms2 = []
        left_list_ms1 = []
        while (RT_left >= MidRT - inputWinRT * 60) & (tmp_index >= cycle):
            tmp_index -= cycle
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
        while (RT_right <= MidRT + inputWinRT * 60) & (tmp_index + cycle <len(inputDataMS2.INDEX_RT)):
            tmp_index += cycle
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

    def __captainfillEvidence(self, inputIndex: int, inputIndexListMS1: list, inputIndexListMS2: list, inputMid: int,
                              inputdataMS1: CFileMS1, inputdataMS2: CFileMS2, inputSeed: CSeedDIACheck):

        myEvidence = CEvidenceDIACheck()
        op_ININT_CEVIDENCE_DIACHECK(myEvidence)

        rt_start = self.dp.myIDForCheck.ID11_RTSTART[inputIndex]
        rt_end = self.dp.myIDForCheck.ID12_RTEND[inputIndex]

        myEvidence.MID_IDNEX = inputMid
        for i in inputIndexListMS1:
            myEvidence.MATRIX_MS1_PEAK_MOZ.append(inputdataMS1.MATRIX_PEAK_MOZ[inputdataMS1.INDEX_SCAN[i]])
            myEvidence.MATRIX_MS1_PEAK_INT.append(inputdataMS1.MATRIX_PEAK_INT[inputdataMS1.INDEX_SCAN[i]])
        for i in inputIndexListMS2:
            myEvidence.MATRIX_MS2_PEAK_INT.append(inputdataMS2.MATRIX_PEAK_INT[inputdataMS2.INDEX_SCAN[i]])
            myEvidence.MATRIX_MS2_PEAK_MOZ.append(inputdataMS2.MATRIX_PEAK_MOZ[inputdataMS2.INDEX_SCAN[i]])
            myEvidence.LIST_RET_TIME.append(inputdataMS2.INDEX_RT[i])
            myEvidence.LIST_SCAN.append(inputdataMS2.INDEX_SCAN[i])

        myEvidence.RT_END = 0
        myEvidence.RT_START = len(myEvidence.LIST_RET_TIME)
        for tmp_i in range(len(myEvidence.LIST_RET_TIME)):
            tmp_rt = myEvidence.LIST_RET_TIME[tmp_i]
            if tmp_rt >= rt_start * 60.0 and tmp_i < myEvidence.RT_START:
                myEvidence.RT_START = tmp_i
            if tmp_rt <= rt_end * 60.0:
                myEvidence.RT_END += 1

        return myEvidence

    def fillEvidence(self, inputdataMS1: CFileMS1, inputdataMS2: CFileMS2, inputSeed: CSeedDIACheck, inputIDIndex: int):

        winRT = self.dp.myCFG.C1_DIA_RT_HALF_WIN_IN_MIN

        iMid = toolFindNeighborFromSortedList1(inputdataMS2.INDEX_SCAN, inputSeed.MID_SCAN)

        if inputdataMS2.INDEX_SCAN[iMid] != inputSeed.MID_SCAN:
            info = 'MSFunctionEvidenceForDIACheck, MSCFileMS2: Can not find scan' + str(
                inputSeed.MID_SCAN) + " in inputDataMS2.INDEX_SCAN. inputDataMS2.INDEX_SCAN[iMid] is " + str(
                inputdataMS2.INDEX_SCAN[iMid])
            logGetError(info)

        # 得到config保留时间窗口范围内的scan号索引列表和DIA鉴定结果记录的保留时间
        index_list_ms1, index_list_ms2, mid = self.__captainGetIndexListFromMS2(inputdataMS1, inputdataMS2, iMid, winRT)

        return self.__captainfillEvidence(inputIDIndex, index_list_ms1, index_list_ms2, mid, inputdataMS1,
                                          inputdataMS2, inputSeed)

class CFunctionEvidenceDDA:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __captainfillEvidence(self, inputwinRT, iMid, inputdataMS1, inputSeed: CSeed):

        rtMid = inputdataMS1.INDEX_RT[iMid]

        # iLeft和iRight是在CFileMS1索引表INDEX_SCAN和INDEX_RT中的位置
        # INDEX_SCAN[iScan]是在MATRIX_PEAK_MOZ和MATRIX_PEAK_INT中的位置
        iLeft = iMid
        borderRTLeft = max(0, rtMid - inputwinRT * 60)

        while iLeft > 0:

            if inputdataMS1.INDEX_RT[iLeft] < borderRTLeft:
                break
            else:
                iLeft = iLeft - 1

        iLeft = iLeft + 1  # 最后一个位置是不符合的，去除

        iRight = iMid
        borderRTRight = min(rtMid + inputwinRT * 60, inputdataMS1.INDEX_RT[-1])

        while iRight < len(inputdataMS1.INDEX_SCAN) - 1:

            if inputdataMS1.INDEX_RT[iRight] > borderRTRight:
                break
            else:
                iRight = iRight + 1

        iRight = iRight - 1  # 最后一个位置是不符合的，去除
        # init
        precursor_MOZ_list = inputSeed.MOZ_CLC
        nScan = iRight - iLeft + 1
        outputEvidence = CEvidence()
        outputEvidence.MATRIX_PROFILE = np.zeros(shape=[len(precursor_MOZ_list), nScan])
        outputEvidence.LIST_RET_TIME = [0] * nScan
        outputEvidence.LIST_SCAN = [0] * nScan

        accuracy = self.dp.myCFG.C2_DIA_ACCURACY_HALF_WIN_PEAK

        for iScan in range(iLeft, iRight + 1):

            tmp_MOZ_list = inputdataMS1.MATRIX_PEAK_MOZ[inputdataMS1.INDEX_SCAN[iScan]]
            tmp_INT_list = inputdataMS1.MATRIX_PEAK_INT[inputdataMS1.INDEX_SCAN[iScan]]

            outputEvidence.LIST_SCAN[iScan - iLeft] = inputdataMS1.INDEX_SCAN[iScan]
            outputEvidence.LIST_RET_TIME[iScan - iLeft] = inputdataMS1.INDEX_RT[iScan]
            if self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['PPM']:
                # print(tmp_MOZ_list)
                # print(precursor_MOZ)
                for i in range(len(precursor_MOZ_list)):

                    indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, precursor_MOZ_list[i])

                    if abs(tmp_MOZ_list[indexINT] - precursor_MOZ_list[i]) / precursor_MOZ_list[i] * 1e6 < accuracy:
                        outputEvidence.MATRIX_PROFILE[i, iScan - iLeft] = tmp_INT_list[indexINT]

            elif self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['Da']:

                for i in range(len(precursor_MOZ_list)):

                    indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, precursor_MOZ_list[i])

                    if abs(tmp_MOZ_list[indexINT] - precursor_MOZ_list[i]) < accuracy:
                        outputEvidence.MATRIX_PROFILE[i, iScan - iLeft] = tmp_INT_list[indexINT]

            else:
                logGetError('MSFunctionEvidence,TYPE_ACCURACY:' + str(
                    self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK) + 'is all right?')

        outputEvidence.I_START = 0
        outputEvidence.I_END = -1

        return outputEvidence

    def fillEvidence(self, inputdataMS1, inputSeed: CSeed):
        winRT = self.dp.myCFG.C1_DIA_RT_HALF_WIN_IN_MIN

        iMid = toolFindIndexFromSortedList1(inputdataMS1.INDEX_SCAN, inputSeed.MID_SCAN)

        if inputdataMS1.INDEX_SCAN[iMid] != inputSeed.MID_SCAN:
            logGetError("MSFunctionEvidence, MK142: Can not find scan " + str(
                inputSeed.MID_SCAN) + " in inputDataMS1.INDEX_SCAN. inputDataMS1.INDEX_SCAN[iMid] is " + str(
                inputdataMS1.INDEX_SCAN[iMid]))

        # print('iMid', iMid)
        return self.__captainfillEvidence(winRT, iMid, inputdataMS1, inputSeed)
