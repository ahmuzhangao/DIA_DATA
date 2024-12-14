from MSLogging import logGetError
from MSData import CFileMS1, CFileMS2, CSeed, CDataPack, CEvidenceDIACheck
from MSTool import toolFindNeighborFromSortedList1
from MSOperator import op_ININT_CEVIDENCE_DIACHECK
from MSSystem import VALUE_ILLEGAL


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
                              inputdataMS1: CFileMS1, inputdataMS2: CFileMS2):

        myEvidence = CEvidenceDIACheck()
        op_ININT_CEVIDENCE_DIACHECK(myEvidence)

        if self.dp.myIDForDIACHeck.ID14_RT_BEGIN:
            rt_start = self.dp.myIDForDIACHeck.ID14_RT_BEGIN[inputIndex]
            rt_end = self.dp.myIDForDIACHeck.ID15_RT_END[inputIndex]
        else:
            rt_start = self.dp.myIDForDIACHeck.ID3_RT_APEX[inputIndex] - self.dp.myCFG.C1_DIA_RT_HALF_WIN_IN_MIN * 30
            rt_end = self.dp.myIDForDIACHeck.ID3_RT_APEX[inputIndex] + self.dp.myCFG.C1_DIA_RT_HALF_WIN_IN_MIN * 30


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
            if tmp_rt >= rt_start  and tmp_i < myEvidence.RT_START:
                myEvidence.RT_START = tmp_i
            if tmp_rt <= rt_end:
                myEvidence.RT_END += 1

        return myEvidence

    def fillEvidence(self, inputdataMS1: CFileMS1, inputdataMS2: CFileMS2, inputSeed: CSeed, inputIDIndex: int):

        winRT = self.dp.myCFG.C1_DIA_RT_HALF_WIN_IN_MIN

        iMid = toolFindNeighborFromSortedList1(inputdataMS2.INDEX_SCAN, inputSeed.MID_SCAN)

        if inputdataMS2.INDEX_SCAN[iMid] != inputSeed.MID_SCAN:
            info = 'MSFunctionEvidenceForDIACheck, MSCFileMS2: Can not find scan' + str(inputSeed.MID_SCAN) + " in inputDataMS2.INDEX_SCAN. inputDataMS2.INDEX_SCAN[iMid] is " + str(inputdataMS2.INDEX_SCAN[iMid])
            logGetError(info)

        # 得到config保留时间窗口范围内的scan号索引列表和DIA鉴定结果记录的保留时间
        index_list_ms1, index_list_ms2, mid = self.__captainGetIndexListFromMS2(inputdataMS1, inputdataMS2, iMid, winRT)

        return self.__captainfillEvidence(inputIDIndex, index_list_ms1, index_list_ms2, mid, inputdataMS1,inputdataMS2)

