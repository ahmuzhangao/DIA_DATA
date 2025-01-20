import math
import pickle
from  tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pc
from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import MultipleLocator
from MSData import CDataPack, CSeed,CSeedRerank, CFileMS2, CEvidenceDIACheck
from MSEmass import CEmass
from MSFunctionComposition import CFunctionComposition
from MSTool import toolFindNeighborFromSortedList1, toolGetNameFromPath, toolMaxContinueBY, toolCosineForList, \
    toolPearsonForNumpy, toolSpearmanForNumpy
from MSLogging import logGetError
from MSOperator import op_CAL_FRAGMENT_MOZ, op_Data_fill_plot, op_INIT_CSEED_DIACHECK, op_INIT_CSEED_RERANK_DIACHECK
from MSFunctionEvidence import CFunctionEvidenceForDIACheck
from MSFunction import CFunctionParseMS1, CFunctionParseMS2
from MSSystem import VALUE_ILLEGAL, CFG_TYPE_ACCURACY_HALF_WIN_PEAK, DRAW_COLOR, TOP_N
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.utils.data as torchData
import torch.utils.model_zoo as model_zoo
from torch.autograd import Variable

class CFunctionDIARerankFeature:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __soliderGetPSMComposition(self, inputSEQ, inputMOD, inputCharge):

        # 得到鉴定信息
        tmpSEQ = inputSEQ
        tmpMOD = inputMOD
        tmpGLC = ''
        tmpLIK = ''

        # 得到最基本的元素组成
        fComposition = CFunctionComposition(self.dp)
        composition = fComposition.getStrComposition(tmpSEQ, tmpMOD, tmpGLC, tmpLIK, self.dp.myINI)
        DICT_COMPOSITION = fComposition.getDictComposition(composition)
        emass = CEmass(self.dp)
        tmpIsoDis = emass.getCalculatedIsotopicPeaks(DICT_COMPOSITION, inputCharge)  # 这一行是做什么？？？
        return tmpIsoDis[0]

    def __soliderCalFragMoz(self, iID):
        ionPrecursor = self.dp.myIDForDIACHeck.ID4_SEQ_WITH_MOD[iID]

        ionTypeList, ionMozList, ionMozAddList, ionMozSubList = op_CAL_FRAGMENT_MOZ(self.dp.myINI, ionPrecursor)

        return ionTypeList, ionMozList, ionMozAddList, ionMozSubList

    def __soliderFillSeed(self, inputIndex, inputdataMS2: CFileMS2):

        # 初始化seed
        myseed = CSeedRerank()
        op_INIT_CSEED_RERANK_DIACHECK(myseed)

        if self.dp.myIDForDIACHeck.ID3_RT_APEX[inputIndex] != VALUE_ILLEGAL:
            indexMS2 = toolFindNeighborFromSortedList1(inputdataMS2.INDEX_RT,self.dp.myIDForDIACHeck.ID3_RT_APEX[inputIndex])
        else:
            info = "MSFunctionDIACheck CFunctionDrawDIACheck in Flow5, MK56, RT:{:.2f} " "or Scan{:d} is right ?".format(self.dp.myIDForDIACHeck.ID3_RT_APEX[inputIndex],self.dp.myIDForDIACHeck.ID2_SCAN_ID[inputIndex])
            logGetError(info)

        tmpSEQ = self.dp.myIDForDIACHeck.ID0_SEQ[inputIndex]
        tmpMOD = self.dp.myIDForDIACHeck.ID16_MOD[inputIndex]
        tmpCharge = self.dp.myIDForDIACHeck.ID7_CHARGE[inputIndex]
        if len(self.__soliderGetPSMComposition(tmpSEQ, tmpMOD, tmpCharge)) < 4:
            myseed.DIS_ISO_INT_CLC = self.__soliderGetPSMComposition(tmpSEQ, tmpMOD, tmpCharge) + [1]
        else:
            myseed.DIS_ISO_INT_CLC = self.__soliderGetPSMComposition(tmpSEQ, tmpMOD, tmpCharge)[0:4]
        myseed.DIS_ISO_INT_CLC = [i / max(myseed.DIS_ISO_INT_CLC ) for i in myseed.DIS_ISO_INT_CLC ]

        Charge = self.dp.myIDForDIACHeck.ID7_CHARGE[inputIndex]
        Precursor_moz = self.dp.myIDForDIACHeck.ID5_PRECURSOR_MOZ_EXP[inputIndex]
        # 计算母离子同位素峰质荷比，在这里算下理论强度，如果是重排序
        for i in [0, 1, 2, 3]:
            isotope_moz = Precursor_moz + i * self.dp.myINI.MASS_PROTON_MONO / Charge
            myseed.DIS_ISO_MOZ_CLC.append(isotope_moz)

        myseed.DIS_FRA_TYPE_CLC, myseed.DIS_FRA_MOZ_CLC, myseed.DIS_FRA_MOZ_ADD_CLC, myseed.DIS_FRA_MOZ_SUB_CLC = self.__soliderCalFragMoz(inputIndex)

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
        num_spectra = []

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
                    num_spectra += [len(i) for i in dataMS2.MATRIX_PEAK_MOZ if len(i) != 0]
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
        num_spectra = np.mean(num_spectra)
        return num_spectra

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


    def __soliderGetDIAFeatures(self, inputSeed: CSeed, inputEvidence: CEvidenceDIACheck, inputIndex: int, spectraNum):

            # 根据seed和evidence计算每一个鉴定结果Lib中的色谱曲线
            moz_frag = inputSeed.DIS_FRA_MOZ_CLC
            type_frag = inputSeed.DIS_FRA_TYPE_CLC
            moz_prec = inputSeed.DIS_ISO_MOZ_CLC
            n_fragment = len(type_frag)
            n_iostope = len(moz_prec)
            n_rt = len(inputEvidence.LIST_RET_TIME)
            accuracy = self.dp.myCFG.C2_DIA_ACCURACY_HALF_WIN_PEAK

            MATRIX_PROFILE_FRAG = np.zeros(shape=(n_fragment, n_rt))  # 碎片离子色谱曲线矩阵
            MATRIX_PROFILE_PREC = np.zeros(shape=(n_iostope, n_rt))  # 母离子色谱曲线矩阵
            MATRIX_MASS_FRAG = np.ones(shape=(n_fragment, n_rt))
            MATRIX_MASS_PREC = np.ones(shape=(n_iostope, n_rt))

            SeqLength = len(self.dp.myIDForDIACHeck.ID0_SEQ[inputIndex])
            ScoreMatch = []

            LdList = np.array([len(i) for i in inputEvidence.MATRIX_MS2_PEAK_MOZ]).reshape(-1,len(inputEvidence.MATRIX_MS2_PEAK_MOZ)
                                                                                           )

            for i in range(n_rt):

                tmp_MOZ_list = inputEvidence.MATRIX_MS2_PEAK_MOZ[i]
                tmp_INT_list = inputEvidence.MATRIX_MS2_PEAK_INT[i]
                tmp_MOZ_list_ms1 = inputEvidence.MATRIX_MS1_PEAK_MOZ[i]
                tmp_INT_list_ms1 = inputEvidence.MATRIX_MS1_PEAK_INT[i]

                match_fragment_num = 0
                match_b_list = [0] * SeqLength
                match_y_list = [0] * SeqLength
                mass_sum = 0

                if self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['PPM']:

                    for j, tmp_isotope_moz in enumerate(moz_prec):
                        indexINT_ms1 = toolFindNeighborFromSortedList1(tmp_MOZ_list_ms1, tmp_isotope_moz)
                        mass_dev = abs(tmp_MOZ_list_ms1[indexINT_ms1] - tmp_isotope_moz) / tmp_isotope_moz * 1e6
                        if mass_dev < accuracy:
                            MATRIX_PROFILE_PREC[j, i] = tmp_INT_list_ms1[indexINT_ms1]

                            mass_bias = mass_dev / accuracy
                            MATRIX_MASS_PREC[j, i] = mass_bias


                    for j, fragment_moz in enumerate(moz_frag):

                        indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, fragment_moz)

                        mass_dev = abs(tmp_MOZ_list[indexINT] - fragment_moz) / fragment_moz * 1e6
                        if mass_dev < accuracy:
                            MATRIX_PROFILE_FRAG[j, i] = tmp_INT_list[indexINT]

                            match_fragment_num += 1

                            mass_bias = (mass_dev / accuracy)
                            mass_sum += mass_bias
                            MATRIX_MASS_FRAG[j, i] = mass_bias

                            if type_frag[j][0] == 'b':
                                match_b_list[int(type_frag[j][1])] = 1
                            elif type_frag[j][0] == 'y':
                                match_y_list[int(type_frag[j][1])] = 1


                elif self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK == CFG_TYPE_ACCURACY_HALF_WIN_PEAK['Da']:

                    for j, tmp_isotope_moz in enumerate(moz_prec):
                        indexINT_ms1 = toolFindNeighborFromSortedList1(tmp_MOZ_list_ms1, tmp_isotope_moz)
                        mass_dev = abs(tmp_MOZ_list_ms1[indexINT_ms1] - tmp_isotope_moz)
                        if mass_dev < accuracy:
                            MATRIX_PROFILE_PREC[j, i] = tmp_INT_list_ms1[indexINT_ms1]

                            mass_bias = mass_dev / accuracy
                            MATRIX_MASS_PREC[j, i] = mass_bias

                    for j, fragment_moz in enumerate(moz_frag):

                        indexINT = toolFindNeighborFromSortedList1(tmp_MOZ_list, fragment_moz)

                        mass_dev = abs(tmp_MOZ_list[indexINT] - fragment_moz)
                        if mass_dev < accuracy:
                            MATRIX_PROFILE_FRAG[j, i] = tmp_INT_list[indexINT]

                            match_fragment_num += 1

                            mass_bias = (mass_dev / accuracy)
                            mass_sum += mass_bias
                            MATRIX_MASS_FRAG[j, i] = mass_bias

                            if type_frag[j][0] == 'b':
                                match_b_list[int(type_frag[j][1])] = 1
                            elif type_frag[j][0] == 'y':
                                match_y_list[int(type_frag[j][1])] = 1


                else:
                    logGetError('MSFunctionEvidence,TYPE_ACCURACY:' + str(
                        self.dp.myCFG.C3_DIA_TYPE_ACCURACY_HALF_WIN_PEAK) + 'is all right?')
            #  根据peak rt的强度对b/y离子和isotope 母离子进行排序

                score_match_num = match_fragment_num / SeqLength
                score_maxMatchLenB = toolMaxContinueBY(match_b_list) / SeqLength
                score_maxMatchLenY = toolMaxContinueBY(match_y_list) / SeqLength
                score_mass_dev = mass_sum / SeqLength
                ScoreMatch.append([score_match_num,score_maxMatchLenB,score_maxMatchLenY,score_mass_dev])

            iMid = inputEvidence.MID_IDNEX  # 中心谱图位次
            score_match_mid = ScoreMatch[iMid] # 中心谱图的matchscore
            ScoreMatch[iMid:iMid+1] = []
            ScoreMatch = np.array(ScoreMatch)
            score_rt = [float(len(inputEvidence.MATRIX_MS2_PEAK_MOZ))]
            score_match_rest = list(np.mean(abs(ScoreMatch - np.array(score_match_mid)), axis=0)) # 除中心谱图以外的matchscore均值


            score_mass_frag_dev_bm = np.cos((MATRIX_MASS_FRAG - np.mean(MATRIX_MASS_FRAG, axis=1).reshape(MATRIX_MASS_FRAG.shape[0], -1)) * np.pi / 2)
            # 质量偏差矩阵均值归一化，求cos
            score_mass_frag_dev = np.cos(MATRIX_MASS_FRAG * np.pi / 2)
            # 不做均值归一的质量偏差矩阵
            # 这里处理一下rest的问题
            frag_rank_list = np.argsort(-1 * MATRIX_PROFILE_FRAG[:, inputEvidence.MID_IDNEX])

            tmp_sub = []
            tmp_add = []
            for i, index in enumerate(range(frag_rank_list.shape[0])):
                if i < TOP_N:
                    inputSeed.DIS_TOP6_TYPE_CLC.append(inputSeed.DIS_FRA_TYPE_CLC[index])
                    inputSeed.DIS_TOP6_MOZ_CLC.append(inputSeed.DIS_FRA_MOZ_CLC[index])
                    tmp_add.append(inputSeed.DIS_FRA_MOZ_ADD_CLC[index])
                    tmp_sub.append(inputSeed.DIS_FRA_MOZ_SUB_CLC[index])
                else:
                    inputSeed.DIS_RES_TYPE_CLC.append(inputSeed.DIS_FRA_TYPE_CLC[index])
                    inputSeed.DIS_RES_MOZ_CLC.append(inputSeed.DIS_FRA_MOZ_CLC[index])

            inputSeed.DIS_FRA_MOZ_ADD_CLC = tmp_add
            inputSeed.DIS_FRA_MOZ_SUB_CLC = tmp_sub
            # inputSeed.DIS_TOP6_TYPE_CLC = [inputSeed.DIS_FRA_TYPE_CLC[i] for i in frag_rank_list[:TOP_N]]
            # inputSeed.DIS_RES_TYPE_CLC = [inputSeed.DIS_FRA_TYPE_CLC[i] for i in frag_rank_list[TOP_N:]]

            # 这里要处理一下sub和add的质荷比
            MATRIX_PROFILE_FRAG = MATRIX_PROFILE_FRAG[frag_rank_list]
            MATRIX_MASS_FRAG = MATRIX_MASS_FRAG[frag_rank_list]
            # FRAG_MOZ = [moz_frag[i] for i in frag_rank_list][:12]
            # FRAG_TYPE = [type_frag[i] for i in frag_rank_list][:12]

            score_mid_bm = []
            score_mid = []
            score_rest_bm = []
            score_rest = []

            k1List = [1.2]
            b = 0.75
            # 超参数
            MATRIX_PROFILE_FRAG_Log10 = np.log10(MATRIX_PROFILE_FRAG + (MATRIX_PROFILE_FRAG == 0.))
            # 强度矩阵取log10，把为0的强度置为1，log后为0
            for k1 in k1List:
                BM_TF = (k1 + 1) * MATRIX_PROFILE_FRAG_Log10 / (MATRIX_PROFILE_FRAG_Log10 + k1 * (1 - b + b * LdList / spectraNum))
                BM_Score_Mass_Dev_bm = BM_TF * score_mass_frag_dev_bm
                BM_Score_Mass_Dev = BM_TF * score_mass_frag_dev
                # 归一化与不归一化的质量偏差矩阵

                BM_Score_RT_bm = np.sum(BM_Score_Mass_Dev_bm, axis=0)
                BM_Score_RT = np.sum(BM_Score_Mass_Dev, axis=0)

                BM_Score_Mid_bm = BM_Score_RT_bm[inputEvidence.MID_IDNEX]
                BM_Score_Mid = BM_Score_RT[inputEvidence.MID_IDNEX]

                BM_Score_Mid_Rest_bm = (np.sum(BM_Score_RT_bm) - 2 * BM_Score_Mid_bm) / (len(BM_Score_RT_bm) - 1)
                BM_Score_Mid_Rest = (np.sum(BM_Score_RT) - 2 * BM_Score_Mid) / (len(BM_Score_RT) - 1)

                score_mid_bm.append(BM_Score_Mid_bm)
                score_mid.append(BM_Score_Mid)
                score_rest_bm.append(BM_Score_Mid_Rest_bm)
                score_rest.append(BM_Score_Mid_Rest)

            match_mass_mid = score_mass_frag_dev[:, iMid]
            match_mass_left = score_mass_frag_dev[:, iMid - 1]
            match_mass_right = score_mass_frag_dev[:, iMid + 1]

            match_mass_mid_frag = (match_mass_mid > 0.)
            match_mass_left_frag = (match_mass_left > 0.)
            match_mass_right_frag = (match_mass_right > 0.)

            match_mass_mid_and_left = match_mass_mid_frag * match_mass_left_frag
            match_mass_mid_and_right = match_mass_mid_frag * match_mass_right_frag

            match_mass_mid_or_left = match_mass_mid_frag + match_mass_left_frag
            match_mass_mid_or_right = match_mass_mid_frag + match_mass_right_frag

            score_match_left_and = match_mass_left * match_mass_mid_and_left # 左谱图质量偏差得分
            score_match_right_and = match_mass_right * match_mass_mid_and_right # 右谱图质量偏差得分
            score_match_mid_left_and = match_mass_mid * match_mass_mid_and_left # 中心谱图质量偏差得分
            score_match_mid_right_and = match_mass_mid * match_mass_mid_and_right # 中心谱图质量偏差得分

            score_match_left_or = match_mass_left * match_mass_mid_or_left # 左谱图质量偏差得分
            score_match_right_or = match_mass_right * match_mass_mid_or_right # 右谱图质量偏差得分
            score_match_mid_left_or = match_mass_mid * match_mass_mid_or_left # 中心谱图质量偏差得分
            score_match_mid_right_or = match_mass_mid * match_mass_mid_or_right

            score_cos_similarity_mid_and_left = toolCosineForList(score_match_left_and,score_match_mid_left_and)
            score_cos_similarity_mid_and_right = toolCosineForList(score_match_right_and,score_match_mid_right_and)
            score_pearson_similarity_mid_and_left = toolPearsonForNumpy(score_match_left_and,score_match_mid_left_and)
            score_pearson_similarity_mid_and_right = toolPearsonForNumpy(score_match_right_and,score_match_mid_right_and)
            score_spearman_similarity_mid_and_left = toolSpearmanForNumpy(score_match_left_and,score_match_mid_left_and)
            score_spearman_similarity_mid_and_right = toolSpearmanForNumpy(score_match_right_and,score_match_mid_right_and)

            score_similarity_match_and = [score_cos_similarity_mid_and_left,score_cos_similarity_mid_and_right,score_pearson_similarity_mid_and_left,
                                      score_pearson_similarity_mid_and_right,score_spearman_similarity_mid_and_left,score_spearman_similarity_mid_and_right]

            score_cos_similarity_mid_or_left = toolCosineForList(score_match_left_or,score_match_mid_left_or)
            score_cos_similarity_mid_or_right = toolCosineForList(score_match_right_or,score_match_mid_right_or)
            score_pearson_similarity_or_and_left = toolPearsonForNumpy(score_match_left_or,score_match_mid_left_or)
            score_pearson_similarity_or_and_right = toolPearsonForNumpy(score_match_right_or,score_match_mid_right_or)
            score_spearman_similarity_or_and_left = toolSpearmanForNumpy(score_match_left_or,score_match_mid_left_or)
            score_spearman_similarity_or_and_right = toolSpearmanForNumpy(score_match_right_or,score_match_mid_right_or)

            score_similarity_match_or = [score_cos_similarity_mid_or_left,score_cos_similarity_mid_or_right,score_pearson_similarity_or_and_left,
                                         score_pearson_similarity_or_and_right,score_spearman_similarity_or_and_left,score_spearman_similarity_or_and_right]

            score_best_frag = self.__soliderBestFrgScore(inputSeed,inputEvidence,MATRIX_PROFILE_PREC,MATRIX_MASS_PREC,MATRIX_PROFILE_FRAG,MATRIX_MASS_FRAG)

            score_profile = [(len(inputSeed.DIS_TOP6_TYPE_CLC) + len(inputSeed.DIS_RES_TYPE_CLC)) / SeqLength]

            return score_mid_bm + score_mid + score_rest_bm + score_rest + score_match_mid + score_rt + score_match_rest + score_similarity_match_and + score_similarity_match_or + score_best_frag + score_profile

    def __soliderBestFrgScore(self, inputSeed: CSeed, inputEvidence: CEvidenceDIACheck,intensity_precursor , mass_dev_precursor, intensity_frag : np.ndarray , mass_dev_frag):
        # 最佳碎片离子打分
        score_top12_frag_similarity = [0.] * 12 # 色谱曲线之间pearson系数
        score_top6_similarity = [0.] * 3 # 碎片离子在1，0.45，0.2倍质量精度下，色谱曲线相似度和
        score_rest_similarity = [0.] * 2
        score_pre_similarity_ms2 = 0
        score_pre_similarity = [0.] * 3
        score_pre_isotope_similarity = [0.] * 3
        score_top6_isotope_add_similarity = 0.
        score_top6_isotope_sub_similarity = 0.
        score_sum_isotope = 0

        #把输入的矩阵进行一个处理，取top6和剩下的rest,这里还没有排序，仅仅是考虑了数量上的问题，在取top色谱曲线哪里把类型也处理好，最好是在seed里就直接排序好
        intensity_top6 = intensity_frag[:6]
        mass_dev_top6 = mass_dev_frag[:6]
        intensity_rest = intensity_frag[6:]
        mass_dev_rest = mass_dev_frag[6:]

        mass_dev_pre_frag = np.ones(shape=(1, len(inputEvidence.MATRIX_MS2_PEAK_MOZ)))
        intensity_pre_frag = np.zeros(shape=(1,len(inputEvidence.MATRIX_MS2_PEAK_MOZ)))
        mass_dev_frag_add = np.ones(shape=(len(inputSeed.DIS_TOP6_MOZ_CLC), len(inputEvidence.MATRIX_MS2_PEAK_MOZ)))
        intensity_frag_add = np.zeros(shape=(len(inputSeed.DIS_TOP6_MOZ_CLC), len(inputEvidence.MATRIX_MS2_PEAK_MOZ)))
        mass_dev_frag_sub = np.ones(shape=(len(inputSeed.DIS_TOP6_MOZ_CLC), len(inputEvidence.MATRIX_MS2_PEAK_MOZ)))
        intensity_frag_sub = np.zeros(shape=(len(inputSeed.DIS_TOP6_MOZ_CLC), len(inputEvidence.MATRIX_MS2_PEAK_MOZ)))

        # 母离子色谱曲线，进行一个起止的处理

        # 求一些没求的色谱曲线
        for tmp_ms2 in range(len(inputEvidence.MATRIX_MS2_PEAK_MOZ)):
            ms2MozList = inputEvidence.MATRIX_MS2_PEAK_MOZ[tmp_ms2]
            ms2IntList = inputEvidence.MATRIX_MS2_PEAK_INT[tmp_ms2]
            mozListFragAdd = inputSeed.DIS_FRA_MOZ_ADD_CLC
            mozListFragSub = inputSeed.DIS_FRA_MOZ_SUB_CLC
            mozListPreMs2 = inputSeed.DIS_ISO_MOZ_CLC[0:1]

            for i, tmp_moz in enumerate(mozListPreMs2):
                index = toolFindNeighborFromSortedList1(ms2MozList, tmp_moz)
                moz_match = ms2MozList[index]
                massDev = abs(tmp_moz - moz_match) / moz_match * 1e6
                if massDev <= 20:
                    mass_bias = massDev / 20
                    mass_dev_pre_frag[i][tmp_ms2] = mass_bias
                    intensity_pre_frag[i][tmp_ms2] = ms2IntList[index]

            for i, tmp_moz in enumerate(mozListFragAdd):
                index = toolFindNeighborFromSortedList1(ms2MozList, tmp_moz)
                moz_match = ms2MozList[index]
                massDev = abs(tmp_moz - moz_match) / moz_match * 1e6
                if massDev <= 20:
                    mass_bias = massDev / 20
                    mass_dev_frag_add[i][tmp_ms2] = mass_bias
                    intensity_frag_add[i][tmp_ms2] = ms2IntList[index]

            for i, tmp_moz in enumerate(mozListFragSub):
                index = toolFindNeighborFromSortedList1(ms2MozList, tmp_moz)
                moz_match = ms2MozList[index]
                massDev = abs(tmp_moz - moz_match) / moz_match * 1e6
                if massDev <= 20:
                    mass_bias = massDev / 20
                    mass_dev_frag_sub[i][tmp_ms2] = mass_bias
                    intensity_frag_sub[i][tmp_ms2] = ms2IntList[index]

        # 求最佳碎片离子,这里可能是个速度瓶颈
        max_similarity = -1
        max_type = 0
        score_similarity_frag_b = 0
        score_similarity_frag_y = 0

        for i in range(intensity_top6.shape[0]):
            sum_similarity = 0
            tmp_similarity_b = 0
            tmp_similarity_y = 0
            tmp_top12 = [0.] * 12
            for j in range(intensity_top6.shape[0]):
                tmp_similarity = toolPearsonForNumpy(intensity_top6[i], intensity_top6[j])
                sum_similarity += tmp_similarity
                if inputSeed.DIS_FRA_TYPE_CLC[i][0] == 'b':
                    tmp_similarity_b += tmp_similarity
                else:
                    tmp_similarity_y += tmp_similarity
                tmp_top12[j] = tmp_similarity
            if sum_similarity > max_similarity:
                max_type = i
                max_similarity = sum_similarity
                score_similarity_frag_y = tmp_similarity_y
                score_similarity_frag_b = tmp_similarity_b
                score_top12_frag_similarity = tmp_top12
        # max_type就是最佳碎片离子的索引

        # b,y各有几个
        # 这里最好处理一下rest碎片
        rest_sum_similarity = 0
        for i in range(intensity_rest.shape[0]):
            tmp_similarity = toolPearsonForNumpy(intensity_rest[i], intensity_top6[max_type])
            if inputSeed.DIS_FRA_TYPE_CLC[i][0] == 'b':
                score_similarity_frag_b += tmp_similarity
            else:
                score_similarity_frag_y += tmp_similarity
            rest_sum_similarity += tmp_similarity
            if i <= 5:
                score_top12_frag_similarity[6 + i] = tmp_similarity

        score_rest_similarity[0] = rest_sum_similarity
        score_rest_similarity[1] = rest_sum_similarity / intensity_rest.shape[0]

        Num_frag_b = inputSeed.DIS_FRA_TYPE_CLC.count('b')
        Num_frag_y = len(inputSeed.DIS_FRA_TYPE_CLC) - Num_frag_b
        Num_frag_b = np.math.factorial(Num_frag_b)
        Num_frag_y = np.math.factorial(Num_frag_y)
        score_Fragger = np.log2(float((Num_frag_b * Num_frag_y * max(score_similarity_frag_b, 1) * max(score_similarity_frag_b, 1))))

        score_top6_similarity[0] = max_similarity

        score_mass_0_45 = 0
        intensity_top6_0_45 = intensity_top6 * ((mass_dev_top6 <= 0.45) + 0)
        for i in range(intensity_top6_0_45.shape[0]):
            score_mass_0_45 += toolPearsonForNumpy(intensity_top6_0_45[i],intensity_top6[max_type])
        score_top6_similarity[1] = score_mass_0_45

        score_mass_0_2 = 0
        intensity_top6_0_2 = intensity_top6 * ((mass_dev_top6 <= 0.2) + 0)
        for i in range(intensity_top6_0_2.shape[0]):
            score_mass_0_2 += toolPearsonForNumpy(intensity_top6_0_2[i],intensity_top6[max_type])
        score_top6_similarity[2] = score_mass_0_2

        score_pre_ms2 = 0
        for i in range(intensity_pre_frag.shape[0]):
            score_pre_ms2 += toolPearsonForNumpy(intensity_pre_frag[i],intensity_top6[max_type])
        score_pre_similarity_ms2 = score_pre_ms2

        pre_intensity_top6 = toolPearsonForNumpy(intensity_precursor[0],intensity_top6[max_type])
        score_pre_similarity[0] = pre_intensity_top6

        pre_intensity_0_45 = intensity_precursor[0] * ((mass_dev_precursor[0] <= 0.45) + 0)
        score_pre_similarity[1] = toolPearsonForNumpy(pre_intensity_0_45, intensity_top6[max_type])

        pre_intensity_0_2 = intensity_precursor[0] * ((mass_dev_precursor[0] <= 0.2) + 0)
        score_pre_similarity[2] = toolPearsonForNumpy(pre_intensity_0_2, intensity_top6[max_type])

        # 同位素母离子共洗脱
        sum_iostope = 0
        for i in range(intensity_precursor.shape[0] -1):
            tmp_similarity = toolPearsonForNumpy(intensity_precursor[i + 1],intensity_top6[max_type])
            sum_iostope += tmp_similarity
            score_pre_isotope_similarity[i] = tmp_similarity
        # 含有1，2，3个同位素的母离子和best碎片离子的皮尔逊

        score_similarity_frag_add = 0
        for i in range(intensity_frag_add.shape[0]):
            tmp_similarity = toolPearsonForNumpy(intensity_frag_add[i], intensity_top6[max_type])
            score_similarity_frag_add += tmp_similarity
        score_top6_isotope_add_similarity = score_similarity_frag_add
        sum_iostope += score_similarity_frag_add

        score_similarity_frag_sub = 0
        for i in range(intensity_frag_sub.shape[0]):
            tmp_similarity = toolPearsonForNumpy(intensity_frag_sub[i], intensity_top6[max_type])
            score_similarity_frag_sub += tmp_similarity
        score_top6_isotope_sub_similarity = score_similarity_frag_sub
        sum_iostope += score_similarity_frag_sub

        score_sum_isotope = sum_iostope

        intensity_top6_new = np.zeros(shape=[6, intensity_top6.shape[1]])
        mass_dev_frag_new = np.ones(shape=[6, intensity_top6.shape[1]])
        intensity_rest_new = np.zeros(shape=[6, intensity_top6.shape[1]])
        intensity_frag_add_new = np.zeros(shape=[6,intensity_frag_add.shape[1]])
        mass_dev_frag_add_new = np.ones(shape=[6,intensity_frag_add.shape[1]])
        intensity_rest_add_new = np.zeros(shape=[6,intensity_frag_add.shape[1]])
        intensity_frag_sub_new = np.zeros(shape=[6,intensity_frag_sub.shape[1]])
        mass_dev_frag_sub_new = np.ones(shape=[6,intensity_frag_sub.shape[1]])
        intensity_rest_sub_new = np.zeros(shape=[6,intensity_frag_sub.shape[1]])

        # 处理母离子
        # list_ratio_frag = np.array(intensity_top6).reshape((-1, 1))
        list_ratio_frag = np.ones(shape=(intensity_top6.shape[0], 1))
        list_ratio_pre = np.array(inputSeed.DIS_ISO_INT_CLC).reshape((-1, 1))

        for i in range(intensity_top6.shape[0]):
            intensity_top6_new[i,:] = intensity_top6[i,:]
            mass_dev_frag_new[i,:] = mass_dev_frag[i,:]
            intensity_frag_add_new[i,:] = intensity_frag_add[i,:]
            mass_dev_frag_add_new[i:] = mass_dev_frag_add[i,:]
            intensity_frag_sub_new[i,:] = intensity_frag_sub[i,:]
            mass_dev_frag_sub_new[i,:] = mass_dev_frag_sub[i,:]
            intensity_rest_new[i,:] = (intensity_top6[i,:] / list_ratio_frag[i]) if list_ratio_frag.any() else intensity_top6[i,:]
            intensity_rest_add_new[i, :] = (intensity_frag_add[i, :] / list_ratio_frag[i]) if list_ratio_frag.any() else intensity_top6[i,:]
            intensity_rest_sub_new[i, :] = (intensity_frag_sub[i, :] / list_ratio_frag[i]) if list_ratio_frag.any() else intensity_top6[i,:]

        intensity_precursor_ratio = intensity_precursor / list_ratio_pre

        mass_dev_precursor = np.log2(1 - mass_dev_precursor + 1)
        mass_dev_frag_new = np.log2(1 - mass_dev_frag_new + 1)
        mass_dev_frag_add_new = np.log2(1 - mass_dev_frag_add_new + 1)
        intensity_top6_new = intensity_top6_new / np.max(intensity_top6_new + 1) * np.max(np.log2(intensity_top6_new + 1))
        intensity_precursor = intensity_precursor / np.max(intensity_precursor + 1) * np.max(np.log2(intensity_precursor + 1))
        intensity_precursor_ratio = intensity_precursor_ratio / np.max(intensity_precursor_ratio + 1) * np.max(np.log2(intensity_precursor_ratio + 1))
        intensity_frag_add_new = intensity_frag_add_new / np.max(intensity_frag_add_new + 1) * np.max(np.log2(intensity_frag_add_new + 1))
        intensity_rest_add_new = intensity_rest_add_new / np.max(intensity_rest_add_new + 1) * np.max(np.log2(intensity_rest_add_new + 1))
        intensity_rest_sub_new = intensity_rest_sub_new / np.max(intensity_rest_sub_new + 1) * np.max(np.log2(intensity_rest_sub_new + 1))


        return [score_Fragger] + score_top12_frag_similarity + score_top6_similarity + score_rest_similarity + \
               [score_pre_similarity_ms2] + score_pre_similarity + score_pre_isotope_similarity + \
               [score_top6_isotope_add_similarity] + [score_top6_isotope_sub_similarity] + [score_sum_isotope]

    def __captainDrawDIACheck(self, inputListSeed, inputListEvidence, inputListIndex, Spectra_Num):

        # 遍历结果列表中所有鉴定结果，对每一个结果分别画结果的library色谱曲线图和在结果rt_start到rt_end保留时间范围内的PSM标图
        out_score = []
        out_index = []
        for tmp_i in tqdm(range(len(inputListIndex)),desc="Get DIA Identification Features"):
            if inputListSeed[tmp_i]:  # 有的鉴定结果可能是无效的（如不存在文件名）
                tmp_index = inputListIndex[tmp_i]
                tmp_seed = inputListSeed[tmp_i]
                tmp_evidence = inputListEvidence[tmp_i]
                tmp_target = self.dp.myIDForDIACHeck.ID22_TARGET[tmp_i]
                score = self.__soliderGetDIAFeatures(tmp_seed, tmp_evidence, tmp_index, Spectra_Num)
                # 提取打分,把打分写到硬盘里
                score.insert(0,self.dp.myIDForDIACHeck.ID23_SCORE[tmp_i])
                score.append(tmp_target)
                out_score.append(score)
                out_index.append(tmp_i)
            else:
                pass

        with open(self.dp.myCFG.D1_PATH_EXPORT + 'Features_Pre_Rerank.tsv','w') as f:
            for i in range(len(out_score)):
                tmp_feature = [str(i) for i in out_score[i]]
                f.write('\t'.join(tmp_feature) + '\n')

        out_score = np.array(out_score)
        out_index = np.array(out_index)
        np.save(self.dp.myCFG.D1_PATH_EXPORT + 'Feature_Pre_Rerank.npy', out_score)
        np.save(self.dp.myCFG.D1_PATH_EXPORT + 'Index_Pre_Rerank.npy', out_index)

        del out_score

    def feature(self):

        # 获取要画图的结果在ID数据结构中的位置
        listIndex = list(range(self.dp.myIDForDIACHeck.N_ID))

        listEvidence = [[] * 1] * len(listIndex)
        listSeed = [[] * 1] * len(listIndex)

        Spectra_Num = self.__captainGetEvidence(listEvidence, listSeed, listIndex)

        self.__captainDrawDIACheck(listSeed, listEvidence, listIndex, Spectra_Num)


class CFunctionDIARerank:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def rerank(self):
        self.dp.myFeature.FEATURE_MATRIX = np.load(self.dp.myCFG.D1_PATH_EXPORT + 'Feature_Pre_Rerank.npy')
        self.dp.myFeature.ID_INDEX_LIST = np.load(self.dp.myCFG.D1_PATH_EXPORT + 'Index_Pre_Rerank.npy')

        FDR = 0.01
        Iteration = 5
        Iteration_target = 2
        Batch_Size = 5000
        net = Net()
        criterion = nn.BCELoss()
        optimizer = optim.SGD(net.parameters(), lr=0.15, momentum=0.9)

        features = self.dp.myFeature.FEATURE_MATRIX
        NoTD_index = features[:,-1] == -1

        Target_Decoy_index = features[:,-1] != -1
        features = features[Target_Decoy_index]

        Target_index = []
        for i in range(len(Target_Decoy_index)):
            if Target_Decoy_index[i] == 1:
                Target_index.append(i)
        # 正例的所有索引

        features_copy = features.copy()
        features_write = features.copy()

        features_decoy = []
        for i in range(features.shape[0]):
            if features[i][-1] == 0:
                features_decoy.append(features[i][:-1])
        # 拿到所有负例
        features_decoy = np.array(features_decoy)
        label_decoy = 0 * np.ones(shape=[features_decoy.shape[0],1])

        features = features[np.lexsort(features[:,::-1].T)]
        features = features[::-1, ]
        # 特征矩阵从大到小

        for i in range(Iteration):
            feature_target = self.__getCalMatirxByFDR(features, FDR)
            label_target = np.ones(shape=[feature_target.shape[0], 1])
            input = np.row_stack((features_decoy,feature_target)).astype(float)
            input = torch.tensor(input, dtype=torch.float32)
            label = np.row_stack((label_decoy, label_target)).astype(float)
            label = torch.tensor(label, dtype=torch.float32)
            dataset = torchData.TensorDataset(input, label)
            loader = torchData.DataLoader(
                dataset=dataset,
                batch_size=Batch_Size,
                shuffle=True,
                num_workers=0,
            )
            for k in range(Iteration_target):
                for j, (batch_x, batch_y) in enumerate(loader):
                    optimizer.zero_grad()
                    output = net(batch_x)
                    loss = criterion(output, batch_y)
                    loss.backward()
                    optimizer.step()
            output = net(input)
            loss = criterion(output,label)
            loss = loss.item() / input.shape[0] * 10000
            # 到这里，一轮训练就算结束了

            features = features.astype(float)
            features = torch.tensor(features, dtype=torch.float32)
            output = net(features[:,:-1])
            predict = output.detach().numpy()
            features = np.column_stack((features,predict))
            features = features[np.lexsort(-features.T)]
            features = features[:,:-1]

        features_copy = features_copy.astype(float)
        features_copy = torch.tensor(features_copy,dtype=torch.float32)
        output = net(features_copy[:, :-1])
        # pre就是预测的分值
        predict = output.detach().numpy()
        features = np.column_stack((features_copy, predict))
        features = features[np.lexsort(-features.T)]
        features = features[:, :-1]
        feature_target = self.__getCalMatirxByFDR(features, FDR)
        index = np.array([[self.dp.myFeature.ID_INDEX_LIST[i], predict[i][0]] for i in range(len(predict))])
        index = index[np.lexsort(-index.T)]

        for j,i in enumerate(index):
            n = int(i[0])
            self.dp.myResult.N_RESULT += 1
            self.dp.myResult.RE0_RAW_NAME.append(self.dp.myIDForDIACHeck.ID1_RAW_NAME[n])
            self.dp.myResult.RE1_SEQ.append(self.dp.myIDForDIACHeck.ID0_SEQ[n])
            self.dp.myResult.RE2_SCAN_NO.append(VALUE_ILLEGAL)
            # self.dp.myResult.RE2_SCAN_NO.append(self.dp.myIDForDIACHeck.ID2_SCAN_ID[n])
            self.dp.myResult.RE3_RT.append(self.dp.myIDForDIACHeck.ID3_RT_APEX[n])
            self.dp.myResult.RE4_MOD.append(self.dp.myIDForDIACHeck.ID16_MOD[n])
            self.dp.myResult.RE5_PRECURSOR_MOZ_CLC.append(self.dp.myIDForDIACHeck.ID5_PRECURSOR_MOZ_EXP[n])
            self.dp.myResult.RE6_CHARGE.append(self.dp.myIDForDIACHeck.ID7_CHARGE[n])
            self.dp.myResult.RE7_PRECURSOR_ID.append(self.dp.myIDForDIACHeck.ID4_SEQ_WITH_MOD[n])
            self.dp.myResult.RE8_SCORE0.append(self.dp.myIDForDIACHeck.ID23_SCORE[n])
            self.dp.myResult.RE9_SCORE1.append('0')
            self.dp.myResult.RE10_SCORE_RERANK.append(i[1])
            self.dp.myResult.RE11_TARGET.append(self.dp.myIDForDIACHeck.ID22_TARGET[n])
            self.dp.myResult.RE12_FEATURE.append(features_write[n,:])

        torch.save(net.state_dict(), self.dp.myCFG.D1_PATH_EXPORT + str(feature_target.shape[0]) + '_' + 'model.pkl')

    def __getCalMatirxByFDR(self, matrix, t, decoyNum=0):
        target_number = 0
        decoy_number = 0
        FDR_LIST = []
        for i in range(matrix.shape[0]):
            if matrix[i][-1] == 1:
                target_number += 1
            else:
                decoy_number += 1
            if target_number == 0:
                FDR_LIST.append(0)
            else:
                FDR_LIST.append(decoy_number / target_number)
        min_FDR = FDR_LIST[-1]
        Q_LIST = []
        for i in reversed(range(matrix.shape[0])):
            if min_FDR > FDR_LIST[i]:
                min_FDR = FDR_LIST[i]
            Q_LIST.append(min_FDR)
        Q_LIST.reverse()
        num = 0
        for i in range(len(Q_LIST) - 1, -1, -1):
            if Q_LIST[i] <= t:
                num = i + 1
                break
        matrix_target = []
        for i in range(num):
            if matrix[i][-1] == 1:
                matrix_target.append(matrix[i][:-1])
        matrix_target = np.array(matrix_target)
        # if matrix_target.shape[0] > decoyNum:
        #     matrix_target = matrix_target[:decoyNum, :]
        return matrix_target


class Net(nn.Module):

    def __init__(self):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(55, 128)
        self.fc2 = nn.Linear(128, 64)
        self.fc3 = nn.Linear(64, 32)
        self.fc4 = nn.Linear(32, 8)
        self.fc5 = nn.Linear(8, 1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        x = F.relu(self.fc4(x))
        x = self.sigmoid(self.fc5(x))
        return x


class CFunctionDIAReport:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def __captainCalqvalue(self, inputScoreList: list, inputTargetList: list):

        # 输出结果文件
        with open(self.dp.myCFG.D1_PATH_EXPORT + 'rerank_report.tsv', 'w') as f:

            tableReport = ['raw_name', 'seq', 'scanNo', 'RT', 'Mod', 'precursor m/z', 'charge', 'precursorID',
                           'cScore', 'Evidence', 'xScore', 'target']

            f.write('\t'.join(tableReport) + '\n')

            for i in range(self.dp.myResult.N_RESULT):
                rawName = self.dp.myResult.RE0_RAW_NAME[i]
                seq = self.dp.myResult.RE1_SEQ[i]
                scanNo = self.dp.myResult.RE2_SCAN_NO[i]
                RT = self.dp.myResult.RE3_RT[i]
                Mod = self.dp.myResult.RE4_MOD[i]
                precursorMoz = self.dp.myResult.RE5_PRECURSOR_MOZ_CLC[i]
                charge = self.dp.myResult.RE6_CHARGE[i]
                precursorID = self.dp.myResult.RE7_PRECURSOR_ID[i]
                score_rerank = self.dp.myResult.RE8_SCORE0[i]
                score_beforeRerank = self.dp.myResult.RE9_SCORE1[i]
                score_xh = self.dp.myResult.RE10_SCORE_RERANK[i]
                target = self.dp.myResult.RE11_TARGET[i]
                feature = list(self.dp.myResult.RE12_FEATURE[i])

                contentReport = [rawName, seq, scanNo, RT, Mod, precursorMoz, charge, precursorID, score_rerank,
                                 score_beforeRerank, score_xh, target]
                contentReport.extend(feature)
                contentReport = [str(i) for i in contentReport]

                f.write('\t'.join(contentReport) + '\n')

        score_list = np.array(inputScoreList)
        target_list = np.array(inputTargetList)
        argsort_list = np.argsort(-1 * score_list)
        score_list = score_list[argsort_list]
        target_list = target_list[argsort_list]
        FDR_list = []
        target_count, decoy_count = 0, 0
        # 生成FDR列表
        for score, target in zip(score_list, target_list):
            if target == 0:
                decoy_count += 1
            elif target == 1:
                target_count += 1
                FDR_list.append(decoy_count / target_count)
        # 生成 q value列表
        qValue_list = []
        min_FDR = FDR_list[-1]
        for i in reversed(FDR_list):
            if min_FDR > i:
                min_FDR = i
            qValue_list.append(min_FDR)
        y_list = []
        x_list = []
        for i in reversed(range(len(qValue_list))):
            q = qValue_list[i]
            x_list.append(q)
            y_list.append(len(qValue_list) - i)

        target_score_list = score_list[target_list == 1]
        decoy_score_list = score_list[target_list == 0]

        return qValue_list, x_list, y_list, target_score_list, decoy_score_list

    def report(self):
       plt.figure(figsize=(8,8))
       gs = GridSpec(8, 8)
       plt.subplot(gs[:,:])
       _, x_listRe, y_listRe, target_score_listRe, decoy_score_listRe = self.__captainCalqvalue(self.dp.myResult.RE10_SCORE_RERANK, self.dp.myResult.RE11_TARGET)
       plt.plot(x_listRe, y_listRe, color='#0000FF')
       plt.plot(x_listRe, y_listRe, color='#00FF00')
       plt.legend(['Rerank', 'myRerank', 'PreRerank'])
       plt.xlabel('Estimated FDR')
       plt.ylabel('#identified precursors')
       plt.xlim(0, 0.01)
       plt.savefig(self.dp.myCFG.D1_PATH_EXPORT + 'rerank_qValue.png')
       plt.clf()
