import os
from MSFunctionDraw import CFunctionDrawDIACheck
from MSFunctionRerank import CFunctionDIARerankFeature, CFunctionDIARerank, CFunctionDIAReport
from MSLogging import logGetError, INFO_TO_USER_TaskReadMS1, INFO_TO_USER_TaskReadID, INFO_TO_USER_TaskReadMS2
from MSSystem import CFG_TYPE_MS1, CFG_TYPE_IDENTIFICATION_RESULT, CFG_TYPE_MS2, UNIMOID_TO_STANDARD_MOD
from MSOperator import op_FILL_LIST_PATH_MS, op_FILL_LIST_PATH_ID, op_FILL_DIC_ADDITIONAL_MOD
from MSFunction import CFunctionParseMS1, CFunctionINI, CFunctionParseMS2, CFunctionReadApuQuantIDForCheck, \
    CFunctionReadDIANNIDForCheck, CFunctionReadDIANNLibIDForCheck, CFunctionReadMSFraggerPinIDForCheck, \
    CFunctionReadMSFraggerIDForCheck
from MSData import CDataPack
from MSTool import toolGetNameFromPath


class CTaskCheck:

    def __init__(self, inputDP):
        self.dp = inputDP

    def work(self, fileList):

        #  输出文件有时候用excel打开了，需要检查下文件是否正常
        for nameFile in fileList:

            path = self.dp.myCFG.D1_PATH_EXPORT + nameFile

            if os.access(self.dp.myCFG.D1_PATH_EXPORT, os.F_OK):

                try:

                    fid1 = open(path, 'w', -1)
                    fid1.close()

                except IOError:

                    logGetError(path + ' is opened! Please close it and run the program again!')

            else:

                os.makedirs(self.dp.myCFG.D1_PATH_EXPORT)

class CTaskReadINI:

    def __init__(self, inputDP):
        self.dp = inputDP

    def work(self):

        # 加载INI文件到全局数据结构self.dp.myINI中，后面用于计算母离子碎片离子质量
        functionLoadINI = CFunctionINI(self.dp)
        functionLoadINI.file2ini()

class CTaskReadMS1:

    def __init__(self, inputDP):

        self.dp = inputDP

    def work(self):

        if self.dp.myCFG.A2_TYPE_MS1 == CFG_TYPE_MS1['MS1']:

            op_FILL_LIST_PATH_MS(self.dp.myCFG.A1_PATH_MS1, self.dp.LIST_PATH_MS1, [".ms1", ".MS1"])

            for path in self.dp.LIST_PATH_MS1:

                print(INFO_TO_USER_TaskReadMS1[0] + path)

                functionMS1 = CFunctionParseMS1(self.dp)
                functionMS1.ms1TOpkl(path)

        self.dp.LIST_PATH_MS = [toolGetNameFromPath(i) for i in self.dp.LIST_PATH_MS1] # 存一下RAW文件名称，主要用于在读ID的时候只读list中的RAW

class CTaskReadMS2:

    def __init__(self,inputDP:CDataPack):

        self.dp = inputDP

    def work(self):

        if self.dp.myCFG.A4_TYPE_MS2 == CFG_TYPE_MS2['MS2']:

            op_FILL_LIST_PATH_MS(self.dp.myCFG.A3_PATH_MS2, self.dp.LIST_PATH_MS2, [".ms2", '.MS2'])

            for path in self.dp.LIST_PATH_MS2:
                print(INFO_TO_USER_TaskReadMS2[0] + path)
                functionParseMS2 = CFunctionParseMS2(self.dp)
                functionParseMS2.ms2topkl(path)

class CTaskReadIDForDIACheck:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def work(self):

        op_FILL_LIST_PATH_ID(self.dp.myCFG.B1_PATH_IDENTIFICATION_RESULT, self.dp.LIST_PATH_ID)
        op_FILL_DIC_ADDITIONAL_MOD(self.dp.myCFG.I3_INI_PATH_ADI_MOD,UNIMOID_TO_STANDARD_MOD)

        if self.dp.myCFG.B2_TYPE_IDENTIFICATION_RESULT == CFG_TYPE_IDENTIFICATION_RESULT['DIA-NN']:
            N_ID = 0
            for i, path in enumerate(self.dp.LIST_PATH_ID):
                print(INFO_TO_USER_TaskReadID[0] + path)
                functionRead = CFunctionReadDIANNIDForCheck(self.dp)
                functionRead.read(path)
                print(INFO_TO_USER_TaskReadID[1] + str(self.dp.myIDForDIACHeck.N_ID - N_ID))
                N_ID = self.dp.myIDForDIACHeck.N_ID

        elif self.dp.myCFG.B2_TYPE_IDENTIFICATION_RESULT == CFG_TYPE_IDENTIFICATION_RESULT['DIA-NN_Lib']:
            N_ID = 0
            for i, path in enumerate(self.dp.LIST_PATH_ID):
                print(INFO_TO_USER_TaskReadID[0] + path)
                functionRead = CFunctionReadDIANNLibIDForCheck(self.dp)
                functionRead.read(path)
                print(INFO_TO_USER_TaskReadID[1] + str(self.dp.myIDForDIACHeck.N_ID - N_ID))
                N_ID = self.dp.myIDForDIACHeck.N_ID

        elif self.dp.myCFG.B2_TYPE_IDENTIFICATION_RESULT == CFG_TYPE_IDENTIFICATION_RESULT['MSFragger']:
            N_ID = 0
            for i, path in enumerate(self.dp.LIST_PATH_ID):
                print(INFO_TO_USER_TaskReadID[0] + path)
                functionRead = CFunctionReadMSFraggerIDForCheck(self.dp)
                functionRead.read(path)
                print(INFO_TO_USER_TaskReadID[1] + str(self.dp.myIDForDIACHeck.N_ID - N_ID))
                N_ID = self.dp.myIDForDIACHeck.N_ID

        elif self.dp.myCFG.B2_TYPE_IDENTIFICATION_RESULT == CFG_TYPE_IDENTIFICATION_RESULT['MSFragger-pin']:
            N_ID = 0
            for i, path in enumerate(self.dp.LIST_PATH_ID):
                print(INFO_TO_USER_TaskReadID[0] + path)
                functionRead = CFunctionReadMSFraggerPinIDForCheck(self.dp)
                functionRead.read(path)
                print(INFO_TO_USER_TaskReadID[1] + str(self.dp.myIDForDIACHeck.N_ID - N_ID))
                N_ID = self.dp.myIDForDIACHeck.N_ID

        elif self.dp.myCFG.B2_TYPE_IDENTIFICATION_RESULT == CFG_TYPE_IDENTIFICATION_RESULT['ApuQuant']:
            N_ID = 0
            for i, path in enumerate(self.dp.LIST_PATH_ID):
                print(INFO_TO_USER_TaskReadID[0] + path)
                functionRead = CFunctionReadApuQuantIDForCheck(self.dp)
                functionRead.read(path)
                print(INFO_TO_USER_TaskReadID[1] + str(self.dp.myIDForDIACHeck.N_ID - N_ID))
                N_ID = self.dp.myIDForDIACHeck.N_ID

        else:
            logGetError('TaskReadIDForCheck Error, config type_identification_result is not legal!')

class CTaskDrawForDIACheck:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def work(self):

        functionDrawDIACheck = CFunctionDrawDIACheck(self.dp)
        functionDrawDIACheck.draw()

class CTaskForDIARerankFeature:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def work(self):

        functionDIARerankFeature = CFunctionDIARerankFeature(self.dp)
        functionDIARerankFeature.feature()


class CTaskForDIARerank:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def work(self):

        functionDIARerank = CFunctionDIARerank(self.dp)
        functionDIARerank.rerank()

class CTaskForDIAReport:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def work(self):

        functionDIAReport = CFunctionDIAReport(self.dp)
        functionDIAReport.report()