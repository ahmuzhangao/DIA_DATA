import os
from MSLogging import logGetError, INFO_TO_USER_TaskReadID, INFO_TO_USER_TaskReadMS1, INFO_TO_USER_TaskReadMS2
from MSSystem import CFG_TYPE_IDENTIFICATION_RESULT, CFG_TYPE_MS1, CFG_TYPE_MS2
from MSOperator import op_FILL_LIST_PATH_ID, op_FILL_LIST_PATH_MS
from MSData import CDataPack
from MSFunction import CFunctionParseIDForDIANN, CFunctionParseMS1, CFunctionParseMS2, \
    CFunctionINI,CFunctionParseIDForMaxDIA, CFunctionParseIDForMSFragger, CFunctionParseIDForMSFraggerPin
from MSFunctionDraw import CFunctionDrawDIAXIC


class CTaskReadINI:

    def __init__(self, inputDP):
        self.dp = inputDP

    def work(self):
        # 加载INI文件到全局数据结构self.dp.myINI中，后面用于计算母离子碎片离子质量
        functionLoadINI = CFunctionINI(self.dp)
        functionLoadINI.file2ini()


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


class CTaskReadIDForDraw:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def work(self):

        op_FILL_LIST_PATH_ID(self.dp.myCFG.B1_PATH_IDENTIFICATION_RESULT, self.dp.LIST_PATH_ID)
        N_ID = 0
        path_library = []
        op_FILL_LIST_PATH_ID(self.dp.myCFG.B4_PATH_LIBRARY_RESULT, path_library)
        # if self.dp.myCFG.B2_TYPE_IDENTIFICATION_RESULT == CFG_TYPE_IDENTIFICATION_RESULT['DIA-NN']:
        #
        #     for path in self.dp.LIST_PATH_ID:
        #         print(INFO_TO_USER_TaskReadID[0] + path)
        #         functionRead = CFunctionParseIDForDIANN(self.dp)
        #         functionRead.read(path)
        #         print(INFO_TO_USER_TaskReadID[1] + str(self.dp.myID.N_ID - N_ID))
        #         N_ID = self.dp.myID.N_ID

        if self.dp.myCFG.B2_TYPE_IDENTIFICATION_RESULT == CFG_TYPE_IDENTIFICATION_RESULT['MSFragger']:

            for path in self.dp.LIST_PATH_ID:
                print(INFO_TO_USER_TaskReadID[0] + path)
                functionRead = CFunctionParseIDForMSFragger(self.dp)
                functionRead.read(path)
                print(INFO_TO_USER_TaskReadID[1] + str(self.dp.myID.N_ID - N_ID))
                N_ID = self.dp.myID.N_ID

        elif self.dp.myCFG.B2_TYPE_IDENTIFICATION_RESULT == CFG_TYPE_IDENTIFICATION_RESULT['MSFragger-pin']:

            for path in self.dp.LIST_PATH_ID:
                print(INFO_TO_USER_TaskReadID[0] + path)
                functionRead = CFunctionParseIDForMSFraggerPin(self.dp)
                functionRead.read(path)
                print(INFO_TO_USER_TaskReadID[1] + str(self.dp.myID.N_ID - N_ID))
                N_ID = self.dp.myID.N_ID

        elif self.dp.myCFG.B2_TYPE_IDENTIFICATION_RESULT == CFG_TYPE_IDENTIFICATION_RESULT['MaxDIA']:

            for i, path in enumerate(self.dp.LIST_PATH_ID):
                print(INFO_TO_USER_TaskReadID[0] + path)
                functionRead = CFunctionParseIDForMaxDIA(self.dp)
                functionRead.read(path, path_library[i])
                print(INFO_TO_USER_TaskReadID[1] + str(self.dp.myDDAID.N_PSM - N_ID))
                N_ID = self.dp.myDDAID.N_PSM



        if N_ID > 1000000:  # 检查要画图的数目，如果大于1000000，输出信息，退出

            InfoError = 'Number of plot is too many, exit...'
            logGetError(InfoError)


class CTaskReadMS1:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def work(self):

        if self.dp.myCFG.A2_TYPE_MS1 == CFG_TYPE_MS1['MS1']:

            op_FILL_LIST_PATH_MS(self.dp.myCFG.A1_PATH_MS1, self.dp.LIST_PATH_MS1, [".ms1", ".MS1"])

            for path in self.dp.LIST_PATH_MS1:
                print(INFO_TO_USER_TaskReadMS1[0] + path)
                functionParseMS1 = CFunctionParseMS1(self.dp)
                functionParseMS1.ms1Topkl(path)


class CTaskReadMS2:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def work(self):

        if self.dp.myCFG.A4_TYPE_MS2 == CFG_TYPE_MS2['MS2']:

            op_FILL_LIST_PATH_MS(self.dp.myCFG.A3_PATH_MS2, self.dp.LIST_PATH_MS2, [".ms2", '.MS2'])

            for path in self.dp.LIST_PATH_MS2:
                print(INFO_TO_USER_TaskReadMS2[0] + path)
                functionParseMS2 = CFunctionParseMS2(self.dp)
                functionParseMS2.ms2topkl(path)


class CTaskDrawDIAXIC:

    def __init__(self, inputDP: CDataPack):
        self.dp = inputDP

    def work(self):
        functionDrawPrecursor = CFunctionDrawDIAXIC(self.dp)
        functionDrawPrecursor.draw()



