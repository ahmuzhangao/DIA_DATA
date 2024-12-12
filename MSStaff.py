import datetime
from MSLogging import logGetError, INFO_TO_USER_Staff, logToUser
from MSSystem import CFG_TYPE_VISUAL
from MSData import CDataPack
from MSFunction import CFunctionConfig
from MSFlow import CFlow0, CFlow1


class CStaff:

    def __init__(self, inputArgv):

        self.dp = CDataPack()
        self.argv = inputArgv

    def start(self):

        dateNow = datetime.datetime.now()
        logToUser(str(dateNow))

        # run flow
        self.__captainRunFlow()

        # finish
        logToUser(INFO_TO_USER_Staff[1])

    def __captainRunFlow(self):

        n = len(self.argv)

        if n == 1:

            logToUser(INFO_TO_USER_Staff[0])
            flow0 = CFlow0()
            flow0.run()

        elif n == 2:

            functionConfig = CFunctionConfig()
            functionConfig.file2config(self.argv[1], self.dp.myCFG)

            if self.dp.myCFG.C0_TYPE_FLOW == CFG_TYPE_VISUAL['DIA_XIC']:
                flow1 = CFlow1(self.dp)
                flow1.run()

            elif self.dp.myCFG.C0_TYPE_FLOW == CFG_TYPE_VISUAL['Fragment']:
                flow2 = CFlow2(self.dp)
                flow2.run()

            elif self.dp.myCFG.C0_TYPE_FLOW == CFG_TYPE_VISUAL['DDA_Precursor']:
                flow3 = CFlow3(self.dp)
                flow3.run()

            elif self.dp.myCFG.C0_TYPE_FLOW == CFG_TYPE_VISUAL['Rerank']:
                flow4 = CFlow4(self.dp)
                flow4.run()

            elif self.dp.myCFG.C0_TYPE_FLOW == CFG_TYPE_VISUAL['DIA_Check']:
                flow5 = CFlow5(self.dp)
                flow5.run()

            elif self.dp.myCFG.C0_TYPE_FLOW == CFG_TYPE_VISUAL['DIA_Phenomenon']:
                flow6 = CFlow6(self.dp)
                flow6.run()

            else:

                logGetError("Get wrong TYPE_VISUAL: " + str(self.dp.myCFG.C1_TYPE_QUANT))

