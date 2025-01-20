import datetime
from MSData import CDataPack
from MSLogging import INFO_TO_USER_Staff, logToUser,logGetError
from MSFlow import CFlow0, CFlow2, CFlow3
from MSFunction import CFunctionConfig
from MSSystem import CFG_TYPE_FLOW, TIME_EXPIRATION

class CStaff:

    def __init__(self, inputArgv):

        self.dp = CDataPack()
        self.argv = inputArgv

    def start(self):

        dateNow = datetime.datetime.now()
        logToUser(str(dateNow))
        logToUser(INFO_TO_USER_Staff[0])

        self.__captainCheckTime()

        self.__captainRunFlow()

        logToUser(INFO_TO_USER_Staff[4])
        dateNow = datetime.datetime.now()
        logToUser(str(dateNow))

    def __captainRunFlow(self):

        n = len(self.argv)

        if n == 1:

            logToUser(INFO_TO_USER_Staff[3])
            flow0 = CFlow0()
            flow0.run()

        elif n == 2:

            functionConfig = CFunctionConfig()
            functionConfig.file2config(self.argv[1], self.dp.myCFG)

            if self.dp.myCFG.C0_TYPE_FLOW ==CFG_TYPE_FLOW['DIA_Evidence_Check']:
                flow2 = CFlow2(self.dp)
                flow2.run()

            elif self.dp.myCFG.C0_TYPE_FLOW ==CFG_TYPE_FLOW['DIA_Rerank']:
                flow3 = CFlow3(self.dp)
                flow3.run()

            else:
                logGetError("Get wrong TYPE_QUANT: " + str(self.dp.myCFG.C0_TYPE_FLOW))


    def __captainCheckTime(self):

        dateNow = datetime.datetime.now()
        dateDead = datetime.datetime(TIME_EXPIRATION['Year'], TIME_EXPIRATION['Month'], TIME_EXPIRATION['Day'], 23, 59)

        n_day = (dateDead - dateNow).days

        if n_day < 0:

            logGetError(INFO_TO_USER_Staff[1])

        elif n_day < 7:

            logToUser(INFO_TO_USER_Staff[2])