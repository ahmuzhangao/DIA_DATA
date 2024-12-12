from MSData import Config
from MSSystem import IO_NAME_FILE_CONFIG, IO_NAME_FILE_EXPORT_Flow1
from MSFunction import CFunctionConfig
from MSLogging import INFO_TO_USER_Flow1,logToUser
from MSTask import CTaskCheck, CTaskReadINI, CTaskReadIDForDraw, CTaskReadMS1, CTaskDrawDIAXIC, CTaskReadMS2
from MSData import CDataPack



class CFlow0:

    def __init__(self):
        pass

    def run(self):

        config = Config()
        functionConfig = CFunctionConfig()
        functionConfig.config2file(IO_NAME_FILE_CONFIG[0], config)


class CFlow1:

    def __init__(self, inputDP: CDataPack):
        self.dp = inputDP

    def run(self):
        logToUser(INFO_TO_USER_Flow1[0])
        taskCheck = CTaskCheck(self.dp)
        taskCheck.work(IO_NAME_FILE_EXPORT_Flow1)

        logToUser(INFO_TO_USER_Flow1[1])
        taskReadINI = CTaskReadINI(self.dp)
        taskReadINI.work()

        logToUser(INFO_TO_USER_Flow1[2])
        taskReadID = CTaskReadIDForDraw(self.dp)
        taskReadID.work()

        logToUser(INFO_TO_USER_Flow1[3])
        taskReadMS1 = CTaskReadMS1(self.dp)
        taskReadMS1.work()

        logToUser(INFO_TO_USER_Flow1[4])
        taskReadMS2 = CTaskReadMS2(self.dp)
        taskReadMS2.work()

        logToUser(INFO_TO_USER_Flow1[5])
        taskDrawPrecursor = CTaskDrawDIAXIC(self.dp)
        taskDrawPrecursor.work()

        logToUser(INFO_TO_USER_Flow1[6])

