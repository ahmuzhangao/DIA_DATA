from MSData import Config
from MSFunction import CFunctionConfig
from MSSystem import IO_NAME_FILE_CONFIG, IO_NAME_FILE_EXPORT_Flow2
from MSLogging import logToUser, INFO_TO_USER_Flow2
from MSTask import CTaskCheck, CTaskReadINI, CTaskReadMS1, CTaskReadMS2, \
	CTaskDrawForDIACheck, CTaskReadIDForDIACheck

class CFlow0:

	def run(self):

		config = Config()

		functionConfig = CFunctionConfig()
		functionConfig.config2file(IO_NAME_FILE_CONFIG[0], config)


class CFlow2:

	def __init__(self, inputDP):

		self.dp = inputDP

	def run(self):

		logToUser(INFO_TO_USER_Flow2[0])
		taskCheck = CTaskCheck(self.dp)
		taskCheck.work(IO_NAME_FILE_EXPORT_Flow2)

		logToUser(INFO_TO_USER_Flow2[1])
		taskReadINI = CTaskReadINI(self.dp)
		taskReadINI.work()

		logToUser(INFO_TO_USER_Flow2[2])
		taskReadID = CTaskReadIDForDIACheck(self.dp)
		taskReadID.work()

		logToUser(INFO_TO_USER_Flow2[3])
		taskFileMS1 = CTaskReadMS1(self.dp)
		taskFileMS1.work()

		logToUser(INFO_TO_USER_Flow2[4])
		taskFileMS2 = CTaskReadMS2(self.dp)
		taskFileMS2.work()

		logToUser(INFO_TO_USER_Flow2[5])
		taskQuant = CTaskDrawForDIACheck(self.dp)
		taskQuant.work()

