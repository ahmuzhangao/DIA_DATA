import time

from MSData import Config
from MSFunction import CFunctionConfig
from MSSystem import IO_NAME_FILE_CONFIG, IO_NAME_FILE_EXPORT_Flow2
from MSLogging import logToUser, INFO_TO_USER_Flow2, INFO_TO_USER_Flow3
from MSTask import CTaskCheck, CTaskReadINI, CTaskReadMS1, CTaskReadMS2, \
	CTaskDrawForDIACheck, CTaskReadIDForDIACheck, CTaskForDIARerankFeature, CTaskForDIARerank, CTaskForDIAReport


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


class CFlow3:

	def __init__(self, inputDP):

		self.dp = inputDP

	def run(self):

		logToUser(INFO_TO_USER_Flow3[0])
		taskCheck = CTaskCheck(self.dp)
		taskCheck.work(IO_NAME_FILE_EXPORT_Flow2)

		logToUser(INFO_TO_USER_Flow3[1])
		taskReadINI = CTaskReadINI(self.dp)
		taskReadINI.work()

		logToUser(INFO_TO_USER_Flow3[2])
		taskReadID = CTaskReadIDForDIACheck(self.dp)
		taskReadID.work()

		logToUser(INFO_TO_USER_Flow3[3])
		taskFileMS1 = CTaskReadMS1(self.dp)
		taskFileMS1.work()

		logToUser(INFO_TO_USER_Flow3[4])
		taskFileMS2 = CTaskReadMS2(self.dp)
		taskFileMS2.work()

		logToUser(INFO_TO_USER_Flow3[5])
		taskRerankFeature = CTaskForDIARerankFeature(self.dp)
		taskRerankFeature.work()

		logToUser(INFO_TO_USER_Flow3[6])
		taskRerank = CTaskForDIARerank(self.dp)
		taskRerank.work()

		logToUser(INFO_TO_USER_Flow3[7])
		taskReport = CTaskForDIAReport(self.dp)
		taskReport.work()


