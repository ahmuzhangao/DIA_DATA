from MSLogging import logGetError


class CEmass:

	EMASS_INTENSITY_MAX = 100.0
	EMASS_DUMMY_MASS = -10000000
	EMASS_DUMMY_PRO = 0.0  # 这个是用于填充S33等情况
	EMASS_CUTOFF_PRO = 0.00001  # 因为很多同位素的丰度就是0.01，两个相乘就是0.00001了。
	EMASS_CUTOFF_PRO_OUTPUT = 0.01  # 最后输出时，控制一下

	def __init__(self, inputDP):
		self.dp = inputDP

	def __soldierGenerateVerginList(self, inputElementName):

		# print("MSEmass55: "+str(inputElementID))

		SuperAtomList = []
		numberAtomType = len(self.dp.myINI.DICT0_ELEMENT_MASS[inputElementName])  # 如果元素不在element.ini中，就会报错

		for i in range(numberAtomType):  # 对于同一个元素中的多种同位素处理, 变成[丰度, 质量]的list

			tmpMass = self.dp.myINI.DICT0_ELEMENT_MASS[inputElementName][i]
			tmpPro = self.dp.myINI.DICT0_ELEMENT_ABDC[inputElementName][i]
			tmpAtomNode = [tmpPro, tmpMass]  # 先丰度再质量

			SuperAtomList.append(tmpAtomNode)  # 有几个同位素就有几个list元素

		return SuperAtomList

	def __soldierGetSuperAtomList(self, dictComposition):  # 这个函数很关键，是java中的admiralGetSuperAtomList

		# init
		salForResult = []
		salForResult.append([1.0, 0.0])

		for name in dictComposition:  # 对组成中的每一个元素处理

			numAtom = dictComposition[name]  # 该元素的数量

			if numAtom > 0:  #  这句话很重要，因为Gln->pyro-Glu[AnyN-termQ]这个修饰，会产生N(-1)。如果是是14N(14)N(-1)，N(-1)是不起作用的
				salForDouble = self.__soldierGenerateVerginList(name)  # element.ini中该元素的同位素和丰度列表

				while numAtom > 0:

					if 1 == numAtom % 2:  # 如果不是偶数 ???
						salForResult = self.__captainMergerSAL(salForResult, salForDouble)
						salForResult = self.__soldierPrune(salForResult, self.EMASS_CUTOFF_PRO)

					salForDouble = self.__captainMergerSAL(salForDouble, salForDouble)
					salForDouble = self.__soldierPrune(salForDouble, self.EMASS_CUTOFF_PRO)  # 这句貌似很重要，不然传进去的都是很大的数组，运算奇慢无比；

					numAtom = int(numAtom / 2)

		salForResult = self.__soldierPrune(salForResult, self.EMASS_CUTOFF_PRO_OUTPUT)  # EMASS_CUTOFF_PRO太小了，同位素会很多
		return salForResult

	def __soldierPrune(self, inputListNode, inputLimit):

		outputList = []

		for i in range(len(inputListNode)):

			tmpNode = inputListNode[i]

			if tmpNode[0] > inputLimit:
				outputList.append(tmpNode)

		return outputList

	def __captainMergerSAL(self, f, g):  # 之所以取名为f、g和h，是为了和论文中公式统一

		h = []

		length_f = len(f)
		length_g = len(g)

		if 0 == length_f or 0 == length_g:

			return h

		for k in range(length_f+length_g):

			sumweight = 0.0
			summass = 0.0

			start = 0
			if k < (length_f - 1):
				start = 0
			else:
				start = k - length_f + 1

			end = k
			if k < length_g - 1:
				end = k
			else:
				end = length_g - 1

			# time20 = time.perf_counter()  # 测试用
			for i in range(start, end+1):

				weight = g[i][0] * f[k-i][0]
				mass = g[i][1] + f[k-i][1]

				sumweight = sumweight + weight
				summass = summass + weight * mass

			if sumweight == 0:
				newMass = self.EMASS_DUMMY_MASS
			else:
				newMass = summass / sumweight

			newPro = sumweight
			newAtomNode = [newPro, newMass]

			h.append(newAtomNode)

		return h

	def __soldierMultiListAndNumber(self, inputList, inputNum):

		for i in range(len(inputList)):
			inputList[i] = inputList[i] * inputNum

	def getCalculatedIsotopicPeaks(self, inputDictComposition, inputCharge):

		sal = self.__soldierGetSuperAtomList(inputDictComposition)  # sal super atom list

		nPeaks = len(sal)

		result = [[0]*nPeaks, [0]*nPeaks]  # 两位数组，一行有nPeaks个数

		maxPro = 0.0
		i = 0

		for tmpNode in sal:

			result[0][i] = tmpNode[0]

			if tmpNode[0] > maxPro:
				maxPro = tmpNode[0]

			if inputCharge > 0:

				result[1][i] = (tmpNode[1] + inputCharge * self.dp.myINI.MASS_PROTON_MONO) / inputCharge

			else:

				logGetError("MSEmass, MK172: " + str(result))

			i = i + 1

		if maxPro == 0:
			logGetError("MSEmass, MK177: " + str(inputDictComposition))
		else:
			factor = self.EMASS_INTENSITY_MAX / maxPro
			self.__soldierMultiListAndNumber(result[0], factor)

		return result

	def getCalculatedMonoMZ(self, inputDictCom, inputCharge):

		mass = 0.0

		for name in inputDictCom:

			tmpNum = inputDictCom[name]
			list_abdc = self.dp.myINI.DICT0_ELEMENT_ABDC[name]
			list_mass = self.dp.myINI.DICT0_ELEMENT_MASS[name]
			tmpMass = list_mass[list_abdc.index(max(list_abdc))]

			mass += tmpMass * tmpNum

		return (mass + self.dp.myINI.MASS_PROTON_MONO * inputCharge) / inputCharge

