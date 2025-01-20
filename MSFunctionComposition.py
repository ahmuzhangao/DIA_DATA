from MSTool import toolCountCharInString, toolGetWord
from MSLogging import logGetError, logToUser, INFO_TO_USER_FunctionComposition, logGetWarning

class CFunctionComposition:

    def __init__(self, inputDP):
        self.dp = inputDP

    def getStrComposition(self, inputSeq, inputMod,  inputGLC,inputLIK,inputINI):

        result = 'H(2)O(1)'  # 写死的

        for aa in inputSeq:

            if inputINI.DICT1_AA_COM.__contains__(aa):

                result = result + inputINI.DICT1_AA_COM[aa][0]

            else:

                continue

        # modification
        nMOD = toolCountCharInString(inputMod, ';')# 标准形式为6,Carbamidomethyl[C];14,Carbamidomethyl[C];
        nGLC = toolCountCharInString(inputGLC, ';')
        nLINK = toolCountCharInString(inputLIK, ';')
        try:
            for iMOD in range(nMOD):
                nameMOD = toolGetWord(toolGetWord(inputMod, iMOD, ';'), 1, ',')
                result = result + inputINI.DICT2_MOD_COM[nameMOD][0]

            for iGLC in range(nGLC):
                nameGLC = toolGetWord(toolGetWord(inputGLC, iGLC, ';'), 1, ',')
                if inputINI.DICT3_GLYCO_COM.__contains__(nameGLC):
                    result = result + inputINI.DICT3_GLYCO_COM[nameGLC]
                else:
                    logGetWarning(inputGLC)

            for iLINK in range(nLINK):
                nameLINK = toolGetWord(toolGetWord(inputLIK, iLINK, ';'), 1, ',')
                result = result + inputINI.DICT4_LINKER_COM[nameLINK]

        except:
            logToUser(INFO_TO_USER_FunctionComposition[0] + inputMod)
            logGetError(INFO_TO_USER_FunctionComposition[1])

        return result

    def getDictComposition(self, inputComposition):

        outputDictComposition = {}  # 这是一个以名字为key，以个数为number的dict

        LastIndexLeftBracket = 0
        LastIndexRightBracket = -1  # 第一开始，+1后从0开始

        lengthInputString = len(inputComposition)
        for i in range(lengthInputString):

            if inputComposition[i] == '(':
                eleName = inputComposition[LastIndexRightBracket + 1:i]
                LastIndexLeftBracket = i

            if inputComposition[i] == ')':

                eleAtomNum = int(inputComposition[LastIndexLeftBracket + 1:i])

                if outputDictComposition.__contains__(eleName):

                    outputDictComposition[eleName] = outputDictComposition[eleName] + eleAtomNum  # eleAtomNum可能是负数

                else:

                    outputDictComposition[eleName] = eleAtomNum

                LastIndexRightBracket = i

        return outputDictComposition