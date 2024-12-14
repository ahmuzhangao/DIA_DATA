from MSLogging import logGetError
import bisect
import numpy
from MSSystem import VALUE_ILLEGAL

def soldierBinarySearch(inputList, start, end, number):
    #二分法搜索，在start到end的inputList里，找到number的索引
    # check
    if end >= len(inputList):
        logGetError("MSTool, MK311: The length of inputList is " + str(len(inputList)) + ", but you want to get " + str(end) + "?")

    if start < 0:
        logGetError(
            "MSTool, MK315: start is " + str(start) + "?")

    if number == inputList[end]:  # 最后一个有问题

        return end

    # find
    while start < end:
        mid = (start + end) // 2
        if number < inputList[mid]:
            end = mid
        elif number > inputList[mid]:
            start = mid + 1
        else:
            return mid

    start = start - 1

    return start

def toolFindNeighborListFromSortedList1(inputSortedList, numberList):

    flag = 0
    pointer = 0
    pointerNext = len(numberList) - 1
    list_border = [0, len(inputSortedList) - 1]
    list_result = [0] * len(numberList)

    for _ in range(len(numberList)):
        start, end = list_border
        number = numberList[pointer]

        if number <= inputSortedList[start]:
            neighbor = start

        elif number >= inputSortedList[end]:
            neighbor = end

        else:
            neighbor = bisect.bisect_left(inputSortedList, number, start, end)

            neighbor = neighbor if abs(inputSortedList[neighbor] - number) < \
                                   abs(inputSortedList[neighbor - 1] - number) else neighbor - 1

        list_border[flag] = neighbor
        list_result[pointer] = neighbor
        # 找下一个numberList的点
        next = pointer + 1 if pointer < pointerNext else pointer - 1
        pointer = pointerNext
        pointerNext = next

        flag = 0 if flag == 1 else 1

    return list_result

def toolGetIndexByWord(line, word, delimiter):
    # 找到word在line里，是以delimiter分开的第几个（即word在line里的索引）

    index = line.find(word)

    if -1 == index:

        logGetError('can not find '+word+' from: '+line)

    else:

        return toolCountCharInString(line[0:index], delimiter)


def toolGetWord(inputString, index, d):
    # 把字符串inputString按照d分隔开，取第index的索引，输出对应的字符串
    if inputString[0] != d:
        inputString = d + inputString

    if inputString[-1] != d:
        inputString = inputString + d

    p_d = []

    i = 0
    for c in inputString:

        if c == d:
            p_d.append(i)

        i = i + 1

    result = inputString[p_d[index] + 1:p_d[index + 1]]

    return result



def toolGetWord1(inputString, d1, d2):
    # 输入的字符串inputString，取d1和d2中间的部分

    start = 0
    end = len(inputString)

    for i in range(len(inputString)):

        if inputString[i] == d1:
            start = i + 1

        if inputString[i] == d2:
            end = i

    return inputString[start:end]

def toolCountCharInString(inputStr, inputChar):
    # 计算输入字符串inputStr，有几个字符串inputChar

    result = 0

    for c in inputStr:
        if c == inputChar:
            result = result + 1

    return result

def toolStr2List(inputStr, inputSeparator):
    #将字符串inputStr，以inputSeparator为间隔分开，转化为包含多个元素list

    outputList = []

    word = ''

    if inputStr[-1] != inputSeparator:
        inputStr = inputStr + inputSeparator

    for c in inputStr:

        if c == inputSeparator:

            number = float(word)
            outputList.append(number)
            word = ''

        else:

            word = word + c

    return outputList

def toolGetNameFromPath(path):
    #从一个文件路径中提取出去除后缀的文件名
    lenStr = len(path)
    iStart = 0
    iEnd = -1

    for i in range(lenStr):

        j = lenStr - 1 - i

        if path[j] == '.':
            iEnd = j  # 亲测必须这么写，不用减一

        if path[j] == '\\' or path[j] == '/':
            iStart = j + 1
            break

    if iStart == 0 and iEnd == -1:
        return path
    else:
        return path[iStart:iEnd]

def toolFindNeighborFromSortedList1(inputSortedList, number):

    #从inputSortedList里找到离number最近的一个数的索引

    start, end = 0, len(inputSortedList) - 1

    if number <= inputSortedList[start]:
        return start

    if number >= inputSortedList[end]:
        return end

    neighbor = bisect.bisect_left(inputSortedList, number)
    if abs(inputSortedList[neighbor] - number) < abs(inputSortedList[neighbor-1] - number):
        return neighbor
    else:
        return neighbor - 1

def toolFindIndexFromSortedList0(inputList, start, end, number):
    # 从inputList里找到number的索引

    if number < inputList[start]:
        return VALUE_ILLEGAL  # 没找到必须返回一个非法值

    if number > inputList[end]:
        return VALUE_ILLEGAL  # 没找到必须返回一个非法值

    index = soldierBinarySearch(inputList, start, end, number)

    if inputList[index] != number:

        return VALUE_ILLEGAL  # 没找到必须返回一个非法值

    else:

        return index

def toolFindIndexFromSortedList1(inputList, number):
    # 从inputList里找到number的索引
    return toolFindIndexFromSortedList0(inputList, 0, len(inputList)-1, number)

def toolSumList0(inputList, start, end):
    #把inputList从start到end加和
    result = 0.0

    for i in range(start, end+1):

        result = result + inputList[i]

    return result

def toolCosineForMatrixAndList(matrix_A, list_b):
    '''
    计算矩阵和列表之间的余弦相似度，返回的是一个相似度向量/列表
    :param matrix_A: 矩阵[M,N]
    :param list_b: 列表[1, N]
    返回shape=[M,]的相似度列表
    '''
    num = numpy.dot(matrix_A, list_b.T).reshape(-1,)
    denom = numpy.sqrt(numpy.multiply(matrix_A, matrix_A).sum(axis=1)) * numpy.sqrt(numpy.dot(list_b, list_b.T).sum(axis=1))

    zeroList = (denom == 0)
    denom[zeroList] = 1.
    cos = num / denom
    sim = 0.5 + 0.5 * cos
    sim[zeroList] = 0

    return sim


def toolCosineForListAndList(arr1, arr2):
    '''
    计算两个形状为 (1, n) 的数组之间的余弦相似度，并将结果映射到 [0, 1] 之间
    :param arr1: 数组1，形状为 (1, n)
    :param arr2: 数组2，形状为 (1, n)
    :return: 映射后的余弦相似度（标量）
    '''

    # 计算点积
    dot_product = numpy.dot(arr1, arr2.T)

    # 计算范数
    norm_arr1 = numpy.linalg.norm(arr1)
    norm_arr2 = numpy.linalg.norm(arr2)

    # 计算余弦相似度
    cosine_sim = dot_product / (norm_arr1 * norm_arr2)

    # 将余弦相似度映射到 [0, 1] 范围
    mapped_similarity = 0.5 + 0.5 * cosine_sim

    return mapped_similarity

def toolSumList1(inputList):
    #把inputList中所有值加和
    result = 0.0

    for num in inputList:

        result = result + num

    return result

def toolCopyList(list1):

    n = len(list1)
    list2 = []

    for i in range(n):
        list2.append(list1[i])

    return list2

def toolStringList2FloatList(inputList):

    result = []

    for word in inputList:
        result.append(float(word))

    return result

def toolPreprocessForModel0(inputList):
    #对数据进行预处理
    inputList = numpy.log2(inputList + 1)
    mean, std = inputList.mean(), inputList.std()
    inputList = (inputList - mean) / (std + 1e-8)

    return inputList

def toolCosineForList(list_a, list_b):

    vector_a = numpy.mat(list_a)
    vector_b = numpy.mat(list_b)
    num = float(vector_a * vector_b.T)
    denom = numpy.linalg.norm(vector_a) * numpy.linalg.norm(vector_b)

    if denom == 0:
        return 0
    else:
        cos = num / denom
        sim = 0.5 + 0.5 * cos
    return sim