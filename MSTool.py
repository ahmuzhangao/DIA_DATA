import numpy as np

from MSLogging import logGetError
from MSSystem import VALUE_ILLEGAL


def toolGetIndexByWord(line, word, delimiter):

    index = line.find(word)

    if -1 == index:

        logGetError('can not find '+word+' from: '+line)

    else:

        return toolCountCharInString(line[0:index], delimiter)


def toolGetNextLine(input_str, input_index):

    line = ''

    len_str = len(input_str)

    i_str = input_index

    while i_str < len_str:

        tmpChar = input_str[i_str]
        i_str = i_str + 1

        if '\n' == tmpChar:
            break
        else:
            line = line + tmpChar

    return line


def toolGetNameFromPath(path):

    lenStr = len(path)
    iStart = 0
    iEnd = -1

    for i in range(lenStr):

        j = lenStr - 1 - i

        if path[j] == '.':
            iEnd = j  # 亲测必须这么写，不用减一
            break

    for i in range(lenStr):

        j = lenStr - 1 - i

        if path[j] == '\\' or path[j] == '/':
            iStart = j + 1
            break

    return path[iStart:iEnd]


def toolFindNeighborFromSortedList1_buzaishiyongban(input_list, start, end, num):
    midIndex = int((start + end) / 2)

    li = -1
    ri = -1

    if start >= end:
        return start

    mid = input_list[midIndex]

    if (midIndex - 1) < start:
        left = input_list[start]
    else:
        left = input_list[midIndex - 1]

    if (midIndex + 1) > end:
        right = input_list[end]
    else:
        right = input_list[midIndex + 1]

    sm = abs(num - mid)
    sl = abs(num - left)
    sr = abs(num - right)

    if sm < sl and sm < sr:
        return midIndex
    else:

        li = toolFindNeighborFromSortedList1(input_list, start, midIndex - 1, num)
        ri = toolFindNeighborFromSortedList1(input_list, midIndex + 1, end, num)

        if abs(num - input_list[li]) < abs(num - input_list[ri]):
            return li

        return ri


def toolFindNeighborFromSortedList2_buzaishiyongban(A, target):
    if len(A) == 0 or target < A[0]:
        return 0

    if target > A[len(A) - 1]:
        return len(A) - 1

    l = 0
    r = len(A) - 1
    while (l < r - 1):
        mid = int(l + (r - l) / 2)
        if A[mid] == target:
            return mid
        elif A[mid] > target:
            r = mid
        else:
            l = mid

    if abs(target - A[l]) < abs(target - A[r]):
        return l
    else:
        return r


# def toolFindNeighborFromSortedList0(input_list, data):
#
#     len_input_list = len(input_list)
#
#     distance = float("inf")
#     result = len_input_list + 1
#     isFirst = True
#
#     for i in range(len_input_list):
#
#         tmpDis = abs(data - input_list[i])
#
#         if tmpDis < distance:
#             distance = tmpDis
#             result = i
#         else:
#             if bool(1-isFirst):
#                 return result
#
#     return result


def toolMass2MZ(inputMH, input_charge):

    MASS_PROTON_MONO = 1.00727645224  # 1.00782503214-0.0005485799
    return (inputMH + (input_charge - 1) * MASS_PROTON_MONO) / input_charge


def toolCountInList(input_list, data):

    n = 0

    for e in input_list:

        if data == e:
            n = n + 1

    return n


def toolSumList0(inputList, start, end):

    result = 0.0

    for i in range(start, end+1):

        result = result + inputList[i]

    return result


def toolSumList1(inputList):

    result = 0.0

    for num in inputList:

        result = result + num

    return result


def toolFindFromListByKey(inputList, key):

    n = len(inputList)

    for i in range(len(inputList)):

        tmpELE = inputList[i]

        if -1 == tmpELE.find(key):
            result = -1
        else:
            result = i
            return result

    return n+1


def toolFindFromSortedList1_buzaishiyongban(list1, data):

    n = len(list1)
    first = 0
    last = n - 1
    while first <= last:
        mid = (last + first) // 2
        if list1[mid] > data:
            last = mid - 1
        elif list1[mid] < data:
            first = mid + 1
        else:
            return mid
    return n+1  # 不能返回-1


def toolCopyList(list1):

    n = len(list1)
    list2 = []

    for i in range(n):
        list2.append(list1[i])

    return list2


def toolGetWord1(inputString, d1, d2):

    start = 0
    end = len(inputString)

    for i in range(len(inputString)):

        if inputString[i] == d1:
            start = i + 1

        if inputString[i] == d2:
            end = i

    return inputString[start:end]

    # a = '0123456'
    # a[1:2] is '1'


def toolGetWord(inputString, index, d):

    if inputString[0] != d:
        inputString = d + inputString

    if inputString[-1] != d:
        inputString = inputString + d

    p_d = []
    
    i = 0
    for c in inputString:
        
        if c == d:
            
            p_d.append(i)
        
        i = i+1
        
    result = inputString[p_d[index]+1:p_d[index+1]]
        
    return result

def toolCountCharInString(inputStr, inputChar):

    result = 0

    for c in inputStr:
        if c == inputChar:
            result = result + 1

    return result


def toolList2Str(inputList, inputSeparator):

    result = ''

    for ele in inputList:

        result = result + str(ele) + inputSeparator

    return result


def toolStr2List(inputStr, inputSeparator):

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


def soldierBinarySearch(inputList, start, end, number):

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


def toolFindIndexFromSortedList0(inputList, start, end, number):

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

    return toolFindIndexFromSortedList0(inputList, 0, len(inputList)-1, number)


def toolFindNeighborFromSortedList0(inputList, start, end, number,):

    if number <= inputList[start]:
        return start

    if number >= inputList[end]:
        return end

    neighbor = soldierBinarySearch(inputList, start, end, number)  # 通常是落在了左边那个索引上

    disLeft = abs(inputList[neighbor] - number)
    disRight = abs(inputList[neighbor+1] - number)

    if disLeft < disRight:
        return neighbor
    else:
        return neighbor + 1


def toolFindNeighborFromSortedList1(inputList, number):

    return toolFindNeighborFromSortedList0(inputList, 0, len(inputList)-1, number)  # 长度必须-1，保证是个数字


def tooFindNeighborListFromSortedList(inputList, numList):

    list_array = np.array(inputList)
    c = np.searchsorted(list_array, numList)
    c[c > len(inputList) - 1] = len(inputList) - 1
    right_index = c
    left_index = c - 1
    left = np.abs(inputList[left_index] - numList)
    right = np.abs(inputList[right_index] - numList)
    min_ = left - right
    min_index = left_index * (min_ < 0) + right_index * (min_ >= 0)
    # 下面用列表计算时间大约是前两行的10倍，这里没有改变算法复杂度，只是用了numpy的广播操作来加速
    #  min_index = [left_index[i] if left[i] < right[i] else right_index[i] for i in range(len(right_index))]

    return min_index


def toolCosineForMatrix(_matrixA, _matrixB):

    # 按行求和，生成一个列向量
    # 即各行向量的模

    _matrixA_matrixB = _matrixA * _matrixB.transpose()
    _matrixA_norm = np.sqrt(np.multiply(_matrixA, _matrixA).sum(axis=1))
    _matrixB_norm = np.sqrt(np.multiply(_matrixB, _matrixB).sum(axis=1))
    return np.divide(_matrixA_matrixB, _matrixA_norm * _matrixB_norm.transpose())


def toolCosineForMatrixAndList(matrix_A, list_b, axis=1):
    '''
    计算矩阵和列表之间的余弦相似度，返回的是一个相似度向量/列表
    :param matrix_A: 矩阵[M,N]
    :param list_b: 列表(长度为L)
    :param axis: 计算相似度的维数，默认为1。1表示对列表和矩阵每一列计算相似度，0表示对列表和矩阵每一行计算相似度
           当axis=1时，list_b的维度为[M,1];
           当axis=0是，list_b的维度为[1,N];
    :return:
    '''
    if axis == 1:
        num = np.dot(list_b.T, matrix_A).reshape(-1,)
        denom = np.sqrt(np.multiply(matrix_A, matrix_A).sum(axis=0)) * np.sqrt(np.dot(list_b.T, list_b).sum(axis=0))
        if sum(denom != 0) < len(denom):
            cos = np.zeros(shape=(len(denom), ))
            cos[denom != 0] = num[denom != 0] / denom[denom != 0]
            sim = 0.5 + 0.5 * cos
        else:
            cos = num / denom
            sim = 0.5 + 0.5 * cos
    else:
        pass

    return sim


def toolCosineForList(list_a, list_b):

    vector_a = np.mat(list_a)
    vector_b = np.mat(list_b)
    num = float(vector_a * vector_b.T)
    denom = np.linalg.norm(vector_a) * np.linalg.norm(vector_b)

    if denom == 0:
        return 0
    else:
        cos = num / denom
        sim = 0.5 + 0.5 * cos
    return sim


def toolPearsonForNumpy(X, Y):

    Array_len = X.shape[0]
    X_mean = X.sum() / Array_len
    Y_mean = Y.sum() / Array_len
    pearson = np.dot(X - X_mean, Y - Y_mean) / (np.sqrt(np.dot(X - X_mean, X - X_mean) * np.dot(Y - Y_mean, Y - Y_mean)) + 1e-6)
    return pearson


def toolSpearmanForNumpy(X, Y):

    Array_len = X.shape[0]
    x_argsort = np.argsort(X)
    y_argsort = np.argsort(Y)

    fenzi = np.dot((x_argsort - y_argsort), (x_argsort - y_argsort).T)

    return 1 - (6 * fenzi) / (Array_len * (Array_len ** 2 - 1))


def toolMaxContinueBY(flagList: list):

    maxLen = 0
    tmpLen = 0
    for i in flagList:
        if i == 1:
            tmpLen += 1
            if tmpLen > maxLen:
                maxLen = tmpLen
        else:
            tmpLen = 0
    return maxLen


