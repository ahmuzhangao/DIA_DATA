

import logging
import sys


INFO_TO_USER_Staff = (

    '\n[2024DIA_XLC]\tWriting config file in the folder...',
    '\n[2024DIA_XLC]\tFinished!'
)

INFO_TO_USER_Flow1 = (

    '\n[2024DIA_XLC]\tChecking the environment...',
    '\n[2024DIA_XLC]\tReading ini files...',
    '\n[2024DIA_XLC]\tReading identification files...',
    '\n[2024DIA_XLC]\tReading MS1 files (The first run will take a few minutes)...',
    '\n[2024DIA_XLC]\tReading MS2 files (The first run will take a few minutes)...',
    '\n[2024DIA_XLC]\tDrawing Curve...',
    '\n[2024DIA_XLC]\tEnd'
)

INFO_TO_USER_Flow2 = (
    '\n[2024DIA_XLC]\tChecking the environment...',
    '\n[2024DIA_XLC]\tReading ini files...',
    '\n[2024DIA_XLC]\tReading identification files...',
    '\n[2024DIA_XLC]\tReading MS1 files (The first run will take a few minutes)...',
    '\n[2024DIA_XLC]\tReading MS2 files (The first run will take a few minutes)...',
    '\n[2024DIA_XLC]\tDrawing Fragment Curve...',
    '\n[2024DIA_XLC]\tEnd'
)

INFO_TO_USER_Flow3 = (
    '\n[2024DIA_XLC]\tChecking the environment...',
    '\n[2024DIA_XLC]\tReading ini files...',
    '\n[2024DIA_XLC]\tReading identification files...',
    '\n[2024DIA_XLC]\tReading MS1 files (The first run will take a few minutes)...',
    '\n[2024DIA_XLC]\tDrawing DDA Precursor Curve...',
    '\n[2024DIA_XLC]\tEnd'
)

INFO_TO_USER_Flow4 = (
    '\n[2024DIA_XLC]\tChecking the environment...',
    '\n[2024DIA_XLC]\tReading ini files...',
    '\n[2024DIA_XLC]\tReading MS1 files (The first run will take a few minutes)...',
    '\n[2024DIA_XLC]\tReading MS2 files (The first run will take a few minutes)...',
    '\n[2024DIA_XLC]\tReading identification files...',
    '\n[2024DIA_XLC]\tExtracting Feature from ID...',
    '\n[2024DIA_XLC]\tReranking...',
    '\n[2024DIA_XLC]\tWriting Result...'
)

INFO_TO_USER_Flow5 = (

    '\n[2024DIA_XLC]\tChecking the environment...',
    '\n[2024DIA_XLC]\tReading ini files...',
    '\n[2024DIA_XLC]\tReading identification files...',
    '\n[2024DIA_XLC]\tReading MS1 files (The first run will take a few minutes)...',
    '\n[2024DIA_XLC]\tReading MS2 files (The first run will take a few minutes)...',
    '\n[2024DIA_XLC]\tDrawing DIA Curve And PSM labeling...',
    '\n[2024DIA_XLC]\tEnd'
)

INFO_TO_USER_FLOW6 = (

    '\n[2024DIA_XLC]\tChecking the environment...',
    '\n[2024DIA_XLC]\tReading ini files...',
    '\n[2024DIA_XLC]\tReading identification files...',
    '\n[2024DIA_XLC]\tReading MS1 files (The first run will take a few minutes)...',
    '\n[2024DIA_XLC]\tReading MS2 files (The first run will take a few minutes)...',
    '\n[2024DIA_XLC]\tDrawing MS1 spectra for ID...',
    '\n[2024DIA_XLC]\tDrawing MS2 spectra for ID...'
    '\n[2024DIA_XLC]\tEnd'
)

INFO_TO_USER_TaskReadMS1 = (

    '\n[2024DIA_XLC]\tCreating index file for: ',

)

INFO_TO_USER_TaskReadMS2 = (

    '\n[2024DIA_XLC]\tCreating index file for: ',

)

INFO_TO_USER_TaskReadID = (

    '\n[2024DIA_XLC]\tReading ',
    '\n[2024DIA_XLC]\tTotal number of IDs: ',

)


INFO_TO_USER_FunctionQuant = (

    '\n[2024DIA_XLC] <Function Quantitation> Getting evidences for references...',
    '\n[2024DIA_XLC] <Function Quantitation> Getting evidences for samples...',
)

myLogPath = '2024DIA_XLC.log'


# 也许向某个日志文件里写东西
def logToUser(strInfo):

    # if os.access(myLogPath, os.W_OK):  # 当文件被excel打开时，这个东东没啥用

    try:
        print(strInfo)
        f_w = open(myLogPath, 'a', encoding='utf8')
        f_w.write(strInfo + '\n')
        f_w.close()
    except IOError:
        print("2024DIA_XLC.log is opened! Please close it and run the program again!")
        sys.exit(0)


def logGetError(strInfo):

    print(strInfo)
    
    logging.basicConfig(filename=myLogPath,
                        filemode='a',
                        format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                        )
    logging.error(strInfo)
    sys.exit(0)
    

def logGetWarning(strInfo):
    
    print(strInfo)
    
    logging.basicConfig(filename=myLogPath,
                        filemode='a',
                        format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                        )
    logging.warning(strInfo)





