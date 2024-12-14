import logging
import sys
from MSSystem import SOFTWARE_NAME

INFO_TO_USER_CTaskCreateLibrary = (

    '\n[{}] \tLibrary PSMs: '.format(SOFTWARE_NAME),

)

INFO_TO_USER_Staff = (

    '\n[{}] Copyright 2024 Beihang University by Andrew Gu. All rights reserved. Version 2024.4.'.format(SOFTWARE_NAME),
    '\n[{}] LICENSE is invalid! Please send e-mail to pQuant-DIA@126.com for the valid version.'.format(SOFTWARE_NAME),
    '\n[{}] Warning! The current license will expired in 7 days. Please send e-mail to pQuant-DIA@126.com for the new version.'.format(SOFTWARE_NAME),
    '\n[{}] Writing config file in the folder...'.format(SOFTWARE_NAME),
    '\n[{}] Finished!'.format(SOFTWARE_NAME),

)

INFO_TO_USER_Flow1 = (

    '\n[{}] <Flow DIA Label Free> Checking the environment...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Label Free> Reading ini files...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Label Free> Reading identification files...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Label Free> Reading MS1 files (The first run will take a few minutes)...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Label Free> Reading MS2 files (The first run will take a few minutes)...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Label Free> Quantifying...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Label Free> Infering...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Label Free> RT Calibrated Information from ApuRT has been read...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Label Free> Predicting MS2 from identifications...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Label Free> Creating dynamic evidence library from IDs...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Label Free> Reading MS files (The first run will take a few minutes)...'.format(SOFTWARE_NAME),
)

INFO_TO_USER_Flow2 = (

    '\n[{}] <Flow DIA Evidence Check> Checking the environment...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Reading ini files...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Reading identification files for evidence check...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Reading MS1 files (The first run will take a few minutes)...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Reading MS2 files (The first run will take a few minutes)...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Drawing DIA Curve...'.format(SOFTWARE_NAME),
)

INFO_TO_USER_Flow3 = (

    '\n[{}] <Flow DIA Evidence Check> Checking the environment...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Reading ini files...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Reading identification files for generating training Data...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Classifying PSM according to precursor ID...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Reading MS1 files (The first run will take a few minutes)...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Reading MS2 files (The first run will take a few minutes)...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Generating training Data...'.format(SOFTWARE_NAME),
    '\n[{}] <Flow DIA Evidence Check> Predicting MS2 from identifications...'.format(SOFTWARE_NAME),
)

INFO_TO_USER_TaskReadID = (

    '\n[{}]\tReading '.format(SOFTWARE_NAME),
    '\n[{}]\tTotal number of IDs: '.format(SOFTWARE_NAME),

)

INFO_TO_USER_TaskReadMS1 = (

    '\n[{}] \tCreating MS1 index file for: '.format(SOFTWARE_NAME),

)

INFO_TO_USER_TaskReadMS2 = (

    '\n[{}] \tCreating MS2 index file for: '.format(SOFTWARE_NAME),

)

INFO_TO_USER_TaskReadRAW = (

    '\n[{}] \tCreating RAW index file for: '.format(SOFTWARE_NAME),

)

myLogPath = '{}.log'.format(SOFTWARE_NAME)

INFO_TO_USER_FunctionQuant = (

    '\n[{}] <Function Quantitation> Getting evidences for references...'.format(SOFTWARE_NAME),
    '\n[{}] <Function Quantitation> Getting evidences for samples...'.format(SOFTWARE_NAME),
)

INFO_TO_USER_FunctionQuant_FOR_DECOY = (

    '\n[{}] <Function Quantitation> Getting decoy evidences for references...'.format(SOFTWARE_NAME),
    '\n[{}] <Function Quantitation> Getting decoy evidences for samples...'.format(SOFTWARE_NAME),
)

INFO_TO_USER_FunctionComposition = (

    '\n[{}] <Function Composition> Wrong format: '.format(SOFTWARE_NAME),
    '\n[{}] <Function Composition> Make sure the format of this modifications is correct!'.format(SOFTWARE_NAME),
    '\n[{}] <Function Composition> Make sure the format of this glyco is correct!'.format(SOFTWARE_NAME),
    '\n[{}] <Function Composition> Make sure the format of this linker is correct!'.format(SOFTWARE_NAME),
)


def logGetWarning(strInfo):
    print(strInfo)

    logging.basicConfig(filename=myLogPath,
                        filemode='a',
                        format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                        )
    logging.warning(strInfo)

def logToUser(strInfo):

    try:
        print(strInfo)
        f_w = open(myLogPath, 'a', encoding='utf8')
        f_w.write(strInfo + '\n')
        f_w.close()

    except IOError:
        print("{}.log is opened! Please close it and run the program again!".format(SOFTWARE_NAME))
        sys.exit(0)


def logGetError(strInfo):
    print(strInfo)

    logging.basicConfig(filename=myLogPath,
                        filemode='a',
                        format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                        )
    logging.error(strInfo)
    sys.exit(0)


INFO_TO_USER_CTaskCreateSamplePairs = (

    '\n[{}] \t Number of Positive sample pairs: '.format(SOFTWARE_NAME),
    '\n[{}] \t Number of Negative sample pairs: '.format(SOFTWARE_NAME),
    '\n[{}] \t Getting evidences from ID files ... '.format(SOFTWARE_NAME),
    '\n[{}] \t Saving Training Pairs from evidences ... '.format(SOFTWARE_NAME),

)

INFO_TO_USER_FunctionFDR = (

    '\n[{}] <Function FDR> Calculate FQR thresholds and filter quantitative results...'.format(SOFTWARE_NAME),
    '\n[{}] <Function FDR> Too few results of match between runs, do not perform FQR filter...'.format(SOFTWARE_NAME),
)