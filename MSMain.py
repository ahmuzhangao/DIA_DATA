import sys
import time
from MSStaff import CStaff


# 打包使用
# if __name__ == "__main__":
#     staff = CStaff(sys.argv)
#     staff.start()
# 
config_name = r'D:\Project\2024DIA_XLC\data\DIANN\DIA_XLC_cfg_DIANN.txt'
# config_name = r'D:\Project\2024DIA_XLC\data\DP\DIA_XLC_cfg_new.txt'
# config_name = r'D:\Project\2024DIA_XLC\DIA_XLC_cfg_msfagger.txt'
# config_name = r'D:\Project\2024DIA_XLC\DIA_XLC_cfg_new.txt'
# config_name = r'D:\Project\2024DIA_XLC\DIA_XLC_cfg_DIANN.txt'

if __name__ == '__main__':
    staff = CStaff(['2', config_name])
    time_start = time.time()
    staff.start()
    time_end = time.time()
    print('all time', time_end - time_start)
