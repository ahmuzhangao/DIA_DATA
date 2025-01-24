import matplotlib
matplotlib.use('TkAgg')
from MSStaff import CStaff
import time

# if __name__ == "__main__":
#     multiprocessing.freeze_support()
#     staff = CStaff(sys.argv)
#     staff.start()

# config_name = r'D:\Project\2025DIA_XLC\Data\Rerank\Main_cfg.txt'
config_name = r'D:\Project\2025DIA_XLC\Data\MSFraggerLib\Main_cfg.txt'
# config_name = r'D:\Project\2025DIA_XLC\Data\MSFraggerPin\Main_cfg.txt'
# config_name = r'D:\Project\2025DIA_XLC\Data\DIANN\Main_cfg.txt'
if __name__ == "__main__":
    staff = CStaff(['2', config_name])
    time_start = time.time()
    staff.start()
    time_end = time.time()