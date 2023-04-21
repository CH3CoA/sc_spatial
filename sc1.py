# 载入包
import numpy as np
import pandas as pd
import scanpy as sc
import os





# 设置
sc.settings.verbosity = 3  # 设置日志等级: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(figsize=(8, 8), facecolor='white')
os.chdir('E:\\Python\\Spatial_Transcriptomics\\data\\spatial')  # 修改路径

# 用于存储分析结果文件的路径
results_file = 'spatial.h5ad'

# 载入文件
adata = sc.read_text(
    '\GSE197064_RAW\GSM5907965_121313_final_matrix.csv.gz'  # mtx 文件目录
    )  # 写入缓存 可以更快的读取文件
