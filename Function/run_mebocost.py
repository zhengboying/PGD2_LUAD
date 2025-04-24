import os,sys
import pandas as pd
import pickle
from scipy.sparse import csc_matrix
# 添加 main.py 所在的目录到 sys.path
main_path = '/project_res/byzheng/biosoft/MEBOCOST'
sys.path.append(main_path)
from main import mebocost_piplines 
import argparse


def parameters():
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description='Process some mebocost parameters.')

    # 添加命令行参数
    parser.add_argument('--input_file', type=str, help='csv.gz or h5ad file, example in data/demo/raw_scRNA/demo_HNSC_200cell/')
    parser.add_argument('--input_type', type=str, help='expr_mat(must csv file), adata (must h5ad file),scEFA,compass')
    parser.add_argument('--cell_anno_csvfile', type=str, help='if input_type is expr_mat, must give a anno_csv as GSM7157674_cell_ann.csv.')
    parser.add_argument('--groups', type=str, help='a list, specify the column names in cell_ann for grouping cells, by default cell_type or cluster will be detected and used')
    parser.add_argument('--species', type=str, help='human or mouse, this determines which database will be used in our collection')
    parser.add_argument('--thread', type=int, help='Number of cores used for running job, default 8')
    parser.add_argument('--output_file', type=str, help='eg. output_file as./commu.pk') 
    return parser.parse_args()
    # 解析命令行参数
    mebocost_para = parser.parse_args()

    # print参数
    # 打印参数
    print(f"Input file: {mebocost_para.input_file}")
    print(f"Input type: {mebocost_para.input_type}")
    print(f"Cell annotation CSV file: {mebocost_para.cell_anno_csvfile}")
    print(f"Groups: {mebocost_para.groups}")
    print(f"Species: {mebocost_para.species}")
    print(f"Thread (number of cores): {mebocost_para.thread}")
    print(f"Output file: {mebocost_para.output_file}")

if __name__ == "__main__":
    mebocost_para = parameters()
    if mebocost_para.input_type == 'adata':
        print('The inputfile type is adata.') 
        mebo_obj = mebocost_piplines.run_mebocost_obj_from_adata(
                            input = mebocost_para.input_file,
                            groups = mebocost_para.groups,##提供的案例是celltype
                            species = mebocost_para.species,
                            config_path = '/project_res/byzheng/biosoft/MEBOCOST/mebocost.conf',
                            thread = mebocost_para.thread,
                            outputfile = mebocost_para.output_file
                            )   

    elif mebocost_para.input_type == 'expr_mat':
        print('The inputfile type is expr_mat.')
        mebo_obj = mebocost_piplines.run_mebocost_obj_from_expr_mat(
                            input = mebocost_para.input_file,
                            cell_ann = mebocost_para.cell_anno_csvfile,
                            groups = mebocost_para.groups,##other colname as groups in GSE229706_two_examples_cell_ann.csv
                            species = mebocost_para.species,
                            config_path = '/project_res/byzheng/biosoft/MEBOCOST/mebocost.conf',
                            thread = mebocost_para.thread,
                            outputfile = mebocost_para.output_file
                            )
    elif mebocost_para.input_type == 'scEFA':
        print('The inputfile type is scEFA output file.')
    elif mebocost_para.input_type == 'compass':
        print('The inputfile type is compass output file.')
    else :
        print('this inputfile type cannot be supported by MEBOCOST.')

##结果提取

with open('../03_mebocost/commu.pk', 'rb') as file:
    data = pickle.load(file)


#dict_keys(['exp_mat', 'adata', 'group_col', 'cell_ann', 'species', 'met_est', 'met_pred', 'config_path', 
#               'met_enzyme', 'met_sensor', 'met_ann', 'scFEA_ann', 'compass_met_ann', 'compass_rxn_ann', 'cutoff_exp', 
#               'cutoff_met', 'cutoff_prop', 'sensor_type', 'thread', 'exp_mat_indexer', 'exp_mat_columns', 
#               'commu_time_stamp', 'avg_exp', 'avg_exp_indexer', 'avg_exp_columns', 'met_mat', 'met_mat_indexer', 
#               'met_mat_columns', 'avg_met', 'avg_met_indexer', 'avg_met_columns', 'commu_res', 'exp_prop', 
#               'met_prop', 'original_result'])
#


#细胞代谢通讯结果
res = 'commu_res'      
df = data[res]
df.to_csv(f'{res}.csv', index=False)





##估算的代谢矩阵
res = 'met_mat'
df = data[res]
dense_df = df.todense()

# 然后使用pandas创建DataFrame
df_pd = pd.DataFrame(dense_df)

res = 'met_mat_indexer'      
indexer = data[res]
res = 'met_mat_columns'      
columns = data[res]
df_pd.index = indexer
df_pd.columns = columns
df_pd.to_csv('met_mat.csv', index=True)





'''
comm_res中包含
res = 'exp_prop'      
df = data[res]
df.to_csv(f'{res}.csv', index=True)


res = 'met_prop'      
df = data[res]
df.to_csv(f'{res}.csv', index=True)
'''


'''
原始结果，大部分不显著
res = 'original_result'      
df = data[res]
'''


#input_type = 'expr_mat'
#input_file = 'data/demo/raw_scRNA/GSM7157674_exp_mat.csv.gz'#'data/demo/raw_scRNA/demo_HNSC_200cell.h5ad'
#cell_anno_csvfile = 'data/demo/raw_scRNA/GSM7157674_cell_ann.csv' ##if input_type = 'adata',cell_anno_csvfile = None
#groups = 'cell_type'
#species = 'human'
#thread = 8
#output_file='./commu.pk'


   





