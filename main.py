import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mebocost import mebocost
import h5py


# ================================
# @auther: BYZheng
# @email: boying.zheng@repugene.com
# @date: Oct 2024
# ================================



class mebocost_piplines:
    # all Params in MEBOCOST projects cited from ./ src/metocost/mebocost.py
    """
    MEBOCOST for predicting metabolite-based cell-cell communication. The modules of the package include communication inference and visualization.

    Params
    -------
    exp_mat
        python pandas data frame, single cell expression matrix, rows are genes, columns are cells
        'exp_mat' is a exclusive parameter to 'adata'
    adata
        scanpy adata object, the expression will be extracted, 'adata' is an exclusive parameter to 'exp_mat'
    cell_ann
        data frame, cell annotation information, cells in row names
    group_col
        a list, specify the column names in 'cell_ann' for grouping cells, by default 'cell_type' or 'cluster' will be detected and used
    species
        human or mouse, this determines which database will be used in our collection

    met_est
        the method for estimating metabolite level in cell, should be one of:
        mebocost: estimated by the enzyme network related to the metabolite
        scFEA-flux: flux result of published software scFEA (https://pubmed.ncbi.nlm.nih.gov/34301623/)
        scFEA-balance: balance result of published software scFEA (https://pubmed.ncbi.nlm.nih.gov/34301623/)
        compass-reaction: reaction result of published software Compass (https://pubmed.ncbi.nlm.nih.gov/34216539/)
        compass-uptake: uptake result of published software Compass (https://pubmed.ncbi.nlm.nih.gov/34216539/)
        compass-secretion: secretion result of published software Compass (https://pubmed.ncbi.nlm.nih.gov/34216539/)
    met_pred
        data frame, if scFEA or Compass is used to impute the metabolite level in cells, please provide the original result from scFEA or Compass, cells in row names, metabolite/reaction/module in column names, 
        Noted that this parameter will be ignored if 'met_est' was set as mebocost.

    config_path
        str, the path for a config file containing the path of files for metabolite annotation, enzyme, sensor, scFEA annotation, compass annotation. These can also be specified separately by paramters as following:

        if config_path not given, please set:
    met_enzyme from ./mebocost.conf
        data frame, metabolite and gene (enzyme) relationships, required columns include HMDB_ID, gene, direction, for instance:
        
        HMDB_ID     gene                                                direction
        HMDB0003375 Cyp2c54[Unknown]; Cyp2c38[Unknown]; Cyp2c50[Un...   substrate
        HMDB0003375 Cyp2c54[Unknown]; Cyp2c38[Unknown]; Cyp2c50[Un...   substrate
        HMDB0003375 Cyp2c54[Unknown]; Cyp2c38[Unknown]; Cyp2c50[Un...   substrate
        HMDB0003450 Cyp2c54[Unknown]; Cyp2c38[Unknown]; Cyp2c50[Un...   product
        HMDB0003948 Tuba8[Unknown]; Ehhadh[Unknown]; Echs1[Enzyme]...   product

    met_sensor ./mebocost.conf
        data frame, metabolite sensor information, each row is a pair of metabolite and sensor, must include columns  HMDB_ID, Gene_name, Annotation, for instance:
        
        HMDB_ID Gene_name   Annotation
        HMDB0006247 Abca1   Transporter
        HMDB0000517 Slc7a1  Transporter
        HMDB0000030 Slc5a6  Transporter
        HMDB0000067 Cd36    Transporter
        
    met_ann 
        ./mebocost.conf：hmdb_info_path
        data frame, the annotation of metabolite collected from HMDB website, these are basic annotation info including HMDB_ID, Kegg_ID, metabolite, etc

    scFEA_ann
        data frame, module annotation of metabolite flux in scFEA, usually is the file at https://github.com/changwn/scFEA/blob/master/data/Human_M168_information.symbols.csv

    compass_met_ann
        data frame, the metabolite annotation used in Compass software, usually is the file at https://github.com/YosefLab/Compass/blob/master/compass/Resources/Recon2_export/met_md.csv

    compass_rxn_ann
        data frame, the reaction annotation used in Compass software, usually is the file at https://github.com/YosefLab/Compass/blob/master/compass/Resources/Recon2_export/rxn_md.csv

    cutoff_exp
        auto or float, used to filter out cells which are lowly expressed for the given gene, by default is auto, meaning that automatically decide cutoffs for sensor expression to exclude the lowly 25% non-zeros across all sensor or metabolites in all cells in addition to zeros 

    cutoff_met
        auto or float, used to filter out cells which are lowly abundant of the given metabolite, by default is auto, meaning that automatically decide cutoffs for metabolite aggregated enzyme to exclude the lowly 25% non-zeros across all sensor or metabolites in all cells in addition to zeros 

    cutoff_prop
        float from 0 to 1, used to filter out metabolite or genes if the proportion of their abundant cells less than the cutoff

    sensor_type
        a list, provide a list of sensor type that will be used in the communication modeling, must be one or more from ['Receptor', 'Transporter', 'Nuclear Receptor'], default is all the three

    thread
        int, number of cores used for running job, default 1
        
    """

    def run_mebocost_obj_from_adata(
                            input='data/demo/raw_scRNA/demo_HNSC_200cell.h5ad',
                            groups='celltype',
                            species='human',
                            config_path = '/project_res/byzheng/biosoft/MEBOCOST/mebocost.conf',
                            thread=8,
                            outputfile="./commu.pk"
                            ):
        adata = sc.read_h5ad(input)
        ## check adata (cells, genes)
        print(adata.shape)
        ### make sure your adata cover all genes and all cells
        ### if not, you may want to load all genes saved in raw by doing:
        # adata = adata.raw.copy()
        ##adata.var 基因名列表
        ##adata.obs 等价于annotation.data
        ##adata.X 表达稀疏矩阵.行为cell,列为基因
        ## initiate the mebocost object
        ### import expression data by scanpy adata object
        mebo_obj = mebocost.create_obj(
                        adata = adata,
                        met_est = 'mebocost',
                        sensor_type=['Receptor', 'Transporter', 'Nuclear Receptor'],##或者选择一部分
                        group_col = [groups], #a list, specify the column names in 'cell_ann' for grouping cells, by default 'cell_type' or 'cluster' will be detected and used
                        species=species,
                        config_path = config_path,
                        #参考数据文件
                        met_enzyme=None, #默认./mebocost.conf提供
                        met_sensor=None, #默认./mebocost.conf提供
                        met_ann=None, #hmdb_info_path
                        ##阈值设置
                        cutoff_exp=0, ## automated cutoff to exclude lowly ranked 25% sensors across all cells
                        cutoff_met=0, ## automated cutoff to exclude lowly ranked 25% metabolites across all cells
                        cutoff_prop=0.25, ## at lease 25% of cells should be expressed the sensor or present the metabolite in the cell group (specified by group_col)
                        thread=thread,

                        ##默认即可，很多参数在以adata为输入时，不执行
                        exp_mat=None,
                        cell_ann=None,##配合exp_mat使用，提供cell注释
                        met_pred=None,#data frame, if scFEA or Compass is used to impute the metabolite level in cells
                        scFEA_ann=None,#scfea_info_path
                        compass_met_ann=None, #compass_rxt_ann_path
                        compass_rxn_ann=None #compass_met_ann_path
                        )
        print('config file path in the object:', mebo_obj.config_path)         
        ##communication inference
        commu_res = mebo_obj.infer_commu(
                                n_shuffle=1000,
                                seed=12345, 
                                Return=True, 
                                thread=None,
                                save_permuation=False,
                                min_cell_number = 1
                            )
        
        print('Number of mCCC(Metabolic Cell-Cell Communication) detected by enzyme and sensor co-expression: ', commu_res.shape[0])
        ## check cutoff
        print('sensor_exp cutoff: %s'%mebo_obj.cutoff_exp)
        print('metabolite_agg_enzyme cutoff: %s'%mebo_obj.cutoff_met)
        print('saving the meto_obj to outputfile:',outputfile,'...')
        mebocost.save_obj(obj = mebo_obj, path = outputfile)
        return mebo_obj


    def run_mebocost_obj_from_expr_mat(
                            input = 'data/demo/raw_scRNA/GSM7157674_exp_mat.csv.gz',
                            cell_ann ='data/demo/raw_scRNA/GSM7157674_cell_ann.csv',
                            groups = 'cell_type',##other colname as groups in GSE229706_two_examples_cell_ann.csv
                            species = 'human',
                            config_path = '/project_res/byzheng/biosoft/MEBOCOST/mebocost.conf',
                            thread = 8,
                            outputfile="./commu.pk"
                            ):
        expr_mat_df = pd.read_csv(input, index_col = 0)
        expr_mat_df = expr_mat_df.T
        cell_ann_df = pd.read_csv(cell_ann, index_col = 0)
        ## check expr_mat (cells, genes)
        print(expr_mat_df.shape)
        print(cell_ann_df)
        ### make sure your exp_mat cover all genes and all cells
        ### if not, you may want to load all genes saved in raw by doing:
        ## initiate the mebocost object
        ### import expression data by exp_mat and cell_ann
        mebo_obj = mebocost.create_obj(
                        exp_mat = expr_mat_df,
                        cell_ann = cell_ann_df,##配合exp_mat使用，提供cell注释

                        met_est = 'mebocost',
                        sensor_type = ['Receptor', 'Transporter', 'Nuclear Receptor'],##或者选择一部分
                        group_col = [groups], #a list, specify the column names in 'cell_ann' for grouping cells, by default 'cell_type' or 'cluster' will be detected and used
                        species = species,
                        config_path = config_path,
                        #参考数据文件
                        met_enzyme=None, #默认./mebocost.conf提供
                        met_sensor=None, #默认./mebocost.conf提供
                        met_ann=None, #hmdb_info_path
                        ##阈值设置
                        cutoff_exp=0, ## automated cutoff to exclude lowly ranked 25% sensors across all cells
                        cutoff_met=0, ## automated cutoff to exclude lowly ranked 25% metabolites across all cells
                        cutoff_prop=0.01, ## at lease 1% of cells should be expressed the sensor or present the metabolite in the cell group (specified by group_col)
                        thread=thread,

                        ##默认即可，很多参数在以adata为输入时，不执行
                        adata = None,
                        met_pred=None,#data frame, if scFEA or Compass is used to impute the metabolite level in cells
                        scFEA_ann=None,#scfea_info_path
                        compass_met_ann=None, #compass_rxt_ann_path
                        compass_rxn_ann=None #compass_met_ann_path
                        )
        print('config file path in the object:', mebo_obj.config_path)         
        ##communication inference
        commu_res = mebo_obj.infer_commu(
                                n_shuffle=1000,
                                seed=12345, 
                                Return=True, 
                                thread=None,
                                save_permuation=False,
                                min_cell_number = 1,
                                pval_method='permutation_test_fdr', 
                                pval_cutoff = 0.9
                            )
        
        print('Number of mCCC(Metabolic Cell-Cell Communication) detected by enzyme and sensor co-expression: ', commu_res.shape[0])
        ## check cutoff
        print('sensor_exp cutoff: %s'%mebo_obj.cutoff_exp)
        print('metabolite_agg_enzyme cutoff: %s'%mebo_obj.cutoff_met)
        print('saving the meto_obj to outputfile:',outputfile,'...')
        mebocost.save_obj(obj = mebo_obj, path = outputfile)
        return mebo_obj
    


  