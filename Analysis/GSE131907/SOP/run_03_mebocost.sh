#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o ../03_mebocost/run_03_mebocost_expr_mat.log
#$ -e ../03_mebocost/run_03_mebocost_expr_mat.error.log
# 指定资源要求
#$ -l mem_free=50G   # 确保有至少10GB的空闲内存
#$ -l h_vmem=50G     # 确保有至少4GB的虚拟内存
#$ -l h_rt=24:00:00 # 作业运行时间限制为24小时
#$ -q new.q  # 集群中包含三个队列，new.q old.q test.q
#$ -pe smp 16




## 结果地址
project_scripts='/project_res/byzheng/project/DMPJ0001-PGD2-byzheng/Function'
workdir='/project_res/byzheng/project/DMPJ0001-PGD2-byzheng/Analysis/GSE131907/'
cd $workdir/03_mebocost




##激活指定环境
source /project_res/byzheng/miniconda3/bin/activate
source activate mebocost

#参数

input_type='expr_mat'
input_file="${workdir}03_mebocost/exp_mat.csv.gz"
cell_anno_csvfile="${workdir}03_mebocost/cell_ann.csv"
groups='celltype'
species='human'
thread=16
output_file="${workdir}03_mebocost/commu.pk"
# 运行mebocost

python ${project_scripts}/run_mebocost.py \
  --input_type $input_type \
  --input_file $input_file \
  --cell_anno_csvfile $cell_anno_csvfile \
  --groups $groups \
  --species $species \
  --thread $thread \
  --output_file $output_file 








## the first step outfile as input for ploting ...
#6种软件提供的可视化结果，根据项目需求调整，例如关注个别细胞类型作为分泌细胞等情况，最好手动画图。
workpath="${workdir}03_mebocost/"  
python ${project_scripts}/mebocost_plot.py $workpath

