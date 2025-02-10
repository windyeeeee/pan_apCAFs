

ubuntu 18:
conda create -y -n pyscenic python=3.10
source activate pyscenic  
pip install pyscenic
pip uninstall numpy
pip install numpy==1.22.4
vim /home/windye/miniconda3/envs/pyscenic/lib/python3.10/site-packages/pyscenic/cli/utils.py          347行
for name, threshold in auc_thresholds.iteritems() 
revised to 修改为
for name, threshold in auc_thresholds.items(）

vim /home/windye/miniconda3/lib/python3.11/site-packages/loompy/utils.py
#解决办法:
#手动把/path/to/loompy/utils.py （见报错的log文件就可以知道这个文件在哪）第27行改成如下形式，然后重新跑一遍第三步AUCell的计算就可以搞定
vf = int("".join(get_loom_spec_version(f)[0].split(".")))



# 1.grn
pyscenic grn \
--num_workers 8 \
--output grn.tsv \
--method grnboost2 \
scenic_CAF.loom  allTFs_hg38.txt

#quiker ubuntu 18
pyscenic grn --num_workers 8 \
  --sparse \
  --method grnboost2 \
  --output sce.adj.csv \
 scenic_CAF2.loom  allTFs_hg38.txt

#or
#slow 
arboreto_with_multiprocessing.py \
    scenic_CAF2.loom \
    allTFs_hg38.txt \
    --method grnboost2 \
    --output grn.tsv \
    --num_workers 8 \
    --seed 777

#1w 8core 12h	随机森林算法GRNBoost2
	
# 2.cistarget
pyscenic ctx \
grn.tsv hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
--annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname scenic_CAF2.loom\
--mode "dask_multiprocessing" \
--output ctx.csv \
--num_workers 8\
--mask_dropouts

pyscenic ctx \
sce.adj.csv hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
--annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname scenic_CAF2.loom \
--mode "dask_multiprocessing" \
--output ctx.csv \
--num_workers 8 \
--mask_dropouts




# 3.AUCell
pyscenic aucell \
scenic_CAF2.loom \
ctx.csv \
--output aucell.loom \
--num_workers 8


# 4.scope loom  不能用
wget https://raw.githubusercontent.com/vib-singlecell-nf/vsn-pipelines/master/src/scenic/bin/add_visualization.py
wget https://raw.githubusercontent.com/vib-singlecell-nf/vsn-pipelines/master/src/scenic/bin/export_to_loom.py
python add_visualization.py \
    --loom_input aucell.loom \
    --loom_output scenic_visualize.loom \
    --num_workers 8
	
	