import os,sys
from mebocost import mebocost
# 添加 main.py 所在的目录到 sys.path
main_path = '/project_res/byzheng/biosoft/MEBOCOST'
sys.path.append(main_path)



#workpath = sys.argv[1]
workpath = "./"
## re-load the previous object if needed
inputfile=os.path.join(workpath, 'commu.pk')
mebo_obj = mebocost.load_obj(inputfile)




## sender and receiver event number
fig1 = mebo_obj.eventnum_bar(
                    sender_focus=[],
                    metabolite_focus=[],
                    sensor_focus=[],
                    receiver_focus=[],
                    xorder=[],
                    and_or='and',
                    pval_method='permutation_test_fdr',
                    pval_cutoff=0.01,
                    comm_score_col='Commu_Score',
                    comm_score_cutoff = 1,
                    cutoff_prop = 0.1,
                    figsize='auto',
                    save=None,
                    show_plot=True,
                    show_num = True,
                    include=['sender-receiver'],
                    group_by_cell=True,
                    colorcmap='tab20',
                    return_fig=True
                )
fig1.savefig(os.path.join(workpath, 'details/mebocost_eventnum.pdf'), format='pdf', dpi=300, bbox_inches='tight')     
fig1.savefig(os.path.join(workpath, 'details/mebocost_eventnum.png'), format='png', dpi=300, bbox_inches='tight')   



## circle plot to show communications between cell groups
## 报错 没过滤prob
fig2 = mebo_obj.commu_network_plot(
                    sender_focus=[],
                    metabolite_focus=[],
                    sensor_focus=[],
                    receiver_focus=[],
                    and_or='and',
                    pval_method='permutation_test_fdr',
                    pval_cutoff=0.01,
                    node_cmap='tab20',
                    figsize='auto',
                    line_cmap='bwr',
                    line_color_vmin=None,
                    line_color_vmax=None,
                    linewidth_norm=(0.2, 1),
                    linewidth_value_range = None,
                    node_size_norm=(50, 200),
                    node_value_range = None,
                    adjust_text_pos_node=True,
                    node_text_hidden = False,
                    node_text_font=8,
                    save=None,
                    #cutoff_prop = 0.25,
                    show_plot=True,
                    comm_score_col='Commu_Score',
                    comm_score_cutoff=1,
                    text_outline=True,
                    return_fig=True
                )
## save figure
#fig2.savefig(os.path.join(workpath, 'details/commu_network_plot.pdf'), format='pdf')
fig2.savefig(os.path.join(workpath, 'details/commu_network_plot.png'), format='png', dpi=300, bbox_inches='tight') 


fig2_1 = mebo_obj.commu_network_plot(
                    sender_focus=["Mast","Fibro","pDC"],
                    metabolite_focus=["Prostaglandin D2"],
                    sensor_focus=[],
                    receiver_focus=["Endo"],
                    and_or='and',
                    pval_method='permutation_test_fdr',
                    pval_cutoff=0.01,
                    node_cmap='tab20',
                    figsize='auto',
                    line_cmap='bwr',
                    line_color_vmin=None,
                    line_color_vmax=None,
                    linewidth_norm=(0.2, 1),
                    linewidth_value_range = None,
                    node_size_norm=(50, 200),
                    node_value_range = None,
                    adjust_text_pos_node=True,
                    node_text_hidden = False,
                    node_text_font=10,
                    save=None,
                    show_plot=True,
                    comm_score_col='Commu_Score',
                    comm_score_cutoff=1,
                    text_outline=True,
                    return_fig=True
                )

#fig2_1.savefig(os.path.join(workpath, 'commu_network_plot.pdf'), format='pdf', dpi=300, bbox_inches='tight')
fig2_1.savefig(os.path.join(workpath, 'details/commu_network_plot_PGD2.png'), format='png', dpi=300, bbox_inches='tight') 


### dot plot to show the number of communications between cells

fig3 = mebo_obj.count_dot_plot(
                        pval_method='permutation_test_fdr',
                        pval_cutoff=0.01,
                        cmap='bwr',
                        figsize='auto',
                        save=None,
                        dot_size_norm =(20, 200),
                        dot_value_range = None,
                        dot_color_vmin=None,
                        dot_color_vmax=None,
                        show_plot=True,
                        comm_score_col='Commu_Score',
                        comm_score_cutoff=1,
                        dendrogram_cluster=True,
                        sender_order=[],
                        cutoff_prop = 0.1,
                        receiver_order=[],
                        return_fig = True
                    )
fig3.savefig(os.path.join(workpath, 'details/count_dot_plot.pdf'), format='pdf', dpi=300, bbox_inches='tight')     
fig3.savefig(os.path.join(workpath, 'details/count_dot_plot.png'), format='png', dpi=300, bbox_inches='tight')  

## Malignant cell was focused, use receiver_focus=[] to include all cell types
fig4 = mebo_obj.commu_dotmap(
                sender_focus=[],
                metabolite_focus=[],
                sensor_focus=[],
                receiver_focus=[],
                and_or='and',
                pval_method='permutation_test_fdr',
                pval_cutoff=0.01,
                figsize='auto',
                cmap='bwr',
                cmap_vmin = None,
                cmap_vmax = None,
                cellpair_order=[],
                met_sensor_order=[],
                dot_size_norm=(10, 150),
                save=None,
                show_plot=True,
                comm_score_col='Commu_Score',
                comm_score_range = None,
                comm_score_cutoff=1,
                cutoff_prop = 0.1,
                swap_axis = False,
                return_fig = True
                )

fig4.savefig(os.path.join(workpath, 'details/commu_dotmap.pdf'), format='pdf', dpi=300, bbox_inches='tight')     
fig4.savefig(os.path.join(workpath, 'details/commu_dotmap.png'), format='png', dpi=300, bbox_inches='tight')  




fig4_1 = mebo_obj.commu_dotmap(
                sender_focus=[],
                metabolite_focus=["Cholesterol","Choline","D-Mannose","Riboflavin","L-Cysteine","Prostaglandin D2"],
                sensor_focus=[],
                receiver_focus=[],
                and_or='and',
                pval_method='permutation_test_fdr',
                pval_cutoff=0.01,
                figsize=(14, 4),#”auto“
                cmap='bwr',
                cmap_vmin = None,
                cmap_vmax = None,
                cellpair_order=[],
                met_sensor_order=[],
                dot_size_norm=(50, 150),
                save=None,
                show_plot=True,
                cutoff_prop = 0.1,
                comm_score_col='Commu_Score',
                comm_score_range = None,
                comm_score_cutoff=1,
                swap_axis = False,
                return_fig = True
                )

fig4_1.savefig(os.path.join(workpath, 'details/commu_dotmap_6mets.pdf'), format='pdf', dpi=300, bbox_inches='tight')     
fig4_1.savefig(os.path.join(workpath, 'details/commu_dotmap_6mets.png'), format='png', dpi=300, bbox_inches='tight')  

## Malignant cell was focused, use receiver_focus=[] to include all cell types
fig5 = mebo_obj.FlowPlot(
                pval_method='permutation_test_fdr',
                pval_cutoff=0.01,
                sender_focus=[],
                metabolite_focus=[],
                sensor_focus=[],
                receiver_focus=[],
                remove_unrelevant = False,
                and_or='and',
                node_label_size=8,
                node_alpha=0.6,
                figsize='auto',
                node_cmap='Set1',
                line_cmap='bwr',
                cutoff_prop = 0.1,
                line_cmap_vmin = None,
                line_cmap_vmax = 15.5,
                node_size_norm=(20, 150),
                node_value_range = None,
                linewidth_norm=(0.5, 5),
                linewidth_value_range = None,
                save=None,
                show_plot=True,
                comm_score_col='Commu_Score',
                comm_score_cutoff=1,
                text_outline=False,
                return_fig = True
            )
fig5.savefig(os.path.join(workpath, 'details/FlowPlot.pdf'), format='pdf', dpi=300, bbox_inches='tight')     
fig5.savefig(os.path.join(workpath, 'details/FlowPlot.png'), format='png', dpi=300, bbox_inches='tight') 



fig5_1 = mebo_obj.FlowPlot(
                pval_method='permutation_test_fdr',
                pval_cutoff=0.01,
                sender_focus=[],
                metabolite_focus=["Prostaglandin D2"],
                sensor_focus=[],
                receiver_focus=[],
                remove_unrelevant = False,
                and_or='and',
                node_label_size=8,
                node_alpha=0.6,
                figsize='auto',
                node_cmap='Set1',
                line_cmap='bwr',
                line_cmap_vmin = None,
                line_cmap_vmax = 15.5,
                node_size_norm=(20, 150),
                node_value_range = None,
                linewidth_norm=(0.5, 5),
                linewidth_value_range = None,
                save=None,
                cutoff_prop = 0.1,
                show_plot=True,
                comm_score_col='Commu_Score',
                comm_score_cutoff=1,
                text_outline=False,
                return_fig = True
            )
fig5_1.savefig(os.path.join(workpath, 'details/FlowPlot_PGD2.pdf'), format='pdf', dpi=300, bbox_inches='tight')     
fig5_1.savefig(os.path.join(workpath, 'details/FlowPlot_PGD2.png'), format='png', dpi=300, bbox_inches='tight') 


fig5_2 = mebo_obj.FlowPlot(
                pval_method='permutation_test_fdr',
                pval_cutoff=0.01,
                sender_focus=[],
                metabolite_focus=["Prostaglandin D2"],
                sensor_focus=[],
                receiver_focus=[],
                remove_unrelevant = False,
                and_or='and',
                node_label_size=8,
                node_alpha=0.6,
                figsize='auto',
                node_cmap='Set1',
                line_cmap='bwr',
                line_cmap_vmin = None,
                line_cmap_vmax = 15.5,
                node_size_norm=(20, 150),
                node_value_range = None,
                linewidth_norm=(0.5, 5),
                linewidth_value_range = None,
                save=None,
                #cutoff_prop = 0,
                show_plot=True,
                comm_score_col='Commu_Score',
                comm_score_cutoff=1,
                text_outline=False,
                return_fig = True
            )
fig5_2.savefig(os.path.join(workpath, 'details/FlowPlot_PGD2-1.pdf'), format='pdf', dpi=300, bbox_inches='tight')     
fig5_2.savefig(os.path.join(workpath, 'details/FlowPlot_PGD2-1.png'), format='png', dpi=300, bbox_inches='tight') 

## violin plot to show the aggregated metabolite enzymes of informative metabolties in communication
### here we show five significant metabolites,
### users can pass several metabolites of interest by provide a list
commu_df = mebo_obj.commu_res.copy()
good_met = commu_df[(commu_df['permutation_test_fdr']<=0.05)]['Metabolite_Name'].sort_values().unique()
print(good_met)
fig6 = mebo_obj.violin_plot(
                    sensor_or_met=["Prostaglandin D2","SLCO2A1"],##good_met[:5], ## only top 5 as example
                    cell_focus=[],
                    cell_order = [],
                    row_zscore = False,
                    cmap=None,
                    vmin=None,
                    vmax=None,
                    figsize='auto',
                    cbar_title='',
                    save=None,
                    show_plot=True
                    )
fig6.savefig(os.path.join(workpath, 'details/enzymes_expr_violin_plot_PGD2.pdf'), format='pdf', dpi=300, bbox_inches='tight')     
fig6.savefig(os.path.join(workpath, 'details/enzymes_expr_violin_plot_PGD2.png'), format='png', dpi=300, bbox_inches='tight') 


## violin plot to show the expression of informative sensors in communication

good_sensor = commu_df[(commu_df['permutation_test_fdr']<=0.05)]['Sensor'].sort_values().unique()
print(good_sensor)
fig7 = mebo_obj.violin_plot(
                    sensor_or_met=good_sensor[:5],## only top 5 as example
                    cell_focus=[],
                    cell_order = [],
                    row_zscore = False,
                    cmap=None,
                    vmin=None,
                    vmax=None,
                    figsize='auto',
                    cbar_title='',
                    save=None,
                    show_plot=True
                    )

fig7.savefig(os.path.join(workpath, 'sensors_expr_violin_plot.pdf'), format='pdf', dpi=300, bbox_inches='tight')     
fig7.savefig(os.path.join(workpath, 'sensors_expr_violin_plot.png'), format='png', dpi=300, bbox_inches='tight') 
                    