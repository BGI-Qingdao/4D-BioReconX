import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'


'''
python cluster_change.py pcgs.6
'''

def draw_cluster_change(c:str,colors:list,dotsize,fntsz,width,height):
    cm = 1/2.54
    draw_data = target_anno[target_anno.classes==c]
    color_number = len(set(draw_data.cluster))
    plt.figure(figsize=(width*cm,height*cm))
    plt.box(False)
    bp = sns.scatterplot(data=draw_data, size='size', sizes=(dotsize,dotsize),  x='time', y='symbol', marker='o', palette=sns.color_palette(colors,n_colors=color_number), hue='cluster', legend=False)
    bp.get_figure().gca().set_xlabel("")
    bp.get_figure().gca().set_ylabel("")
    #plt.legend(frameon=False, loc=(1.04,0))
    plt.tick_params(axis='both', length = 0)#, labelsize=6)
    plt.yticks(fontsize=6)
    plt.xticks(ticks = np.arange(len(time_list)), fontsize=fntsz, labels=time_list, rotation=45)
    plt.savefig(f'{c}_scatters_cmp_dotplot.pdf', bbox_inches='tight', format='pdf')



if __name__ == '__main__':
    time_list = ['WT', '0hpa','12hpa','36hpa','3dpa','5dpa','7dpa','10dpa','14dpa']
    anno = pd.read_csv('/dellfsqd2/ST_OCEAN/USER/liyao1/scripts/pcg_ap_pattern/stereo_pcg/clustering_projection/project_ap/3.meta_data.csv')
    pcgs = pd.read_csv(sys.argv[1], sep='\t', header=None)
    pcgs.columns=['symbol','gene','classes']

    target_genes = list(pcgs.gene)
    target_anno = anno[anno.gene.isin(target_genes)]
    target_anno = pcgs.merge(target_anno, on='gene')
    target_anno['size'] = [6]*target_anno.shape[0]
    target_anno.to_csv('pattern_id_data.csv',index=False)

    colors = ['#5b8ff9','#5ad8a6','#6dc8ec','#269a99','#9270CA','#5d7092','#f6bd16','#e8684a']
    colors2 = ['#5ad8a6', '#9270CA', '#5d7092','#f6bd16','#e8684a']
    fntsz = 4.5
    dotsize = 40
    draw_cluster_change('Anterior',colors,dotsize,fntsz,3,5.8)
    draw_cluster_change('Posterior',colors2,dotsize,fntsz,3,4)


