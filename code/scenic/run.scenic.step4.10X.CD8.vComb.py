#!/usr/bin/env python3

import os
import matplotlib as mpl

if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize
#from scanpy.plotting._tools.scatterplots import plot_scatter

# set a working directory
#wdir = '/path/to/scenic/OUT.scenic/grn.comb.pos'
wdir = os.path.dirname(os.path.realpath(__file__))
os.chdir( wdir )
# path to loom output, generated from a combination of Scanpy and pySCENIC results:
out_prefix = 'zhangLab10X.CD8.pyScenic'
f_final_loom = 'zhangLab10X.CD8.pyScenic.loom'
f_nrow = 4
f_ncol = 5

sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=300)

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

# scenic output
lf = lp.connect( f_final_loom, mode='r', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)

# create a dictionary of regulons:
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)

# cell annotations from the loom column attributes:
cellAnnot = pd.concat(
    [
        pd.DataFrame(lf.ca['cancerType'], index=lf.ca.CellID ),
        pd.DataFrame(lf.ca['dataset'], index=lf.ca.CellID ),
        pd.DataFrame(lf.ca['dataset.old'], index=lf.ca.CellID ),
        pd.DataFrame(lf.ca['meta.cluster'], index=lf.ca.CellID ),
        pd.DataFrame(lf.ca['meta.cluster.coarse'], index=lf.ca.CellID ),
        pd.DataFrame(lf.ca['miniCluster.size'], index=lf.ca.CellID ),
        pd.DataFrame(lf.ca['nGene'], index=lf.ca.CellID ),
        pd.DataFrame(lf.ca['stype'], index=lf.ca.CellID ),
    ],
    axis=1
)
cellAnnot.columns = [ 'cancerType', 'dataset', 'dataset.old',
                     'meta.cluster', 'meta.cluster.coarse',
                     'miniCluster.size', 'nGene', 'stype']
# capture embeddings:
dr = [
    pd.DataFrame( lf.ca.Embedding, index=lf.ca.CellID )
]
dr_names = [
    meta['embeddings'][0]['name'].replace(" ","_")
]
# rename columns:
for i,x in enumerate( dr ):
    x.columns = ['X','Y']
lf.close()

################################

rss_cellType = regulon_specificity_scores( auc_mtx, cellAnnot['meta.cluster'] )
#rss_cellType
rss_cellType.to_csv("%s.RSS.meta.cluster.csv" % out_prefix)

cats = sorted(list(set(cellAnnot['meta.cluster'])))

fig = plt.figure(figsize=(15, 12))
for c,num in zip(cats, range(1,len(cats)+1)):
    x=rss_cellType.T[c]
    ax = fig.add_subplot(f_nrow,f_ncol,num)
    plot_rss(rss_cellType, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )

fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large' ,
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
plt.savefig("%s-RSS-top5.pdf" % out_prefix,
            dpi=600, bbox_inches = "tight")
plt.show()

topreg = []
for i,c in enumerate(cats):
    topreg.extend(
        list(rss_cellType.T[c].sort_values(ascending=False)[:5].index)
    )
topreg = list(set(topreg))

auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
    auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
####auc_mtx_Z.sort_index(inplace=True)
#auc_mtx_Z.to_csv("%s.auc_mtx_Z.csv" % out_prefix)
#####
binary_mtx, auc_thresholds = binarize( auc_mtx, num_workers=20 )
binary_mtx.head()
binary_mtx.to_csv("%s.binary_mtx.csv" % out_prefix)



######
def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f

colors = sns.color_palette('bright',n_colors=len(cats) )
colorsd = dict( zip( cats, colors ))
colormap = [ colorsd[x] for x in cellAnnot['meta.cluster'] ]

#sns.set()
#sns.set(font_scale=0.8)
#fig = palplot( colors, cats, size=1.0)
#plt.savefig("%s.meta.cluster-heatmap-legend-top5.pdf" % out_prefix,
#            dpi=600, bbox_inches = "tight")
#
#sns.set(font_scale=1.2)
#g = sns.clustermap(auc_mtx_Z[topreg], annot=False,  square=False,
#                   linecolor='gray',
#                   yticklabels=False, xticklabels=True,
#                   vmin=-2, vmax=6, row_colors=colormap,
#                   cmap="YlGnBu", figsize=(21,16) )
#g.cax.set_visible(True)
#g.ax_heatmap.set_ylabel('')
#g.ax_heatmap.set_xlabel('')
#plt.savefig("%s.meta.cluster-heatmap-top5.pdf" % out_prefix,
#            dpi=600, bbox_inches = "tight")

# select regulons:
#r = [ 'TP73(+)', 'TRPS1(+)', 'VDR(+)' ]
#
#fig, axs = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=False)
#for i,ax in enumerate(axs):
#    sns.distplot(auc_mtx[ r[i] ], ax=ax, norm_hist=True, bins=100)
#    ax.plot( [ auc_thresholds[ r[i] ] ]*2, ax.get_ylim(), 'r:')
#    ax.title.set_text( r[i] )
#    ax.set_xlabel('')
#
#fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='large')
#fig.text(0.5, -0.01, 'AUC', ha='center', va='center', rotation='horizontal', size='large')
#
#fig.tight_layout()
#fig.savefig('zhangLab10X.CD8.meta.cluster-binaryPlot2.pdf', dpi=600, bbox_inches='tight')





#
#################################
#
## helper functions (not yet integrated into pySCENIC):
#
#from pyscenic.utils import load_motifs
#import operator as op
#from IPython.display import HTML, display
#
#BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
#COLUMN_NAME_LOGO = "MotifLogo"
#COLUMN_NAME_MOTIF_ID = "MotifID"
#COLUMN_NAME_TARGETS = "TargetGenes"
#
#def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
#    """
#    :param df:
#    :param base_url:
#    """
#    # Make sure the original dataframe is not altered.
#    df = df.copy()
#
#    # Add column with URLs to sequence logo.
#    def create_url(motif_id):
#        return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
#    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
#
#    # Truncate TargetGenes.
#    def truncate(col_val):
#        return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
#    df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
#
#    MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
#    pd.set_option('display.max_colwidth', 200)
#    display(HTML(df.head().to_html(escape=False)))
#    pd.set_option('display.max_colwidth', MAX_COL_WIDTH)
#
#
#df_motifs = load_motifs('zhangLab10X.CD8.reg.csv')
#
#selected_motifs = ['BATF','ETV1']
#df_motifs_sel = df_motifs.iloc[ [ True if x in selected_motifs else False for x in df_motifs.index.get_level_values('TF') ] ,:]
#
##display_logos(df_motifs.head())
#display_logos( df_motifs_sel.sort_values([('Enrichment','NES')], ascending=False).head(9))
#
#
#######
#def colorMap( x, palette='bright' ):
#    import natsort
#    from collections import OrderedDict
#    #
#    n=len(set(x))
#    cpalette = sns.color_palette(palette,n_colors=n )
#    cdict = dict( zip( list(set(x)), cpalette ))
#    cmap = [ cdict[i] for i in x ]
#    cdict = OrderedDict( natsort.natsorted(cdict.items()) )
#    return cmap,cdict
#
#def drplot( dr, colorlab, ax, palette='bright', title=None, **kwargs ):
#    cmap,cdict = colorMap( colorlab, palette )
#    for lab,col in cdict.items():
#        ix = colorlab.loc[colorlab==lab].index
#        ax.scatter( dr['X'][ix], dr['Y'][ix], c=[col]*len(ix), alpha=0.7, label=lab, edgecolors='none')
#    if( title is not None ):
#        ax.set_title(title, fontsize='x-large');
#    #
#    ax.set_xticks([])
#    ax.set_yticks([])
#    ax.spines['top'].set_visible(False)
#    ax.spines['right'].set_visible(False)
#    ax.spines['bottom'].set_visible(False)
#    ax.spines['left'].set_visible(False)
#
#plt.rcParams.update({'font.size':12})
#fig, (ax1,ax2) = plt.subplots(1,2, figsize=(20,10), dpi=150 )
#drplot( dr[0], colorlab=cellAnnot['meta.cluster'], ax=ax1, palette='bright', s=2, title='Highly variable genes - UMAP' )
##drplot( dr[4], colorlab=cellAnnot['meta.cluster'], ax=ax2, palette='bright', s=2, title='SCENIC AUC - UMAP' )
#ax1.legend(loc='right', bbox_to_anchor=(1.15, 0.5), ncol=1, markerscale=2, fontsize='x-large', frameon=False, title="meta.cluster")
#plt.tight_layout()
#plt.savefig("zhangLab10X.CD8_umap.meta.cluster.pdf", dpi=300, bbox_inches = "tight")
#
#
######################
#
#
#
#
#
#
#
