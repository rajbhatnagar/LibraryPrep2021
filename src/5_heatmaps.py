import os, sys, time, glob
import pandas as pd
import numpy as np

# display
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import seaborn as sns
sns.set(style="ticks", color_codes=True, font_scale=.52)
plt.rc('font', family='Helvetica')
np.set_printoptions(linewidth=120)


#------------------------------------------------------------------------------
# FILE I/O

# Where to store results
FIG_DIR  =  'fig/expression_sets/'
if not os.path.exists(FIG_DIR):
	os.makedirs(FIG_DIR)
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# COLORS

color_map  = {'Illumina_truseq': '#BD9020', 
				'Illumina Truseq': '#BD9020', 
				'Swift':'#315498', 
				'Swift_Rapid': '#A22382',
				'Swift Rapid': '#A22382'
				}

greens_map = {'10 ng': '#edf8e9',
				'010 ng': '#edf8e9',
				'50 ng': '#bae4b3',
				'050 ng': '#bae4b3',
				'100 ng': '#74c476',
				'200 ng': '#31a354',
				'500 ng': '#006d2c'
				}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# LOAD DATA

# expression data
filename = 'results/normalized_counts/LibraryPrep_rlog.csv'
df = pd.read_csv(filename)
df = df.set_index('Ensembl').drop('Unnamed: 0', axis=1)
df.columns = df.columns.str.replace('X','')


# Get metadata for expression set
pt0 = pd.read_csv('data/metadata/LibraryPrep_metadata.csv')
pt = pt0.set_index('samplename')

pt['Library_prep'].replace({'Illumina_truseq': 'Illumina Truseq', 
			'Swift_Rapid': 'Swift Rapid'},
			inplace=True
			)
pt['Input_concentration'] = pt['Input_concentration'].replace('10_ng', '010_ng').replace('50_ng', '050_ng')
pt['Input_concentration'] = pt['Input_concentration'].str.replace('_',' ')
pt = pt[['Library_prep','Input_concentration']]
pt = pt.sort_values(by=['Input_concentration','Library_prep'])
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
# DEFINE GENE SETS

# Gene set: housekeeping
eisenberg = [
	('C1orf43','ENSG00000143612'),
	('CHMP2A','ENSG00000130724'),
	('EMC7','ENSG00000134153'),
	('GPI','ENSG00000105220'),
	('PSMB2','ENSG00000126067'),
	('PSMB4','ENSG00000159377'),
	('RAB7A','ENSG00000075785'),
	('REEP5','ENSG00000129625'),
	('SNRPD3','ENSG00000100028'),
	('VCP','ENSG00000165280'),
	('VPS29','ENSG00000111237')
]

selection = [x[1] for x in eisenberg]
new_index = {x[1]: x[1] + ' ' + x[0] for x in eisenberg}
dfs1 = df[df.index.isin(selection)]
dfs1 = dfs1.rename(new_index)


# Gene set: oncogenic
sanchez = [
	('ENSG00000105810', 'CDK6'),
	('ENSG00000136997', 'MYC'),
	('ENSG00000110092', 'CCND1'),
	('ENSG00000168036', 'CTNNB1'),
	('ENSG00000141027', 'NCOR1'),
	('ENSG00000134250', 'NOTCH2'),
	('ENSG00000083857', 'FAT1'),
	('ENSG00000142208', 'AKT1'),
	('ENSG00000081059', 'TCF7'),
	('ENSG00000135446', 'CDK4'),
	('ENSG00000135679', 'MDM2'),
	('ENSG00000149311', 'ATM'),
	('ENSG00000137693', 'YAP1'),
	('ENSG00000139687', 'RB1'),
	('ENSG00000105221', 'AKT2'),
	('ENSG00000005339', 'CREBBP'),
	('ENSG00000148400', 'NOTCH1'),
	('ENSG00000166949', 'SMAD3'),
	('ENSG00000141510', 'TP53'),
	('ENSG00000079999', ''),
	('ENSG00000163513', 'TGFBR2'),
	('ENSG00000106799', 'TGFBR1'),
	('ENSG00000198625', 'MDM4'),
	('ENSG00000036257', 'CUL3'),
	('ENSG00000175387', 'SMAD2'),
	('ENSG00000116044', 'NFE2L2'),
	('ENSG00000145675', 'PIK3R1'),
	('ENSG00000171862', 'PTEN'),
	('ENSG00000174197', 'MGA'),
	('ENSG00000105173', 'CCNE1'),
	('ENSG00000109670', 'FBXW7'),
	('ENSG00000150457', 'LATS2'),
	('ENSG00000074181', 'NOTCH3'),
	('ENSG00000134982', 'APC'),
	('ENSG00000141646', 'SMAD4'),
	('ENSG00000125952', 'MAX'),
	('ENSG00000186575', 'NF2'),
	('ENSG00000147889', 'CDKN2A'),
	('ENSG00000121879', 'PIK3CA'),
	('ENSG00000117020', 'AKT3'),
	('ENSG00000131023', 'LATS1')
]
selection = [x[0] for x in sanchez]
new_index = {x[0]: x[0] + ' ' + x[1] for x in sanchez}
dfs2 = df[df.index.isin(selection)]
dfs2 = dfs2.rename(new_index)
#------------------------------------------------------------------------------






#------------------------------------------------------------------------------
# Plot

def make_cluster_map(data, phenotable, zscore, title, figsize=[6,3]):
	annots = pd.concat([phenotable['Library_prep'].map(color_map), 
						phenotable['Input_concentration'].map(greens_map) ], 
						axis=1
	)

	if zscore == 0:
		(vmin, vmax) = (-2,2)
	elif zscore == None:
		(vmin, vmax) = (7,11)

	cm = sns.clustermap(data=data, 
			method='average', 
			metric='euclidean', 
			z_score=zscore, 
			standard_scale=None, 
			figsize=figsize,
			row_cluster=True, 
			col_cluster=True,
			cmap="vlag",
			linewidths=.25,
			col_colors=annots,
			xticklabels=False,
			vmin=vmin, vmax=vmax
	)

	for label in sorted(np.unique(phenotable['Library_prep'])):
		cm.ax_col_dendrogram.bar(0, 0, 
			color=color_map[label],
			label=label, 
			linewidth=0
		)
	l1 = cm.ax_col_dendrogram.legend(title='', loc="center", 
			ncol=3, bbox_to_anchor=(0.47, 0.89), 
			bbox_transform=plt.gcf().transFigure,
			fontsize=5
	)


	for label in sorted(np.unique(phenotable['Input_concentration'])):
		cm.ax_col_dendrogram.bar(0, 0, 
			color=greens_map[label],
			label=label, 
			linewidth=0
		)
	l2 = cm.ax_col_dendrogram.legend(title='', loc="center", 
			ncol=3, bbox_to_anchor=(0.87, 0.89), 
			bbox_transform=plt.gcf().transFigure,
			fontsize=5
	)

	cm.ax_row_dendrogram.set_visible(False)
	# cm.ax_col_dendrogram.set_visible(False)

	filename = FIG_DIR+title+'_'+str(zscore)
	for ext in ['.pdf', '.svg']:
		cm.savefig(filename+ext)
		print('\n* Wrote', filename+ext)

make_cluster_map(data=dfs1, phenotable=pt, zscore=0, title='eisenberg', figsize=[8,2.5])
make_cluster_map(data=dfs2, phenotable=pt, zscore=0, title='sanchez', figsize=[8,7])

make_cluster_map(data=dfs1, phenotable=pt, zscore=None, title='eisenberg', figsize=[8,2.5])
make_cluster_map(data=dfs2, phenotable=pt, zscore=None, title='sanchez', figsize=[8,7])
#------------------------------------------------------------------------------