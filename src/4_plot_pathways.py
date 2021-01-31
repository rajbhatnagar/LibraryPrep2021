import os, sys, time, glob
import pandas as pd
import numpy as np

# display
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import seaborn as sns
sns.set(style="ticks", color_codes=True, font_scale=.7)
plt.rc('font', family='Helvetica')
np.set_printoptions(linewidth=120)


#------------------------------------------------------------------------------
# Color definitions

# col1 = '#08519c'
# col2 = '#3182bd'
# col3 = '#6baed6'
# col4 = '#bdd7e7'

col1 = '#003B46'
col2 = '#07575B'
col3 = '#66A5AD'
col4 = '#C4DFE6'

#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# FILE I/O

# Where to store results?
FIG_DIR  =  'fig/pathways/'
if not os.path.exists(FIG_DIR):
	os.makedirs(FIG_DIR)

# Where are the pathway results stored? (from dex_analysis.R)
PATHWAY_DIR = 'results/pathways/'

# Parameters
THRESHOLD = 1e-2
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
print('\n* CLEANING PATHWAY DATA')

# get pathway files 
pathway_files = glob.glob(PATHWAY_DIR+'/*/*.csv')

# Loop over each dex analysis
all_entries = []
for csv in pathway_files:
	print('\n\t+ Reading file', csv)

	# What is case and control?
	reference = csv.split('/')[-2]
	if reference not in ['Illumina_truseq__500_ng', 'Swift__100_ng', 'Swift_Rapid__200_ng']:
		continue

	case = csv.split('/')[-1].replace('.csv', '')
	print('\t\t+ Control:', reference)
	print('\t\t+ Case:', case)

	(reference_library, reference_conc) = reference.split('__')
	(case_library, case_conc) = case.split('__')

	# load data
	dp = pd.read_csv(csv)
	criterion = dp['adj.p.val'] < THRESHOLD
	print('\t\t+ Pathways expressed with padj < threshold:', np.sum(criterion), '\n')

	for (j,row) in dp.iterrows():
		entry = [
				reference_library, reference_conc, 
				case_library, case_conc,
				row['Unnamed: 0'],
				row['fold.change'],
				row['adj.p.val']
				]
		all_entries.append(entry)


header = ['reference_library', 'reference_concentration', 
			'case_library', 'case_concentration',
			'pathway', 'fold_change', 'adj_pval']
df0 = pd.DataFrame(all_entries, columns=header)

df = df0.loc[ df0.case_library==df0.reference_library ]
df = df.drop(['reference_concentration', 'case_library'], axis=1)
print('\t Summary dataframe dimensions:', df.shape)
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
# PLOTS

print('\n* MAKING PLOTS')

for reference in ['Illumina_truseq', 'Swift_Rapid', 'Swift']:
	# subset to entries with the reference
	df1 = df[df.reference_library==reference]

	# Get p values first, pick pathways for plotting
	df1_pval = df1.pivot(index='pathway', 
					columns='case_concentration', 
					values='adj_pval')
	pathways_to_keep = df1_pval.index[ np.sum(df1_pval < THRESHOLD, axis=1)>0 ]

	# Get fold changes, subset to chosen pathways
	df1_fold = df1.pivot(index='pathway', 
					columns='case_concentration', 
					values='fold_change')
	df1_fold = df1_fold.loc[pathways_to_keep]
	df1_fold.sort_values(by='50_ng', inplace=True)

	# Make plot
	yrange = range(1, len(df1_fold.index)+1)
	height = 2 + .15*len(df1_fold.index)
	f  = plt.figure(figsize=[9,height])
	ax = plt.subplot(111)

	plt.hlines(y=yrange, 
				xmin=np.min(df1_fold, axis=1), 
				xmax=np.max(df1_fold, axis=1), 
				color='grey', alpha=0.8)

	try:
		plt.scatter(df1_fold['200_ng'], yrange, color=col1, alpha=1.0, label='200_ng')
	except:
		pass
	try:
		plt.scatter(df1_fold['100_ng'], yrange, color=col2, alpha=1.0, label='100_ng')
	except:
		pass
	try:
		plt.scatter(df1_fold['50_ng'],  yrange, color=col3, alpha=1.0, label='50_ng')
	except:
		pass
	try:
		plt.scatter(df1_fold['10_ng'],  yrange, color=col4, alpha=1.0, label='10_ng')
	except:
		pass

	plt.legend(frameon=False)
	plt.yticks(yrange, [x[0:37]+'...'*(len(x)>37) for x in df1_fold.index])
	plt.title(reference, loc='left')
	plt.xlabel('fold change')
	plt.ylabel('pathway')
	ax.set_xlabel('log2 fold change')
	ax.set_ylabel('')
	sns.despine(ax=ax, offset=6)
	f.tight_layout()


	for ext in ['pdf', 'svg']:
		outf  =  FIG_DIR + reference + '_pathways.' + ext
		f.savefig(outf)
		print('\t+ Wrote:', outf)
#------------------------------------------------------------------------------





print('\n*** Finished ***\n')
