import os, sys, time, glob
import pandas as pd
import numpy as np

# display
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import seaborn as sns
sns.set(style="ticks", color_codes=True, font_scale=.8)
plt.rc('font', family='Helvetica')
np.set_printoptions(linewidth=120)


#------------------------------------------------------------------------------
# UTILITY FUNCTIONS

def load_dex_results(filename, eps=1e-8):
	dr = pd.read_csv(filename)
	dr = dr.sort_values(by='padj')
	dr = dr.dropna(subset=['padj'], axis=0)
	dr['rank'] = dr['padj'].rank()
	dr['logp'] = -np.log10(dr['padj']+eps)
	return dr

def get_cdf(a, bins):
	(hist, bin_edges) = np.histogram(a, bins)
	return np.cumsum(hist), bin_edges
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
# FILE I/O

# Where to store results?
FIG_DIR  =  'fig/dex/'
if not os.path.exists(FIG_DIR):
	os.makedirs(FIG_DIR)

# Where are the DEX results stored? (from dex_analysis.R)
DEX_RESULTS_DIR = 'results/dex/'

# Parameters
eps = 1e-8
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
print('\n\n======== PLOTTING COMBINED DEX CDFS ========')

# Plot options
color_map  = {'Illumina_truseq': '#BD9020', 'Swift':'#315498', 'Swift_Rapid': '#A22382'}
marker_map = {'10_ng': '+', '50_ng': 'x', '100_ng': '.', '200_ng': '^', '500_ng': None}

# The path to data is like: <DEX_DIR>/<CONTROL GROUP>/<CASE GROUP>.csv
control_group_paths = [DEX_RESULTS_DIR+'Illumina_truseq__500_ng',
						DEX_RESULTS_DIR+'Swift__100_ng',
						DEX_RESULTS_DIR+'Swift_Rapid__200_ng'
						]

(fig,ax) = plt.subplots(figsize=[9,3], sharex=True, sharey=True, nrows=1, ncols=3)

for cg_path in control_group_paths:
	# the reference group
	control_group = cg_path.split('/')[-1]
	(control_lib, control_conc) = control_group.split('__')
	print('\nControl group:', control_group)

	dex_csvs = glob.glob(cg_path+'/*csv')

	for (k,dex_csv) in enumerate(dex_csvs):
		# string parsing
		label = os.path.basename(dex_csv).split('.')[0]
		(this_lib, this_conc) = label.split('__')
		if this_lib != control_lib:
			continue

		label = label.replace('__','_').replace('_',' ')
		print('\t['+str(k)+']', label, ':\t', dex_csv)

		dr = load_dex_results(dex_csv, eps)
		bins  = np.linspace(0,1,201)
		(y,x) = get_cdf(dr['padj'], bins)
		y     = y/1e3

		print('\t+ Genes diff expressed at .05 level:', np.sum(dr['padj']<.05), '\n')


		if this_lib=='Illumina_truseq':
			ax[0].plot(x[1:], y, color=color_map[this_lib], 
						linestyle='-', 
						marker=marker_map[this_conc],
						linewidth=.8,
						markersize=4,
						markevery=10,
						label=label)
		if this_lib=='Swift':
			ax[1].plot(x[1:], y, color=color_map[this_lib], 
						linestyle='-', 
						marker=marker_map[this_conc],
						linewidth=.8,
						markersize=4,
						markevery=10,
						label=label)
		if this_lib=='Swift_Rapid':
			ax[2].plot(x[1:], y, color=color_map[this_lib], 
						linestyle='-', 
						marker=marker_map[this_conc],
						linewidth=.8,
						markersize=4,
						markevery=10,
						label=label)


for a in ax:
	sns.despine(ax=a, offset=5)
	a.set_title('')

ax[1].set_xlabel('p adjusted')
ax[0].set_ylabel('Cumulative number\nof differentially expressed genes\nthousands')

ax[0].legend(loc='upper left', fontsize=7, frameon=False)
ax[1].legend(loc='upper left', fontsize=7, frameon=False)
ax[2].legend(loc='upper left', fontsize=7, frameon=False)

fig.tight_layout()
for ext in ['pdf', 'svg']:
	outf  =  FIG_DIR + 'combined_dex_cdf.' + ext
	fig.savefig(outf)
	print('\t+ Wrote:', outf)
#------------------------------------------------------------------------------





#------------------------------------------------------------------------------
print('\n\n======== PLOTTING MA PLOTS ========')

def dummy(x):
	if (x['padj'] < 5e-3) and (x['log2FoldChange']<0):
		return 'red'
	if (x['padj'] < 5e-3) and (x['log2FoldChange']>0):
		return 'blue'
	return 'grey'

(fig,ax) = plt.subplots(figsize=[9,8], sharex=True, sharey=True, nrows=3, ncols=3)

# Where to place each subplot:
dmap =  {'Illumina_truseq__500_ng': 
			{'Illumina truseq 100 ng': (0,1),
			 'Illumina truseq 200 ng': (0,2),
			 'Illumina truseq 50 ng': (0,0)
			 },
		'Swift__100_ng': 
			{'Swift 10 ng': (1,0),
			 'Swift 50 ng': (1,1)
			},
		'Swift_Rapid__200_ng': 
			{'Swift Rapid 50 ng': (2,0),
			 'Swift Rapid 100 ng': (2,1)
			}
		}

# Make each subplot
for (j,cg_path) in enumerate(control_group_paths):
	# the reference group
	control_group = cg_path.split('/')[-1]
	(control_lib, control_conc) = control_group.split('__')
	print('\nControl group:', control_group)

	dex_csvs = glob.glob(cg_path+'/*csv')

	for (k,dex_csv) in enumerate(dex_csvs):
		# string parsing
		label = os.path.basename(dex_csv).split('.')[0]
		(this_lib, this_conc) = label.split('__')
		if this_lib != control_lib:
			continue

		label = label.replace('__','_').replace('_',' ')
		print('\t['+str(k)+']', label, ':\t', dex_csv)

		dr = load_dex_results(dex_csv, eps)
		bins  = np.linspace(0,1,201)
		(y,x) = get_cdf(dr['padj'], bins)
		y     = y/1e3
		dr['color'] = dr.apply(dummy, axis=1)

		(rowk, colk) = dmap[control_group][label]
		ax[rowk,colk].scatter(np.log10(dr.baseMean), dr.log2FoldChange, 
						c=dr.color, marker='.', s=4, linewidths=0,
						alpha=0.5, edgecolors=None)
		ax[rowk,colk].set_title(label)


ax[0,0].set_ylim([-1,1])
for j in range(3):
	for k in range(3):
		sns.despine(ax=ax[j,k], offset=5)

fig.tight_layout()
for ext in ['pdf']:
	outf  =  FIG_DIR + 'dex_MA.' +  ext
	fig.savefig(outf)
	print('\t+ Wrote:', outf)

plt.close('all')
#------------------------------------------------------------------------------



print('\n*** Finished ***\n')
