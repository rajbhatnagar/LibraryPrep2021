import os, sys
import pandas as pd
import numpy as np
import scipy
from scipy.stats import pearsonr

# display
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
plt.rc('font', family='Helvetica')
import seaborn as sns
np.set_printoptions(linewidth=120)


#------------------------------------------------------------------------------
# UTILITY FUNCTIONS

def jaccard_similarity(list1, list2):
	s1 = set(list1)
	s2 = set(list2)
	return float(len(s1.intersection(s2)) / len(s1.union(s2)))

def pwcorr(a):
	""" Return pairwise correlations and pearson p-values"""
	assert type(a) == pd.DataFrame
	n  = len(a.columns)

	# correlations and pvalues
	o  = np.zeros([n,n])
	oo = np.zeros([n,n])
	for i in range (0,n):
		for j in range(i,n):
			(rho, p) = pearsonr(a.iloc[:,i], a.iloc[:,j])
			o[i][j]  = rho
			o[j][i]  = o[i][j]
			oo[i][j] = p
			oo[j][i] = oo[i][j]
	return pd.DataFrame(o,columns=a.columns,index=a.columns), pd.DataFrame(oo,columns=a.columns,index=a.columns)
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# FILE I/O

# Where to store results
FIG_DIR  =  'fig/expression_sets/'
if not os.path.exists(FIG_DIR):
	os.makedirs(FIG_DIR)

# get expression data
df_rcf  = pd.read_csv('data/counts/LibraryPrep_downsampled_filtered_counts.csv')
df_rlog = pd.read_csv('results/normalized_counts/LibraryPrep_rlog.csv')

# Get metadata for expression set
pt0 = pd.read_csv('data/metadata/LibraryPrep_metadata.csv')
pt = pt0.set_index('samplename')
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
# EXPRESSION SET DATA

print('\n======== LOADING DATA ========')

# get experimental groups in expressions sets
d = {}
groups = list(pt.groupby(['Input_concentration', 'Library_prep']))
for (key, group) in groups:
	print('\nLoading:', key)
	print('\tSamples:', group.index.values)
	d[key[0]+'__'+key[1]] = group.index.values

# average samples over groups
dmean = {}
for key in d.keys():
	try:
		members    = list(d[key])
		members += ['Ensembl']
		dmean[key] = df_rlog[members].set_index('Ensembl').mean(axis=1)
	except:
		members = list(d[key])
		members = ['X'+m for m in members]
		members += ['Ensembl']
		dmean[key] = df_rlog[members].set_index('Ensembl').mean(axis=1)
dfmean = pd.DataFrame(dmean)

# pairwise correlation
(corrX, pval) = pwcorr(dfmean)
print('\nAll correlations > 0.97?', np.all(corrX > 0.97))
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# EXPRESSION SET PLOTS
print('\n\n======== EXPRESSION SET PLOTS ========')

sns.set(style="ticks", color_codes=True)
plt.rc('font', family='Helvetica')

# Plot #1
sns.set(style="ticks", color_codes=True)
g = sns.pairplot(dfmean.sample(300), 
	plot_kws=dict(color='k', s=5, linewidth=0),
	corner=True,
	)
g.fig.suptitle("Pairwise Correlations of Regularized Log Values Averaged in Groups")

xyloc = (.1,.85)
for i in range(10):
	for j in range(10):
		if i <= j:
			continue
		label =  r'$\rho$ = ' + str(round(corrX.iloc[i,j], 2))
		g.axes[i,j].annotate(label, xy=xyloc, size=18, xycoords='axes fraction')

filename = FIG_DIR+'pw_rlog1'
for ext in ['.pdf']:
	g.savefig(filename+ext, width=8, height=8)
	print('\n* Wrote', filename+ext)


# Plot #2 -- just a subset of the first plot
cols_to_keep = ['500_ng__Illumina_truseq', '200_ng__Swift_Rapid', '100_ng__Swift']
dfmean1 = dfmean[cols_to_keep]
dfmean1 = dfmean1.rename(lambda x: x.replace('_',' '), axis=1)
(corrX, pval) = pwcorr(dfmean1)

g = sns.pairplot(dfmean1.sample(300), 
	plot_kws=dict(color='k', s=5, linewidth=0),
	corner=True,
	height=2
	)
xyloc = (.1,.85)
label = r'$\rho$ = ' + str(round(corrX.iloc[2,1], 2))
g.axes[2,1].annotate(label, xy=xyloc, size=18, xycoords='axes fraction')

label = r'$\rho$ = ' + str(round(corrX.iloc[2,0], 2))
g.axes[2,0].annotate(label, xy=xyloc, size=18, xycoords='axes fraction')

label = r'$\rho$ = ' + str(round(corrX.iloc[1,0], 2))
g.axes[1,0].annotate(label, xy=xyloc, size=18, xycoords='axes fraction')

#g.fig.suptitle("Pairwise Correlations of Regularized Log Values Averaged in Groups")
filename = FIG_DIR+'pw_rlog2'
for ext in ['.pdf', '.svg']:
	g.savefig(filename+ext)
	print('\n* Wrote', filename+ext)
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
# JACCARD PLOTS

lowest = {}
highest = {}
T = 250
for column in dfmean.columns:
	dfmeanX = dfmean[column].copy()
	dfmeanX = dfmeanX.sort_values()
	lowest[column]  = set(dfmeanX[0:T].index)
	highest[column] = set(dfmeanX[-T:].index)

ngroups = len(dfmean.columns)
jaccard_low = pd.DataFrame(np.zeros([ngroups, ngroups]), 
				index=dfmean.columns, 
				columns=dfmean.columns)
jaccard_high = pd.DataFrame(np.zeros([ngroups, ngroups]), 
				index=dfmean.columns, 
				columns=dfmean.columns)

for groupA in dfmean.columns:
	for groupB in dfmean.columns:
		jsl = jaccard_similarity(lowest[groupA], lowest[groupB])
		jsh = jaccard_similarity(highest[groupA], highest[groupB])
		jaccard_low.loc[groupA, groupB] = jsl
		jaccard_high.loc[groupA, groupB] = jsh


g = sns.clustermap(jaccard_low, figsize=[7,7], cmap="Blues", vmin=0, vmax=1, annot=True, annot_kws={'size': 6})
g.ax_row_dendrogram.set_visible(False)
g.ax_col_dendrogram.set_visible(False)
g.fig.suptitle('Pairwise Jaccard Similarity\nBottom 250 genes') 
filename = FIG_DIR+'jaccard_bottom'
for ext in ['.pdf', '.svg']:
	g.savefig(filename+ext)
	print('\n* Wrote', filename+ext)


g = sns.clustermap(jaccard_high, figsize=[7,7], cmap="Blues", vmin=0, vmax=1, annot=True, annot_kws={'size': 6})
g.ax_row_dendrogram.set_visible(False)
g.ax_col_dendrogram.set_visible(False)
g.fig.suptitle('Pairwise Jaccard Similarity\nTop 250 genes') 
filename = FIG_DIR+'jaccard_top'
for ext in ['.pdf', '.svg']:
	g.savefig(filename+ext)
	print('\n* Wrote', filename+ext)
#------------------------------------------------------------------------------





#------------------------------------------------------------------------------
# CLEAN METADATA
pt = pt0[['samplename', 'Fastq_file_prefixes', 'Library_prep', 'Input_concentration', 'FRAC_RIBOSOMAL_BASES']]
pt['Fastq_file_prefixes'] = pt['Fastq_file_prefixes'].str.replace('_R1_001','')

pt['group'] = pt['Library_prep'] + ' ' + pt['Input_concentration']
pt['Library_prep'].replace({'Illumina_truseq': 'Illumina Truseq', 
			'Swift_Rapid': 'Swift Rapid'},
			inplace=True
			)
pt['group'] = pt.group.str.replace('_', ' ')

pt['Percentage Ribosomal Bases'] = pt['FRAC_RIBOSOMAL_BASES']*100
pt['group'] = pt['group'].str.replace('10 ng', '010 ng')
pt['group'] = pt['group'].str.replace('50 ng', '050 ng')
pt = pt.sort_values(by=['group'])
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# PICARD METRICS PLOTS
sns.set(style="ticks", color_codes=True, font_scale=.7)
plt.rc('font', family='Helvetica')

color_map  = {'Illumina_truseq': '#BD9020', 
				'Illumina Truseq': '#BD9020', 
				'Swift':'#315498', 
				'Swift_Rapid': '#A22382',
				'Swift Rapid': '#A22382'
				}

# Statistics
group_IL = pt.where(pt.Library_prep=='Illumina Truseq').dropna()['Percentage Ribosomal Bases']
group_S  = pt.where(pt.Library_prep=='Swift').dropna()['Percentage Ribosomal Bases']
group_SR = pt.where(pt.Library_prep=='Swift Rapid').dropna()['Percentage Ribosomal Bases']

print('\tIL vs S:', scipy.stats.ttest_ind(group_IL, group_S))
print('\tIL vs SR:', scipy.stats.ttest_ind(group_IL, group_SR))
print('\tSR vs S:', scipy.stats.ttest_ind(group_SR, group_S))


# Plot
f  = plt.figure(figsize=[5,3])
ax = plt.subplot(111)
sns.stripplot(y='group', x='Percentage Ribosomal Bases', 
			data=pt,
			hue='Library_prep', 
			palette=color_map,
			dodge=False,
			ax=ax)
ax.set_xlim([0.4,1.2])
sns.despine(ax=ax, offset=6)
ax.set_title('Percentage of reads mapped to ribosomal genes')

# Legend
ax.legend(frameon=False)
new_title = 'Library Prep Method'
ax.legend_.set_title(new_title)
new_labels = ['Illumina','Swift','Swift Rapid']
for t, l in zip(ax.legend_.texts, new_labels): 
	t.set_text(l)

f.tight_layout()

filename = FIG_DIR+'ribosomal_percentage'
for ext in ['.pdf', '.svg']:
	f.savefig(filename+ext)
	print('\n* Wrote', filename+ext)
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
# MAPPING EFFICIENCY WITHIN GROUPS

# Construct the STAR summary file in bash
# grep '%' *_Log.final.out | tr -s ' ' | tr -d "\t" | sed  's/_Log.final.out:/\t/' | tr "|" "\t" > summary_log_final_out.tsv

df = pd.read_csv('data/metrics/STAR_log_summary.tsv', 
			delimiter='\t', 
			names=['Fastq_file_prefixes','variable','percentage'])
df = df.sort_values(by='Fastq_file_prefixes')
df.variable   = df.variable.str.strip()

# Clean
df = df[df.variable != 'Mismatch rate per base, %']
df = df[df.variable != 'Deletion rate per base']
df = df[df.variable != 'Insertion rate per base']
df = df[df.variable != '% of chimeric reads']
df = df[df.variable != '% of reads unmapped: too many mismatches']


df['variable'] = df['variable'].str.replace('Uniquely mapped reads %', 'Uniquely mapped')
df['variable'] = df['variable'].str.replace('% of reads unmapped: too short',     'Unmapped')
df['variable'] = df['variable'].str.replace('% of reads unmapped: other',         'Unmapped')
df['variable'] = df['variable'].str.replace('% of reads mapped to too many loci', 'Multimapped')
df['variable'] = df['variable'].str.replace('% of reads mapped to multiple loci', 'Multimapped')


# Join
df = df.set_index('Fastq_file_prefixes').join( pt[['group','Fastq_file_prefixes']].set_index('Fastq_file_prefixes') )
df = df.sort_values(by=['group', 'variable'])

# Type cast
df.percentage = df.percentage.str.replace('%','').astype(float)

# Average by group
dfg  =  df.groupby(['variable', 'group']).mean().reset_index()

dfg = dfg.pivot(index='group', columns='variable')
dfg.columns = dfg.columns.droplevel()
dfg = dfg[['Uniquely mapped','Multimapped','Unmapped']]
dfg['No feature'] = 100 - dfg['Uniquely mapped'] - dfg['Multimapped'] - dfg['Unmapped']


# Plot
f = plt.figure(figsize=[5,3])
ax = plt.subplot(111)
dfg.plot(kind='barh', stacked=True, alpha=0.7, sort_columns=True, ax=ax)
plt.legend(loc="lower center", bbox_to_anchor=(0.5, -0.3), ncol=2, frameon=False)
plt.title('STAR: Gene Counts')
plt.xticks(rotation=0, ha='center')

sns.despine(ax=ax, offset=6)
f.tight_layout()

filename = FIG_DIR+'star_gene_counts'
for ext in ['.pdf', '.svg']:
	f.savefig(filename+ext)
	print('\n* Wrote', filename+ext)
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# COVERAGE AS FUNCTION OF GENE LENGTH

filename = 'data/metrics/picard_rna_coverage.tsv'
ptg = pd.read_csv(filename, delimiter='\t')

ptg = ptg.melt(id_vars='Percent through gene')
ptg = ptg.rename({'variable':'Fastq_file_prefixes', 'value':'Coverage'}, axis=1)
ptg = ptg.set_index('Fastq_file_prefixes').join( pt.set_index('Fastq_file_prefixes') )

# Plot
f = plt.figure(figsize=[5,3])
ax = plt.subplot(111)
sns.lineplot(x="Percent through gene", y="Coverage",
			hue='Library_prep', 
			data=ptg,
			palette=color_map,
			ax=ax)
ax.set_title('')
ax.set_xlabel('Percent through gene')
ax.set_ylabel('Coverage')
ax.get_legend().remove()
sns.despine(ax=ax, offset=6)

f.tight_layout()

filename = FIG_DIR+'picard_coverage'
for ext in ['.pdf', '.svg']:
	f.savefig(filename+ext)
	print('\n* Wrote', filename+ext)
#------------------------------------------------------------------------------





plt.close('all')
print('\n*** Finished ***\n')