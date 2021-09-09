library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
source('src/dex_lib.R')
source('src/KEGG_pathway_enrichment.R')
options(width=160)


#----------------------------------------------------------------
# SCRIPT OPTIONS

# Where are the expected input files
DATA_DIR = 'data'

# Reads matrix: genes by samples
READS_FILE = file.path(DATA_DIR, 'counts', 'LibraryPrep_downsampled_filtered_counts.csv')

# Sample metadata file: samples by attributes
METADATA_FILE = file.path(DATA_DIR, 'metadata', 'LibraryPrep_metadata.csv')

# Use gene symbols or Entrez IDs?
USE_ENTREZID = TRUE

# Pathway definitions (see e.g. MSigDB: http://www.broadinstitute.org/gsea/downloads.jsp)
if (USE_ENTREZID) {
	PATHWAY_FILE  =  file.path(DATA_DIR, 'kegg', 'c2.cp.kegg.v7.1.entrez.gmt')
} else {
	PATHWAY_FILE  =  file.path(DATA_DIR, 'kegg', 'c2.cp.kegg.v7.1.symbols.gmt')
}

# Where to store results
OUT_DIR = 'results/'
if (!dir.exists(OUT_DIR)) {dir.create(OUT_DIR)}

OUT_COUNT_DIR = 'results/normalized_counts'
if (!dir.exists(OUT_COUNT_DIR)) {dir.create(OUT_COUNT_DIR)}

OUT_DEX_DIR = 'results/dex'
if (!dir.exists(OUT_DEX_DIR)) {dir.create(OUT_DEX_DIR)}

OUT_PATHWAYS_DIR = 'results/pathways'
if (!dir.exists(OUT_PATHWAYS_DIR)) {dir.create(OUT_PATHWAYS_DIR)}

# Parameters for pathway enrichment
MIN_NUM_GENES   =  10
NUM_BOOTSTRAPS  =  50000
ABSOLUTE        =  FALSE
#----------------------------------------------------------------



#----------------------------------------------------------------
# LOAD READ COUNTS AND CLEAN
df = read.csv(READS_FILE)
df = df[ !duplicated(df['ensemblname']),]
df = as.data.frame(df)

if ('genesymbol' %in% colnames(df)) {
	df$genesymbol = NULL
}

# Load Ensembl to symbol/entrez id mapping using cluster profiler
if (USE_ENTREZID) {
	mapping = clusterProfiler::bitr(df$ensemblname, 
				fromType='ENSEMBL', 
				toType=c('ENTREZID'),
				OrgDb = org.Hs.eg.db
				)
	mapping = mapping[ !duplicated(mapping['ENTREZID']),]
	mapping = mapping[ !duplicated(mapping['ENSEMBL']),]
} else {
	mapping = clusterProfiler::bitr(df$ensemblname, 
			fromType='ENSEMBL', 
			toType=c('SYMBOL'),
			OrgDb = org.Hs.eg.db
			)
	mapping$SYMBOL = toupper(mapping$SYMBOL)
	mapping = mapping[ !duplicated(mapping['SYMBOL']),]
	mapping = mapping[ !duplicated(mapping['ENSEMBL']),]
}


# Add annotations, clean read counts data frame
criteria = df$ensemblname %in% mapping$ENSEMBL
df       = df[criteria,]
colnames(df)[1] = 'ENSEMBL'
df  =  merge(x=df, y=mapping, by="ENSEMBL")
if (USE_ENTREZID) {
	rownames(df) = df$ENTREZID
	df$ENTREZID  = NULL
} else {
	rownames(df) = df$SYMBOL
	df$SYMBOL    = NULL
}
ensembl_buffer = df$ENSEMBL
df$ENSEMBL  = NULL


# Normalize the raw counts with the rlog transformation from DESeq2
df_rlog  =  rlog(as.matrix(df))
df_rlog  =  as.data.frame(df_rlog)
rownames(df_rlog) = rownames(df)
df_rlog$Ensembl   = ensembl_buffer

# Write out normalized read counts
filename_rlog = file.path(OUT_COUNT_DIR, 'LibraryPrep_rlog.csv')
write.csv(df_rlog, filename_rlog)
cat(paste("\n* Wrote out regularized counts to:", filename_rlog, "\n"))
#----------------------------------------------------------------



#----------------------------------------------------------------
# LOAD METADATA AND CLEAN
pt = read.csv(METADATA_FILE)
pt$samplename = make.names(pt$samplename)

# This defines the experimental groups
pt$Library_input = paste0(pt$Library_prep, '__', pt$Input_concentration)
pt$Library_input = as.factor(pt$Library_input)
cat(paste('\n\n* Frequency of library prep and input concentrations:\n'))
print(table(pt[, c('Input_concentration', 'Library_prep')]))

# Align the dataframe and phenotable
rownames(pt) = pt$samplename
df = df[,rownames(pt)]

# Error checking
check  =  all( colnames(df) == rownames(pt) )
if (!check) {
	stop('\n *** ERROR 128: mismatch between reads dataframe and phenotype table ***\n')
}
#----------------------------------------------------------------




#----------------------------------------------------------------
# DEX PARIWISE COMPARISONS


# All experimental groups are defined by Library_input
groups   = sort(unique(pt$Library_input))
n_groups = length(groups)

# Set the control group:
control_groups   = c('Illumina_truseq__500_ng', 'Swift__100_ng', 'Swift_Rapid__200_ng')
n_control_groups = length(control_groups)

# Remember which library and concentration correspond to each group
library_conc = strsplit(as.character(groups), '__')
library_conc_control = strsplit(as.character(control_groups), '__')


# Loop over control groups
for (jj in 1:n_control_groups) 
{
	control_group = as.character(control_groups[jj])

	# Set up directory for DEX outputs, use control group as the base path, case group as filename
	results_dex_dir = file.path(OUT_DEX_DIR, control_group)
	if (!dir.exists(results_dex_dir)) {
		dir.create(file.path(OUT_DEX_DIR, control_group), recursive=T)
	}

	# Set up directory for pathway outputs, use control group as the base path, case group as filename
	results_pathway_dir = file.path(OUT_PATHWAYS_DIR, control_group)
	if (!dir.exists(results_pathway_dir)) {
		dir.create(file.path(OUT_PATHWAYS_DIR, control_group), recursive=T)
	}

	# Loop over case groups
	for (kk in 1:n_groups) {
		case_group = as.character(groups[kk])
		if (control_group == case_group) {next}

		# Comment this line out to perform ALL pairwise comparisons 
		if (library_conc[[kk]][1] != library_conc_control[[jj]][1]) {next}

		cat(paste('\n\n----------', case_group, 'vs', control_group, '----------\n'))

		# Isolate the samples of interest 
		ix  = (case_group == pt$Library_input) | (control_group == pt$Library_input)
		dfx = df[,ix]
		ptx = pt[ix,]
		
		print(ptx[,c('samplename','Library_input','Input_concentration')])
		check  =  all(colnames(dfx) == rownames(ptx))
		if (!check) {stop('ERROR 176: Expression set and phenotype table are not matched correctly.')}

		ptx$Library_input = droplevels(ptx$Library_input)
		ptx$Library_input = relevel(ptx$Library_input, ref=control_group)

		# Step 1: Get differential expression results with DESeq2 (see dex_lib.R)
		dd = run_dex_with_formula(
				datalist         = dfx,
				coldata          = ptx,
				design_formula   = as.formula('~ Library_input'),
				traitval_control = control_group,
				traitval_case    = case_group
			)

		outfilename = file.path(results_dex_dir, paste0(case_group,'.csv'))
		write.csv(as.data.frame(dd$res), file=outfilename)
		cat(paste("\n* Wrote out to", outfilename, "\n"))

		res  =  dd$res
		res  =  res[order(res$padj),]
		cat(paste('\n* Number of differentially expressed genes:', sum(na.omit(res$padj) < .05)))

		# Step 2: Pathway enrichment (see KEGG_pathway_enrichment.R)
		cat(paste('\n* Running pathway enrichment...\n'))
		these_pathway_results = geneset_enrichment_based_on_statistic(
			gene.sets.file = PATHWAY_FILE,
			gene.sets.list = NULL,
			num.genes      = MIN_NUM_GENES,
			test.results   = res,
			statistic      = 'stat',
			fold.change    = 'log2FoldChange',
			num.samples    = NUM_BOOTSTRAPS,
			abs            = ABSOLUTE
			)
		geneset_enrichments = these_pathway_results$pathway.results
		n_sig_pathways      = sum(these_pathway_results$pathway.results$adj.p.val < 0.05)
		cat(paste("\n* Number of significant pathways:", n_sig_pathways, '\n'))

		outfilename = file.path(results_pathway_dir, paste0(case_group,'.csv') )
		write.csv(as.data.frame(geneset_enrichments), file=outfilename)
		cat(paste("\n* Wrote out to", outfilename, "\n\n"))

	}
}
#----------------------------------------------------------------



cat('\n\n*** Finished ***\n\n')