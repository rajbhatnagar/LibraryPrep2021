run_dex_with_formula <- function(
	datalist = datalist,
	coldata = coldata,
	design_formula = design_formula,
	traitval_control = traitval_control,
	traitval_case = traitval_case)
{
	coldata <- as.data.frame(coldata)

	# Get all columns/variables of interest
	columns <- strsplit(as.character(design_formula)[2],
						" + ",
						fixed = TRUE)[[1]]

	# In case some trait names in the formula have spaces which were converted to '\"'
	columns <- gsub('\"', "", columns)

	# Check if all variables specified in the formula exist in the sample traits:
	columns_present <- is.element(columns, colnames(coldata))
	names(columns_present) <- columns
	if(any(!columns_present))
	{
		cat(paste0("\nERROR: The following variables in the design_formula could not be found: ", 
			names(columns_present)[!columns_present], "."))
		cat(paste0("The design formula is this: ", design_formula))
		return(FALSE)
	}

	# Make all trait names are syntactically valid names (no spaces etc)
	columns <- make.names(columns)
	colnames(coldata) <- make.names(colnames(coldata))
	design_formula <- as.formula(paste("~ ", paste(columns, collapse= " + ")))

	# Subset coldata to the variables of interest:
	coldata <- coldata[, columns, drop = FALSE]

	# Are there NAs in the columns of interest?
	coldata_NA <- apply(
		coldata[, columns, drop = FALSE],
		2,
		function(x) any(is.na(x)))
	if (any(coldata_NA)) {
		print(paste0("\nERROR: The following variables in the design_formula contain NAs: ", 
			names(coldata_NA)[coldata_NA]))
		return(FALSE)
	}

	# Make sure the first variable does only have two levels
	if (length(unique(coldata[, columns[1]])) != 2) {
		print(paste0("\nERROR: The trait of interest should have exactly two levels."))
		return(FALSE)
	}

	# Drop unused factor levels for the traits that aren't numerics
	for (i in 1:ncol(coldata)) {
		if (!is.numeric(coldata[, i])) {
			coldata[, i] = as.factor(as.character(coldata[, i]))
		}
	}

	# Relevel first variable:
	coldata[, columns[1]] <- relevel(coldata[, columns[1]], ref=traitval_control)

	# Build a DESeq Object
	if (any(is.na(datalist))) {
		print('ERROR: NAs in the datalist')
		print(sum(is.na(datalist)))
		head(datalist)
	}
	if (any(is.na(coldata))) {
		print('ERROR: NAs in the coldata')
	}

	dds <- DESeqDataSetFromMatrix(
		countData = datalist,
		colData = coldata,
		design = as.formula(paste("~", design_formula[2])))
	dds <- DESeq(dds)

	# Take the first variable and compare case vs. control
	res <- results(dds, contrast = c(columns[1], traitval_case, traitval_control))

	cat(paste("Computed DESeq2", res@elementMetadata@listData$description[2]), "\n\n")
	return(list(res=res, dds=dds))
}
