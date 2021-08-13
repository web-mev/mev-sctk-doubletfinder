suppressMessages(suppressWarnings(library("singleCellTK")))
suppressMessages(suppressWarnings(library("Matrix")))
suppressMessages(suppressWarnings(library("optparse")))

# args from command line:
args<-commandArgs(TRUE)

option_list <- list(
    make_option(
        c('-f','--input_file'),
        help='Path to the count matrix input.'
    ),
    make_option(
        c('-o','--output_file_prefix'),
        help='The prefix for the output file'
    ),
    make_option(
        c('-d', '--doublet_file_prefix'),
        help = 'The prefix for the doublet class file'
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

# Check that the file was provided:
if (is.null(opt$input_file)){
    message('Need to provide a count matrix with the -f/--input_file arg.')
    quit(status=1)
}

if (is.null(opt$output_file_prefix)) {
    message('Need to provide the prefix for the output file with the -o arg.')
    quit(status=1)
}

# Define the appended filename
method_str <- "DoubletFinder_filtered"


# Import counts as a data.frame
cnts <- read.table(
    file = opt$input_file,
    sep = "\t",
    row.names = 1,
    header=T,
    check.names = FALSE
)

# above, we set check.names=F to prevent the mangling of the sample names.
# Now, we stash those original sample names and run make.names, so that any downstream
# functions, etc. don't run into trouble. In the end, we convert back to the original names
orig_cols = colnames(cnts)
new_colnames = make.names(orig_cols)
colnames(cnts) = new_colnames

colname_mapping = data.frame(
    orig_names = orig_cols,
    row.names=new_colnames,
    stringsAsFactors=F
)

# change to a sparse matrix representation, necessary for SCE
cnts <- as(as.matrix(cnts), "sparseMatrix")

# Create an SCE object from the counts
sce <- SingleCellExperiment(
    assays=list(counts=cnts)
)

# Run DoubletFinder on the sce object.
# NOTE: if you want to parametize the seuratRes, it will affect the name of
#   doubletFinder_doublet_label_resolution_XX. This will need coding around
#   if we want to allow adjusting parameters from the front end.
sce <- runDoubletFinder(
    sce,
    seuratNfeatures = 2000,
    seuratPcs = 1:15,
    seuratRes = 1.5,
    formationRate = 0.075
)


# Create a data frame from the DoubletFinder factors and the SCE object
# cell barcodes
# Additionally, normalize DoubletFinder output to DoubletCells output
df.doublets <- data.frame(
    cell_barcode = as.vector(colnames(cnts)),
    doublet_class = as.vector(
        as.factor(tolower(sce$doubletFinder_doublet_label_resolution_1.5))
    )
)

# remap the barcodes back to the originals. Do this by merging with the mapping
# dataframe, dropping the possibly mangled cell_barcode column, and then renaming
# the columns
df.doublets = merge(df.doublets,colname_mapping,by.x='cell_barcode', by.y=0)
df.doublets = df.doublets[,c('orig_names','doublet_class')]
colnames(df.doublets) = c('cell_barcode','doublet_class')

# Export doublet classification to file
doublet_filename <- paste(
    opt$doublet_file_prefix, 
    method_str,
    'tsv', 
    sep='.'
)
write.table(
    df.doublets,
    doublet_filename,
    sep = "\t",
    quote = F,
    row.names = F
)

# Create a subset SCE object by filtering for singlets
# Again hard coded name, but it's technically dependent on the parameters
# Also, "Singlet" is different from DoubletCells output.
sce.sub <- subsetSCECols(
    sce,
    colData = "doubletFinder_doublet_label_resolution_1.5 == 'Singlet'"
)

# Export count matrix to file, remapping the cell barcodes
output_cnts <- data.frame(as.matrix(counts(sce.sub)))
colnames(output_cnts) = colname_mapping[colnames(output_cnts), 'orig_names']
counts_filename <- paste(
    opt$output_file_prefix, 
    method_str,
    'tsv', 
    sep='.'
)
write.table(
    output_cnts,
    counts_filename,
    sep = "\t",
    quote = F,
    row.names = T
)
