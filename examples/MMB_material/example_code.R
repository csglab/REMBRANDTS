setwd('/home/alkallr/Documents/book_chapter/REMBRANDTS/')

# ### RUNNING REMBRANDTS ####
# library(DESeq2); library(gplots); library(tibble);
# library(magrittr); library(data.table);
#
# ### Prepare sample metadata ####
# metadata <- data.table(
#     Label = NA_character_, File = list.files(path = './examples/MMB_material/htseq_count', pattern = 'htseq_ixs_(exon|intron)\\.txt', full.names = T),
#     ReadType = NA_character_, Batch = NA_integer_
# )
#
# metadata[ , Label := gsub('.*/|\\.bam.*', '', File), ]
# metadata[ , ReadType := paste0(gsub('.*_|\\.txt', '', File), 'ic'), ]
# metadata[ , Batch := 1L, ]
#
# write.table(metadata, './examples/MMB_material/sample_metadata.tsv', sep = '\t', row.names = F, quote = F)

# cd /home/alkallr/Documents/book_chapter/REMBRANDTS/
# bash ./REMBRANDTS.sh AD_analysis \
# ./examples/MMB_material/sample_metadata.tsv ./ 0.7 linear

### DOWNSTREAM ANALYSIS ####
setwd('/home/alkallr/Documents/book_chapter/REMBRANDTS/')
library(data.table)
library(coin)

### Import differential stability measurements ####
ds <- fread(input = './out/AD_analysis/stability.filtered.mx.txt')
ds <- as.matrix(x = ds, rownames = 'GeneID')

str(ds, strict.width = 'cut')

### Import SRA sample metadata ####
sra.metadata <- fread('./examples/MMB_material/SraRunTable_PRJNA232669.txt')
sra.metadata <- sra.metadata[ , c('Run', 'disease_status'), ]
sra.metadata[ , disease_status := factor(
    x = gsub(x = disease_status, pattern = "advanced Alzheimer's Disease", replacement = 'AD'),
    levels = c('control', 'AD')
), ]

str(sra.metadata, strict.width = 'cut')

### Match metadata rows and stability matrix columns ####
sample.intersect <- intersect(colnames(ds), sra.metadata$Run)

ds <- ds[ , sample.intersect]
sra.metadata <- sra.metadata[match(sample.intersect, sra.metadata$Run), ]

all(colnames(ds) == sra.metadata$Run) # check that sample names match

### Fit linear models relating dStability to disease status
linear.models <- lm(formula = t(ds) ~ disease_status, data = sra.metadata)

### Compute p-values of model coefficients
linear.models.summaries <- summary(object = linear.models)
linear.models.coef <- coef(object = linear.models.summaries)

linear.models.coef.table <- lapply(
    X = linear.models.coef,
    FUN = function(x) as.data.table(x['disease_statusAD', , drop = F])
)

linear.models.coef.table <- rbindlist(l = linear.models.coef.table, idcol = 'GeneID')

linear.models.coef.table[ , GeneID :=
    gsub(x = GeneID, pattern = '^Response ', replacement = ''), ]

linear.models.coef.table[ , FDR := p.adjust(p = `Pr(>|t|)`, method = 'fdr'), ]

linear.models.coef.table[ , GeneID_noVersion :=
    gsub(x = `GeneID`, pattern = '\\..*$', replacement = ''), ]

str(linear.models.coef.table, strict.width = 'cut')

### Load predicted RNA targets of microRNA families from targetScan ####
tgtScan <- fread('./examples/MMB_material/Predicted_Targets_Info.default_predictions.txt.gz')
tgtScan[ , GeneID_noVersion :=
    gsub(x = `Gene ID`, pattern = '\\..*$', replacement = ''), ]

str(tgtScan, strict.width = 'cut')

microRNA_bm <- table(tgtScan$GeneID_noVersion, tgtScan$`miR Family`)
microRNA_bm <- as.matrix(microRNA_bm)

str(microRNA_bm, strict.width = 'cut')

microRNA_bm[microRNA_bm > 0] <- 1 # binarize matrix

# keep microRNAs with at least 5 targets
microRNA_bm <- microRNA_bm[ , colSums(microRNA_bm) >= 5]

### Match gene names in microRNA_bm and linear.models.coef.table ####
gene.intersect <- intersect(
    rownames(microRNA_bm),
    linear.models.coef.table$GeneID_noVersion
)

microRNA_bm <- microRNA_bm[gene.intersect, ]
gene_dStability <- setNames(
    object = linear.models.coef.table$Estimate,
    nm = linear.models.coef.table$GeneID_noVersion
)[gene.intersect]

all( # check that gene names match
    names(gene_dStability) ==
    rownames(microRNA_bm)
)

### Use Mann-Whitney U test to test whether the stability of each
### microRNA is altered in AD patients
diff_miRNA_activity <- apply(
    X = microRNA_bm, MARGIN = 2,
    FUN = function(gene_binding_status) {
        out <- coin::wilcox_test(
            formula = gene_dStability ~ factor(gene_binding_status, levels = 0:1),
            alternative = 'two.sided'
        )

        data.table(
            test_statistic = statistic(object = out, type = 'test'),
            p.value = pvalue(out)
        )
    }
)

diff_miRNA_activity <- rbindlist(diff_miRNA_activity, idcol = 'miRNA_family')
diff_miRNA_activity[ , FDR := p.adjust(p = p.value, method = 'fdr'), ]

head(diff_miRNA_activity[order(p.value), ], n = 15)
