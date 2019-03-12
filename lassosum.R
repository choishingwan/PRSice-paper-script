library(lassosum)
library(data.table)
library(methods)
library(optparse)
library(dplyr)

option_list <- list(
    make_option(c("-b", "--base"), type = "character", help = "Summary Statistic File"),
    make_option(c("-r", "--ref"), type = "character", help = "Reference file prefix"),
    make_option(c("-t", "--target"), type = "character", help = "Target file prefix"),
    make_option(c("-v", "--valid"), type = "character", help = "Validate file prefix"),
    make_option(c("-g", "--genotype"), type = "character", help = "Genotype file prefix"),
    make_option(c("-l", "--ld"), type = "character", help = "LD region BED file"),
    make_option(c("-s", "--size"), type = "numeric", help = "Number of samples in summary statistic file"),
    make_option(c("-o", "--out"), type = "character", help = "Output prefix")
)

# we don't worry about QCed SNPs as they are not included in the GWAS (in my pipeline)
argv <- parse_args(OptionParser(option_list = option_list))

sum.stat <- argv$base

ref.file <- argv$ref
bfile <- argv$genotype
ld.file <- argv$ld
size <- argv$size
prefix <- argv$out
valid.file <- argv$valid
if (is.null(ld.file) |
    is.null(sum.stat) |
    is.null(ref.file)  |
    is.null(prefix) |
    is.null(size) |
    is.null(argv$target) |
    is.null(bfile)) {
    stop("Please provide all required arguments")
}
ref.bfile <- ref.file
if (is.null(argv$valid)) {
    # We don't want to time lassosum for running the validation step as that isn't used in the PRSice-2 paper
    # So we separate out the steps
    fread(argv$target, data.table = F) %>%
        select(FID, IID, Pheno) -> target.pheno
    target.pheno %>%
        select(FID, IID) -> target.keep
    ss <- fread(sum.stat)
    ss <- ss[!P == 0, ]
    ld <- fread(ld.file)
    cor <- p2cor(p = ss$P,
                 n = size,
                 sign = (ss$BETA))
    out <- lassosum.pipeline(
        cor = cor,
        chr = ss$CHR,
        pos = ss$BP,
        A1 = ss$A1,
        A2 = ss$A2,
        ref.bfile = ref.bfile,
        keep.ref = target.keep,
        test.bfile = bfile,
        keep.test = target.keep,
        LDblocks = ld,
        max.ref.bfile.n = 200000,
        trace = 2
    )
    target.res <- validate(out, pheno = target.pheno)
    r2 <- max(target.res$validation.table$value)
    out2 <-
        subset(out, s = target.res$best.s, lambda = target.res$best.lambda)
    # output the object
    save(out2, r2, file = paste0(prefix, ".Rdata"))
    
    data.frame(Program = "lassosum",
               R2 = r2 ^ 2,
               Validate.R2 = NA) %>%
        write.table(paste0(prefix, ".lassosum"),
                    quote = F,
                    row.names = F)
} else{
    # Validation
    # We on purposely separate out the validation step and the optimization step so that we only time the optimization step
    load(argv$target)
    fread(argv$valid, data.table = F) %>%
        select(FID, IID, Pheno) -> valid.pheno
    valid.pheno %>%
        select(FID, IID) -> valid.keep
    validate.res <-
        validate(out2,
                 test.bfile = bfile,
                 keep.test = valid.keep,
                 pheno = valid.pheno)
    validate.r2 <- max(validate.res$validation.table$value)
    data.frame(Program = "lassosum",
               R2 = r2 ^ 2,
               Validate.R2 = validate.r2 ^ 2) %>%
        write.table(paste0(prefix, ".lassosum"),
                    quote = F,
                    row.names = F)
}
