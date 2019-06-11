library(data.table)
library(dplyr)
library(tidyr)
library(optparse)
library(codetools)
options(scipen = 999)



option_list <- list(
    make_option(c("--bfile", "-b"), type = "character", help = "PLINK prefix"),
    make_option(c("--fam", "-f"), type = "character", help = "QCed sample list"),
    make_option(c("--snp", "-s"), type = "character", help = "QCed SNP list"),
    make_option(c("--bbs", "-w"), type = "character", help = "Biobank simulation software"),
    make_option(
        c("--out", "-o"),
        type = "character",
        help = "Output prefix",
        default = "Out"
    ),
    make_option(c("--batch"), type = "character", help = "Fam file containing batch information"),
    make_option(c("--pc"), type = "character", help = "File containing PCA information"),
    make_option(c("--centre"), type = "character", help = "File containing centre information"),
    make_option(
        c("--quick"),
        action = "store_true",
        default = F,
        help = "Skip pehnotype generation"
    )
)

argv <- parse_args(OptionParser(option_list = option_list))
if (is.na(argv$bfile) |
    is.na(argv$fam) | is.na(argv$snp) | is.na(argv$bbs)  |
    is.na(argv$batch) | is.na(argv$pc)) {
    stop("Please provide all required arguments")
}
fam.read <- function(..., data.table = FALSE) {
    fread(
        ...,
        data.table = data.table,
        col.names = c("FID", "IID", "Dad", "Mum", "Sex", "Pheno")
    )
}


# Parameters  -------------------------------------------------------------

num.causal <- c(100, 1000, 10000, 100000, 560173)
# Update for paper, but should run full if time allows
#target.size <- 10000
target.size <- c(100, 1000, 10000, 100000)
heritability <- c(0.2, 0.6)
base.size <- c(50000, 200000)
max.batch <- 80
# Biobank simulator
bbs <- argv$bbs

# read in covariates ------------------------------------------------------
print("Read in covariates")
pcs <- fread(argv$pc, data.table = F)
batch <- fam.read(argv$batch) %>%
    mutate(Batch = as.factor(Pheno)) %>%
    select(FID, IID, Batch)
covariate <- batch %>%
    inner_join(pcs) %>%
    drop_na()


# Read in sample IDs ------------------------------------------------------
print("Read in Sample IDs")
samples <- fam.read(argv$fam) %>%
    filter(FID > 0 & FID %in% covariate$FID)
write.table(
    samples,
    paste0(argv$out, ".valid.samples"),
    quote = F,
    row.names = F
)


# Generate Phenotype ------------------------------------------------------

print("Start generating phenotypes")
if (!argv$quick) {
    for (causal in num.causal) {
        print(paste("Generating trait with", causal, "causal SNPs"))
        name <- paste(argv$out, causal, sep = "-")
        system2(
            bbs,
            args = c(
                "-i",
                argv$bfile,
                # use the biobank genotype file
                "-o",
                name,
                "-n",
                causal,
                "-e 2",
                # Effect sizes are simulated under normal distribution
                "-E",
                argv$snp,
                "-k",
                paste0(argv$out, ".valid.samples"),
                # samples to be included in the simulation
                "-H",
                paste(heritability, collapse = ",") # Heritability of phenotype to be simulated
            )
        )
    }
}

# LDpred require the radius, which is recommended to be total number of SNPs / 3000
total.num.snp <- nrow(fread(argv$snp, header = F))
radius <- ceiling(total.num.snp / 3000)
# Generate the file for running GWAS
# Prefix = output prefix
# Heritability = Heritability of trait
# Num.Causal = number of causal SNPs
# Pheno.File = phenotype file
# Sample.Size = number of sample in base
gwas.guide <- data.frame(
    Prefix = character(),
    Heritability = numeric(),
    Num.Causal = numeric(),
    Pheno.File = character(),
    Sample.Size = numeric()
)

# Generate teh file for running PRS
# Sample.File = file containing the target samples
# Prefix = output prefix
# Pheno.File = phenotype file
# Num.Causal = number of causal SNPs
# Target.Size = Size of target sample
# Heritability = heritability
# Base.Size = Size of base data
# Radius = parameter required by LDPred. Should be num SNP / 3000
prs.guide <- data.frame(
    Base.Assoc = character(),
    Base.Size = numeric(),
    Prefix = character(),
    Target.Pheno = character(),
    Target.Size = numeric(),
    Validate.Pheno = character(),
    Num.Causal = numeric(),
    Heritability = numeric(),
    Radius = numeric(),
    Target.Plink = character(),
    Valid.Plink = character()
)
plink.guide <- data.frame(Prefix = character())
for (causal in num.causal) {
    name <- paste(argv$out, causal, sep = "-")
    pheno <- fread(name, header = T, data.table = F)
    info <- inner_join(pheno, covariate) %>%
        drop_na()
    phenotypes <- data.frame(FID = info$FID, IID = info$IID)
    # select maximum number of required base
    base.samples <- info %>%
        sample_n(max(base.size))
    # All target samples should be sampled from this pool
    # And we can do the residualization within these sample in one go
    non.base.samples <- info %>%
        filter(!FID %in% base.samples$FID)
    
    print("Residualizing phenotypes")
    generated <- F
    for (target in target.size) {
        # Other loop is target, because we want to reuse genotype files
        # required by LDpred, which doesn't support remove or keeping samples
        target.samples <- non.base.samples
        cur.target.samples <- target.samples %>%
            sample_n(target)
        target.filter.name <-
            paste(name, target, "target", sep = "-")
        cur.target.samples %>%
            select(FID, IID) %>%
            write.table(target.filter.name,
                        quote = F,
                        row.names = F)
        plink.guide <-
            rbind(plink.guide,
                  data.frame(Prefix = target.filter.name))
        base.name <- paste0(name, "-base")
        for (h in heritability) {
            h.name <- paste0("h_", h)
            cur.name <- paste0(name, "-", h.name)
            target.name <-
                paste(cur.name, target, "target", sep = "-")
            equation <-
                as.formula(paste0(
                    h.name,
                    "~Batch+",
                    paste("PC",
                          1:40,
                          sep = "",
                          collapse = "+")
                ))
            # valid.name <- paste0(cur.name, "-valid")
            if (!generated) {
                base.pheno <-
                    rstandard(lm(equation, data = base.samples))
                base.samples[, h.name] <- base.pheno
                target.pheno <-
                    rstandard(lm(equation , data = target.samples))
                target.samples[,h.name] <- target.pheno
                # Now generate the base files
                for (b in base.size) {
                    base.res <- base.samples %>%
                        select(FID, IID) %>%
                        mutate(Pheno = base.samples[, h.name]) %>%
                        sample_n(b)
                     fwrite(
                        base.res,
                        file = paste(base.name,h, b, sep = "-"),
                        sep = "\t",
                        quote = F,
                        row.names = F
                    )
                    gwas.guide <- rbind(
                        gwas.guide,
                        data.frame(
                            Prefix = paste(base.name,h,b, sep = "-"),
                            Heritability = h,
                            Num.Causal = causal,
                            Pheno.File = paste(base.name,h, b, sep = "-"),
                            Sample.Size = b
                        )
                    )
                    
                }
            }
            target.res <- target.samples %>%
                select(FID, IID) %>%
                mutate(Pheno = target.samples[,h.name]) %>%
                filter(FID %in% cur.target.samples$FID)
            # target.res <- target.data[target.select,c("FID", "IID")]
            #target.res$Pheno <- target.data[target.select, h.name]

            fwrite(
                target.res,
                file = target.name,
                sep = "\t",
                quote = F,
                row.names = F
            )
            for (b in base.size) {
                prs.guide <-
                    rbind(
                        prs.guide,
                        data.frame(
                            Base.Assoc = paste0(base.name, "-",h,"-", b, ".assoc.gz"),
                            Base.Size = b,
                            Prefix = paste0(target.name,"-",b),
                            Target.Pheno = target.name,
                            Target.Size = target,
                            Validate.Pheno = NA,
                            Num.Causal = causal,
                            Heritability = h,
                            Radius = radius,
                            Target.Plink = target.filter.name,
                            Valid.Plink = NA
                        )
                    )
            }
        }
        generated <- T
    }
}



total.num <- nrow(gwas.guide)
chunk2 <-
    function(x, n)
        split(x, cut(seq_along(x), n, labels = FALSE))
chunks <- chunk2(1:total.num, max.batch)
chunk.id <- 1

for (i in chunks) {
    gwas.guide[i,] %>%
        write.table(
            paste(argv$out, "gwas", chunk.id, sep = "."),
            quote = F,
            row.names = F
        )
    chunk.id <- chunk.id + 1
}


total.prs <- nrow(prs.guide)
chunks <- chunk2(1:total.prs, max.batch)
chunk.id <- 1
for (i in chunks) {
    prs.guide[i,] %>%
        write.table(
            paste(argv$out, "prs", chunk.id, sep = "."),
            quote = F,
            row.names = F
        )
    chunk.id <- chunk.id + 1
}
write.table(gwas.guide,
            paste0(argv$out, ".full.gwas"),
            quote = F,
            row.names = F)
write.table(prs.guide,
            paste0(argv$out, ".full.prs"),
            quote = F,
            row.names = F)



total.num <- nrow(plink.guide)
chunk2 <-
    function(x, n)
        split(x, cut(seq_along(x), n, labels = FALSE))
chunks <- chunk2(1:total.num, max.batch)
chunk.id <- 1

for (i in chunks) {
    plink.guide[i,] %>%
        write.table(
            paste(argv$out, "plink", chunk.id, sep = "."),
            quote = F,
            row.names = F
        )
    chunk.id <- chunk.id + 1
}

write.table(
    plink.guide,
    paste(argv$out, "plink", "full", sep = "."),
    quote = F,
    row.names = F ,
    col.names = F
)
