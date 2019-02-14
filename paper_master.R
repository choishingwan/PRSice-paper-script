library(data.table)
library(dplyr)
library(tidyr)
library(optparse)
library(codetools)
options(scipen = 999)



#qsub -v bfile=/mnt/lustre/groups/ukbiobank/usr/sam/project/PRS/software/raw/ukb_binary_v2 -v file=sim3 -v prsice=/mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/PRSice-cpp_development/PRSice ../scripts/run_prsice.sh

#qsub -v bfile=/mnt/lustre/groups/ukbiobank/usr/sam/project/PRS/software/raw/ukb_binary_v2 -v file=sim3 -v lassosum=/mnt/lustre/groups/ukbiobank/usr/sam/project/PRS/software/scripts/lassosum.R -v LD=/mnt/lustre/groups/ukbiobank/usr/sam/speed_test/raw/Berisa.EUR.hg19.bed ../scripts/run_lassosum.sh

#qsub -v bfile=/mnt/lustre/groups/ukbiobank/usr/sam/project/PRS/software/raw/ukb_binary_v2 -v file=sim3 -v ldpred=/mnt/lustre/groups/ukbiobank/usr/sam/project/PRS/software/scripts/ldpred/LDpred.py ../scripts/run_ldpred.sh 
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
    is.na(argv$batch) | is.na(argv$pc) | is.na(argv$centre)) {
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

num.causal <- c(100, 1000, 10000, 100000)
target.size <- c(100, 1000, 10000, 100000)
validate.size <- 10000
heritability <- c(0, 0.2, 0.4, 0.6, 0.8)
base.size <- 200000
max.batch <- 80
# Biobank simulator
bbs <- argv$bbs

# read in covariates ------------------------------------------------------
print("Read in covariates")
pcs <- fread(argv$pc, data.table = F)
centre <- fread(argv$centre, header = T, data.table = F) %>%
    mutate(Centre = as.factor(Centre))
batch <- fam.read(argv$batch) %>%
    mutate(Batch = as.factor(Pheno)) %>%
    select(FID, IID, Batch)
covariate <- batch %>%
    inner_join(pcs) %>%
    inner_join(centre) %>%
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

# Note that the phenotype was calculated without standardizing the genotype (which require the -x flag)
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

# For each analysis we will need to run PRSice 3 times
# 1. Generate the inflated R2 by running with the original summary stat
# 2. Generate the adjusted R2 by running with the adjusted summary stat
# 3. Generate the reference R2 by running with the summary stat generated from removing overlapped samples

total.num.snp <- nrow(fread(argv$snp, header = F))
radius <- ceiling(total.num.snp / 3000)
# now we will generate the overlapped samples
# and also generate the file for running GWAS
# Prefix = output prefix
# Heritability = Heritability of trait
# Num.Causal = number of causal SNPs
# Pheno.File = phenotype file
gwas.guide <- data.frame(
    Prefix = character(),
    Heritability = numeric(),
    Num.Causal = numeric(),
    Pheno.File = character()
)

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
    
    base.validate <-
        sample.int(n = nrow(phenotypes), size = base.size + validate.size)
    base.samples <- sort(head(base.validate, n = base.size))
    validate.samples <-
        sort(base.validate[!base.validate %in% base.samples])
    # All target samples should be sampled from this pool
    # And we can do the residualization within these sample in one go s
    target.pool <-
        c(1:nrow(phenotypes))[!c(1:nrow(phenotypes)) %in% base.validate]
    print("Residualizing phenotypes")
    base.data <- phenotypes[base.samples,]
    validate.data <- phenotypes[validate.samples,]
    valid.filter.name <- paste(name, "valid", sep = "-")
    write.table(validate.data,
                valid.filter.name,
                quote = F,
                row.names = F)
    plink.guide <-
        rbind(plink.guide, data.frame(Prefix = valid.filter.name))
    generated <- F
    target.data <- info[target.pool,]
    for (target in target.size) {
        # We do target before heritability such that we can use the same subset of target
        # for different heritability. This help us to reduce the number bninary file we
        # need to generate for LDPred as it doesn't support sample filtering or
        # SNP filtering
        target.select <- sort(sample.int(n=length(target.pool), size=target))
        target.filter.name <- paste(name, target, "target", sep = "-")
        write.table(target.data[target.select,c("FID","IID")],
                    target.filter.name,
                    quote = F,
                    row.names = F)
        plink.guide <-
            rbind(plink.guide,
                  data.frame(Prefix = target.filter.name))
        for (h in heritability) {
            h.name <- paste0("h_", h)
            cur.name <- paste0(name, "-", h.name)
            target.name <-
                paste(cur.name, target, "target", sep = "-")
            equation <-
                as.formula(paste0(
                    h.name,
                    "~Centre+Batch+",
                    paste("PC", 1:40, sep = "", collapse = "+")
                ))
            # disadvantage of this procedure is that we will repeatly do the
            # regression multiple time.
            base.name <- paste0(cur.name, "-base")
            valid.name <- paste0(cur.name, "-valid")
            if (!generated) {
                base.pheno <-
                    rstandard(lm(equation, data = info[base.samples, ]))
                validate.pheno <-
                    rstandard(lm(equation , data = info[validate.samples, ]))
                base.data$Pheno <- base.pheno
                validate.data$Pheno <- validate.pheno
                fwrite(
                    base.data,
                    file = base.name,
                    sep = "\t",
                    row.names = F,
                    quote = F
                )
                fwrite(
                    validate.data,
                    file = valid.name,
                    sep = "\t",
                    row.names = F,
                    quote = F
                )
                gwas.guide <-
                    rbind(
                        gwas.guide,
                        data.frame(
                            Prefix = base.name,
                            Heritability = h,
                            Num.Causal = causal,
                            Pheno.File = base.name
                        )
                    )
                target.pheno <-
                    rstandard(lm(equation , data = target.data))
                target.data[,h.name] <- target.pheno
            }
            target.res <- target.data[target.select,c("FID", "IID")]
            target.res$Pheno <- target.data[target.select, h.name]
            
            fwrite(
                target.res,
                file = target.name,
                sep = "\t",
                quote = F,
                row.names = F
            )
            prs.guide <-
                rbind(
                    prs.guide,
                    data.frame(
                        Base.Assoc = paste0(base.name, ".assoc.gz"),
                        Base.Size = base.size,
                        Prefix = target.name,
                        Target.Pheno = target.name,
                        Target.Size = target,
                        Validate.Pheno = valid.name,
                        Num.Causal = causal,
                        Heritability = h,
                        Radius = radius,
                        Target.Plink = target.filter.name,
                        Valid.Plink = valid.filter.name
                    )
                )
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
    gwas.guide[i, ] %>%
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
    prs.guide[i, ] %>%
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
    plink.guide[i, ] %>%
        write.table(
            paste(argv$out, "plink", chunk.id, sep = "."),
            quote = F,
            row.names = F
        )
    chunk.id <- chunk.id + 1
}

write.table(plink.guide, paste(argv$out, "plink","full", sep = "."), quote=F, row.names=F , col.names = F)
