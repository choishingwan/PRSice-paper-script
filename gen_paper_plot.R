library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
files <- list.files(pattern = "result.csv$", recursive=T)
res <- NULL
idx <- 1
for (i in files) {
    tmp <- read.csv(i)
    tmp$Perm <- idx
    idx <- idx + 1
    res <- rbind(res, tmp)
}
levels(res$Program)[2] <- "LDpred"
# Generate Time plot
# Heritability=0 is only done for PRSice-1.25 vs PRSice-2 comparison
# While it doesn't really matter if we merge that into our final 
# results, it might be best for me to separate it out such that 
# PRSice-2 doesn't have the benefit of having more run (thus 
# more "stable" run-time)
res %>%
    filter(!is.na(Target.R2)) %>%
    filter(Heritability!=0) %>% 
    group_by(Target.Size, Program) %>%
    summarize(
        Avg.Time = mean(Real / 60, na.rm = T),
        SE.Time = sd(Real / 60, na.rm = T) / sqrt(n()),
        SD.Time = sd(Real / 60, na.rm = T) ,
        Avg.Mem = mean(MaxKB / 1024 ^ 2),
        SE.Mem = sd(MaxKB / 1024 ^ 2) / sqrt(n()),
        SD.Mem = sd(MaxKB / 1024 ^ 2) ,
        N = n()
    ) -> time.res
time.res$Program <-
    factor(time.res$Program,
           levels = c("PRSice-1.25", "LDpred","lassosum", "PRSice-2"))


theme_sam <-
    theme_classic() + theme(
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size =
                                        18),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust =
                                       1),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()
    )
line.size <- 1
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

res %>%
    filter(!is.na(Target.R2)) %>%
    filter(Target.Size< 100000) %>%
    filter(Heritability==0) %>%
    group_by(Target.Size, Program) %>%
    summarize(
        Avg.Time = mean(Real / 60, na.rm = T),
        SE.Time = sd(Real / 60, na.rm = T) / sqrt(n()),
        SD.Time = sd(Real / 60, na.rm = T) ,
        Avg.Mem = mean(MaxKB / 1024 ^ 2),
        SE.Mem = sd(MaxKB / 1024 ^ 2) / sqrt(n()),
        SD.Mem = sd(MaxKB / 1024 ^ 2) ,
        N = n()
    ) -> prsice.time.res

print(prsice.time.res)
# We want to keep the coloring for each program consistent across plots
# Green = PRSice-2 
# Red = lassosum
# Blue = LDpred
default.color <- gg_color_hue(4)
d <- default.color[1]
default.color[1] <- default.color[3]
default.color[3] <- d
prsice.time.graph <- prsice.time.res %>% 
	ggplot(
		aes(x = log10(Target.Size),
			y = Avg.Time,
			ymin = Avg.Time - SE.Time,
			ymax = Avg.Time + SE.Time,
			color = Program
		)
	) + 
	geom_point() +
    geom_pointrange() +
    geom_line(size = line.size) +
    scale_x_continuous(breaks = c(2, 3, 4, 5), labels = comma(c(100, 1000, 10000, 100000))) +
    labs(x = "Number of Samples", y = "Average Time (in Minute)") +
    guides(color = guide_legend(title = "Program")) +
    theme_sam+
	scale_color_manual(values=c(default.color[4],default.color[3]))

ggsave("PRSice.time.png", prsice.time.graph, height = 7, width = 7)		

prsice.mem.graph <- prsice.time.res %>% 
	ggplot(aes(
        x = log10(Target.Size),
        y = Avg.Mem,
        ymin = Avg.Mem - SE.Mem,
        ymax = Avg.Mem + SE.Mem,
        color = Program
    )) +
    geom_point() +
    geom_pointrange() +
    geom_line(size = line.size) +
    scale_x_continuous(breaks = c(2, 3, 4, 5), labels = comma(c(100, 1000, 10000, 100000))) +
    guides(color = guide_legend(title = "Program")) +
    labs(x = "Number of Samples", y = "Average Memory (in GB)") +
    theme_sam+
	scale_color_manual(values=c(default.color[4],default.color[3]))

ggsave("PRSice.Memory.png", width = 7, height = 7, prsice.mem.graph)

time.graph <- time.res %>%
    filter(Program != "PRSice-1.25") %>%
    ggplot(
        aes(
            x = log10(Target.Size),
            y = Avg.Time,
            ymin = Avg.Time - SE.Time,
            ymax = Avg.Time + SE.Time,
            color = Program
        )
    ) +
    geom_point() +
    geom_pointrange() +
    geom_line(size = line.size) +
    scale_x_continuous(breaks = c(2, 3, 4, 5), labels = comma(c(100, 1000, 10000, 100000))) +
    labs(x = "Number of Samples", y = "Average Time (in Minute)") +
    guides(color = guide_legend(title = "Program")) +
    scale_color_manual(values=c(default.color[2], default.color[1],default.color[3]))+
    theme_sam

ggsave("Time.png", width = 7, height = 7, time.graph)

mem.graph <- time.res %>%
    filter(Program != "PRSice-1.25") %>%
    ggplot(aes(
        x = log10(Target.Size),
        y = Avg.Mem,
        ymin = Avg.Mem - SE.Mem,
        ymax = Avg.Mem + SE.Mem,
        color = Program
    )) +
    geom_point() +
    geom_pointrange() +
    geom_line(size = line.size) +
    scale_x_continuous(breaks = c(2, 3, 4, 5), labels = comma(c(100, 1000, 10000, 100000))) +
    guides(color = guide_legend(title = "Program")) +
    labs(x = "Number of Samples", y = "Average Memory (in GB)") +
    scale_color_manual(values=c(default.color[2], default.color[1],default.color[3]))+
	theme_sam
	
ggsave("Memory.png", width = 7, height = 7, mem.graph)

# Now the overall performance graphs
performance <- subset(res, Program != "PRSice-1.25")
performance$Program <-
    factor(performance$Program,
           levels = c("lassosum", "PRSice-2", "LDpred"))
performance %>%
    filter(Heritability == 0.6) %>%
    filter(Target.Size == 10000) %>%
    ggplot(aes(
        x = log10(Num.Causal),
        y = Target.R2,
        color = Program,
        fill = Program,
        group = interaction(Program, Num.Causal)
    )) +
    geom_boxplot(position = position_dodge()) +
    guides(color = guide_legend(title = "Program")) +
	scale_color_manual(values=c(default.color[1],default.color[3], default.color[2]))+
    scale_fill_manual(values=c(default.color[1],default.color[3], default.color[2]))+
    labs(x = "Number of Causal SNPs",
         y = expression(paste("Trait variance explained (", R ^ 2, ")", sep = " "))) +
    theme_sam +
    scale_x_continuous(breaks = c(2, 3, 4, 5), labels = comma(c(100, 1000, 10000, 100000))) -> perform.plot
ggsave("Performance.png",
       width = 7,
       height = 7,
       perform.plot)

performance %>%
    filter(Heritability != 0) %>%
    ggplot(aes(
        x = log10(Num.Causal),
        y = Target.R2,
        color = Program,
        fill = Program,
        group = interaction(Program, Num.Causal)
    )) +
    geom_boxplot(position = position_dodge()) +
    guides(color = guide_legend(title = "Program")) +
    labs(x = "Number of Causal SNPs",
         y = expression(paste("Trait variance explained (", R ^ 2, ")", sep = " "))) +
    theme_sam +
	scale_color_manual(values=c(default.color[1],default.color[3], default.color[2]))+
	scale_fill_manual(values=c(default.color[1],default.color[3], default.color[2]))+
    scale_x_continuous(breaks = c(2, 3, 4, 5), labels = comma(c(100, 1000, 10000, 100000))) +
    facet_grid(Heritability ~ Target.Size, scales = "free") -> perform.plot.full

ggsave(
    "Performance.supplementary.png",
    width = 7,
    height = 7,
    perform.plot.full
)


t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
    if( equal.variance==FALSE ) 
    {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else{
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
        df <- n1+n2-2
    }      
    t <- (m1-m2-m0)/se 
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat) 
}
# We want to compare the average performance of all program 
# And know the fold change between PRSice-2 when compared to others at the 
# "extreme" end
prsice <- subset(time.res, Target.Size==100000 & Program=="PRSice-2")
lassosum <- subset(time.res, Target.Size==100000 & Program=="lassosum")
ldpred <- subset(time.res, Target.Size==100000 & Program=="LDpred")
print(paste("lassosum / PRSice-2 =  ", lassosum$Avg.Time/prsice$Avg.Time))
print(paste("LDpred / PRSice-2 =  ", ldpred$Avg.Time/prsice$Avg.Time))
prsice <- subset(prsice.time.res, Target.Size==10000 & Program=="PRSice-2")
prsice.old <- subset(prsice.time.res, Target.Size==10000 & Program=="PRSice-1.25")
print(paste("PRSice-1.25 / PRSice-2 =  ", prsice.old$Avg.Time/prsice$Avg.Time))


print("Avg.Mem")
time.res %>%
    select(Program, Target.Size, Avg.Mem) %>%
    spread(Target.Size, Avg.Mem) %>%
    print
print("SE.Mem")
time.res %>%
    select(Program, Target.Size, SE.Mem) %>%
    spread(Target.Size, SE.Mem) %>%
    print
print("Avg.Time")
time.res %>%
    select(Program, Target.Size, Avg.Time) %>%
    spread(Target.Size, Avg.Time) %>%
    print
print("Time.SE")
time.res %>%
    select(Program, Target.Size, SE.Time) %>%
    spread(Target.Size, SE.Time) %>%
    print



time.res %>%
    select(Program, Target.Size, N) %>%
    spread(Target.Size, N) %>%
    print
levels(res$Program)[3:4] <- c("PRSice", "PRSice2")
res %>% 
    filter(Program=="PRSice2" | Program=="PRSice") %>% 
    filter(Heritability==0 & Target.Size < 100000) %>% 
    select(Program, Real, Perm, Target.Size, Num.Causal)%>%
    spread(Program, Real) -> old.vs.new
prsice.comp.res <- t.test(old.vs.new$PRSice2, old.vs.new$PRSice, alternative="less")
fold <- mean(old.vs.new$PRSice/old.vs.new$PRSice2, na.rm=T)
print(paste(
    "PRSice-2 vs PRSice-1.25: p =",
    format(prsice.comp.res$p.value, digits = 3),
    "t =",
    format(prsice.comp.res$statistic, digits = 3),
    "fold =",
    format(fold, digits = 5)
))

res %>% 
    filter(Program=="PRSice2" | Program=="lassosum") %>% 
    filter(Heritability!=0) %>% 
    select(Program, Real, Perm, Target.Size, Num.Causal, Heritability)%>%
    spread(Program, Real) -> prsice.vs.lassosum
fold <- mean(prsice.vs.lassosum$lassosum/prsice.vs.lassosum$PRSice2, na.rm=T)
prsice.vs.lassosum.res <- t.test(prsice.vs.lassosum$PRSice2, prsice.vs.lassosum$lassosum, alternative="less")

print(paste(
    "Time PRSice-2 vs lassosum: p =",
    format(prsice.vs.lassosum.res$p.value, digits = 3),
    "t =",
    format(prsice.vs.lassosum.res$statistic, digits = 3),
    "fold =",
    format(fold, digits = 5)
))

res %>% 
    filter(Program=="PRSice2" | Program=="lassosum") %>% 
    filter(Heritability!=0) %>% 
    select(Program, MaxKB, Perm, Target.Size, Num.Causal, Heritability)%>%
    spread(Program, MaxKB) -> prsice.vs.lassosum
fold <- mean(prsice.vs.lassosum$lassosum/prsice.vs.lassosum$PRSice2, na.rm=T)
prsice.vs.lassosum.res <- t.test(prsice.vs.lassosum$PRSice2, prsice.vs.lassosum$lassosum, alternative="less")

print(paste(
    "Mem PRSice-2 vs lassosum: p =",
    format(prsice.vs.lassosum.res$p.value, digits = 3),
    "t =",
    format(prsice.vs.lassosum.res$statistic, digits = 3),
    "fold =",
    format(fold, digits = 5)
))

res %>% 
    filter(Program=="PRSice2" | Program=="LDpred") %>% 
    filter(Heritability!=0) %>% 
    select(Program, Real, Perm, Target.Size, Num.Causal, Heritability)%>%
    spread(Program, Real) -> prsice.vs.ldpred
fold <- mean(prsice.vs.ldpred$LDpred/prsice.vs.ldpred$PRSice2, na.rm=T)
prsice.vs.ldpred.res <- t.test(prsice.vs.ldpred$PRSice2, prsice.vs.ldpred$LDpred, alternative="less")

print(paste(
    "Time PRSice-2 vs LDpred: p =",
    format(prsice.vs.ldpred.res$p.value, digits = 3),
    "t =",
    format(prsice.vs.ldpred.res$statistic, digits = 3),
    "fold =",
    format(fold, digits = 5)
))


res %>% 
    filter(Program=="PRSice2" | Program=="LDpred") %>% 
    filter(Heritability!=0) %>% 
    select(Program, MaxKB, Perm, Target.Size, Num.Causal, Heritability)%>%
    spread(Program, MaxKB) -> prsice.vs.ldpred
fold <- mean(prsice.vs.ldpred$LDpred/prsice.vs.ldpred$PRSice2, na.rm=T)
prsice.vs.ldpred.res <- t.test(prsice.vs.ldpred$PRSice2, prsice.vs.ldpred$LDpred, alternative="less")

print(paste(
    "Mem PRSice-2 vs LDpred: p =",
    format(prsice.vs.ldpred.res$p.value, digits = 3),
    "t =",
    format(prsice.vs.ldpred.res$statistic, digits = 3),
    "fold =",
    format(fold, digits = 5)
))
