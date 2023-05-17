library(GEOquery)
library(limma)
library(umap)


gset <- getGEO("GSE57338", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL11532", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


fvarLabels(gset) <- make.names(fvarLabels(gset))

gsms <- paste0("00000000000000000011111111111111111111111111111111",
               "11111111111111111011111111111111111111111111111111",
               "11110011111111110000000000000000000000000000000001",
               "10011111111111111111111111111100000000000000001011",
               "11110000010000010000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "0000000000111")
sml <- strsplit(gsms, split="")[[1]]

ex <- exprs(gset)
ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) # log2 transform


gs <- factor(sml)
groups <- make.names(c("HF","no HF"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] 


v <- vooma(gset, design, plot=T)
v$genes <- fData(gset) 
fit  <- lmFit(v)

cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")


dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=1.5)

#venn diagram
vennDiagram(dT, circle.col=palette())
#Q-Q plot
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
# volcano plot
colnames(fit2) 
ct <- 1        
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE57338", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE57338", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)





DE_genes <- subset(tT2, P.Value < 0.05 & abs(logFC) > 0.2, select=c("ID"))
DE_genes$ID <- as.character(DE_genes$ID)
head(DE_genes$ID)
write.table(DE_genes$ID, file = "DE_genespbtrial.txt", sep = "\t", quote = FALSE, row.names = FALSE)

