setwd("D:/research/RA/figures/figure_3/")
library(corrplot)
library(ggplot2)
library(reshape2)
library(factoextra)
library(cluster)
library(tidyverse)
options(stringsAsFactors = F)

source("./scripts/STIA_helpers.R")
############################################
####################### Figure 3 #######################################
### load processed matrix files ##########
cpm <- read.csv("./objects/cpm.csv",row.names = 1)
log.cpm <- read.csv("./objects/log_cpm.csv",row.names = 1)
cpm.mean <- read.csv("./objects/cpm_mean.csv",row.names = 1)
log.cpm.mean <- read.csv("./objects/log_cpm_mean.csv",row.names = 1)
genes <- read.csv("./objects/genes.csv",row.names = 1)
stia.sample.table <- read.csv("./objects/sample.table.csv",row.names = 1)

######## Figure 3D
#### k means clustering

##### Set Basement expression
# Set basement for average groups
cpm.mean.base <- cpm.mean 
cpm.mean.base[cpm.mean.base<16] <- 16
# Set basement for all Samples
cpm.base <- cpm
cpm.base[cpm.base<16] <- 16

# calculate foldchanges to D0 for each condition
log.cpm.mean.base = log2(cpm.mean.base+1)
stia.fc.change <- cbind(log.cpm.mean.base[,1:5]-log.cpm.mean.base[,1],
                        log.cpm.mean.base[,6:10]-log.cpm.mean.base[,6],
                        log.cpm.mean.base[,11:15]-log.cpm.mean.base[,11],
                        log.cpm.mean.base[,16:20]-log.cpm.mean.base[,16])

# select genes that change from D0
stia.fc.change <- subset(fc.change, rowSums(abs(fc.change)>=1)>=1) # 1772 genes
write.csv(fc.change, "./objects/fc_change_full.csv", row.names = T)

####### K-means clustering in morpheus###########

########### Figure 3E ############
##### Expression line plots of select genes
#### read in clustering results
fc.change.cluster <- read.csv("./output/k_means_clusters.csv",row.names = 1) %>% dplyr::select(ensembl=gene,everything())
fc.change.cluster <- left_join(fc.change.cluster,genes,by="ensembl")
write.csv(fc.change.cluster,"./output/k_means.csv")


################# individual gene expressions #######
stia.cpm <- cpm
lineplot(stia.cpm, "Nr4a3", stia.sample.table) 
lineplot(stia.cpm, "Ccl3", stia.sample.table) 
lineplot(stia.cpm, "Trem2", stia.sample.table) 
lineplot(stia.cpm, "C1qb", stia.sample.table) 

############# Cluster trend line ############
draw_trend_line(fc.change.cluster,"c1",TRUE)
draw_trend_line(fc.change.cluster,"c2",TRUE)
draw_trend_line(fc.change.cluster,"c3",TRUE)
draw_trend_line(fc.change.cluster,"c4",TRUE)
draw_trend_line(fc.change.cluster,"c5",TRUE)
draw_trend_line(fc.change.cluster,"c6",TRUE)

########## Figure 3F ###################
##### monocyte genes over the course of STIA
# get control expressions of monocytes
raw.cpm.d0 <- read.table("./data/CPM_Nurr_MACS_RNASeq_final_count.txt",
                         header = T, row.names = 1)

# remove monocyte samples and gene symbols
genes <- raw.cpm.d0 %>% mutate(ensembl=rownames(.)) %>% 
  dplyr::select(ensembl,symbol=Gene_Symbol) # save gene names

raw.cpm.d0 <- raw.cpm.d0 %>% dplyr::select(starts_with("B6Anna_")) %>% 
  dplyr::select(-ends_with("Ly6clo")) %>%   # remove ankle non-classical samples
  rename_all(funs(str_replace(., "B6Anna_3Mo_", ""))) %>% # remove prefix
  dplyr::select(str_order(colnames(.))) %>% # reorder samples
  dplyr::select(ends_with("MA"),ends_with("MB"),ends_with("MC"),ends_with("MD"),
                ends_with("Class"),ends_with("NC"))

# remove weird gene
raw.cpm.d0 <- raw.cpm.d0[!rownames(raw.cpm.d0)=="ENSMUSG00000098178",]
raw.cpm.log.d0 <- log2(raw.cpm.d0+1)
raw.cpm.log.mean.d0 <- data.frame(MA = rowMeans(raw.cpm.log.d0[,1:5]),
                                  MB = rowMeans(raw.cpm.log.d0[,6:9]),
                                  MC = rowMeans(raw.cpm.log.d0[,10:14]),
                                  MD = rowMeans(raw.cpm.log.d0[,15:19]),
                                  CM = rowMeans(raw.cpm.log.d0[,20:24]),
                                  NC = rowMeans(raw.cpm.log.d0[,25:27]))

# at least 1 macrophage pop. have avg log cpm >=4
log.cpm.mean.d0 <- raw.cpm.log.mean.d0[rowSums(raw.cpm.log.mean.d0[,1:6]>=4)>= 1,]
nrow(log.cpm.mean.d0) # 8135

# count up regulation of steady state K-means clusters in timecourse dataset
clusters <- read.csv("../figure_2/output/subset_signatures/k_means_latest.csv")
CM.genes <- clusters %>% dplyr::filter(cluster==1) %>% pull(ensembl)
CM.genes.cpm <- subset(cpm.base, rownames(cpm.base) %in% CM.genes)
CM.genes.cpm.mean <- subset(cpm.mean.base, rownames(cpm.mean.base) %in% CM.genes)

sample.table <- sample.table %>% mutate(condition=paste(day,celltype,sep="_"))
CM.DE <- DE_timepoints_barplot(CM.genes.cpm,CM.genes.cpm.mean,0.2)
CM.DE[[1]] #up regulated classical monocyte genes
CM.DE[[2]] #down regulated classical monocyte genes

############## overlap of STIA and CIA clusters###########
#### STIA
stia.clusters <- read.csv("../figure_3/output//k_means_clusters.csv") %>%
  select(-1)
#### CIA
cia.clusters <- read.csv("./output/k_means_clusters.csv")

# calculate overlap of two dataset clusters
shared.genes <- intersect(fc.change.cluster$ensembl, clusters$ensembl)
length(shared.genes) # 1223 genes overlap

ss.indexes <- sapply(shared.genes, function(x) match(x, clusters$ensembl))
stia.indexes <- sapply(shared.genes, function(x) match(x, fc.change.cluster$ensembl))

shared.clusters <- data.frame(ensembl = shared.genes,
                              stia.cluster = fc.change.cluster[stia.indexes,2],
                              ss.cluster = clusters[ss.indexes,2])

shared.clusters <- shared.clusters %>% 
  mutate(stia.cluster=paste0("STIA_",stia.cluster),
         ss.cluster=paste0("SS_",ss.cluster)) %>% arrange(stia.cluster,ss.cluster)

cluster.numbers <- as.data.frame(table(clusters$cluster)) %>% pull(Freq)

# table of number of shared genes between each stia and cia cluster pairs
cont.table <- as.data.frame(table(shared.clusters$ss.cluster, 
                                  shared.clusters$stia.cluster)) %>%
  select(ss.cluster=Var1,stia.cluster=Var2,everything()) %>%
  mutate(ss_num = rep(as.data.frame(table(clusters$cluster)) %>% pull(Freq), 6)) %>%
  arrange(stia.cluster) %>% 
  mutate(stia_num = rep(as.data.frame(table(fc.change.cluster$cluster)) %>% pull(Freq), 6),
         freq_ss = Freq/ss_num, freq_stia = Freq/stia_num)

ggplot(cont.table, aes(x=stia.cluster,y=fct_rev(ss.cluster), fill=freq_stia))+
  geom_tile(color="black")+scale_fill_distiller(direction=1)+
  theme_minimal(base_size = 16)


####### Figure 3I
##### function for calculating foldchanges between each timepoint
MA_fc <- between_timepoint_fc(cpm.mean.base,"MA")
colnames(MA_fc) <- c("D0-3_MA", "D3-7_MA", "D7-13_MA", "D13-21_MA")
MB_fc <- between_timepoint_fc(cpm.mean.base,"MB")
colnames(MB_fc) <- c("D0-3_MB", "D3-7_MB", "D7-13_MB", "D13-21_MB")
MC_fc <- between_timepoint_fc(cpm.mean.base,"MC")
colnames(MC_fc) <- c("D0-3_MC", "D3-7_MC", "D7-13_MC", "D13-21_MC")
MD_fc <- between_timepoint_fc(cpm.mean.base,"MD")
colnames(MD_fc) <- c("D0-3_MD", "D3-7_MD", "D7-13_MD", "D13-21_MD")

preceding_fc <- cbind(MA_fc, MB_fc, MC_fc, MD_fc)
rownames(preceding_fc) <- rownames(cpm)

# count up regulated genes between each timepoint
pos <- rbind(MA = colSums(MA_fc>=1),MB = colSums(MB_fc>=1),
             MC = colSums(MC_fc>=1),MD = colSums(MD_fc>=1))

# count down regulated genes between each timepoint
neg <- rbind(MA = colSums(MA_fc<=-1),MB = colSums(MB_fc<=-1),
             MC = colSums(MC_fc<=-1),MD = colSums(MD_fc<=-1))

#################### bar plots #############
# up regulation
df <- data.frame(avg = as.vector(pos),
                 pop = rep(c("MA","MB","MC","MD"),4),
                 day = c(rep("0-3",4),rep("3-7",4),rep("7-13",4),rep("13-21",4)))
df$day <- factor(df$day,levels=c("0-3","3-7","7-13","13-21"))

ggplot(df, aes(x=day, y=avg, fill=pop))+
  geom_bar(stat="identity", position=position_dodge(),width=0.8)+
  scale_fill_manual(values = c("#07693A","#834698","#8CC544","#BC7EAF"))+
  coord_cartesian(ylim = c(0,350))+  theme_classic(base_size = 16)+
  theme(axis.title=element_blank(),
        legend.position="none",
        axis.text=element_text(face="bold"))

# down regulation
df <- data.frame(avg = as.vector(neg),
                 pop = rep(c("MA","MB","MC","MD"),4),
                 day = c(rep("0-3",4),rep("3-7",4),rep("7-13",4),rep("13-21",4)))
df$day <- factor(df$day,levels=c("0-3","3-7","7-13","13-21"))

ggplot(df, aes(x=day, y=avg, fill=pop))+
  geom_bar(stat="identity", position=position_dodge(),width=0.8)+
  scale_y_reverse()+
  coord_cartesian(ylim = c(350,0))+
  scale_fill_manual(values = c("#07693A","#834698","#8CC544","#BC7EAF"))+
  theme_classic(base_size = 16)+
  theme(axis.title=element_blank(),
        legend.position="none",
        axis.text=element_text(face="bold"))

############ Figure S3A #############
### PCA plot of conditions
sample.table <- read.csv("./objects/sample.table.csv")
PCA <- prcomp(t(cpm), scale. = T)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

dat.gg <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                     sample = sample.table$names,
                     day = factor(sample.table$day),
                     celltype=sample.table$celltype)
# pca plot
ggplot(dat.gg,aes(PC1,PC2,color=celltype,label=day))+
  geom_text(size=I(5))+
  scale_color_manual(values = c("#07693A","#834698","#8CC544","#BC7EAF"))+
  labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
       y = paste0("PC2, VarExp:", round(percentVar[2],4)))+
  theme_bw(base_size = 16)+
  theme(axis.text=element_text(face="bold"))

########### Figure S3B ################
### correlation dot plot with monocytes
# find genes shared between two datasets
shared.genes <- intersect(rownames(cpm),rownames(log.cpm.mean.d0))
cpm.subset <- cpm[shared.genes,]
mono.mean.d0 <- data.frame(CM = rowMeans(raw.cpm.d0[,20:24]+1),
                           NC = rowMeans(raw.cpm.d0[,25:27]+1))
mono.mean.d0.subset <- mono.mean.d0[shared.genes,]

cor.with.mono <- data.frame(CM=cor(cpm.subset, mono.mean.d0.subset$CM,method = "spearman"), 
                            NC=cor(cpm.subset, mono.mean.d0.subset$NC,method = "spearman"),
                            day=factor(sample.table$day),
                            celltype=sample.table$celltype)

ggplot(cor.with.mono, aes(x=day, y=CM,color=celltype))+
  geom_point(size=3)+coord_cartesian(ylim=c(0.525,0.725))+
  scale_color_manual(values = c("#07693A","#834698","#8CC544","#BC7EAF"))+
  theme_classic(base_size = 16)+
  theme(axis.text=element_text(face="bold"))

ggplot(cor.with.mono, aes(x=day, y=NC,color=celltype))+
  geom_point(size=3)+coord_cartesian(ylim=c(0.525,0.725))+
  scale_color_manual(values = c("#07693A","#834698","#8CC544","#BC7EAF"))+
  theme_classic(base_size = 16)+
  theme(axis.text=element_text(face="bold"))


########## Figure S3C ###############
### pathway bubble plots
## cluster 1
leukocyte_activation <- read.table("./objects/GO_genesets/leukocye_activation.txt")
leukocyte_activation <- genes %>% filter(symbol %in% leukocyte_activation$V2) %>% pull(ensembl)
leukocyte_activation <- GO.bubble.setup(cpm.mean,leukocyte_activation,"Leukocyte activation",fc.change.cluster,"c1")

cell_chemotaxis <- read.table("./objects/GO_genesets/cell_chemotaxis.txt")
cell_chemotaxis <- genes %>% filter(symbol %in% cell_chemotaxis$V2) %>% pull(ensembl)
cell_chemotaxis <- GO.bubble.setup(cpm.mean,cell_chemotaxis,"Cell chemotaxis",fc.change.cluster,"c1")

collagen_catabolism <- read.table("./objects/GO_genesets/collagen_catabolism.txt")
collagen_catabolism  <- genes %>% filter(symbol %in% collagen_catabolism $V2) %>% pull(ensembl)
collagen_catabolism <- GO.bubble.setup(cpm.mean,collagen_catabolism,"Collegen catabolism",fc.change.cluster,"c1")


df <- rbind(leukocyte_activation, cell_chemotaxis,collagen_catabolism)
df$condition <- factor(df$condition, levels = colnames(log.cpm.mean))
ggplot(df, aes(x=condition,y=fct_rev(pathway),color=e, size=n))+
  geom_point()+theme_classic(base_size = 16)+ 
  labs(color="Z-score",size="Proportion above average")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(face = "bold"))+
  scale_color_gradient2(midpoint=0, low="#B10DC9", mid="grey",
                        high="#FFDC00", space ="Lab",limits = c(-2.2,2.2))

df$exp <- "STIA"
write.csv(df, "./output/STIA_GO_processes.csv",row.names = F)

########## Figure S3D ##############
### scatter plot of D7 vs D21 fold changes
fc.change.modified <- data.frame(D7_MA=log.cpm.mean.base[,3]-log.cpm.mean.base[,1],
                                 D7_MB=log.cpm.mean.base[,8]-log.cpm.mean.base[,6],
                                 D7_MC=log.cpm.mean.base[,13]-log.cpm.mean.base[,11],
                                 D7_MD=log.cpm.mean.base[,18]-log.cpm.mean.base[,16],
                                 D21_MA=log.cpm.mean.base[,5]-log.cpm.mean.base[,1],
                                 D21_MB=log.cpm.mean.base[,10]-log.cpm.mean.base[,6],
                                 D21_MC=log.cpm.mean.base[,15]-log.cpm.mean.base[,11],
                                 D21_MD=log.cpm.mean.base[,20]-log.cpm.mean.base[,16])
g1 <- D7_vs_D21_scatterplot(fc.change.modified,"MA")
g1+geom_abline(slope = 0.3641147, intercept = 0, color="#07693A", size = 1.5, linetype="dashed")
g2 <- D7_vs_D21_scatterplot(fc.change.modified,"MB")
g2+geom_abline(slope = 0.4443093, intercept = 0, color="#834698", size = 1.5, linetype="dashed")
g3 <- D7_vs_D21_scatterplot(fc.change.modified,"MC")
g3+geom_abline(slope = 0.1948378, intercept = 0, color="#8CC544", size = 1.5, linetype="dashed")
g4 <- D7_vs_D21_scatterplot(fc.change.modified,"MD")
g4+geom_abline(slope = 0.1866999, intercept = 0, color="#BC7EAF", size = 1.5, linetype="dashed")
gridExtra::grid.arrange(g1,g2,g3,g4)

########### Figure S3E ############
### barplots of peak, late, and both genes
D7_vs_D21_barplot(fc.change.modified, T)
D7_vs_D21_barplot(fc.change.modified, F)

######## extract upregulated peak and late genes ##
peak.late.lists <- lapply(c("MA","MB","MC","MD"),function(x) extract_peak_late_genes(fc.change, genes, x, "STIA"))
names(peak.late.lists) <- c("MA","MB","MC","MD")

MA.df <- stia.fc.change %>% filter(rownames(.) %in% c(peak.late.lists$MA$peak_up, 
                                                      peak.late.lists$MA$peak_down,))
test <- lm(D21_MA ~ D7_MA, stia.fc.change)
test$coefficients[[2]]
summary(test)$adj.r.squared


########## UpSet plots ##################
library(UpSetR)
up.input <- list(MA=c(peak.late.lists$MA$both_up,peak.late.lists$MA$peak_up),
                 MB=c(peak.late.lists$MB$both_up,peak.late.lists$MB$peak_up),
                 MC=c(peak.late.lists$MC$both_up,peak.late.lists$MC$peak_up),
                 MD=c(peak.late.lists$MD$both_up,peak.late.lists$MD$peak_up))
upset(fromList(up.input),sets=c("MD","MC","MB","MA"),
      sets.bar.color = c("#BC7EAF","#8CC544","#834698","#07693A"),
      order.by="freq",keep.order = T,text.scale = 1.5,point.size = 3,line.size = 1)

library(Vennerable)

down.input <- list(MA=c(peak.late.lists$MA$both_down,peak.late.lists$MA$peak_down),
                   MB=c(peak.late.lists$MB$both_down,peak.late.lists$MB$peak_down),
                   MC=c(peak.late.lists$MC$both_down,peak.late.lists$MC$peak_down),
                   MD=c(peak.late.lists$MD$both_down,peak.late.lists$MD$peak_down))
upset(fromList(down.input),sets=c("MD","MC","MB","MA"),
      sets.bar.color = c("#BC7EAF","#8CC544","#834698","#07693A"),
      order.by="freq",keep.order = T,text.scale = 1.5,point.size = 3,line.size = 1)


################# parse homer results ###########
##### TF motif bar plots
motif.names <- data.frame(homer_name=c(
  "Mef2a(MADS)/HL1-Mef2a.biotin-ChIP-Seq(GSE21529)/Homer",
  "Klf9(Zf)/GBM-Klf9-ChIP-Seq(GSE62211)/Homer",
  "USF1(bHLH)/GM12878-Usf1-ChIP-Seq(GSE32465)/Homer",
  "MITF(bHLH)/MastCells-MITF-ChIP-Seq(GSE48085)/Homer",
  "KLF1(Zf)/HUDEP2-KLF1-CutnRun(GSE136251)/Homer",
  "Sp5(Zf)/mES-Sp5.Flag-ChIP-Seq(GSE72989)/Homer",
  "KLF3(Zf)/MEF-Klf3-ChIP-Seq(GSE44748)/Homer",
  "c-Myc(bHLH)/LNCAP-cMyc-ChIP-Seq(Unpublished)/Homer",
  "AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer",
  "JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer"
),
motif=c("MEF2A (MADS)","KLF9 (Zf)","USF1 (bHLH)","MITF (bHLH)","KLF1 (Zf)",
        "SP5 (Zf)","KLF3 (Zf)","c-MYC (bHLH)", 'AP-1 (bZIP)', 'JUNB(bZIP)')) %>%
  arrange(homer_name)


parse_motif_results <- function(folder_name){
  path <- file.path("./output/motifs/",folder_name)
  motif_df <- read_delim(file.path(path,"knownResults.txt"),delim = "\t") %>%
    filter(`Motif Name` %in% motif.names$homer_name) %>%
    arrange(`Motif Name`) %>%
    mutate(motif=motif.names$motif,
           percent=`% of Target Sequences with Motif`,
           percent=as.numeric(sub("%", "", percent)),
           background=`% of Background Sequences with Motif`,
           background=as.numeric(sub("%", "", background)),
           log_pvalue= -(`Log P-value`),
           pvalue=`P-value`,
           cluster=folder_name) %>%
    dplyr::select(motif:cluster)
  return(motif_df)
}

celltypes <- c(paste0("c",1:6))
cluster.list <- as.list(c(celltypes, paste0(celltypes,"_genome")))
df <- lapply(cluster.list,parse_motif_results) %>% bind_rows()
background <- df %>% filter(!endsWith(cluster,"genome")) %>% group_by(motif) %>% summarise(bg.mean=mean(background))
gn.background <- df %>% filter(endsWith(cluster,"genome")) %>% group_by(motif) %>% summarise(gn.mean=mean(background))

tf_barplot <- function(df,mt,background){
  bg.pc <- background %>% filter(motif==mt) %>% pull(bg.mean)
  gn.pc <- gn.background %>% filter(motif==mt) %>% pull(gn.mean)
  return(ggplot(df %>% filter(motif==mt),aes(x=cluster,y=percent))+
           geom_bar(stat = "identity",width = 0.8)+ #ggtitle(label = mt)+
           #scale_fill_manual(values = c("#07693A","#834698","#8CC544","#BC7EAF"))+
           #scale_fill_manual(values = colors)+
           geom_hline(yintercept = bg.pc, linetype="dashed",size=1.5,color="blue")+
           geom_hline(yintercept = gn.pc, linetype="dashed",size=1.5,color="red")+
           scale_x_discrete(labels= c("I","II","III","IV","V","VI"))+
           theme_classic(base_size = 16)+ scale_fill_brewer(palette = "Set2")+
           theme(plot.title = element_text(hjust = 0.5),
                 legend.position = "none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank()))
}
df <- df %>% filter(!endsWith(cluster,"genome"))
#df$cluster <- factor(df$cluster, levels = c(paste0("c",1:10)))
p1 <- tf_barplot(df,"USF1 (bHLH)",background)
p2 <- tf_barplot(df,"Esrrb(NR)",background)
p3 <- tf_barplot(df,"TFE3(bHLH)",background)
p4 <- tf_barplot(df,"BATF(bZIP)",background)
p5 <- tf_barplot(df,"FOXK2(Forkhead)",background)
p6 <- tf_barplot(df,"IRF4(IRF)",background)
p7 <- tf_barplot(df,"HINFP(Zf)", background)
p8 <- tf_barplot(df,"Srebp1a(bHLH)", background)
gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=4)

#######  dorathea target genes ###############
data(dorothea_mm, package = "dorothea")

regulons = dorothea_mm %>%
  filter(confidence %in% c("A", "B", "C","D","E"))

######### GSVA calculation
library(GSVA)

MEF2A <- regulons %>% filter(tf == "Mef2a", mor==1) %>% pull(target)
KLF3 <- regulons %>% filter(tf =="Klf3", mor==1) %>% pull(target)
NFKB1 <- regulons %>% filter(tf =="Nfkb1", mor==1) %>% pull(target)

gs <- list(MEF2A, KLF3,NFKB1)
names(gs) <- c("MEF2A","KLF3","NFKB1")

gsva.cpm <- log.cpm %>% rownames_to_column("ensembl") %>%
  left_join(genes) %>%
  mutate(symbol = make.names(symbol, unique = TRUE)) %>% 
  select(-ensembl) %>% column_to_rownames("symbol")
gsva.es <- gsva(as.matrix(gsva.cpm), gs, method="gsva", mx.diff=T)

gsva.scores <- as.data.frame(gsva.es) %>% rownames_to_column("tf") #%>% 
pivot_longer(DOII_M3_A_1011:D21_M5_D, names_to = "names") %>% 
  left_join(sample.table) %>% mutate(condition = paste(day, celltype, sep = "_")) %>% 
  group_by(tf, condition, day, celltype) %>% 
  summarise(N=length(value),mean=mean(value),sd=sd(value),se=sd/sqrt(N))
gsva.scores$day <- factor(gsva.scores$day, levels = c("0","3","7","13","21"))

ggplot(gsva.scores %>% filter(tf=="MEF2A"), aes(x=day, y=mean, group=celltype, color=celltype, fill = celltype))+
  geom_ribbon(aes(ymin = mean-se,
                  ymax = mean+se),alpha=0.3) + 
  geom_line(size=2)+geom_point(size=3.5)+ 
  #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2,size=1.2)+
  expand_limits(y=0)+ #facet_wrap(~tf, scales = "free")+
  scale_color_manual(values = c("#07693A","#834698","#8CC544","#BC7EAF"))+
  scale_fill_manual(values = c("#07693A","#834698","#8CC544","#BC7EAF"))+
  theme_classic(base_size = 16)+ #xlab("Day")+ylab("GSVA scores")+
  theme(#axis.title = element_text(size=14,face = "bold"),
    axis.title = element_blank(),
    axis.text = element_text(size=14,face = "bold"),
    plot.title = element_text(hjust = 0.5,face = "bold"),
    legend.position = "none")
