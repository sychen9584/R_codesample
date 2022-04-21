#################### helper functions for the STIA timecourse

#### between time points fold changes quantification ########
between_timepoint_fc <- function(cpm_matrix,celltype){
  fc1 <- log2(cpm_matrix[,paste0("D3_",celltype)]+1)-log2(cpm_matrix[,paste0("D0_",celltype)]+1)
  fc2 <- log2(cpm_matrix[,paste0("D7_",celltype)]+1)-log2(cpm_matrix[,paste0("D3_",celltype)]+1)
  fc3 <- log2(cpm_matrix[,paste0("D13_",celltype)]+1)-log2(cpm_matrix[,paste0("D7_",celltype)]+1)
  fc4 <- log2(cpm_matrix[,paste0("D21_",celltype)]+1)-log2(cpm_matrix[,paste0("D13_",celltype)]+1)
  return(cbind(fc1,fc2,fc3,fc4))
}

######## line plots for individual gene expression #######
lineplot <- function(cpm ,gene, sample.table){
  ensembl <- genes %>% filter(symbol == gene) %>% pull(ensembl)
  df <- cpm[ensembl,] %>% melt() %>%
    mutate(condition=paste(sample.table$day,sample.table$celltype,sep="_"),
           day=as.factor(sample.table$day),celltype=sample.table$celltype)
  
  df_c <- df %>% group_by(condition,day,celltype) %>% 
    summarise(N=length(value),mean=mean(value),sd=sd(value),se=sd/sqrt(N))
  
  ggplot(df_c, aes(x=day, y=mean, group=celltype, color=celltype))+
    geom_line(size=2)+geom_point(size=3.5)+ggtitle(gene)+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2,size=1.2)+
    expand_limits(y=0)+
    scale_color_manual(values = c("#07693A","#834698","#8CC544","#BC7EAF"))+
    theme_classic(base_size = 16)+ xlab("Day")+ylab("CPM")+
    theme(axis.title = element_text(size=14,face = "bold"),
          axis.text = element_text(size=14,face = "bold"),
          plot.title = element_text(hjust = 0.5,face = "bold"))
}

########  compute FC between all samples compared to avg D0  ############
timepoints_d0_fc <- function(cpm_matrix,mean_matrix,pop,sample.info){
  fc <- matrix()
  conditions <- sample.info %>% dplyr::filter(celltype==pop) %>% 
    pull(condition) %>% unique() # compile celltype specific timepoints to compare
  for (i in 1:4){
    fc2=log2(cpm_matrix[,sample.info$condition==conditions[i+1]]+1)-log2(mean_matrix[,paste0("D0_",pop)]+1) 
    fc=cbind(fc,fc2)
  }
  return(fc[,-1])
}
se <- function(x) sd(x)/sqrt(length(x)) # function to calculate standard error
set_up_DE_table <- function(df){
  return(data.frame(num = c(rowMeans(df[,1:4]),rowMeans(df[,5:8]),
                            rowMeans(df[,9:12]),rowMeans(df[,13:16])),
                    se = c(apply(df[,1:4],1,se),apply(df[,5:8],1,se),
                           apply(df[,9:12],1,se),apply(df[,13:16],1,se)),
                    pop = rep(c("MA","MB","MC","MD"),4),
                    day = factor(c(rep(3,4),rep(7,4),rep(13,4),rep(21,4))))
  )
}

###### final wrapper function
DE_timepoints_barplot <- function(cpm.matrix,mean.matrix,ycutoff){
  MA_fc <- timepoints_d0_fc(cpm.matrix,mean.matrix,"MA",sample.table)
  MB_fc <- timepoints_d0_fc(cpm.matrix,mean.matrix,"MB",sample.table)
  MC_fc <- timepoints_d0_fc(cpm.matrix,mean.matrix,"MC",sample.table)
  MD_fc <- timepoints_d0_fc(cpm.matrix,mean.matrix,"MD",sample.table)
  # positive
  pos <- rbind(MA = colSums(MA_fc>=1)/nrow(MA_fc),MB = colSums(MB_fc>=1)/nrow(MA_fc),
               MC = colSums(MC_fc>=1)/nrow(MA_fc),MD = colSums(MD_fc>=1)/nrow(MA_fc))
  
  pos.table <- set_up_DE_table(pos)
  g1 <- ggplot(pos.table, aes(x=day, y=num, ymin=num-se, ymax=num+se, color=pop))+
    geom_bar(stat="identity",fill="white", position=position_dodge(),size=1.2,width=0.8)+
    geom_errorbar(position = position_dodge(0.8),size=0.8, width=0.75)+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
    scale_color_manual(values =c("#07693A","#834698","#8CC544","#BC7EAF"))+
    theme_classic(base_size = 16)+coord_cartesian(ylim = c(0,ycutoff))+
    theme(axis.text=element_text(face="bold"), 
          legend.position = "None", axis.title = element_blank())
  # negative
  neg <- rbind(MA = colSums(MA_fc<=-1)/nrow(MA_fc),MB = colSums(MB_fc<=-1)/nrow(MA_fc),
               MC = colSums(MC_fc<=-1)/nrow(MA_fc),MD = colSums(MD_fc<=-1)/nrow(MA_fc))
  neg.table <- set_up_DE_table(neg)
  g2 <- ggplot(neg.table, aes(x=day, y=num, ymin=num-se, ymax=num+se, color=pop))+
    geom_bar(stat="identity", fill="white",position=position_dodge(),size=1.2,width=0.8)+
    geom_errorbar(position = position_dodge(0.8),size=0.8, width=0.75)+
    scale_color_manual(values =c("#07693A","#834698","#8CC544","#BC7EAF"))+
    coord_flip()+ scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
    theme_classic(base_size = 16)+coord_cartesian(ylim = c(ycutoff,0))+
    theme(axis.text=element_text(face="bold"), 
          legend.position = "None", axis.title = element_blank())
  return(list(g1,g2))
}


############# scatter plot between Day 7 and Day 21 ############
D7_vs_D21_scatterplot <- function(cpm.matrix,celltype){
  df <- data.frame(gene = rownames(cpm.matrix),
                   D7= cpm.matrix[,paste0("D7_",celltype)],
                   D21 = cpm.matrix[,paste0("D21_",celltype)])
  df <- df %>% mutate(label=case_when(
    (df$D21 >= df$D7-1 & df$D21 <= df$D7+1) &
      ((df$D21 >= 1 | df$D7 >=1)|(df$D21<=-1 | df$D7 <=-1)) ~ "Both",
    (df$D21 <= df$D7-1 & df$D7>=1 & abs(df$D7)>abs(df$D21))| 
      (df$D21 >= df$D7+1 & df$D7 <=-1 & abs(df$D7)>abs(df$D21)) ~ "Peak",
    (df$D21 <= df$D7-1 & df$D21<=-1 & abs(df$D21)>abs(df$D7)) | 
      (df$D21 >= df$D7+1 & df$D21 >=1 & abs(df$D21)>abs(df$D7)) ~ "Late",
    TRUE ~ "No Change"))

  
  g <- ggplot(df, aes(x=D7, y=D21, color=label))+
    geom_point()+
    coord_cartesian(xlim = c(-6,6),ylim =c(-6,6))+
    geom_hline(yintercept=0, color="black",size=1)+
    geom_vline(xintercept=0, color="black",size=1)+
    scale_color_manual(values=c("#7851A9","#88C2C9","grey","#DD2A41"))+
    theme_bw(base_size = 16)+theme(legend.position="none")+
    theme(axis.text=element_text(face="bold"))
  
  return(g)
}

#### quantify acute and chronic genes
get_peak_both_genes <- function(cpm.matrix, celltype){
  df <- data.frame(ensembl = rownames(cpm.matrix),
                   D7= cpm.matrix[,paste0("D7_",celltype)],
                   D21 = cpm.matrix[,paste0("D21_",celltype)])
  df.count <- data.frame(celltype = celltype)
  df.count$both.up <- sum((df$D21 >= df$D7-1 & df$D21 <= df$D7+1) & (df$D21>=1 | df$D7>=1))  # both up
  df.count$peak.up <- sum((df$D21 <= df$D7-1 & df$D7 >=1) & (abs(df$D7)>abs(df$D21))) # peak up
  df.count$late.up <- sum((df$D21 >= df$D7+1 & df$D21>=1) & (abs(df$D21)>abs(df$D7))) # late up
  
  df.count$both.down <- sum((df$D21 >= df$D7-1 & df$D21 <= df$D7+1) & (df$D21<=-1 | df$D7 <=-1)) # both down
  df.count$peak.down <- sum((df$D21 >= df$D7+1 & df$D7 <=-1) & (abs(df$D7)>abs(df$D21))) # peak down
  df.count$late.down <- sum((df$D21 <= df$D7-1 & df$D21<=-1) & (abs(df$D21)>abs(df$D7))) # late down
  df.count$celltype <- celltype
  
  return(df.count)
}

D7_vs_D21_barplot <- function(cpm.matrix,is.up){
  
  MA <- get_peak_both_genes(cpm.matrix, "MA")
  MB <- get_peak_both_genes(cpm.matrix, "MB")
  MC <- get_peak_both_genes(cpm.matrix, "MC")
  MD <- get_peak_both_genes(cpm.matrix, "MD")
  df <- rbind(MA,MB,MC,MD)
  df <- melt(df)
  
  if (is.up){
    df.prop <- df %>% filter(str_ends(variable, ".up")) %>% 
      mutate(variable = factor(variable, levels = c("peak.up","both.up","late.up")))
    y_limits <- c(0,300)
  } else {
    df.prop <- df %>% filter(str_ends(variable, ".down")) %>% 
      mutate(value = -value,
             variable = factor(variable, levels = c("peak.down","both.down","late.down")))
    y_limits <- c(-300,0)
  }
  
  g <- ggplot(df.prop, aes(y=value, x=celltype,fill=variable))+
    geom_bar(stat = "identity",width = 0.6, position="dodge")+
    scale_fill_manual(values = c("#DD2A41","#7851A9","#88C2C9"))+
    coord_cartesian(ylim=y_limits)+
    guides(fill=FALSE)+
    theme_classic(base_size = 16)+ 
    theme(axis.title=element_blank(),
          legend.position="none",
          axis.text=element_text(face="bold"))
  
  return(g)
}

extract_peak_late_genes <- function(cpm.matrix, genes, celltype, exp){
  
  if (exp=="STIA"){
    
    df <- data.frame(ensembl = rownames(cpm.matrix),
                     D7= cpm.matrix[,paste0("D7_",celltype)],
                     D21 = cpm.matrix[,paste0("D21_",celltype)])
    
    df$both.up <- (df$D21 >= df$D7-1 & df$D21 <= df$D7+1) & (df$D21>=1 | df$D7>=1)  # both up
    df$peak.up <- (df$D21 <= df$D7-1 & df$D7 >=1) & (abs(df$D7)>abs(df$D21)) # peak up
    df$late.up <- (df$D21 >= df$D7+1 & df$D21>=1) & (abs(df$D21)>abs(df$D7)) # late up
    
    df$both.down <- (df$D21 >= df$D7-1 & df$D21 <= df$D7+1) & (df$D21<=-1 | df$D7 <=-1) # both down
    df$peak.down <- (df$D21 >= df$D7+1 & df$D7 <=-1) & (abs(df$D7)>abs(df$D21)) # peak down
    df$late.down <- (df$D21 <= df$D7-1 & df$D21<=-1) & (abs(df$D21)>abs(df$D7)) # late down
    
  } else if (exp=="CIA"){
    
    df <- data.frame(ensembl = rownames(cpm.matrix),
                     D41= cpm.matrix[,paste0("D41_",celltype)],
                     D62 = cpm.matrix[,paste0("D62_",celltype)])
    
    df$both.up <- (df$D62 >= df$D41-1 & df$D62 <= df$D41+1) & (df$D62>=1 | df$D41>=1)  # both up
    df$peak.up <- (df$D62 <= df$D41-1 & df$D41 >=1) & (abs(df$D41)>abs(df$D62)) # peak up
    df$late.up <- (df$D62 >= df$D41+1 & df$D62>=1) & (abs(df$D62)>abs(df$D41)) # late up
    
    df$both.down <- (df$D62 >= df$D41-1 & df$D62 <= df$D41+1) & (df$D62<=-1 | df$D41 <=-1) # both down
    df$peak.down <- (df$D62 >= df$D41+1 & df$D41 <=-1) & (abs(df$D41)>abs(df$D62)) # peak down
    df$late.down <- (df$D62 <= df$D41-1 & df$D62 <=-1) & (abs(df$D62)>abs(df$D41)) # late down
    
  }
  
  df$celltype <- celltype
  
  df <- left_join(df,genes,by="ensembl")
  both.up <- df %>% filter(both.up) %>% pull(ensembl)
  both.down <- df %>% filter(both.down) %>% pull(ensembl)
  peak.up <- df %>% filter(peak.up) %>% pull(ensembl)
  peak.down <- df %>% filter(peak.down) %>% pull(ensembl)
  late.up <- df %>% filter(late.up) %>% pull(ensembl)
  late.down <- df %>% filter(late.down) %>% pull(ensembl)
  gene.list <- list(both_up=both.up,both_down=both.down,
                    peak_up=peak.up,peak_down=peak.down,
                    late_up=late.up,late_down=late.down)
  #return(t(plyr::ldply(gene.list, rbind)))
  return(gene.list)
}

############# bubble plots ###################
z <- function(x){(x-mean(x))/sd(x)}
bubble.setup <- function(df,input,name){
  x <- subset(df, symbols %in% input)
  n <- colSums(x >= rowMeans(x))/nrow(x)
  e <- rowMeans(apply(x,1,z))
  x <- data.frame(condition = colnames(x),
                  pathway = rep(name,20),
                  e = e,n = n, num = nrow(x))
  return(x)
}

############## cluster trend line ###########
#### draw trend line of each cluster
draw_trend_line <- function(cluster.df,cluster,STIA){
  cluster <- cluster.df[cluster.df$cluster==cluster,-c(1:2)]
  cluster.means <- colMeans(cluster)
  
  if (STIA){
    pop <- c(rep("MA",5),rep("MB",5),rep("MC",5),rep("MD",5))
    day <- factor(rep(c(0,3,7,13,21),4))
  } else {
    pop <- c(rep("MA",4),rep("MB",4),rep("MC",4),rep("MD",4))
    day <- factor(rep(c(0,27,41,62),4))
  }
  
  cluster <- data.frame(days=day,celltype=pop,mean = cluster.means)
  return(ggplot(cluster, aes(x=days,y=mean,group=celltype,colour=celltype)) + 
           geom_line(size=2) + geom_point()+
           scale_color_manual(values = c("#07693A","#834698","#8CC544","#BC7EAF"))+
           coord_cartesian(ylim = c(-1.3, 1.3))+
           geom_hline(yintercept = 0, linetype="dashed",size=1.5,color="grey")+
           theme_classic(base_size = 16)+theme(legend.position="none"))
}

STIA.process.foldchange <- function(df, input){
  temp <- df %>% filter(ensembl %in% input) 
  geneNum <- nrow(temp)
  
  temp <- temp %>% select(starts_with("D")) %>% colMeans()
  
  df <- data.frame(fc = temp, condition = names(temp),
                   exp = "STIA")
  #print(df)
  return(list(df = df, num = geneNum))
}

STIA.fc.plot <- function(stia, name, geneNum){
  colors <- c("#07693A","#834698","#8CC544","#BC7EAF")
  df <- stia %>%
    separate(condition, into=c("day","celltype")) %>% 
    mutate(day = as.numeric(str_extract(.$day, "[[:digit:]]+")))
  
  plot.title <- paste0(name, " (", geneNum, ")")
  
  return(ggplot(df, aes(x=day,y=fc,color=celltype,linetype=exp, shape=exp))+
           geom_line(size=1.3)+
           geom_point(size=2.5)+
           scale_color_manual(values = colors)+
           geom_hline(yintercept = 0,color="grey50",linetype="dashed",size=1)+
           ggtitle(plot.title)+
           theme_classic(base_size = 18)+
           theme(plot.title = element_text(hjust = 0.5),
                 axis.title = element_blank(),
                 legend.position = "none"))
}