##差异分析
library(limma)
dat=read.table("TCGA.txt",sep="\t",header=T,check.names=F)    #读取输入文件
#数据处理，如果一个基因有多行，取均值
dat=as.matrix(dat)
rownames(dat)=dat[,1]
exp=dat[,2:ncol(dat)]
dimnames=list(rownames(exp),colnames(exp))
dat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
dat=avereps(dat)
gene <- read.table("gene.txt",head=F,sep='\t',check.names = F) 
expr<-dat[intersect(gene$V1,row.names(dat)),]
range(expr)
expr<-log2(expr+1)
write.table(t(expr),"Entosis_expr.txt",quote=F,sep='\t')
library(reshape2)
library(ggpubr)
library(ggsci)
df<-read.table("Entosis_expr.txt",head=T,sep='\t',check.names = F,row.names = 1)
data<-melt(df,
           id.vars = c('Group'),
           measure.vars = colnames(df[-1]),
           variable.name='Immune_check',
           value.name='Expression')
p=ggboxplot(data, x="Immune_check", y="Expression", fill = "Group",alpha=0.8,orientation = "horizontal",notch = TRUE,
            ylab="Gene_expression",
            xlab="",
            palette = c("lancet"))
p=p+rotate_x_text(60)                    
p+stat_compare_means(aes(group=Group),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),label = "p.signif")+theme(axis.text = element_text(size = 13, face = "bold"),axis.title = element_text(size = 13, face = "bold"))
##富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
pvalueFilter=0.05        #p值过滤条件
p.adjustFilter=0.05 
rt=read.table("gene.txt", header=F, sep="\t", check.names=F)     #读取输入文件

#提取交集基因的名称, 将基因名称转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$p.adjust<p.adjustFilter),]
#输出显著富集的结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

dt=read.table('GO_input.txt', header = T, check.names = F,sep='\t')
dt$Description <- factor(dt$Description, levels = rev(dt$Description))
p1 = ggplot() +
    geom_bar(data = dt,
             aes(x = -log10(p.adjust), y = Description, fill = ONTOLOGY),
             width=0.8, #柱子宽度调整
             stat='identity') +
    theme_classic()+scale_x_continuous(expand = c(0,0))+theme(axis.text.y = element_blank()) + 
    geom_text(data = dt,
              aes(x = 0.1, #用数值向量控制文本标签起始位置
                  y=Description, 
                  label=Description),
              size=4.5,
              hjust=0)
mytheme<- theme(
    legend.position = 'none',
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13,face = 'bold'),
    axis.text = element_text(size = 11,face = 'bold'),
    axis.ticks.y = element_blank())
p2 <- p1 + mytheme
p2 <- p1+labs(x = '-Log10P', 
         y= 'MF                                                              CC                                                                             BP')
mycol<- c('#6fbedd', '#d9543d',"#00A087FF")
p3<- p2 + 
    scale_fill_manual(values = mycol) +
    scale_color_manual(values = mycol)
p3


#KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)

pvalueFilter=0.05        #p值过滤条件
p.adjustFilter=0.05      #矫正后的p值过滤条件

rt=read.table("gene.txt", header=F, sep="\t", check.names=F)     #读取输入文件

#提取交集基因的名称,将基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#KEGG富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<p.adjustFilter),]
#保存显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#绘制
dt=read.table('KEGG_input.txt', header = T, check.names = F,sep='\t')
dt$Description <- factor(dt$Description, levels = rev(dt$Description))
p=ggplot() +
    geom_bar(data = dt,
             aes(x = -log10(pvalue), y = Description, fill = pvalue),
             width=0.8, #柱子宽度调整
             stat='identity') +scale_fill_continuous(type = "viridis")+
    theme_classic()
p1 <- p+scale_x_continuous(expand = c(0,0))
p2 = p1 +
    theme(axis.text.y = element_blank()) + #去掉y轴标签
    geom_text(data = dt,
              aes(x = 0.1, #用数值向量控制文本标签起始位置
                  y=Description, 
                  label=Description),
              size=6,
              hjust=0)
mytheme<- theme(
    legend.position = 'none',
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13,face = 'bold'),
    axis.text = element_text(size = 15,face = 'bold'),
    axis.ticks.y = element_blank())
p3 <- p2 + mytheme
p4 <- p3 +
    labs(x = '-Log10P', 
         y= 'KEGG')
p4

##CNV 变异
inputFile="CNVmatrix.txt"     #输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
GAIN=rowSums(rt> 0)       #拷贝数增加的样品数目
LOSS=rowSums(rt< 0)       #拷贝数缺失的样品数目
GAIN=GAIN/ncol(rt)*100      #拷贝数增加的百分率
LOSS=LOSS/ncol(rt)*100      #拷贝数缺失的百分率
data=cbind(GAIN, LOSS)
data=data[order(data[,"GAIN"],decreasing = T),]

#绘制图形
data.max = apply(data, 1, max)
pdf(file="CNVfreq.pdf", width=9, height=6)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(data.max, col="#e0f3f8", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(data.max)))
points(bar,data[,"GAIN"], pch=20, col="#E64B35FF", cex=3)
points(bar,data[,"LOSS"], pch=20, col="#00A087FF", cex=3)
legend("top", legend=c('GAIN','LOSS'), col=c("#E64B35FF","#00A087FF"), pch=20, bty="n", cex=2, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(data), adj=1)
dev.off()

#SNP变异
library(maftools)
KIRC<-read.maf("KIRC.maf")
gene<-read.table("gene.txt",head=F,sep='\t',check.names = F)
oncoplot(KIRC,genes = c("MTOR","TP53","TNFSF10","ROCK1","CTNNA1","ATG7","MYH14","PIKFYVE","DIAPH1","MCOLN1","CDH1","KIF2C","FOXO1","CDC42"))

######Circos
library(circlize)             #引用包
#读取基因位置信息文件
genepos=read.table("geneREF.txt", header=T, sep="\t", check.names=F)
colnames(genepos)=c('genename','chr','start','end')
genepos=genepos[,c('chr','start','end','genename')]
row.names(genepos)=genepos[,'genename']

#读取基因列表文件,获取特征基因的位置信息
geneRT=read.table("gene.txt", header=T, sep="\t", check.names=F)
genepos=genepos[as.vector(geneRT[,1]),]
bed0=genepos
bed0 <- bed0[-c(35,27),]
#绘制图形
pdf(file="circlize.pdf", width=6, height=6)
#初始化圈图
circos.clear()
circos.initializeWithIdeogram(species="hg38", plotType=NULL)
#展示每条染色体的注释信息
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col=rand_color(24))
    circos.text(mean(xlim), mean(ylim), chr, cex=0.6, col = "white",
                facing = "inside", niceFacing = TRUE)
}, track.height=0.15, bg.border = NA)
#绘制基因组的图形
circos.genomicIdeogram(species = "hg38", track.height=mm_h(6))
#在染色体相应位置上标注特征基因的名称
circos.genomicLabels(bed0, labels.column=4, side = "inside", cex=0.8)
circos.clear()
dev.off()

##Foresplot
dat<-read.table("uniq.symbol.txt",head=T,sep='\t',check.names = F,row.names = 1)
expr<-dat[intersect(gene$V1,row.names(dat)),]
time<-read.table("time.txt",head=T,sep='\t',check.names = F,row.names = 1)
expr<-expr[,intersect(row.names(time),colnames(expr))]
time<-time[intersect(row.names(time),colnames(expr)),]
rt<-cbind(time,t(expr))
library(survival)
pFilter=0.05  
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    if(coxP<pFilter){
        sigGenes=c(sigGenes,i)
        outTab=rbind(outTab,
                     cbind(id=i,
                           HR=coxSummary$conf.int[,"exp(coef)"],
                           HR.95L=coxSummary$conf.int[,"lower .95"],
                           HR.95H=coxSummary$conf.int[,"upper .95"],
                           pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
        )
    }
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

rt <- read.table("uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#输出图形
pdf(file="forest.pdf", width = 8,height = 6.5)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

#绘制森林图左边的基因信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#4DBBD5FF",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "#E64B35FF", "#00A087FF")
points(as.numeric(hr), n:1, pch = 20, col = boxcolor, cex=1.3)
axis(1)
dev.off()

##亚型构建
data<-read.table("Consensus_expr.txt",head=T,sep='\t',check.names = F,row.names = 1)
data<-as.matrix(expr)
maxK=9
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(data,
                               maxK=maxK,
                               reps=50,
                               pItem=0.8,
                               pFeature=1,
                               title="Cluster",
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot="pdf")
Kvec = 2:9
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
    M = results[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK
PAC <- as.data.frame(PAC)
PAC$K <- 2:9
library(ggplot2)
ggplot(PAC,aes(factor(K),PAC,group=1))+
    geom_line(color="#00A087FF",size=2)+
    theme_bw(base_rect_size = 1.5)+
    geom_point(size=6,shape=21,color="#E64B35FF",fill="#E64B35FF")+
    ggtitle('Proportion of ambiguous clustering')+
    xlab('Cluster number K')+ylab(NULL)+
    theme(axis.text = element_text(size=16,color="black",face="bold"),
          plot.title = element_text(size=16,color="black",face="bold",hjust = 0.5),
          axis.title = element_text(size=16,color="black",face="bold"))
clusterNum=2      
cluster=results[[clusterNum]][["consensusClass"]]
sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('C',sub$Cluster)
table(sub$Cluster)
head(sub)
my <- results[[2]][["ml"]]
library(pheatmap)
rownames(my) <- sub$Sample
colnames(my) <- sub$Sample
library(ggsci)
pheatmap(1-my,show_colnames = F,show_rownames = F,
         treeheight_row = 20,treeheight_col = 20,
         clustering_method = 'complete',
         color = colorRampPalette(c("navy", "white"))(50),
         annotation_names_row = F,annotation_names_col = F,
         annotation_row = data.frame(Cluster=sub$Cluster,row.names = sub$Sample),
         annotation_col = data.frame(Cluster=sub$Cluster,row.names = sub$Sample),
         annotation_colors = list(Cluster=c('C2'="#E64B35FF",'C1'="#4DBBD5FF")),legend_breaks = c(0,0.2,0.4,0.6,0.8,1),legend_labels = c(1,0.8,0.6,0.4,0.2,0))

#Survival analysis
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
    original.p <- p
    if(is.ggplot(original.p)) list.plots <- list(original.p)
    else if(is.list(original.p)) list.plots <- original.p
    else stop("Can't handle an object of class ", class (original.p))
    .set_font <- function(font){
        font <- ggpubr:::.parse_font(font)
        ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
    }
    for(i in 1:length(list.plots)){
        p <- list.plots[[i]]
        if(is.ggplot(p)){
            if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
            if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
            if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
            if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
            if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
            if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
            if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
            list.plots[[i]] <- p
        }
    }
    if(is.ggplot(original.p)) list.plots[[1]]
    else list.plots
}
library(survival)
library("survminer")
rt=read.table("Cluster.txt",header=T,sep="\t");rt$OS.time=rt$OS.time/12
diff=survdiff(Surv(OS.time,OS) ~Cluster,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(OS.time,OS) ~ Cluster, data = rt)
p<-ggsurvplot(fit, 
              data=rt,
              conf.int=T,conf.int.style='step', size=1.5,
              pval=paste0 ("P = ",pValue),
              pval.size=12,
              legend.title="Cluster",
              legend.labs=levels(factor(rt[,"Cluster"])),
              legend = c(0.9, 0.9),
              font.legend=10,
              xlab="Time(years)",
              break.time.by = 2,
              palette = c("jco"),
              surv.median.line = "hv",
              risk.table=T,
              cumevents=F,
              risk.table.height=.30)
p$table <- customize_labels(
    p$table,
    font.title    = c(16, "bold", "darkblue"),         
    font.subtitle = c(15, "bold.italic", "purple"), 
    font.caption  = c(14, "plain", "orange"),        
    font.x        = c(20, "bold", "black"),          
    font.y        = c(20, "bold", "black"),      
    font.xtickslab = c(20, "bold", "black"),font.ytickslab = c(20, "bold", "black")
)
p$plot<-customize_labels(
    p$plot,
    font.title    = c(16, "bold", "darkblue"),         
    font.subtitle = c(15, "bold.italic", "purple"), 
    font.caption  = c(14, "plain", "orange"),        
    font.x        = c(25, "bold", "black"),          
    font.y        = c(25, "bold", "black"),      
    font.xtickslab = c(20, "bold", "black"),
    font.ytickslab = c(20, "bold", "black")
)
#ComplexHeatmap
library(ComplexHeatmap)
ee<-read.table("Cluster_expr.txt",head=T,sep='\t',check.names = F,row.names = 1)
my<-read.table("Clinical.txt",head=T,sep='\t',check.names = F,row.names = 1)
ee<-t(ee)
ee<-ee[,intersect(row.names(my),colnames(ee))]
identical(row.names(my),colnames(ee))
ee <- t(scale(t(ee)))
ee[ee > 2] <- 2 
ee[ee < -2] <- -2 
my$Status <- factor(my$Status)
my$Cluster <- factor(my$Cluster)
my$Stage <- factor(my$Stage)
my$Gender <- factor(my$Gender,levels = c('Female','Male'))
my$Age <- factor(my$Age,levels = c('<=65','>65'))
my$T=factor(my$T)
my$N=factor(my$N)
my$M=factor(my$M,levels = c('M0','M1','NA'))
Cluster <- c("#00a087","#e64b35")
names(Cluster) <- levels(my$Cluster)
Age <- c(pal_nejm(alpha = 0.9)(8)[3],'#CF4E27')
names(Age) <- levels(my$Age)
Gender <- c('#E0864A','rosybrown')
names(Gender) <- levels(my$Gender)
Stage <- c('cornsilk','paleturquoise','goldenrod','firebrick')
names(Stage) <- levels(my$Stage)
T <- c('cornsilk','paleturquoise','goldenrod','firebrick')
names(T) <- levels(my$T)
N <- c('paleturquoise','goldenrod','firebrick')
names(N) <- levels(my$N)
M<- c('paleturquoise','firebrick',"white")
names(M) <- levels(my$M)
Status <- c('lavenderblush','slategray')
names(Status) <- levels(my$Status)
col<-list(Cluster,Age,Gender,Stage,T,N,M,Status)


Top = HeatmapAnnotation(Cluster=my$Cluster,
                        Age=my$Age,
                        Gender=my$Gender,
                        Stage= my$Stage,T=my$T,N=my$N,N=my$N,
                        Status = my$Status,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,title_gp = gpar(fontsize = 10,fontface = "bold"),ncol=1),
                        col=list(Cluster = Cluster,
                                 Age = Age,
                                 Gender = Gender,
                                 Stage= Stage,T=T,N=N,M=M,
                                 Status = Status
                        ),
                        show_annotation_name = TRUE,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))


Heatmap(ee,name='Z-score',
        cluster_rows = TRUE,top_annotation = Top,
        col=colorRamp2(c(-2,0,2),c('#1d6fae','white','#ef7700')),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        column_split = c(rep(1,149),rep(2,380)),
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 13), border = T,
                                  title_gp = gpar(fontsize = 13, fontface = "bold")),
        column_gap = unit(2,'mm'))

#PCA analysis
library(ggplot2)
library(ggh4x)
data<-read.table("Cluster_expr.txt",head=T,sep='\t',check.names = F,row.names = 1)
pca1 <- prcomp(data,center = TRUE,scale. = TRUE)
write.table(pca1$x[,c(1,2)],"PCA.txt",quote=F,sep='\t')
dt<-read.table("PCA.txt",head=T,sep='\t',check.names = F,row.names = 1)
mytheme <- theme(panel.grid = element_blank(),
                 panel.background = element_blank(),
                 legend.key = element_blank(),
                 legend.position = "top",
                 axis.line = element_line(colour = "grey30"),axis.text=element_text(color="black",face="bold",size=16),
                 axis.ticks.length = unit(1.8, "mm"),axis.title=element_text(color="black",face="bold",size=16),legend.text = element_text(color="black",face="bold",size=18),
                 ggh4x.axis.ticks.length.minor = rel(0.6))
ggplot(dt,aes(x=PC1,y=PC2,fill=Cluster))+
    stat_centroid(aes(xend = PC1, yend = PC2, colour = Cluster),
                  geom = "segment", crop_other = F,
                  alpha=0.3,size = 1,show.legend = F)+
    geom_point(size=3,alpha=0.7,
               color="white",shape = 21,show.legend = T)+
    scale_color_manual(name="",
                       values = c("#00a087","#e64b35"))+
    scale_fill_manual(name="",
                      values = c("#00a087","#e64b35"))+
    scale_x_continuous(expand=expansion(add = c(0.7,0.7)))+
    scale_y_continuous(expand=expansion(add = c(0.5,0.5)))+
    guides(x = "axis_truncated",y = "axis_truncated")+mytheme
	
##WGCNA 分析
library(limma)
library(gplots)
library(WGCNA)
data<-read.table("uniq.symbol.txt",head=T,sep='\t',check.names = F,row.names = 1)
datTraits<-read.table("Cluster.txt",head=T,sep='\t',check.names = F,row.names = 1)
selectGenes=names(tail(sort(apply(data,1,sd)), n=round(nrow(data)*0.25)))
data=data[selectGenes,]
datExpr0=t(data)
datExpr0<-datExpr0[intersect(row.names(datTraits),row.names(datExpr0)),]
identical(row.names(datTraits),row.names(datExpr0))
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###样品聚类
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "01.sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
###剪切线
abline(h = 20000, col="red")
dev.off()
###样品聚类,得到样品聚类的热图
sampleTree2 = hclust(dist(datExpr0), method="average")
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="02.sample_heatmap.pdf", width=12, height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5,RsquaredCut = 0.90)
pdf(file="03.scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="#00A087FF");
abline(h=0.90,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="#4DBBD5FF")
dev.off()

###邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate     #最佳power值
adjacency = adjacency(datExpr0, power = softPower)
softPower


###TOM矩阵
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

###基因聚类
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="04.gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


###动态模块的识别
minModuleSize = 100      #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="05.Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


###对模块进行聚类,找出模块之间的相似性

MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="06.Clustering_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()

#绘制模块基因的热图
moduleColors=dynamicColors
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
select = sample(nGenes, size=1000)      #随机选择基因进行可视化
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method="average")
selectColors = moduleColors[select]
#sizeGrWindow(9,9)
plotDiss=selectTOM^softPower
diag(plotDiss)=NA
myheatcol = colorpanel(250, "red", "orange", "lemonchiffon")    #设置热图颜色(白色背景）
pdf(file="07.TOMplot.pdf", width=7, height=7)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col=myheatcol)
dev.off()


###模块与性状数据的热图
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="08.Module_trait.pdf", width=6.5, height=5.5)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3.5, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


###计算MM和GS值
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

###输出每个模块的散点图
trait="Cluster"
traitColumn=match(trait,traitNames)  
for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
        outPdf=paste("09.", trait, "_", module,".pdf",sep="")
        pdf(file=outPdf,width=7,height=7)
        par(mfrow = c(1,1))
        verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                           abs(geneTraitSignificance[moduleGenes, traitColumn]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = paste("Gene significance for ",trait),
                           main = paste("Module membership vs. gene significance\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
        dev.off()
    }
}
###输出GS_MM的表格
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                           GSPvalue[, Tra])
    names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                         names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                           MMPvalue[, mod])
    names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                         names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls",sep="\t",row.names=F)

###亚型功能评估
library(ComplexHeatmap)
library(ggsci)
library(circlize)
ee<-read.table("ssgsea.txt",head=T,sep='\t',check.names = F,row.names = 1)
my<-read.table("Clinical.txt",head=T,sep='\t',check.names = F,row.names = 1)
ee<-ee[,intersect(row.names(my),colnames(ee))]
identical(row.names(my),colnames(ee))
ee <- t(scale(t(ee)))
ee[ee > 2] <- 2 
ee[ee < -2] <- -2 
my$Status <- factor(my$Status)
my$Cluster <- factor(my$Cluster)
my$Stage <- factor(my$Stage)
my$Gender <- factor(my$Gender,levels = c('Female','Male'))
my$Age <- factor(my$Age,levels = c('<=65','>65'))
my$T=factor(my$T)
my$N=factor(my$N)
my$M=factor(my$M,levels = c('M0','M1','NA'))
Cluster <- c("#00a087","#e64b35")
names(Cluster) <- levels(my$Cluster)
Age <- c(pal_nejm(alpha = 0.9)(8)[3],'#CF4E27')
names(Age) <- levels(my$Age)
Gender <- c('#E0864A','rosybrown')
names(Gender) <- levels(my$Gender)
Stage <- c('cornsilk','paleturquoise','goldenrod','firebrick')
names(Stage) <- levels(my$Stage)
T <- c('cornsilk','paleturquoise','goldenrod','firebrick')
names(T) <- levels(my$T)
N <- c('paleturquoise','goldenrod','firebrick')
names(N) <- levels(my$N)
M<- c('paleturquoise','firebrick',"white")
names(M) <- levels(my$M)
Status <- c('lavenderblush','slategray')
names(Status) <- levels(my$Status)
col<-list(Cluster,Age,Gender,Stage,T,N,M,Status)

Top = HeatmapAnnotation(Cluster=my$Cluster,
                        Age=my$Age,
                        Gender=my$Gender,
                        Stage= my$Stage,T=my$T,N=my$N,N=my$N,
                        Status = my$Status,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,title_gp = gpar(fontsize = 10,fontface = "bold"),ncol=1),
                        col=list(Cluster = Cluster,
                                 Age = Age,
                                 Gender = Gender,
                                 Stage= Stage,T=T,N=N,M=M,
                                 Status = Status
                        ),
                        show_annotation_name = TRUE,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))
Heatmap(ee,name='Z-score',
        cluster_rows = TRUE,top_annotation = Top,
        col=colorRamp2(c(-2,0,2),c('#1d6fae','white','#ef7700')),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        column_split = c(rep(1,149),rep(2,380)),
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 13), border = T,
                                  title_gp = gpar(fontsize = 13, fontface = "bold")),
        column_gap = unit(2,'mm'))

####免疫细胞浸润
library(ComplexHeatmap)
library(ggsci)
library(circlize)
ee<-read.table("immune_cell.txt",head=T,sep='\t',check.names = F,row.names = 1)
my<-read.table("Clinical.txt",head=T,sep='\t',check.names = F,row.names = 1)
ee <- t(ee)
ee<-ee[,intersect(row.names(my),colnames(ee))]
identical(row.names(my),colnames(ee))
ee <- t(scale(t(ee)))
ee[ee > 2] <- 2 
ee[ee < -2] <- -2 
my$Status <- factor(my$Status)
my$Cluster <- factor(my$Cluster)
my$Stage <- factor(my$Stage)
my$Gender <- factor(my$Gender,levels = c('Female','Male'))
my$Age <- factor(my$Age,levels = c('<=65','>65'))
my$T=factor(my$T)
my$N=factor(my$N)
my$M=factor(my$M,levels = c('M0','M1','NA'))
Cluster <- c("#00a087","#e64b35")
names(Cluster) <- levels(my$Cluster)
Age <- c(pal_nejm(alpha = 0.9)(8)[3],'#CF4E27')
names(Age) <- levels(my$Age)
Gender <- c('#E0864A','rosybrown')
names(Gender) <- levels(my$Gender)
Stage <- c('cornsilk','paleturquoise','goldenrod','firebrick')
names(Stage) <- levels(my$Stage)
T <- c('cornsilk','paleturquoise','goldenrod','firebrick')
names(T) <- levels(my$T)
N <- c('paleturquoise','goldenrod','firebrick')
names(N) <- levels(my$N)
M<- c('paleturquoise','firebrick',"white")
names(M) <- levels(my$M)
Status <- c('lavenderblush','slategray')
names(Status) <- levels(my$Status)
col<-list(Cluster,Age,Gender,Stage,T,N,M,Status)

Top = HeatmapAnnotation(Cluster=my$Cluster,
                        Age=my$Age,
                        Gender=my$Gender,
                        Stage= my$Stage,T=my$T,N=my$N,N=my$N,
                        Status = my$Status,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,title_gp = gpar(fontsize = 10,fontface = "bold"),ncol=1),
                        col=list(Cluster = Cluster,
                                 Age = Age,
                                 Gender = Gender,
                                 Stage= Stage,T=T,N=N,M=M,
                                 Status = Status
                        ),
                        show_annotation_name = TRUE,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))
Heatmap(ee,name='Z-score',
        cluster_rows = TRUE,top_annotation = Top,
        col=colorRamp2(c(-2,0,2),c('#1d6fae','white','#ef7700')),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        column_split = c(rep(1,149),rep(2,380)),
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 13), border = T,
                                  title_gp = gpar(fontsize = 13, fontface = "bold")),
        column_gap = unit(2,'mm'))

###临床与亚型相关性
library(ggplot2)
library(ggsci)
data<-read.table("Clinical.txt",head=T,sep='\t',check.names = F,row.names = 1)

#绘制柱状图
for(i in colnames(data)[1:(ncol(data)-1)]){
    rt=data[,c("Cluster",i)]
    colnames(rt)=c("Cluster","type")
    #统计检验
    tableStat=table(rt)
    stat=chisq.test(tableStat)
    pvalue=stat$p.value
    if(pvalue<0.001){
        pvalue="p<0.001"
    }else{
        pvalue=paste0("p=",sprintf("%.03f",pvalue))
    }
    #绘制柱状图
    p=ggplot(rt, aes(Cluster)) + 
        geom_bar(aes(fill=type), position="fill")+scale_fill_npg()+
        labs(x = 'Cluster',y = '',title=paste0(i," (",pvalue,")"))+
        guides(fill = guide_legend(title =i))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5))
    pdf(file=paste0("cliCor.",i,".pdf"),width=4.5,height=6)	
    print(p)
    dev.off()
}
#功能富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
pvalueFilter=0.05        #p值过滤条件
p.adjustFilter=0.05 
rt=read.table("gene.txt", header=F, sep="\t", check.names=F)     #读取输入文件

#提取交集基因的名称, 将基因名称转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$p.adjust<p.adjustFilter),]
#输出显著富集的结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

dt=read.table('GO_input.txt', header = T, check.names = F,sep='\t')
dt$Description <- factor(dt$Description, levels = rev(dt$Description))
p1 = ggplot() +
    geom_bar(data = dt,
             aes(x = -log10(p.adjust), y = Description, fill = ONTOLOGY),
             width=0.8, #柱子宽度调整
             stat='identity') +
    theme_classic()+scale_x_continuous(expand = c(0,0))+theme(axis.text.y = element_blank()) + 
    geom_text(data = dt,
              aes(x = 0.1, #用数值向量控制文本标签起始位置
                  y=Description, 
                  label=Description),
              size=4.5,
              hjust=0)
mytheme<- theme(
    legend.position = 'none',
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13,face = 'bold'),
    axis.text = element_text(size = 11,face = 'bold'),
    axis.ticks.y = element_blank())
p2 <- p1 + mytheme
p2 <- p1+labs(x = '-Log10P', 
         y= 'MF                                                              CC                                                                             BP')
mycol<- c('#6fbedd', '#d9543d',"#00A087FF")
p3<- p2 + 
    scale_fill_manual(values = mycol) +
    scale_color_manual(values = mycol)
p3


#KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)

pvalueFilter=0.05        #p值过滤条件
p.adjustFilter=0.05      #矫正后的p值过滤条件

rt=read.table("gene.txt", header=F, sep="\t", check.names=F)     #读取输入文件

#提取交集基因的名称,将基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#KEGG富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<p.adjustFilter),]
#保存显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#绘制
dt=read.table('KEGG_input.txt', header = T, check.names = F,sep='\t')
dt$Description <- factor(dt$Description, levels = rev(dt$Description))
p=ggplot() +
    geom_bar(data = dt,
             aes(x = -log10(pvalue), y = Description, fill = pvalue),
             width=0.8, #柱子宽度调整
             stat='identity') +scale_fill_continuous(type = "viridis")+
    theme_classic()
p1 <- p+scale_x_continuous(expand = c(0,0))
p2 = p1 +
    theme(axis.text.y = element_blank()) + #去掉y轴标签
    geom_text(data = dt,
              aes(x = 0.1, #用数值向量控制文本标签起始位置
                  y=Description, 
                  label=Description),
              size=6,
              hjust=0)
mytheme<- theme(
    legend.position = 'none',
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13,face = 'bold'),
    axis.text = element_text(size = 15,face = 'bold'),
    axis.ticks.y = element_blank())
p3 <- p2 + mytheme
p4 <- p3 +
    labs(x = '-Log10P', 
         y= 'KEGG')
p4


####机器学习构建预后模型
##先在TCGA数据集筛选预后相关基因
dat<-read.table("uniq.symbol.txt",head=T,sep='\t',check.names = F,row.names = 1)
time<-read.table("time.txt",head=T,sep='\t',check.names = F,row.names = 1)
gene<-read.table("Intersect_gene.txt",head=F,sep='\t',check.names = F)
expr<-dat[intersect(gene$V1,row.names(dat)),]
expr<-expr[,intersect(row.names(time),colnames(expr))]
time<-time[intersect(row.names(time),colnames(expr)),]
rt<-cbind(time,t(expr))
library(survival)
pFilter=0.05  
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    if(coxP<pFilter){
        sigGenes=c(sigGenes,i)
        outTab=rbind(outTab,
                     cbind(id=i,
                           HR=coxSummary$conf.int[,"exp(coef)"],
                           HR.95L=coxSummary$conf.int[,"lower .95"],
                           HR.95H=coxSummary$conf.int[,"upper .95"],
                           pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
        )
    }
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
uniSigExp$futime=uniSigExp$futime*30
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)
#EMTB-1980
dat<-read.table("E-MATB-1980.txt",head=T,sep='\t',check.names = F,row.names = 1)
time<-read.table("time.txt",head=T,sep='\t',check.names = F,row.names = 1)
gene<-read.table("Intersect_gene.txt",head=F,sep='\t',check.names = F)
expr<-dat[intersect(gene$V1,row.names(dat)),]
expr<-expr[,intersect(row.names(time),colnames(expr))]
time<-time[intersect(row.names(time),colnames(expr)),]
write.table(expr,"MATB_expr.txt",quote=F,sep='\t')
write.table(time,"MATB_time.txt",quote=F,sep='\t')
#GSE29609
dat<-read.table("GSE29609.txt",head=T,sep='\t',check.names = F,row.names = 1,quote='')
time<-read.table("time.txt",head=T,sep='\t',check.names = F,row.names = 1)
gene<-read.table("Intersect_gene.txt",head=F,sep='\t',check.names = F)
expr<-dat[intersect(gene$V1,row.names(dat)),]
expr<-expr[,intersect(row.names(time),colnames(expr))]
time<-time[intersect(row.names(time),colnames(expr)),]
write.table(expr,"GSE29609_expr.txt",quote=F,sep='\t')
write.table(time,"GSE29609_time.txt",quote=F,sep='\t')


###开始机器学习

work.path <- "I:/肾透明细胞癌-Entosis细胞侵入性死亡/15.机器学习构建预后模型"; setwd(work.path) 
# 设置其他路径
code.path <- file.path(work.path, "Codes") # 存放脚本
data.path <- file.path(work.path, "InputData") # 存在输入数据（需用户修改）
res.path <- file.path(work.path, "Results") # 存放输出结果
fig.path <- file.path(work.path, "Figures") # 存放输出图片

# 如不存在这些路径则创建路径
if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)

source(file.path(code.path, "ML.R"))

# 选择最后生成的模型类型：panML代表生成由不同算法构建的模型； multiCox表示抽取其他模型所用到的变量并建立多变量cox模型
FinalModel <- c("panML", "multiCox")[2]

## Training Cohort ---------------------------------------------------------
# 训练集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与测试集保持相同类型，表达谱需有一定变异性，以免建模过程报错）
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 训练集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]

## Validation Cohort -------------------------------------------------------
# 测试集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与训练集保持相同类型）
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 测试集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

# 提取相同基因
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因

# 按队列对数据分别进行标准化（根据情况调整centerFlags和scaleFlags）
Train_set = scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
names(x = split(as.data.frame(Test_expr), f = Test_surv$Cohort)) # 注意测试集标准化顺序与此一致
Test_expr1 <- Test_expr[-c(102:140),]
Test_surv1 <- Test_surv[-c(102:140),]
Test_set1 = scaleData(data = Test_expr1, cohort = Test_surv1$Cohort, centerFlags = T, scaleFlags = T)
Test_set<-rbind(Test_set1,Test_expr[c(102:140),])
identical(row.names(Test_set2),row.names(Test_expr))
# 目前仅有StepCox和RunEnet支持输入算法参数
methods <- read.xlsx(file.path(code.path, "41467_2022_28421_MOESM4_ESM.xlsx"), startRow = 2)$Model
methods <- gsub("-| ", "", methods)

## Train the model --------------------------------------------------------
min.selected.var <- 5 # 筛选变量数目的最小阈值
timeVar = "OS.time"; statusVar = "OS" # 定义需要考虑的结局事件，必须出现在Train_surv以及Test_surv中

## Pre-training 
Variable = colnames(Train_expr)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))
preTrain.method

set.seed(seed = 123) # 设置建模种子，使得结果可重复
preTrain.var <- list()
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, # 机器学习方法
                                 Train_expr = Train_set, # 训练集有潜在预测价值的变量
                                 Train_surv = Train_surv, # 训练集生存数据
                                 mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                 classVar = classVar) # 用于训练的生存变量，必须出现在Train_surv中
}
preTrain.var[["simple"]] <- colnames(Train_expr)

model <- list() # 初始化模型结果列表
set.seed(seed = 123) # 设置建模种子，使得结果可重复
for (method in methods){ # 循环每一种方法组合
  # method <- "CoxBoost+plsRcox" # [举例]若遇到报错，请勿直接重头运行，可给method赋值为当前报错的算法来debug
  cat(match(method, methods), ":", method, "\n") # 输出当前方法
  method_name = method # 本轮算法名称
  method <- strsplit(method, "\\+")[[1]] # 各步骤算法名称
  
  if (length(method) == 1) method <- c("simple", method)
  
  selected.var = preTrain.var[[method[1]]]
  # 如果筛选出的变量小于阈值，则该算法组合无意义，置空（尤其针对以RSF筛选变量的情况，需在ML脚本中尝试调参）
  if (length(selected.var) <= min.selected.var) {
    model[[method_name]] <- NULL
  } else {
    model[[method_name]] <- RunML(method = method[2], # 用于构建最终模型的机器学习方法
                                  Train_expr = Train_expr[, selected.var], # 训练集有潜在预测价值的变量
                                  Train_surv = Train_surv, # 训练集生存数据
                                  mode = "Model",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                  classVar = classVar)  # 用于训练的生存变量，必须出现在Train_surv中
  }
  
  # 如果最终筛选出的变量小于阈值，则该算法组合也无意义，置空
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
saveRDS(model, file.path(res.path, "model.rds")) # 保存所有模型输出

# 当要求最终模型为多变量cox时，对模型进行更新
if (FinalModel == "multiCox"){
  coxmodel <- lapply(model, function(fit){ # 根据各算法最终获得的变量，构建多变量cox模型，从而以cox回归系数和特征表达计算单样本风险得分
    tmp <- coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                 data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
    return(tmp)
  })
}
saveRDS(coxmodel, file.path(res.path, "coxmodel.rds")) # 保存最终以多变量cox拟合所筛选变量的模型

## Evaluate the model -----------------------------------------------------

# 读取已保存的模型列表（请根据需要调整）
model <- readRDS(file.path(res.path, "model.rds")) # 若希望使用各自模型的线性组合函数计算得分，请运行此行
# model <- readRDS(file.path(res.path, "coxmodel.rds")) # 若希望使用多变量cox模型计算得分，请运行此行

methodsValid <- names(model) # 取出有效的模型（变量数目小于阈值的模型视为无效）

# 根据给定表达量计算样本风险评分
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]], 
                                    new_data = rbind.data.frame(Train_set,Test_set), # 4.0更新
                                    type = "lp") # 同原文，使用linear Predictor计算得分

}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) # 输出风险评分文件

# 提取所筛选的变量（列表格式）
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]]) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
}

# 提取所筛选的变量（数据框格式）
fea_df <- lapply(model, function(fit){ data.frame(ExtractVar(fit)) }) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"  # 数据框有两列，包含算法以及算法所筛选出的变量
write.table(fea_df, file.path(res.path, "fea_df.txt"),sep = "\t", row.names = F, col.names = T, quote = F)

# 对各模型计算C-index
Cindexlist <- list()
for (method in methodsValid){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], # 预后模型
                                  Test_expr = Test_set, # 测试集预后变量，应当包含训练集中所有的变量，否则会报错
                                  Test_surv = Test_surv, # 训练集生存数据，应当包含训练集中所有的变量，否则会报错
                                  Train_expr = Train_set, # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
                                  Train_surv = Train_surv, # 若需要同时评估训练集，则给出训练集生存数据，否则置NULL
                                  Train_name = "TCGA", # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
                                  #Train_expr = NULL,
                                  #Train_surv = NULL, 
                                  cohortVar = "Cohort", # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
                                  timeVar = timeVar, # 用于评估的生存时间，必须出现在Test_surv中；这里是OS.time
                                  statusVar = statusVar) # 用于评估的生存状态，必须出现在Test_surv中；这里是OS
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

Cindex_mat <- read.table(file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T) # 计算每种算法在所有队列中平均C-index，并降序排列
Cindex_mat <- Cindex_mat[names(avg_Cindex), ] # 对C-index矩阵排序
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数
fea_sel <- fea_list[[rownames(Cindex_mat)[1]]] # 最优模型（即测试集[或者训练集+测试集]C指数均值最大）所筛选的特征

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") # 设置绘图时的队列颜色
names(CohortCol) <- colnames(Cindex_mat)

# 调用简易绘图函数
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, # 主矩阵
                    avg_Cindex = avg_Cindex, # 侧边柱状图
                    CohortCol = CohortCol, # 列标签颜色
                    barCol = "steelblue", # 右侧柱状图颜色
                    col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 热图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类

pdf(file.path(fig.path, "heatmap of cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 3, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") # 热图注释均放在右侧
invisible(dev.off())

###预后分析
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
    original.p <- p
    if(is.ggplot(original.p)) list.plots <- list(original.p)
    else if(is.list(original.p)) list.plots <- original.p
    else stop("Can't handle an object of class ", class (original.p))
    .set_font <- function(font){
        font <- ggpubr:::.parse_font(font)
        ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
    }
    for(i in 1:length(list.plots)){
        p <- list.plots[[i]]
        if(is.ggplot(p)){
            if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
            if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
            if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
            if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
            if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
            if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
            if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
            list.plots[[i]] <- p
        }
    }
    if(is.ggplot(original.p)) list.plots[[1]]
    else list.plots
}
library(survminer)
library(survival)
rt<-read.table("Train_risk.txt",head=T,sep='\t',check.names = F,row.names = 1)
res.cut <- surv_cutpoint(rt, time = "OS.time", event = "OS", variables = "riskscore",minprop = 0.5)
cutoff <- res.cut$cutpoint[1, "cutpoint"]
rt$risk <- ifelse(rt$riskscore > cutoff, "high", "low")
rt$OS.time= rt$OS.time/365

diff=survdiff(Surv(OS.time,OS) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(OS.time,OS) ~ risk, data = rt)
p<-ggsurvplot(fit, 
              data=rt,
              conf.int=T,conf.int.style='step', size=1.5,
              pval=paste0 ("P = ",pValue),
              pval.size=12,
              legend.title="risk",
              legend.labs=levels(factor(rt[,"risk"])),
              legend = c(0.9, 0.9),
              font.legend=10,
              xlab="Time(years)",
              break.time.by = 2,
              palette = c("jco"),
              surv.median.line = "hv",
              risk.table=T,
              cumevents=F,
              risk.table.height=.30)
p$table <- customize_labels(
    p$table,
    font.title    = c(16, "bold", "darkblue"),         
    font.subtitle = c(15, "bold.italic", "purple"), 
    font.caption  = c(14, "plain", "orange"),        
    font.x        = c(20, "bold", "black"),          
    font.y        = c(20, "bold", "black"),      
    font.xtickslab = c(20, "bold", "black"),font.ytickslab = c(20, "bold", "black")
)
p$plot<-customize_labels(
    p$plot,
    font.title    = c(16, "bold", "darkblue"),         
    font.subtitle = c(15, "bold.italic", "purple"), 
    font.caption  = c(14, "plain", "orange"),        
    font.x        = c(25, "bold", "black"),          
    font.y        = c(25, "bold", "black"),      
    font.xtickslab = c(20, "bold", "black"),
    font.ytickslab = c(20, "bold", "black")
)
#ROC analysis
library(timeROC)
library(survival)
rt<-read.table("Train_risk.txt",head=T,sep='\t',check.names = F,row.names = 1)
time_roc_res <- timeROC(
  T = rt$OS.time,
  delta = rt$OS,
  marker = rt$riskscore,
  cause = 1,
  weighting="marginal",
  times = c(1*365, 3*365, 5*365),
  ROC = TRUE,
  iid = TRUE)
time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3])

library(ggplot2)
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF") +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF") +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF") +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(axis.text = element_text(face = "bold", size = 14, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0))))

###预后独立性评估
library(survival)
rt=read.table("Independence.txt",header=T,sep="\t",check.names=F,row.names=1)            #读取输入文件

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
	 cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
	 coxSummary = summary(cox)
	 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	 outTab=rbind(outTab,
	              cbind(id=i,
	              HR=coxSummary$conf.int[,"exp(coef)"],
	              HR.95L=coxSummary$conf.int[,"lower .95"],
	              HR.95H=coxSummary$conf.int[,"upper .95"],
	              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
	              )
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

######绘制森林图######
#读取输入文件
rt <- read.table("uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#输出图形
pdf(file="Unicox_forest.pdf", width = 7,height = 4)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

#绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00A087FF",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "#E64B35FF", "#4DBBD5FF")
points(as.numeric(hr), n:1, pch = 16, col = boxcolor, cex=1.3)
axis(1)
dev.off()
#多因素cox分析
library(survival)
rt=read.table("Multi_independence.txt",header=T,sep="\t",check.names=F,row.names=1)

multiCox=coxph(Surv(OS.time, OS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

######绘制森林图######
#读取输入文件
rt <- read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#输出图形
pdf(file="Multi_forest.pdf", width = 7,height = 4)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

#绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00A087FF",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "#E64B35FF", "#4DBBD5FF")
points(as.numeric(hr), n:1, pch = 16, col = boxcolor, cex=1.3)
axis(1)
dev.off()

#Nomogram构建
library(survival)
library(regplot)
library(rms)
library(ggsci)
rt<-read.table("Nomogram.txt",head=T,sep='\t',check.names = F,row.names = 1)
res.cox=coxph(Surv(OS.time, OS) ~ . , data = rt)
nom1=regplot(res.cox,
              plots = c("density", "boxes"),
              clickable=F,
              title="",
              points=TRUE,
              droplines=TRUE,
              observation=rt[9,],
              rank="sd",
              failtime = c(365,1095,1825),
              prfail = F)
nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(rt, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

#Calberation
#校准曲线
pdf(file="calibration.pdf", width=5, height=5)
#1年校准曲线
f <- cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=365)
cal <- calibrate(f, cmethod="KM", method="boot", u=365, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="#E64B35FF", sub=F)
#3年校准曲线
f <- cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1095)
cal <- calibrate(f, cmethod="KM", method="boot", u=1095, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="#4DBBD5FF", sub=F, add=T)
#5年校准曲线
f <- cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1825)
cal <- calibrate(f, cmethod="KM", method="boot", u=1825, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="#00A087FF", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("#E64B35FF", "#4DBBD5FF", "#00A087FF"), lwd=1.5, bty = 'n')
dev.off()
#DCA curve
library(survival)
library(ggDCA)
rt<- read.table("OS_Independence.txt",head=T,sep='\t',check.names = F,row.names = 1)
#DCA分析
predictTime=12    #预测时间
Risk<-coxph(Surv(OS.time,OS)~riskscore,rt)
Age<-coxph(Surv(OS.time,OS)~Age,rt)
Gender<-coxph(Surv(OS.time,OS)~Gender,rt)
Grade<-coxph(Surv(OS.time,OS)~Grade,rt)
Stage<-coxph(Surv(OS.time,OS)~Stage,rt)

#绘制决策曲线
pdf(file="1_DCA.pdf", width=6.5, height=5.2)
d_train=dca(Risk,Age,Gender,Grade,Stage, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()

#3年DCA
predictTime=36    #预测时间
Risk<-coxph(Surv(OS.time,OS)~riskscore,rt)
Age<-coxph(Surv(OS.time,OS)~Age,rt)
Gender<-coxph(Surv(OS.time,OS)~Gender,rt)
Grade<-coxph(Surv(OS.time,OS)~Grade,rt)
Stage<-coxph(Surv(OS.time,OS)~Stage,rt)

#绘制决策曲线
pdf(file="3_DCA.pdf", width=6.5, height=5.2)
d_train=dca(Risk,Age,Gender,Grade,Stage, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()

#5年DCA
predictTime=60    #预测时间
Risk<-coxph(Surv(OS.time,OS)~riskscore,rt)
Age<-coxph(Surv(OS.time,OS)~Age,rt)
Gender<-coxph(Surv(OS.time,OS)~Gender,rt)
Grade<-coxph(Surv(OS.time,OS)~Grade,rt)
Stage<-coxph(Surv(OS.time,OS)~Stage,rt)

#绘制决策曲线
pdf(file="5_DCA.pdf", width=6.5, height=5.2)
d_train=dca(Risk,Age,Gender,Grade,Stage, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()

#靶向药物筛选
library(tidyverse) # 用于读取MAF文件
library(ISOpureR) # 用于纯化表达谱
library(impute) # 用于KNN填补药敏数据
library(pRRophetic) # 用于药敏预测
library(SimDesign) # 用于禁止药敏预测过程输出的信息
library(ggplot2) # 绘图
library(cowplot) # 合并图像

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}  

dat<-read.table("merge.txt",head=T,sep='\t',check.names = F,row.names = 1)
Sinfo<-read.table("TCGA_risk.txt",head=T,sep='\t',check.names = F,row.names = 1)
sample<-intersect(row.names(Sinfo),colnames(dat))
dat<-dat[,sample]
tumoexpr =dat[rowMeans(dat)>0.5,]
runpure <- F # 如果想运行就把这个改为T
if(runpure) {
    set.seed(123)
    # Run ISOpureR Step 1 - Cancer Profile Estimation
    ISOpureS1model <- ISOpure.step1.CPE(tumoexpr, normexpr)
    # For reproducible results, set the random seed
    set.seed(456);
    # Run ISOpureR Step 2 - Patient Profile Estimation
    ISOpureS2model <- ISOpure.step2.PPE(tumoexpr,normexpr,ISOpureS1model)
    pure.tumoexpr <- ISOpureS2model$cc_cancerprofiles
}

if(!runpure) {
    pure.tumoexpr <- tumoexpr
}
auc <- read.table("CTRP_AUC_raw.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 3
auc$comb <- paste(auc$master_cpd_id,auc$master_ccl_id,sep = "-")
auc <- apply(auc[,"area_under_curve",drop = F], 2, function(x) tapply(x, INDEX=factor(auc$comb), FUN=max, na.rm=TRUE)) # 重复项取最大AUC
auc <- as.data.frame(auc)
auc$master_cpd_id <- sapply(strsplit(rownames(auc),"-",fixed = T),"[",1)
auc$master_ccl_id <- sapply(strsplit(rownames(auc),"-",fixed = T),"[",2)
auc <- reshape(auc, 
               direction = "wide",
               timevar = "master_cpd_id",
               idvar = "master_ccl_id")
colnames(auc) <- gsub("area_under_curve.","",colnames(auc),fixed = T)
ctrp.ccl.anno <- read.table("CTRP_ccl_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
ctrp.cpd.anno <- read.table("CTRP_cpd_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T,quote="")
write.table(auc,"CTRP_AUC.txt",sep = "\t",row.names = F,col.names = T,quote = F)

# 2.加载药敏AUC矩阵并进行数据处理
ctrp.auc <- read.table("CTRP_AUC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
prism.auc <- read.delim("PRISM_AUC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) # 数据来自https://depmap.org/portal/download/ Drug sensitivity AUC (PRISM Repurposing Secondary Screen) 19Q4
prism.ccl.anno <- prism.auc[,1:5] # 前5列为细胞系注释信息
prism.auc <- prism.auc[,-c(1:5)]

## a. 移除缺失值大于20%的药物
ctrp.auc <- ctrp.auc[,apply(ctrp.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(ctrp.auc)]
prism.auc <- prism.auc[,apply(prism.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(prism.auc)]

## b. 移除CTRP数据里源自haematopoietic_and_lymphoid_tissue的细胞系

rmccl <- paste0("CCL", na.omit(ctrp.ccl.anno[which(ctrp.ccl.anno$ccle_primary_site != "bone"), "master_ccl_id"]))
rownames(ctrp.auc) <- paste0("CCL",rownames(ctrp.auc))
ctrp.auc <- ctrp.auc[setdiff(rownames(ctrp.auc),rmccl),]

## c. KNN填补缺失值
ctrp.auc.knn <- impute.knn(as.matrix(ctrp.auc))$data

prism.auc.knn <- impute.knn(as.matrix(prism.auc))$data

## d. 数据量级修正（与作者沟通得知）
ctrp.auc.knn <- ctrp.auc.knn/ceiling(max(ctrp.auc.knn)) # 参考Expression Levels of Therapeutic Targets as Indicators of Sensitivity to Targeted Therapeutics (2019, Molecular Cancer Therapeutics)
prism.auc.knn <- prism.auc.knn/ceiling(max(prism.auc.knn))

#药敏预测
# 加载CCLE细胞系的表达谱，作为训练集
ccl.expr <- read.table("CCLE_RNAseq_rsem_genes_tpm_20180929.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 

# 加载基因注释文件，用于基因ID转换
Ginfo <- read.table("overlapTable27.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) # 参考FigureYa34count2FPKMv2制作的基因注释文件

# 把基因的ensembl ID转换为gene symbol
ccl.expr <- ccl.expr[,-1]; rownames(ccl.expr) <- sapply(strsplit(rownames(ccl.expr),".",fixed = T),"[",1)
comgene <- intersect(rownames(ccl.expr),rownames(Ginfo))
ccl.expr <- ccl.expr[comgene,]
ccl.expr$gene <- Ginfo[comgene,"genename"]; ccl.expr <- ccl.expr[!duplicated(ccl.expr$gene),]; rownames(ccl.expr) <- ccl.expr$gene; ccl.expr <- ccl.expr[,-ncol(ccl.expr)]
keepgene <- apply(ccl.expr, 1, mad) > 0.5 # 保留表达值有效的基因
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1) # 重置细胞系名
trainPtype <- as.data.frame(ctrp.auc.knn)
ccl.name <- ccl.miss <- c() # 替换细胞系名
for (i in rownames(trainPtype)) {
    if(!is.element(gsub("CCL","",i),ctrp.ccl.anno$master_ccl_id)) {
        cat(i,"\n")
        ccl.miss <- c(ccl.miss, i) # 没有匹配到的细胞系
        ccl.name <- c(ccl.name, i) # 插入未匹配的细胞系
    } else {
        ccl.name <- c(ccl.name,  ctrp.ccl.anno[which(ctrp.ccl.anno$master_ccl_id == gsub("CCL","",i)),"ccl_name"]) # 插入匹配的细胞系
    }
}

cpd.name <- cpd.miss <- c() # 替换药物名
for (i in colnames(trainPtype)) {
    if(!is.element(i,ctrp.cpd.anno$master_cpd_id)) {
        cat(i,"\n")
        cpd.miss <- c(cpd.miss, i) # 没有匹配到的药物
        cpd.name <- c(cpd.name, i) # 插入未匹配的药物
    } else {
        cpd.name <- c(cpd.name,  ctrp.cpd.anno[which(ctrp.cpd.anno$master_cpd_id == i),"cpd_name"]) # 插入匹配的药物
    }
}

rownames(trainPtype) <- ccl.name
trainPtype <- trainPtype[setdiff(rownames(trainPtype),ccl.miss),] # 去除未匹配的细胞系
colnames(trainPtype) <- cpd.name
trainPtype <- trainPtype[,setdiff(colnames(trainPtype),cpd.miss)] # 去除未匹配的药物
comccl <- intersect(rownames(trainPtype),colnames(trainExpr)) # 提取有表达且有药敏的细胞系
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5 # 纯化的测试集取表达稳定的基因
testExpr <- pure.tumoexpr[keepgene,]
# 取训练集和测试集共有的基因
comgene <- intersect(rownames(trainExpr),rownames(testExpr)) 
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- as.matrix(testExpr[comgene,])
outTab <- NULL
# 循环很慢，请耐心
for (i in 1:ncol(trainPtype)) { 
    display.progress(index = i,totalN = ncol(trainPtype))
    d <- colnames(trainPtype)[i]
    tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # 由于CTRP的AUC可能有0值，因此加一个较小的数值防止报错
    
    # 岭回归预测药物敏感性
    ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                    trainingPtype = tmp,
                                    testExprData = testExpr,
                                    powerTransformPhenotype = F,
                                    selection = 1))
    ptypeOut <- 2^ptypeOut - 0.00001 # 反对数
    outTab <- rbind.data.frame(outTab,ptypeOut)
}

dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc <- outTab
keepgene <- apply(ccl.expr, 1, mad) > 0.5
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1)
trainPtype <- as.data.frame(prism.auc.knn)
rownames(trainPtype) <- prism.ccl.anno[rownames(trainPtype),"cell_line_display_name"]
#colnames(trainPtype) <- sapply(strsplit(colnames(trainPtype)," (",fixed = T), "[",1)
comccl <- intersect(rownames(trainPtype),colnames(trainExpr))
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]

# 测试集
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5
testExpr <- pure.tumoexpr[keepgene,]
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- as.matrix(testExpr[comgene,])


outTab <- NULL

for (i in 1:ncol(trainPtype)) { 
    display.progress(index = i,totalN = ncol(trainPtype))
    d <- colnames(trainPtype)[i]
    tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # 由于PRISM的AUC可能有0值，因此加一个较小的数值防止报错
    ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                    trainingPtype = tmp,
                                    testExprData = testExpr,
                                    powerTransformPhenotype = F,
                                    selection = 1))
    ptypeOut <- 2^ptypeOut - 0.00001 # 反对数
    outTab <- rbind.data.frame(outTab,ptypeOut)
}

dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
prism.pred.auc <- outTab
top.pps <- Sinfo[Sinfo$riskscore >= quantile(Sinfo$riskscore,probs = seq(0,1,0.1))[10],] # 定义上十分位的样本
bot.pps <- Sinfo[Sinfo$riskscore <= quantile(Sinfo$riskscore,probs = seq(0,1,0.1))[2],] # 定义下十分位的样本
darkblue <- "#4DBBD5FF"
lightblue <- "#E64B35FF"
ctrp.log2fc <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
    display.progress(index = i,totalN = nrow(ctrp.pred.auc))
    d <- rownames(ctrp.pred.auc)[i]
    a <- mean(as.numeric(ctrp.pred.auc[d,rownames(top.pps)])) # 上十分位数的AUC均值
    b <- mean(as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])) # 下十分位数的AUC均值
    fc <- b/a
    log2fc <- log2(fc); names(log2fc) <- d
    ctrp.log2fc <- c(ctrp.log2fc,log2fc)
}

candidate.ctrp <- ctrp.log2fc[ctrp.log2fc > 0.1]
prism.log2fc <- c()
for (i in 1:nrow(prism.pred.auc)) {
    display.progress(index = i,totalN = nrow(prism.pred.auc))
    d <- rownames(prism.pred.auc)[i]
    a <- mean(as.numeric(prism.pred.auc[d,rownames(top.pps)])) # 上十分位数的AUC均值
    b <- mean(as.numeric(prism.pred.auc[d,rownames(bot.pps)])) # 下十分位数的AUC均值
    fc <- b/a
    log2fc <- log2(fc); names(log2fc) <- d
    prism.log2fc <- c(prism.log2fc,log2fc)
}
candidate.prism <- prism.log2fc[prism.log2fc > 0.25]

ctrp.cor <- ctrp.cor.p <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
    display.progress(index = i,totalN = nrow(ctrp.pred.auc))
    d <- rownames(ctrp.pred.auc)[i]
    a <- as.numeric(ctrp.pred.auc[d,rownames(Sinfo)]) 
    b <- as.numeric(Sinfo$riskscore)
    r <- cor.test(a,b,method = "pearson")$estimate; names(r) <- d
    p <- cor.test(a,b,method = "pearson")$p.value; names(p) <- d
    ctrp.cor <- c(ctrp.cor,r)
    ctrp.cor.p <- c(ctrp.cor.p,p)
}
candidate.ctrp2 <- ctrp.cor[ctrp.cor < -0.35]  # 这里我调整了阈值，控制结果数目
ctrp.candidate <- intersect(names(candidate.ctrp),names(candidate.ctrp2))

prism.cor <- prism.cor.p <- c()
for (i in 1:nrow(prism.pred.auc)) {
    display.progress(index = i,totalN = nrow(prism.pred.auc))
    d <- rownames(prism.pred.auc)[i]
    a <- as.numeric(prism.pred.auc[d,rownames(Sinfo)]) 
    b <- as.numeric(Sinfo$riskscore)
    r <- cor.test(a,b,method = "pearson")$estimate; names(r) <- d
    p <- cor.test(a,b,method = "pearson")$p.value; names(p) <- d
    prism.cor <- c(prism.cor,r)
    prism.cor.p <- c(prism.cor.p,p)
}
candidate.prism2 <- prism.cor[prism.cor < -0.5]  
prism.candidate <- intersect(names(candidate.prism),names(candidate.prism2))
cor.data <- data.frame(drug = ctrp.candidate,
                       r = ctrp.cor[ctrp.candidate],
                       p = -log10(ctrp.cor.p[ctrp.candidate]))
p1 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
    geom_segment(aes(xend=0,yend=drug),linetype = 2) +
    geom_point(aes(size=p),col = darkblue) +
    scale_size_continuous(range =c(2,8)) +
    scale_x_reverse(breaks = c(0, -0.3, -0.5),
                    expand = expansion(mult = c(0.01,.1))) + #左右留空
    theme_classic() +
    labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
    theme(legend.position = "bottom", 
          axis.line.y = element_blank(),axis.title = element_text(color="black",face="bold",size=16),axis.text = element_text(color="black",face="bold",size=15))

cor.data <- data.frame(drug = prism.candidate,
                       r = prism.cor[prism.candidate],
                       p = -log10(prism.cor.p[prism.candidate]))
cor.data$drug <- sapply(strsplit(cor.data$drug," (",fixed = T), "[",1)

p2 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
    geom_segment(aes(xend=0,yend=drug),linetype = 2) +
    geom_point(aes(size=p),col = darkblue) +
    scale_size_continuous(range =c(2,8)) +
    scale_x_reverse(breaks = c(0, -0.3, -0.5),
                    expand = expansion(mult = c(0.01,.1))) + #左右留空
    theme_classic() +
    labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
    theme(legend.position = "bottom", 
          axis.line.y = element_blank(),axis.title = element_text(color="black",face="bold",size=16),axis.text = element_text(color="black",face="bold",size=15))

ctrp.boxdata <- NULL
for (d in ctrp.candidate) {
    a <- as.numeric(ctrp.pred.auc[d,rownames(top.pps)]) 
    b <- as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])
    p <- wilcox.test(a,b)$p.value
    s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
    ctrp.boxdata <- rbind.data.frame(ctrp.boxdata,
                                     data.frame(drug = d,
                                                auc = c(a,b),
                                                p = p,
                                                s = s,
                                                group = rep(c("High riskscore","Low riskscore"),c(nrow(top.pps),nrow(bot.pps))),
                                                stringsAsFactors = F),
                                     stringsAsFactors = F)
}
p3 <- ggplot(ctrp.boxdata, aes(drug, auc, fill=group)) + 
    geom_boxplot(aes(col = group),outlier.shape = NA) + 
    # geom_text(aes(drug, y=min(auc) * 1.1, 
    #               label=paste("p=",formatC(p,format = "e",digits = 1))),
    #           data=ctrp.boxdata, 
    #           inherit.aes=F) + 
    geom_text(aes(drug, y=max(auc)), 
              label=ctrp.boxdata$s,
              data=ctrp.boxdata, 
              inherit.aes=F) + 
    scale_fill_manual(values = c(darkblue, lightblue)) + 
    scale_color_manual(values = c(darkblue, lightblue)) + 
    xlab(NULL) + ylab("Estimated AUC value") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 15,face="bold",color="black"),axis.text.y = element_text(size = 15,face="bold",color="black"),axis.title.y = element_text(size = 16,face="bold",color="black"),
          legend.position = "bottom",
          legend.title = element_blank()) 
dat <- ggplot_build(p3)$data[[1]]

p3 <- p3 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)

prism.boxdata <- NULL
for (d in prism.candidate) {
    a <- as.numeric(prism.pred.auc[d,rownames(top.pps)]) 
    b <- as.numeric(prism.pred.auc[d,rownames(bot.pps)])
    p <- wilcox.test(a,b)$p.value
    s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
    prism.boxdata <- rbind.data.frame(prism.boxdata,
                                      data.frame(drug = d,
                                                 auc = c(a,b),
                                                 p = p,
                                                 s = s,
                                                 group = rep(c("High riskscore","Low riskscore"),c(nrow(top.pps),nrow(bot.pps))),
                                                 stringsAsFactors = F),
                                      stringsAsFactors = F)
}
prism.boxdata$drug <- sapply(strsplit(prism.boxdata$drug," (",fixed = T), "[",1)

p4 <- ggplot(prism.boxdata, aes(drug, auc, fill=group)) + 
    geom_boxplot(aes(col = group),outlier.shape = NA) + 
    # geom_text(aes(drug, y=min(auc) * 1.1, 
    #               label=paste("p=",formatC(p,format = "e",digits = 1))),
    #           data=prism.boxdata, 
    #           inherit.aes=F) + 
    geom_text(aes(drug, y=max(auc)), 
              label=prism.boxdata$s,
              data=prism.boxdata, 
              inherit.aes=F) + 
    scale_fill_manual(values = c(darkblue, lightblue)) + 
    scale_color_manual(values = c(darkblue, lightblue)) + 
    xlab(NULL) + ylab("Estimated AUC value") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 15,face="bold",color="black"),axis.text.y = element_text(size = 15,face="bold",color="black"),axis.title.y = element_text(size = 16,face="bold",color="black"),
          legend.position = "bottom",
          legend.title = element_blank())
dat <- ggplot_build(p4)$data[[1]]

p4 <- p4 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
plot_grid(p1, p3, p2, p4, labels=c("A", "", "B", ""), 
          ncol=2, 
          rel_widths = c(2, 2))

##免疫细胞浸润
library(ComplexHeatmap)
library(ggsci)
library(circlize)
ee<-read.table("immune_cell.txt",head=T,sep='\t',check.names = F,row.names = 1)
my<-read.table("Clinical.txt",head=T,sep='\t',check.names = F,row.names = 1)
ee <- t(ee)
ee<-ee[,intersect(row.names(my),colnames(ee))]
identical(row.names(my),colnames(ee))
ee <- t(scale(t(ee)))
ee[ee > 2] <- 2 
ee[ee < -2] <- -2 
my$Status <- factor(my$Status)
my$risk <- factor(my$risk)
my$Stage <- factor(my$Stage)
my$Gender <- factor(my$Gender,levels = c('Female','Male'))
my$Age <- factor(my$Age,levels = c('<=65','>65'))
my$T=factor(my$T)
my$N=factor(my$N)
my$M=factor(my$M,levels = c('M0','M1','NA'))
risk <- c("#4DBBD5FF","#E64B35FF")
names(risk) <- levels(my$risk)
Age <- c(pal_nejm(alpha = 0.9)(8)[3],'#CF4E27')
names(Age) <- levels(my$Age)
Gender <- c('#E0864A','rosybrown')
names(Gender) <- levels(my$Gender)
Stage <- c('cornsilk','paleturquoise','goldenrod','firebrick')
names(Stage) <- levels(my$Stage)
T <- c('cornsilk','paleturquoise','goldenrod','firebrick')
names(T) <- levels(my$T)
N <- c('paleturquoise','goldenrod','firebrick')
names(N) <- levels(my$N)
M<- c('paleturquoise','firebrick',"white")
names(M) <- levels(my$M)
Status <- c('lavenderblush','slategray')
names(Status) <- levels(my$Status)
col<-list(risk,Age,Gender,Stage,T,N,M,Status)

Top = HeatmapAnnotation(risk=my$risk,
                        Age=my$Age,
                        Gender=my$Gender,
                        Stage= my$Stage,T=my$T,N=my$N,N=my$N,
                        Status = my$Status,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,title_gp = gpar(fontsize = 10,fontface = "bold"),ncol=1),
                        col=list(risk = risk,
                                 Age = Age,
                                 Gender = Gender,
                                 Stage= Stage,T=T,N=N,M=M,
                                 Status = Status
                        ),
                        show_annotation_name = TRUE,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))
Heatmap(ee,name='Z-score',
        cluster_rows = TRUE,top_annotation = Top,
        col=colorRamp2(c(-2,0,2),c('#1d6fae','white','#ef7700')),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        column_split = c(rep(1,264),rep(2,265)),
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 13), border = T,
                                  title_gp = gpar(fontsize = 13, fontface = "bold")),
        column_gap = unit(2,'mm'))

##细胞周期
library(reshape2)
library(ggpubr)
library(ggsci)
df<-read.table("Entosis_expr.txt",head=T,sep='\t',check.names = F,row.names = 1)
data<-melt(df,
           id.vars = c('Group'),
           measure.vars = colnames(df[-1]),
           variable.name='Immune_check',
           value.name='Expression')
p=ggboxplot(data, x="Immune_check", y="Expression", fill = "Group",alpha=0.8,orientation = "horizontal",notch = TRUE,
            ylab="Gene_expression",
            xlab="",
            palette = c("lancet"))
p=p+rotate_x_text(60)                    
p+stat_compare_means(aes(group=Group),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),label = "p.signif")+theme(axis.text = element_text(size = 13, face = "bold"),axis.title = element_text(size = 13, face = "bold"))
##药物靶点分析
library(readr)
library(ggplot2)
library(ggrepel)
library(cowplot)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
Tinfo <- read.delim("table s10 target.txt",sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
Cinfo <- read.csv("sample_info.csv", row.names = 1,check.names = F,stringsAsFactors = F,header = T)
expr.ccle <- read.csv("CCLE_expression.csv", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
ceres <- read.csv("CRISPR_gene_effect.csv", row.names = 1,check.names = F, stringsAsFactors = F, header = T)
Cinfo.hcc <- Cinfo[which(Cinfo$primary_disease == "Kidney Cancer"),]
expr.ccle<-expr.ccle[intersect(row.names(Cinfo.hcc),row.names(expr.ccle)),]
colnames(expr.ccle) <- sapply(strsplit(colnames(expr.ccle), " (", fixed = T), "[", 1)
expr.ccle <- as.data.frame(t(expr.ccle))
gene<-read.table("gene.txt",head=F,sep='\t',check.names = F)
ccle.set<-expr.ccle[intersect(gene$V1,row.names(expr.ccle)),]
ccle.set<-t(ccle.set)
library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggrepel)
model <- readRDS("model.rds") 
methodsValid <- names(model)

source("ML.R")
ccle.set = scaleData(data = ccle.set, centerFlags = T, scaleFlags = T)
RS <- CalRiskScore(fit = model[["Enet[alpha=0.2]"]], 
                                      new_data = ccle.set, 
                                      type = "lp") 
write.table(RS,"CCLE_risk.txt",quote=F,sep='\t')
#提取共有基因，既来自TCGA，又来自细胞系
TCGA<-read.table("uniq.symbol.txt",head=T,sep='\t',check.names = F,row.names = 1)
comtarget <- intersect(Tinfo$`Target genes`, rownames(TCGA))
comtarget <- intersect(comtarget, rownames(expr.ccle))
NKS<-read.table("Train_risk.txt",head=T,sep='\t',check.names = F,row.names = 1)
corTCGA <- corCERES <- NULL
for (i in comtarget) {
    
    drug <- paste(Tinfo[which(Tinfo$`Target genes` == i), "Agent names"], collapse = " | ") # 确定该靶点所对应的药物
    
    ## TCGA表达和临床样本的PPS相关性分析
    cor <- cor.test(as.numeric(TCGA[i, names(NKS)]),
                    as.numeric(NKS),
                    method = "spearman") # 原文采用spearman相关性
    corTCGA <- rbind.data.frame(corTCGA,
                                data.frame(target = i,
                                           r = cor$estimate,
                                           p = cor$p.value,
                                           drug = drug,
                                           row.names = i,
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
    
    ## CERES矩阵和细胞系的PPS相关性分析
    cor <- cor.test(as.numeric(expr.ccle[i, names(RS)]),
                    as.numeric(RS),
                    method = "spearman")
    corCERES <- rbind.data.frame(corCERES,
                                 data.frame(target = i,
                                            r = cor$estimate,
                                            p = cor$p.value,
                                            drug = drug,
                                            row.names = i,
                                            stringsAsFactors = F),
                                 stringsAsFactors = F)
    
}
write.table(corTCGA, file = "output_correlation_tcga.txt", sep = "\t",row.names = F, col.names = T, quote = F)
write.table(corCERES, file = "output_correlation_ccle.txt", sep = "\t", row.names = F, col.names = T, quote = F)
# 根据原文阈值筛选对应的候选靶点
candidate.TCGA <- corTCGA[which(corTCGA$r > 0.05 & corTCGA$p < 0.05),] # 蛋白谱靶点值与PPS需正相关
candidate.CERES <- corCERES[which(corCERES$r < -0.2 & corCERES$p < 0.05),] # CERES靶点值与PPS需负相关
candidate.target <- intersect(candidate.TCGA$target, candidate.CERES$target) # 匹配上了原文中的3个，也许还存在数据标准化问题导致结果不完全一致
# 把筛选出的药物靶点保存到文件
write.table(candidate.target, "output_candidate.target.txt", quote = F, row.names = F)
# 设置颜色
grey <- "#BFBFBF"
lightred <- "#FDC9B5"
red <- "#E9583A"
lightblue <- "#66C2A4"
blue <- "#006D2C"
corTCGA$color <- "A" # 给颜色做编号，基础为A
corTCGA[which(corTCGA$r > 0.1 & corTCGA$p < 0.05), "color"] <- "B" # 显著的为B
corTCGA[candidate.target,"color"] <- "C" # 感兴趣的为C
corTCGA$size <- "A" # 同理给大小做编号
corTCGA[which(corTCGA$r > 0.1 & corTCGA$p < 0.05), "size"] <- "B"
corTCGA[candidate.target,"size"] <- "C"
corTCGA <- corTCGA[order(corTCGA$color),] # 让想突出的靶点出现在列表的末尾，这样绘图时不会被其他点所遮挡

selecttargets <- corTCGA[candidate.target,]
selecttargets$label <- rownames(selecttargets)
p1 <- ggplot(corTCGA, aes(r, -log10(p))) + 
    geom_point(aes(color = color, size = size)) + 
    scale_color_manual(values = c(grey, lightred, red))+ # 将编号映射到对应颜色
    scale_size_manual(values = c(1,2,3)) +  # 将编号映射到对应大小
    xlab("Spearman's rank correlation coefficient") + 
    ylab("-log10(P-value)") +
    scale_x_continuous(
        breaks = c(-0.6,-0.3,0,0.3,0.6), # x轴的一些修饰
        labels = c(-0.6,-0.3,0,0.3,0.6),
        limits = c(-0.6, 0.6)) + 
    geom_vline(xintercept = 0.1, color="grey70", # 添加垂直相关性阈值线
               linetype = "longdash", lwd = 0.6) + 
    geom_hline(yintercept = -log10(0.05), color = "grey70", # 添加水平P值阈值线
               linetype = "longdash", lwd = 0.6) +
    theme_bw() +
    theme(axis.ticks = element_line(size = 0.2, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.text = element_text(size = 16, color = "black",face="bold"),
          axis.title = element_text(size = 18, color = "black",face="bold"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
p2 <- p1 + 
    geom_text_repel(data = selecttargets,
                    aes(x = r, y = -log10(p), 
                        label = label),
                    colour="black", nudge_x = .15,
                    box.padding = 0.5,
                    nudge_y = 1,
                    segment.curvature = -0.1,
                    segment.ncp = 3,
                    segment.angle = 20,
                    size = 10, min.segment.length = 0,
                    point.padding = unit(1, "lines"))

p.tcga <- p2

#CCLE
corCERES$color <- "A"
corCERES[which(corCERES$r < -0.2&corCERES$p < 0.05), "color"] <- "B"
corCERES[candidate.target, "color"] <- "C"
corCERES$size <- "A"
corCERES[which(corCERES$r < -0.2&corCERES$p < 0.05), "size"] <- "B"
corCERES[candidate.target, "size"] <- "C"
corCERES <- corCERES[order(corCERES$color),] # 让想突出的靶点出现在列表的末尾，这样绘图时不会被其他点所遮挡

selecttargets <- corCERES[candidate.target,]
selecttargets$label <- rownames(selecttargets)

p1 <- ggplot(corCERES, aes(r, -log10(p))) + 
    geom_point(aes(color = color, size = size)) + 
    scale_color_manual(values = c(grey, lightblue, blue))+
    scale_size_manual(values = c(1,2,3)) + 
    xlab("Spearman's rank correlation coefficient") + 
    ylab("-log10(P-value)") +
    scale_x_continuous(
        breaks = c(-0.8,-0.5,0,0.5,0.8), 
        labels = c(-0.8,-0.5,0,0.5,0.8),
        limits = c(-0.8, 0.8)) + 
    geom_vline(xintercept = -0.325, color = "grey70", 
               linetype = "longdash", lwd = 0.6) + 
    geom_hline(yintercept = -log10(0.05), color = "grey70", 
               linetype = "longdash", lwd = 0.6) +
    theme_bw() +
    theme(axis.ticks = element_line(size = 0.2, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.text = element_text(size = 16, color = "black",face="bold"),
          axis.title = element_text(size = 18, color = "black",face="bold"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
p2 <- p1 + 
    geom_text_repel(data = selecttargets,
                    aes(x = r, y = -log10(p), 
                        label = label),
                    colour="black", nudge_x = .15,
                    box.padding = 0.5,
                    nudge_y = 1,
                    segment.curvature = -0.1,
                    segment.ncp = 3,
                    segment.angle = 20,
                    size = 10, min.segment.length = 0,
                    point.padding = unit(1, "lines"))
p.ccle <- p2
p.ccle

scatterProteome <- scatterCERES <- list()
for (i in candidate.target) {
    tmp <- data.frame(var = as.numeric(TCGA[i, names(NKS)]), NKS = as.numeric(NKS))
    cor <- cor.test(tmp$var,
                    tmp$NKS,
                    method = "spearman")
    txt <- paste0("r = ", round(cor$estimate,2), "\n", "P < 0.00001 ") # 构建相关性值的文字标签
    scatterProteome[[i]] <- 
        ggplot(tmp, aes(NKS, var)) + 
        geom_ribbon(stat = "smooth", method = "lm", se = TRUE, # 先画置信区间的彩带以免遮挡散点
                    fill = alpha(lightred, 0.6)) + 
        geom_smooth(span = 2, method = lm, color = red, fill = NA) + # 绘制回归线
        geom_point(color = red, size = 2) + 
        xlab("PARS scores") + 
        ylab(paste0("Protein abundance of ", i)) +
        theme_bw() +
        theme(axis.ticks = element_line(size = 0.2, color = "black"),
              axis.ticks.length = unit(0.2, "cm"),
              axis.text = element_text(size = 16, color = "black"),
              axis.title = element_text(size = 18, color = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank()) +
        annotate("text", # 添加相关性的文字标签
                 x = min(tmp$NKS), 
                 y = max(tmp$var), 
                 hjust = 0, fontface = 4, 
                 label = txt)
    scatterProteome[[i]]
    ggsave(filename = paste0("scatter plot between risk score and TCGA abundance of ", i, ".pdf"), width = 5, height = 5)
    
    tmp <- data.frame(var = as.numeric(expr.ccle[i, names(RS)]), RS = as.numeric(RS),cell=names(RS))
    cor <- cor.test(tmp$var,
                    tmp$RS,
                    method = "spearman")
    txt <- paste0("r = ", round(cor$estimate,2), "\n", "P = ", round(cor$p.value, 4))
    scatterCERES[[i]] <- 
        ggplot(tmp, aes(RS, var)) + 
        geom_ribbon(stat = "smooth", method = "lm", se = TRUE,
                    fill = alpha(lightblue, 0.6)) + 
        geom_smooth(span = 2, method = lm, color = blue, fill = NA) +
        geom_point(color = blue, size = 2) + 
        geom_hline(yintercept = 0, color=red, # 添加CERES的0值水平线
                   linetype="longdash", lwd = 0.6) +
        xlab("risk scores") + 
        ylab(paste0("CERES score of ", i)) +
        theme_bw() + geom_text_repel(aes(label = cell), size = 3,color="#006D2C")+
        theme(axis.ticks = element_line(size = 0.2, color = "black"),
              axis.ticks.length = unit(0.2, "cm"),
              axis.text = element_text(size = 16, color = "black"),
              axis.title = element_text(size = 18, color = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank()) +
        annotate("text", 
                 x = min(tmp$RS), 
                 y = max(tmp$var), 
                 hjust = 0, fontface = 4, 
                 label = txt)
    scatterCERES[[i]]
    ggsave(filename = paste0("scatter plot between risk score and ceres score of ", i, ".pdf"), width = 5, height = 5)
}

#靶点泛癌分析
library(TCGAplot)
gene<-"TNFSF11"
pan_boxplot(gene)
pan_tumor_boxplot(gene)
gene_TMB_radar(gene,method = "pearson")
gene_MSI_radar(gene,method = "pearson")
gene_chemokine_heatmap(gene,method="pearson",lowcol="blue",highcol="red")
gene_receptor_heatmap(gene,method="pearson",lowcol="blue",highcol="red")
gene_immustimulator_heatmap(gene,method="pearson",lowcol="blue",highcol="red")
gene_immuinhibitor_heatmap(gene,method="pearson",lowcol="blue",highcol="red")
gene_immucell_heatmap(gene,method="pearson",lowcol="blue",highcol="red")

##单细胞分析
library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)
library(stringr)
library(harmony)
library(recall)
library(SCP)
library(SCpubr)
library(ggsci)

dir.create("merged_data", showWarnings = FALSE)
#读取GSE156632数据集
read_GSE156632 <- function() {
  cat("正在读取GSE156632数据集...\n")
  
  # 定义文件路径（请根据实际情况修改）
  files <- c(
    "GSE156632/GSM4735365_RCC1n.csv.gz", 
    "GSE156632/GSM4735367_RCC2n.csv.gz",
	"GSE156632/GSM4735369_RCC3n.csv.gz",
	"GSE156632/GSM4735371_RCC4n.csv.gz",
	"GSE156632/GSM4735373_RCC5n.csv.gz",
    "GSE156632/GSM4735364_RCC1t.csv.gz",
	"GSE156632/GSM4735366_RCC2t.csv.gz",
	"GSE156632/GSM4735368_RCC3t.csv.gz",
    "GSE156632/GSM4735370_RCC4t.csv.gz",
    "GSE156632/GSM4735372_RCC5t.csv.gz"
	)
  
  sample_names <- c("GSE156632N1", "GSE156632N2", "GSE156632N3", "GSE156632N4", "GSE156632N5",
					"GSE156632T1", "GSE156632T2", "GSE156632T3", "GSE156632T4", "GSE156632T5")
  
  seurat_list <- list()
  
  for (i in 1:length(files)) {
    cat(paste0("读取文件: ", files[i], "\n"))
    
    # 读取压缩的CSV文件
    # 注意：这里假设CSV文件是基因在行，细胞在列
    counts_data <- fread(files[i], header = TRUE)
    
	# 查看数据结构
    cat(paste0("数据维度: ", nrow(counts_data), "行 × ", ncol(counts_data), "列\n"))
	
	# 提取基因名（第二列SYMBOL）
    gene_symbols <- counts_data[[2]]
	
	  if (any(duplicated(gene_symbols))) {
      cat("警告: 发现重复的基因符号，正在处理...\n")
      
      # 对于重复的基因符号，我们可以：
      # 1. 保留第一个出现的
      # 2. 或者合并表达值（这里选择保留第一个）
      dup_genes <- gene_symbols[duplicated(gene_symbols)]
      cat(paste0("重复基因数量: ", length(unique(dup_genes)), "\n"))
      
      # 保留第一次出现的基因
      keep_rows <- !duplicated(gene_symbols)
      counts_data <- counts_data[keep_rows, ]
      gene_symbols <- gene_symbols[keep_rows]
    }
	
	# 提取计数数据（从第三列开始）
    count_matrix_data <- counts_data[, 3:ncol(counts_data), with = FALSE]
	
	# 转换为矩阵
    counts_matrix <- as.matrix(count_matrix_data)
	 # 设置行名（基因符号）和列名（细胞）
    rownames(counts_matrix) <- gene_symbols
    
    # 如果列名是数字，添加前缀
    if (is.numeric(colnames(counts_matrix)) || 
        all(grepl("^[0-9]", colnames(counts_matrix)))) {
      colnames(counts_matrix) <- paste0(sample_names[i], "_Cell", colnames(counts_matrix))
    } else {
      colnames(counts_matrix) <- paste0(sample_names[i], "_", colnames(counts_matrix))
    }
	
	# 转换为稀疏矩阵
    counts_sparse <- as(counts_matrix, "dgCMatrix")
	
    # 创建Seurat对象
    seurat_obj <- CreateSeuratObject(
      counts = counts_matrix,
      project = "GSE156632",
      min.cells = 3,
      min.features = 200
    )
    
    # 添加元数据
    seurat_obj$Dataset <- "GSE156632"
    seurat_obj$Sample <- sample_names[i]
    seurat_obj$Condition <- ifelse(grepl("t$", sample_names[i]), "Tumor", "Normal")
    
    seurat_list[[sample_names[i]]] <- seurat_obj
    
    rm(counts_data, counts_matrix)  # 清理内存
    gc()
  }
  
  return(seurat_list)
}
##读取GSE159115数据集（H5格式）
read_GSE159115 <- function() {
  cat("正在读取GSE159115数据集...\n")
  
  # 定义文件路径
  files <- c(
	"GSE159115/GSM4819726_SI_18856_filtered_gene_bc_matrices_h5.h5",
	"GSE159115/GSM4819727_SI_18855_filtered_gene_bc_matrices_h5.h5",
    "GSE159115/GSM4819728_SI_19704_filtered_gene_bc_matrices_h5.h5",
	"GSE159115/GSM4819729_SI_19703_filtered_gene_bc_matrices_h5.h5",
    "GSE159115/GSM4819733_SI_22369_filtered_gene_bc_matrices_h5.h5",
	"GSE159115/GSM4819734_SI_22368_filtered_gene_bc_matrices_h5.h5",
	"GSE159115/GSM4819735_SI_22605_filtered_gene_bc_matrices_h5.h5",
	"GSE159115/GSM4819736_SI_22604_filtered_gene_bc_matrices_h5.h5"
  )
  
  sample_names <- c("GSE159115P2T", "GSE159115P2N", "GSE159115P3T", "GSE159115P3N", 
                    "GSE159115P5T", "GSE159115P5N", "GSE159115P6T", "GSE159115P6N")
  
  seurat_list <- list()
  
  for (i in 1:length(files)) {
    cat(paste0("读取文件: ", files[i], "\n"))
    
    # 读取10x Genomics H5文件
    counts_matrix <- Read10X_h5(files[i])
    
    # 创建Seurat对象
    seurat_obj <- CreateSeuratObject(
      counts = counts_matrix,
      project = "GSE159115",
      min.cells = 3,
      min.features = 200
    )
    
    # 添加元数据
    seurat_obj$Dataset <- "GSE159115"
    seurat_obj$Sample <- sample_names[i]
    seurat_obj$Condition <- ifelse(grepl("N", sample_names[i]), "Normal", "Tumor")
    seurat_list[[sample_names[i]]] <- seurat_obj
    
    rm(counts_matrix)  # 清理内存
    gc()
  }
  
  return(seurat_list)
}

##读取GSE178481数据集（CSV格式）
read_GSE178481 <- function() {
  cat("正在读取GSE178481数据集...\n")
  
  # 定义文件路径
  files <- c(
    "GSE178481/GSM5392398_RCC-PR6-Normal.count.csv.gz",
    "GSE178481/GSM5392399_RCC-PR6-PTumor.count.csv.gz",
    "GSE178481/GSM5392405_RCC-PR5-Normal.count.csv.gz",
	"GSE178481/GSM5392406_RCC-PR5-PTumor1.count.csv.gz",
	"GSE178481/GSM5392407_RCC-PR5-PTumor2.count.csv.gz",
	"GSE178481/GSM5392408_RCC-PR5-PTumor3.count.csv.gz",
	"GSE178481/GSM5392409_RCC-BM1-PTumor.count.csv.gz",
	"GSE178481/GSM5392410_RCC-BM1-Normal.count.csv.gz",
	"GSE178481/GSM5392416_RCC-PR2-Normal.count.csv.gz",
	"GSE178481/GSM5392417_RCC-PR2-PTumor.count.csv.gz",
	"GSE178481/GSM5392418_RCC-PR4-Normal.count.csv.gz",
	"GSE178481/GSM5392419_RCC-PR4-PTumor.count.csv.gz",
	"GSE178481/GSM5392420_RCC-PR3-Normal.count.csv.gz",
	"GSE178481/GSM5392421_RCC-PR3-PTumor1.count.csv.gz",
	"GSE178481/GSM5392422_RCC-PR3-PTumor2.count.csv.gz",
	"GSE178481/GSM5392423_RCC-PR3-PTumor3.count.csv.gz"
  )
  
  sample_names <- c("RCC.PR6.Normal", "RCC.PR6.PTumor", "RCC.PR5.Normal", "RCC.PR5.PTumor1", "RCC.PR5.PTumor2",
                    "RCC.PR5.PTumor3", "RCC.BM1.PTumor", "RCC.BM1.Normal", "RCC.PR2.Normal", "RCC.PR2.PTumor",
					"RCC-PR4.Normal", "RCC.PR4.PTumor", "RCC.PR3.Normal", "RCC.PR3.PTumor1", "RCC.PR3.PTumor2",
					"RCC.PR3.PTumor3")
  
  seurat_list <- list()
  
  for (i in 1:length(files)) {
    cat(paste0("读取文件: ", files[i], "\n"))
    
    # 读取压缩的CSV文件
    counts_data <- fread(files[i], header = TRUE)
    
    # 转换数据
    if (class(counts_data)[1] == "data.table") {
      # 检查数据格式：可能是基因在行，也可能在列
      # 假设第一列是基因名，其余列是细胞
      if (any(grepl("^ENSG", counts_data[[1]])) | any(grepl("^MT-", counts_data[[1]]))) {
        # 基因在行
        genes <- counts_data[[1]]
        counts_data[[1]] <- NULL
        counts_matrix <- as(as.matrix(counts_data), "dgCMatrix")
        rownames(counts_matrix) <- genes
        colnames(counts_matrix) <- paste0(sample_names[i], "_", colnames(counts_matrix))
      } else {
        # 可能需要转置
        # 尝试将第一行作为基因名
        genes <- as.character(counts_data[1, ])[-1]
        samples <- as.character(counts_data[-1, 1])
        counts_data <- counts_data[-1, -1]
        counts_matrix <- as(as.matrix(counts_data), "dgCMatrix")
        rownames(counts_matrix) <- samples
        colnames(counts_matrix) <- genes
        counts_matrix <- t(counts_matrix)
        colnames(counts_matrix) <- paste0(sample_names[i], "_", colnames(counts_matrix))
      }
    }
    
    # 创建Seurat对象
    seurat_obj <- CreateSeuratObject(
      counts = counts_matrix,
      project = "GSE178481",
      min.cells = 3,
      min.features = 200
    )
    
    # 添加元数据
    seurat_obj$Dataset <- "GSE178481"
    seurat_obj$Sample <- sample_names[i]
    seurat_obj$Condition <- ifelse(grepl("Normal", sample_names[i]), "Normal", "Tumor")
    
    seurat_list[[sample_names[i]]] <- seurat_obj
    
    rm(counts_data, counts_matrix)
    gc()
  }
  
  return(seurat_list)
}

merge_all_datasets <- function() {
  cat("开始合并所有数据集...\n")
  
  # 读取所有数据集
  gse156632_list <- read_GSE156632()
  gse159115_list <- read_GSE159115()
  gse178481_list <- read_GSE178481()
  
  # 合并所有Seurat对象
  all_objects <- c(gse156632_list, gse159115_list, gse178481_list)
  
  # 方法1: 直接合并（如果内存足够）
  cat("正在合并数据集...\n")
  merged_seurat <- merge(
    x = all_objects[[1]],
    y = all_objects[-1],
    add.cell.ids = names(all_objects),
    project = "Combined_ccRCC"
  )
  
  # 添加详细的元数据
  merged_seurat$Dataset_Sample <- paste(merged_seurat$Dataset, merged_seurat$Sample, sep = "_")
  
  # 查看合并后的基本信息
  cat("\n合并完成！\n")
  cat(paste0("总细胞数: ", ncol(merged_seurat), "\n"))
  cat(paste0("总基因数: ", nrow(merged_seurat), "\n"))
  cat(paste0("数据集数: ", length(unique(merged_seurat$Dataset)), "\n"))
  cat(paste0("样本数: ", length(unique(merged_seurat$Sample)), "\n"))
  
  # 保存合并后的对象
  saveRDS(merged_seurat, file = "merged_data/merged_seurat_object.rds")
  
  return(merged_seurat)
}
#数据质量控制
quality_control <- function(seurat_obj) {
  cat("正在进行质量控制...\n")
  
  # 计算线粒体基因百分比
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # 计算核糖体基因百分比
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  
  # 计算血红蛋白基因百分比
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P)]")
  
  # 可视化QC指标
  pdf("merged_data/QC_metrics.pdf", width = 12, height = 8)
  
  # 按数据集显示QC指标
  v1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                group.by = "Dataset", pt.size = 0.1, ncol = 3)
  print(v1)
  
  # 特征之间的相关性
  p1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Dataset")
  p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Dataset")
  print(p1 + p2)
  
  dev.off()
  
  # 应用QC过滤（可根据数据调整阈值）
  cat("应用QC过滤...\n")
  seurat_obj_filtered <- subset(seurat_obj,
                                subset = nFeature_RNA > 200 & nFeature_RNA < 7500 &
                                nCount_RNA > 1000&
                                percent.mt < 15)
  
  cat(paste0("过滤后细胞数: ", ncol(seurat_obj_filtered), "\n"))
  cat(paste0("过滤掉的细胞数: ", ncol(seurat_obj) - ncol(seurat_obj_filtered), "\n"))
  
  saveRDS(seurat_obj_filtered, file = "merged_data/merged_seurat_filtered.rds")
  
  return(seurat_obj_filtered)
}
##运行代码
merged_data <- merge_all_datasets()
filtered_data <- quality_control(merged_data)
scRNA <- filtered_data
#Harmony数据校正批次效应

#1. 数据标准化
scRNA <- NormalizeData(scRNA, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000)
# 2. 寻找高变基因
scRNA <- FindVariableFeatures(scRNA, 
                                     selection.method = "vst", 
                                     nfeatures = 2000)
# 3. 识别前10个高变基因
top10 <- head(VariableFeatures(scRNA), 10)
cat("前10个高变基因:\n")
print(top10)

# 4. 可视化高变基因
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
pdf("merged_data/variable_features.pdf", width = 12, height = 6)
print(plot2)
dev.off()

# 5. 数据scale标准化
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA))

# 6. 运行PCA
scRNA <- RunPCA(scRNA, npcs = 10,
                    verbose = FALSE)
# 7. Harmnony 校正批次效应
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "Dataset")

# 8.对Harmony降维结果进行可视化
p1 <- DimPlot(scRNA_harmony, reduction = "harmony", group.by = "Dataset")

#9. 寻找最佳阈值
merged_seurat <- scRNA_harmony
merged_seurat <- FindNeighbors(merged_seurat,reduction = "harmony", dims = 1:10)#%>%
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)

#10. 可视化聚类簇
merged_seurat <- RunUMAP(merged_seurat, dims = 1:10,reduction = "harmony", reduction.name = "UMAP",min.dist =0.1)
p2 <- DimPlot(merged_seurat , reduction = "UMAP", group.by = "Sample")
p3 <- DimPlot(merged_seurat , reduction = "UMAP", group.by = "Dataset")
p4 <- DimPlot(merged_seurat, reduction = "UMAP", group.by = "seurat_clusters")
p5 <- DimPlot(merged_seurat, reduction = "UMAP", group.by = "Condition")

##11. 单细胞注释
###手动注释
cell_markers <- list(
  "B cell" = c("CD79A", "CD79B", "MS4A1"),
  "Macrophages" = c("CD163", "S100A9", "LYZ"),
  "Fibroblast cell" = c("DCN", "LOX","LUM"),
  "Monocytes" = c("C1QA", "CD14", "C1QB"),
  "Mast cell" = c("CAP","HPGD", "TPSB2"),
  "Mesenchymal cell" = c("IGFBP5", "TM4SF1"),
  "Endothelial cell" = c("PECAM1", "ADGRL4", "AQP1"),
  "Smooth muscle cell" = c("PLN", "MYH11","KCNMB1"),
  "T cell" = c("CD3D", "CD2","CD3E","CCL5"),
  "Plasmacyte" = c("IGHG4", "IGKC", "JCHIAN"))

p6 <- DotPlot(merged_seurat, 
                          features = unique(unlist(cell_markers)),
                          group.by = "seurat_clusters",
                          cols = c("blue", "red"),
                          dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Cell Type Marker Expression Across Clusters")

scRNA <- merged_seurat
scRNA <- subset(scRNA,RNA_snn_res.0.5!= 5)
scRNA <- subset(scRNA,RNA_snn_res.0.5!= 9)
scRNA <- subset(scRNA,RNA_snn_res.0.5!= 13)
scRNA <- subset(scRNA,RNA_snn_res.0.5!= 14)

scRNA$celltype <- recode(scRNA$RNA_snn_res.0.5,
                         "0" = "T cell", 
                         "1" = "Myeloid", 
                         "2" = "T cell", 
                         "3" = "Fibroblast cell",
                         "4" = "Endothelial cell", 
                         "5" = "unknown",
                         "6" = "Myeloid",
                         "7" = "Epithelial cell", 
                         "8" = "Epithelial cell",
                         "9" = "unknown",
                         "10" = "Epithelial cell",
                         "11" = "Endothelial cell",
                         "12" = "T cell",
                         "13" = "unknown",
                         "14" = "unknown", 
                         "15" = "Mast cell",#
                         "16" = "Endothelial cell",#
                         "17" = "Myeloid",#
                         "18" = "Fibroblast cell",
                         "19" = "T cell",
                         "20" = "Endothelial cell",
						 "21" = "Myeloid",
						 "22" = "Myeloid"
						)

##绘制细胞图
p7 <- DimPlot(scRNA, reduction = "UMAP", group.by='celltype',pt.size = 1.5) +scale_color_aaas()

##细胞-marker气泡图
markers <- c("CD3D", "CD2", "CD3E", "CCL5", "LYZ", "CD14", "CD68", "VCAN", "COL1A2", "ACTA2", 
	         "VWF", "FLT1", "KRT18", "KRT19", "KRT8", "EPCAM", "KIT", "HPGDS" )
p8 <- DotPlot(scRNA, features = markers,group.by = "celltype",)+coord_flip()+
    theme_bw()+
    theme(panel.grid = element_blank(), axis.text =element_text(hjust = 1,vjust = 0.5,face = "bold", size = 12,color = "black"))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

saveRDS(scRNA, file = "scRNA.rds")

##细胞百分比图
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
my_theme <- function(base_size = 14, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # 标题样式
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5, 
                               margin = margin(b = 20)),
      plot.subtitle = element_text(size = base_size + 2, face = "bold", hjust = 0.5,
                                  margin = margin(b = 15)),
      
      # 坐标轴标题
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      
      # 坐标轴文本
      axis.text = element_text(size = base_size, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                margin = margin(t = 5)),
      axis.text.y = element_text(margin = margin(r = 5)),
      
      # 图例
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1, face = "bold"),
      legend.position = "right",
      
      # 分面
      strip.text = element_text(size = base_size, face = "bold"),
      
      # 网格线
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      
      # 绘图区域
      plot.margin = margin(20, 20, 20, 20),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# 设置默认主题
theme_set(my_theme(base_size = 16))

#.计算不同细胞类型在样本间的差异
cell_proportions <- scRNA@meta.data %>%
    group_by(Dataset, Sample, celltype) %>%
    summarise(n_cells = n()) %>%
    group_by(Dataset, Sample) %>%
    mutate(total_cells = sum(n_cells),
           proportion = n_cells / total_cells * 100)
p9 <- ggplot(cell_proportions, aes(x = Sample, y = proportion, fill = celltype)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~Dataset, scales = "free_x") +
    scale_fill_npg() +
    labs(title = "Cell Type Proportions by Dataset",
         x = "Sample", 
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "lightgray"))+my_theme(base_size = 16)+theme(
              legend.title = element_text(size = 14, face = "bold"),
              legend.text = element_text(size = 12, face = "bold"),
              legend.key.size = unit(1, "cm"))

#.计算不同细胞类型在不同组中的差异
table(scRNA$celltype)
prop.table(table(scRNA$celltype))
table(scRNA$celltype, scRNA$Condition)
Cellratio <- prop.table(table(scRNA$celltype, scRNA$Condition), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
p10 <- ggplot(Cellratio) + 
    geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
    scale_fill_manual(values=c(pal_npg()(9)))+
    labs(x='Sample',y = 'Ratio')+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.text = element_text(color="black",face="bold",size=16),axis.title = element_text(color="black",face="bold",size=18),legend.title = element_blank())
ggsave(p10,filename = "p10_barplot.pdf")

####RO/e偏好性分析
library(Startrac)
library(circlize)
library(ComplexHeatmap)
sco_platelets <- scRNA
#构建Roe计算需要的输入表格
data <- sco_platelets@meta.data
data <- data[,c(5,6,13)]
colnames(data) <- c("Sample", "Condition", "celltype")
#核心计算，Roe计算
Roe <- calTissueDist(data,
                     byPatient = F,
                     colname.cluster = "celltype", # 不同细胞亚群
                     colname.patient = "Sample", # 不同样本
                     colname.tissue = "Condition", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0)
Roe
#可视化
col_fun  <-  colorRamp2(c(min(Roe,na.rm  =  TRUE),1,max(Roe,na.rm  =TRUE)),  
                        c("#f6f8e6",  "#f9a33e",  "red"))
pdf("p11.STARTRAC_Roe.pdf",width = 6,height = 8)
Heatmap(as.matrix(Roe),
        show_heatmap_legend = TRUE, 
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = 'right',
        column_names_side = "top",
        show_column_names = TRUE,
        show_row_names = TRUE,
        #row_split = 3,
        col = col_fun,
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        heatmap_legend_param = list(
          title = "Ro/e",
          at = c(0, max(Roe)),
          labels = c("0", "Max.")  # 对应标签
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          value <- Roe[i, j]
          # 符号标记逻辑，+++, Ro/e > 1; ++, 0.8 < Ro/e ≤ 1; +, 0.2 ≤ Ro/e ≤ 0.8; +/−, 0 < Ro/e < 0.2; −, Ro/e = 0
          symbol <- if(value == 0) {
            "−"
          } else if(value > 0 & value < 0.2) {
            "+/−"
          } else if(value >= 0.2 & value <= 0.8) {
            "+"
          } else if(value > 0.8 & value <= 1) {
            "++"
          } else if(value > 1) {
            "+++"
          }
          # 显示符号
          grid.text(symbol, x, y, gp = gpar(
            fontsize = 10, 
            col = "black")
          )
        }
)
dev.off()

#######不同亚型功能富集分析
library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)
library(clustree)
library(Matrix)
library(dplyr)
library(ggsci)
library(ClusterGVis)
library(org.Hs.eg.db)
scRNA <- readRDS("scRNA.rds")
Idents(scRNA) <- scRNA@meta.data$celltype
scRNA <- JoinLayers(scRNA)
markers_genes.all <- FindAllMarkers(scRNA, 
                               logfc.threshold = 0.25, 
                               test.use = "wilcox", 
                               min.pct = 0.25, 
                               min.diff.pct = 0.2, 
                               only.pos = TRUE, 
                               assay = "RNA")
markers <- markers_genes.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

# prepare data from seurat object
st.data <- prepareDataFromscRNA(object = scRNA,
                                diffData = markers,
                                showAverage = TRUE)
str(st.data)

# enrich GO (BP) for clusters
enrich_BP <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)

# add gene name
markGenes = unique(markers$gene)[sample(1:length(unique(markers$gene)),40,
                                             replace = F)]

# line plot
visCluster(object = st.data,
           plotType = "line")
#Heatmap plot
pdf('p12.Enrichmentplot_GO.pdf',height = 10,width = 16,onefile = F)
visCluster(object = st.data,
           plotType = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenesSide = "left",
           annoTermData = enrich_BP,
		   goCol = rep(jjAnno::useMyCol("calm",n = 7),each = 5),
		   lineSide = "left",
           addBar = T)
dev.off()

# enrich KEGG for clusters
enrich_KEGG <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "KEGG",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)

# add gene name
markGenes = unique(markers$gene)[sample(1:length(unique(markers$gene)),40,
                                             replace = F)]
# line plot
visCluster(object = st.data,
           plotType = "line")
#Heatmap plot
pdf('p12.Enrichmentplot_KEGG.pdf',height = 10,width = 16,onefile = F)
visCluster(object = st.data,
           plotType = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenesSide = "left",
           annoKeggData = enrich_KEGG,
		   lineSide = "left",
           addBar = T)
dev.off()

##基因在不同细胞类型中表达情况
library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)
library(clustree)
library(Matrix)
library(dplyr)
library(ggsci)
library(ClusterGVis)
library(org.Hs.eg.db)
scRNA <- readRDS("scRNA.rds")
Idents(scRNA) <- scRNA@meta.data$celltype
scRNA <- JoinLayers(scRNA)
gene <- read.table("gene.txt",head=F,sep='\t',check.names = F)
marker <- gene$V1
p13 <- DotPlot(scRNA, features = marker)+coord_flip()+
    theme_bw()+
    theme(panel.grid = element_blank(), axis.text=element_text(hjust = 1,vjust = 0.5,face = "bold", size = 14))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

##细胞通讯
library(limma)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(CellChat)

# 1. 创建CellChat对象
# 使用标准化数据和细胞注释信息
cellchat <- createCellChat(object = scRNA, 
                           group.by = "celltype", # 替换为你的细胞类型注释列名
                           assay = "RNA") # 默认assay

# 2. 设置配体-受体数据库（使用默认的CellChatDB）
CellChatDB <- CellChatDB.human # 如果是小鼠数据，使用 CellChatDB.mouse
cellchat@DB <- CellChatDB

# 3. 预处理：对表达数据进行过表达分析
cellchat <- subsetData(cellchat) # 可选的子集化，加速计算
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 4. 计算细胞通讯概率
# 使用默认的“truncatedMean”方法，设定值0.1表示截取高表达基因
cellchat <- computeCommunProb(cellchat, 
                              population.size = TRUE, # 考虑细胞群体大小
                              type = "truncatedMean", 
                              trim = 0.1)

# 过滤掉低可靠性的交互（通常由少数细胞贡献）
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 5. 从细胞水平推断到通路水平
cellchat <- computeCommunProbPathway(cellchat)

# 6. 整合计算所有信号通路
cellchat <- aggregateNet(cellchat)

# 7. 整体通讯数量热图（按细胞群）
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# 8. 单个细胞通讯强度 （按细胞群）
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#比较不同条件（肿瘤 vs 正常）
# 根据样本条件拆分对象
# 按条件拆分数据
tumor_seu <- subset(scRNA, subset = condition == "Tumor")
normal_seu <- subset(scRNA, subset = condition == "Normal")

# 打印拆分后的细胞数
cat("Tumor样本细胞数:", ncol(tumor_seu), "\n")
cat("Normal样本细胞数:", ncol(normal_seu), "\n")

# 检查每个条件的细胞类型分布
table(tumor_seu$celltype)
table(normal_seu$celltype)

#分别创建CellChat对象并推断通讯网络

# 创建肿瘤样本的CellChat对象
cellchat_tumor <- createCellChat(object = tumor_seu,
                                 group.by = "celltype",  # 按细胞类型分组
                                 assay = "RNA")  # 使用RNA assay

# 添加细胞通讯数据库（CellChatDB包含人和小鼠的配体-受体对）
CellChatDB <- CellChatDB.human  # 如果是小鼠数据，使用CellChatDB.mouse
cellchat_tumor@DB <- CellChatDB

# 对肿瘤样本进行细胞通讯分析
cellchat_tumor <- subsetData(cellchat_tumor)  # 子集化数据
cellchat_tumor <- identifyOverExpressedGenes(cellchat_tumor)  # 识别过表达基因
cellchat_tumor <- identifyOverExpressedInteractions(cellchat_tumor)  # 识别过表达的互作
cellchat_tumor <- projectData(cellchat_tumor, PPI.human)  # 投影到蛋白互作网络（可选）

# 计算通讯概率（核心步骤）
cellchat_tumor <- computeCommunProb(cellchat_tumor, 
                                    population.size = TRUE,  # 考虑细胞群体大小
                                    type = "truncatedMean",  # 计算方法
                                    trim = 0.1)  # 截断参数

# 过滤弱通讯
cellchat_tumor <- filterCommunication(cellchat_tumor, min.cells = 10)

# 计算整合的通讯网络
cellchat_tumor <- computeCommunProbPathway(cellchat_tumor)

# 计算整合的通讯网络（细胞类型水平）
cellchat_tumor <- aggregateNet(cellchat_tumor)

# 对正常样本重复同样步骤
cellchat_normal <- createCellChat(object = normal_seu,
                                  group.by = "celltype",
                                  assay = "RNA")
cellchat_normal@DB <- CellChatDB
cellchat_normal <- subsetData(cellchat_normal)
cellchat_normal <- identifyOverExpressedGenes(cellchat_normal)
cellchat_normal <- identifyOverExpressedInteractions(cellchat_normal)
cellchat_normal <- computeCommunProb(cellchat_normal, population.size = TRUE)
cellchat_normal <- filterCommunication(cellchat_normal, min.cells = 10)
cellchat_normal <- computeCommunProbPathway(cellchat_normal)
cellchat_normal <- aggregateNet(cellchat_normal)

#肿瘤 vs 正常：系统比较分析
# 合并两个CellChat对象进行比较
object.list <- list(Normal = cellchat_normal, Tumor = cellchat_tumor)
cellchat <- mergeCellChat(object.list, 
                          add.names = names(object.list),
                          cell.prefix = TRUE)

# 4.1 比较总体通讯强度和模式
# 整体通讯数量比较
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), 
                           measure = "count", color.use = c("#5B9BD5", "#ED7D31"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), 
                           measure = "weight", color.use = c("#5B9BD5", "#ED7D31"))
gg1 + gg2

# 4.2 比较通讯网络差异
# 差异通讯数量热图
netVisual_diffInteraction(cellchat, 
                          weight.scale = TRUE,
                          color.use = c("#5B9BD5", "#ED7D31"),
                          title.name = "Differential interaction numbers")

# 差异通讯强度热图
netVisual_diffInteraction(cellchat, 
                          measure = "weight",
                          weight.scale = TRUE,
                          color.use = c("#5B9BD5", "#ED7D31"),
                          title.name = "Differential interaction strength")

# 4.3 比较信号通路水平的通讯
# 识别差异活化的信号通路
gg <- rankNet(cellchat, 
              mode = "comparison", 
              stacked = TRUE, 
              color.use = c("#5B9BD5", "#ED7D31"),
              font.size = 10)
print(gg)

# 4.4 查看特定信号通路的变化
# 查看在肿瘤中显著增强或减弱的信号通路
pathways.show <- cellchat@netP$pathways  # 获取所有通路

# 选择最显著的几个通路进行比较
top_pathways <- pathways.show[1:4]  # 查看前4个通路

for (pathway in top_pathways) {
  # 并行显示肿瘤和正常的该通路通讯网络
  par(mfrow = c(1,2))
  netVisual_aggregate(cellchat_normal, 
                      signaling = pathway, 
                      layout = "circle",
                      color.use = colors_set3,  # 使用之前定义的配色
                      title.name = paste("Normal:", pathway))
  netVisual_aggregate(cellchat_tumor, 
                      signaling = pathway, 
                      layout = "circle",
                      color.use = colors_set3,
                      title.name = paste("Tumor:", pathway))
}
