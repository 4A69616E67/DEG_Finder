library("getopt", quietly = TRUE, warn.conflicts = F)
library("plotly", quietly = TRUE, warn.conflicts = F)
source("Utils.R")

##==========================================================================================================================================================
parameter <- matrix(c(
  "Input", "i", 1, "character", "<file> \t input file, gene experiment count matrix",
  "OutDir", "o", 2, "character", "<dir> \t output dir (default ./)",
  "Prefix", "p", 2, "character", "<string> \t output prefix",
  "Control", "C", 1, "character", "<string> \t control name (must exist in group name)",
  "Treat", "T", 2, "character", "<string> \t treat names (defalut all group names exclude control name)",
  "Group", "g", 2, "character", "<file> \t group setting (include header and sample names must same as input file header)",
  "Species", "s", 2, "character", "<string> \t species database, used to process GO rich and KEGG rich (\"xtr\" or \"xla \" or hsa \")\"",
  "GO", "G", 0, "logical", "\t\t execute GO term",
  "KEGG", "K", 0, "logical", "\t\t execute KEGG term",
  "Genes", "E", 2, "character", "<file> \t related gene list",
  "Config", "c", 2, "character", "<file> \t configure file"
),
byrow = T, ncol = 5
)
##默认参数
# default_value <- list(Input=NULL, OutDir="./", Prefix="DES_OUT", Control=NULL, Treat=NULL, Group=NULL, GeneLens=NULL, Species="xtr", GO=F, KEGG=F, TPM=F, MaxP=0.01, MinLog2FC=1) #默认参数值
##测试参数
default_value <- list(Input="f:/Project/LGH/Heart-RNA/data/xt_featureCount.matrix",
                      OutDir="f:/Project/LGH/Heart-RNA/result/",
                      Prefix="DES_OUT",
                      Control="xt_wt",
                      Treat=NULL,
                      Group="f:/Project/LGH/Heart-RNA/script/sample.info.txt",
                      GeneLens="f:/Project/Xenopus/XENTR_10.0_geneLen.txt",
                      Species="xtr", GO=T, KEGG=F, TPM=T, MaxP=0.01, MinLog2FC=1) #测试参数
## options----------------------------------------------------------------------------------------------------------------------------------------------------
# 从命令行接收参数
opt <- getopt(spec = parameter)
###读取配置文件的信息
if(!is.null(opt$Config)){
  if (file.exists(opt$Config)) {
    opt <- get_configure(opt,config_file = opt$Config)
  }else{
    warning("Incorrect configure file: ", opt$Config)
  }
}
#########参数初始化
opt <- opt_init(opt, default_value)
options(nwarnings = 1000)
##############若必须参数的值为空，则打印帮助信息
if (is.null(opt$Input) || is.null(opt$Control)) {
  cat(getopt(spec = parameter, usage = T))
  stop()
}
show_opt(opt)
## ---------------------------
library("DESeq2", quietly = TRUE, warn.conflicts = F)
library("ggplot2", quietly = TRUE, warn.conflicts = F)
library("pheatmap", quietly = TRUE, warn.conflicts = F)
library("openxlsx", quietly = TRUE, warn.conflicts = F)
## ---------------------------
InFile <- opt$Input
OutPrefix <- opt$Prefix #输出前缀
OutDir <- opt$OutDir #输出目录
Species <- opt$Species  #物种名 xtr xla hsa
max_pvalue <- opt$MaxP # 最大p值
min_log2_fold_change <- opt$MinLog2FC # 最小倍数差异
SignIndicator <- "pvalue" #使用p值做差异分析
# factor的存储方式为，所有非重复的值存在level中，而factor中存的是level的索引。所以factor直接转character时可能会转换成对应的索引值，而不是字符串
RawData <- read.table(file = InFile, header = F, sep = "\t", colClasses = "character") # 默认情况下每一列都会是factor，这里不用factor
names(RawData) <- RawData[1, ] # 不直接读header，防止R改名
RawData <- RawData[-1, ] #去掉第一行
Group <- NULL 
# 判断condition文件是否存在
if (!is.null(opt$Group)) {
  fit <- try({Group <- read.table(file = opt$Group, header = T)}) # 使用try方法，即使读取错误也无妨
}
Pre <- PreProcess(data = RawData, group = Group) # 数据预处理
Group <- Pre$group
GeneLen <- Pre$gene_length
CountData <- Pre$data
#-------------------判断是否使用TPM值---------------------------------------------------------------------------------------------------------------------------
# 读取基因长度
if (!is.null(opt$GeneLens) && file.exists(opt$GeneLens)) {
  GeneLen <- read.table(opt$GeneLens, header = F)
  row.names(GeneLen) <- GeneLen[[1]]
  names(GeneLen) <- c("Gene", "Length")
}
if (opt$TPM && !is.null(GeneLen)) {
  TPMData <- CalculateTPM(data = CountData, geneLen = GeneLen)
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf(paste(OutDir, "/", OutPrefix, ".plot.pdf", sep = ""), height = 10, width = 16) # 暂时注释
pheatmap(cor(CountData,method = "pearson"), display_numbers = T, show_colnames = F)
dds <- DESeqDataSetFromMatrix(CountData, colData = Group, design = ~condition) # 由于data的列必须与group的行相同，所以得去掉group中没有的样本
dds <- DESeq(dds)
write.table(assay(normTransform(dds)), file = paste(OutDir, "/", OutPrefix, ".nor.tsv", sep = ""), quote = F, sep = "\t") # 打印标准化后的序列
# PCA分析
pcaData <- plotPCA(vst(dds, blind = F), intgroup = "condition")
pcaData
write.table(pcaData$data, file = paste(OutDir, "/", OutPrefix, ".pca.tsv", sep = ""), quote = F, sep = "\t", row.names = F) # 打印PCA的结果
# 识别或加载对照组和处理组的名字
Control <- strsplit(x = opt$Control, split = ",", fixed = T)[[1]]
TreatList <- levels(Group$condition)
if (!is.null(opt$Treat)) {
  TreatList <- strsplit(x = opt$Treat, split = ",", fixed = T)[[1]]
}

#------------------------------------------------------------加载数据库-----------------------------------------------------------------------------------------
# 判断物种并加载数据库
data_base <- NULL
library("AnnotationHub", quietly = T, warn.conflicts = F)
if (Species == "xla") {
  data_base <- AnnotationDbi::loadDb("Xenopus_laevis.db")
  # data_base<-AnnotationHub()[["AH75748"]]
  entrez_id <- mapIds(x = data_base, keys = rownames(CountData), keytype = "SYMBOL", column = "ENTREZID")
} else if (Species == "xtr") {
  data_base <- AnnotationDbi::loadDb("Xenopus_tropicalis.db")
  # data_base<-AnnotationHub()[["AH76241"]]
  entrez_id <- mapIds(x = data_base, keys = rownames(CountData), keytype = "SYMBOL", column = "ENTREZID")
} else if (Species == "hsa") {
  library(org.Hs.eg.db)
  data_base <- org.Hs.eg.db
  gene_id_trans_table <- read.table("f:/Project/Xenopus/gene_id_trans_table.txt", sep = "\t", header = T, quote = "")
  entrez_id <- as.character(gene_id_trans_table$Human_Entrez_ID[match(rownames(CountData), gene_id_trans_table$XB_GENEPAGE_NAME)]) ## 转换成人的Entrez id
} else{
  Species=NULL
}
# 加载GO和KEGG分析需要的包
if (!is.null(opt$GO) || !is.null(opt$KEGG)) {
  library("clusterProfiler", quietly = T, warn.conflicts = F)
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
foldChangeList <- data.frame(row.names = rownames(assay(dds)), gene = rownames(assay(dds)))
res_xlsx_list <- list()
Changed_gene_list <- list()
go_term_list <- list()
# 输出配置信息
cat("Max p value: ", max_pvalue, "\n")
cat("Min log2FoldChange: ", min_log2_fold_change, "\n")
#-------------------------进行差异分析----------------------------------------
for (i in 1:length(TreatList)) {
  aSample <- as.character(TreatList[i])
  aControl <- as.character(Control[(i - 1) %% length(Control) + 1])
  title_name <- paste(sep = "_vs_", aSample, aControl)
  if (aSample == aControl) {
    next
  }
  cat("process DESeq2:\t", aSample, "\t", aControl, "\n")
  res <- results(dds, contrast = c("condition", aSample, aControl)) # 获取差异表达结果
  res_dataframe <- cbind(data.frame(genes = rownames(res)), as.data.frame(res), significant = res$log2FoldChange)
  res_dataframe$significant <- factor("normal", levels = c("up", "normal", "down")) # 添加是否是差异基因
  gene_index <- which(res[[SignIndicator]] < max_pvalue & abs(res$log2FoldChange) > min_log2_fold_change)
  up_gene_index <- which(res[[SignIndicator]] < max_pvalue & res$log2FoldChange > min_log2_fold_change)
  res_dataframe$significant[up_gene_index] <- "up"
  down_gene_index <- which(res[[SignIndicator]] < max_pvalue & res$log2FoldChange < -min_log2_fold_change)
  res_dataframe$significant[down_gene_index] <- "down"
  print(ggplot(data = res_dataframe) + theme_classic() + geom_point(mapping = aes(x = log2(baseMean), y = log2FoldChange, colour = significant)) +
          ggtitle(title_name) + theme(plot.title = element_text(hjust = 0.5))) # 绘制MA图
  print(ggplot(data = res_dataframe) + theme_classic() + geom_point(mapping = aes(x = log2FoldChange, y = -log10(pvalue), colour = significant)) +
          ggtitle(title_name) + theme(plot.title = element_text(hjust = 0.5))) # 绘制火山图 
  res_xlsx_list[[paste(aSample, aControl, sep = "-")]] <- res_dataframe
  Changed_gene_list[[title_name]] <- data.frame(gene = c(rownames(res)[up_gene_index], rownames(res)[down_gene_index]), count = c(rep(1, length(up_gene_index)), rep(-1, length(down_gene_index))))
  names(Changed_gene_list[[title_name]]) <- c("gene", title_name)
  if (length(gene_index) > 0) {
    # 有差异表达基因才输出
    pheatmap(assay(normTransform(dds))[gene_index, which(Group$condition %in% c(aSample, aControl))], cluster_rows = T, show_rownames = F, cluster_cols = T, annotation_col = Group, fontsize = 6, main = paste(sep = "-", aSample, aControl))
    if (!is.null(data_base) && !is.null(opt$GO) && opt$GO) {
      cat("process GO enrich:\t", aSample, "\t")
      go_result <- process.go(entrez_id = na.omit(entrez_id[up_gene_index]), data_base = data_base, pvalue = 0.05, qvalue = 0.05, universe = na.omit(entrez_id))
      go_term_list[[paste(title_name,"_Up",sep = "")]] <- go_result$ALL
      if (!is.null(go_result$ALL)) {
        print(dotplot(go_result$ALL, title = paste(aSample, " vs ", aControl, " Up genes\tGO term", sep = ""), showCategory = 30))
        
      }
      go_result <- process.go(entrez_id = na.omit(entrez_id[down_gene_index]), data_base = data_base, pvalue = 0.05, qvalue = 0.05, universe = na.omit(entrez_id))
      go_term_list[[paste(title_name,"_Down",sep = "")]] <- go_result$ALL
      if (!is.null(go_result$ALL)) {
        print(dotplot(go_result$ALL, title = paste(aSample, " vs ", aControl, " Down genes\tGO term", sep = ""), showCategory = 30))
      }
    }
    if (!is.null(opt$KEGG) && !is.null(Species) && opt$KEGG) {
      cat("process KEGG enrich:\t", aSample)
      if (length(na.omit(entrez_id[up_gene_index])) > 0) {
        kegg_term <- enrichKEGG(gene = na.omit(entrez_id[up_gene_index]), organism = Species, qvalueCutoff = 0.05, universe = na.omit(entrez_id))
        if (!is.null(kegg_term)) {
          print(dotplot(kegg_term, title = paste(aSample, " vs ", aControl, " Up genes\tKEGG term", sep = ""), showCategory = 30))
        }
      }
      if (length(na.omit(entrez_id[down_gene_index])) > 0) {
        kegg_term <- enrichKEGG(gene = na.omit(entrez_id[down_gene_index]), organism = Species, qvalueCutoff = 0.05, universe = na.omit(entrez_id))
        if (!is.null(kegg_term)) {
          print(dotplot(kegg_term, title = paste(aSample, " vs ", aControl, " Down genes\tKEGG term", sep = ""), showCategory = 30))
        }
      }
    }
    cat("\n")
  } else {
    warning("Can not find different expression gene")
  }
  foldChangeList[[paste(aSample, aControl, sep = "-")]] <- res$log2FoldChange
}
test<-warnings()
#--------------------------------------------------------------------------------------
Changed_gene <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), Changed_gene_list)
Changed_gene[is.na(Changed_gene)] <- 0
res_xlsx_list[["ChangedGene"]] <- Changed_gene
res_xlsx_list[["FoldChangeList"]] <- foldChangeList
write.xlsx(res_xlsx_list, file = paste(OutDir, "/", OutPrefix, ".res.xlsx", sep = ""))
if (opt$GO) {
  write.xlsx(go_term_list, file = paste(OutDir, "/", OutPrefix, ".go.xlsx", sep = ""))
}
dev.off()
quit()

# changed_name <- Changed_gene$gene
# changed_entrez_id <- as.character(gene_id_trans_table$Human_Entrez_ID[match(changed_name, gene_id_trans_table$XB_GENEPAGE_NAME)])
# go_result <- process.go(entrez_id = na.omit(changed_entrez_id[Changed_gene[[2]]==1 & Changed_gene[[3]]==0 & Changed_gene[[4]]==1]), data_base = data_base, pvalue = 0.05, qvalue = 0.05, universe = na.omit(changed_entrez_id))
# print(dotplot(go_result$ALL, title = "1,0,0", showCategory = 30))


