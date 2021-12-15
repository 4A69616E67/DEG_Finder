library("getopt", quietly = TRUE, warn.conflicts = F)

#region
get_configure <- function(opt, config_file) {
  config_data <- read.table(config_file, sep = "=", header = F, strip.white = T, encoding="UTF-8", colClasses = "character")
  for (i in 1:dim(config_data)[1]) {
    if(is.null(opt[[config_data[i,1]]]) && !config_data[i,2]==""){
      opt[[config_data[i,1]]] <- gsub("\\\\","/",config_data[i,2])
    }
  }
  for(name in names(default_value)){
    if (is.null(opt[[name]])) {
      opt[[name]] <- default_value[[name]]
    }
  }
  return(opt)
}

show_opt <- function(opt){
  for(name in names(opt)){
    cat(name," : ", opt[[name]], "\n")
  }
}

PreProcess <- function(data, group) {
  #-------------将featurecount的结果转换成需要的格式-----------
  GeneName <- data[[1]] # 先储存基因名
  GeneLength <- NULL # 储存基因长度（所有不重合的外显子长度之和）
  if (names(data)[2] == "Chr" && names(data)[6] == "Length") {
    GeneLength <- data.frame(Gene = GeneName, Length = data[["Length"]])
    row.names(GeneLength) <- GeneName # 行名替换成基因名
    data <- data[, -(2:6)]
  }
  data <- data[, -1]
  ColName <- names(data) # 先将列名储存
  data <- as.data.frame(lapply(data, as.numeric)) # 转换成数值，列名可能会变动
  names(data) <- ColName # 再将列名赋值回去
  row.names(data) <- GeneName # 转换成数值后重新赋值一遍行名
  # 输入格式转换完成, 添加分组
  # 分组至少有两列，分别是sample和condition
  # 先初始化分组信息
  # 判断是否为空
  if (!is.null(group)) {
    if (dim(group)[1] <= 0 || (!"sample" %in% names(group))) {
      cat("can not find corrected group content, use file header as sample name\n")
      group <- data.frame(sample = colnames(data), condition = gsub(x = colnames(data), pattern = "_[^_]*$", replacement = "", perl = T))
    }
    # 若无condition，则依据sample值设置
    if (!"condition" %in% names(group)) {
      cat("set condition depend on sample name\n")
      group <- data.frame(group, condition = gsub(x = group[["sample"]], pattern = "_[^_]*$", replacement = "", perl = T))
    }
  } else {
    # 如果group为空则使用输入数据的header作为sample
    cat("Empty group, use file header as sample name\n")
    group <- data.frame(sample = colnames(data), condition = gsub(x = colnames(data), pattern = "_[^_]*$", replacement = "", perl = T))
  }
  # 判断sample名是否一致
  if (!all(group[["sample"]] %in% names(data))) {
    # 不一致，退出
    cat("sample name include unkonow contnet\n")
    quit(status = 1)
  }
  data <- data[, match(Group$sample, names(data))] # 只保留sample中含有的样本,which返回的时前一个数组的索引值
  row.names(group) <- as.character(group$sample)
  if ("name" %in% names(group)) {
    # 如果有name词条，则将样本名改成name
    names(data) <- as.character(group$name)
    row.names(group) <- as.character(group$name)
  }
  data <- data[which(rowSums(data) > dim(data)[2]), ] # 去除全为0的行
  cat("Sample name:\t", as.character(group[["sample"]]), "\n")
  cat("Group name:\t", as.character(group[["condition"]]), "\n")
  res <- list(data = data, group = group, gene_length = GeneLength)
  return(res)
}
#endregion
# =================================计算TPM=======================================
CalculateTPM <- function(data, geneLen) {
  OverlapGene <- intersect(rownames(data), rownames(geneLen)) # 计算重叠的基因
  if (length(OverlapGene) < dim(data)[1]) {
    # 若不能得到全部的基因长度，输出警告信息
    warning("Can not get all gene's length!")
    warning(dim(data)[1] - length(OverlapGene), " genes will be remove!")
  }
  # 只保留有基因长度的基因
  data <- data[OverlapGene, ]
  geneLen <- geneLen[OverlapGene, ]
  cat("remove gene length bias\n")
  TPM_data <- data
  RPK_sum <- list()
  # 计算TPM
  for (j in names(TPM_data)) {
    TPM_data[[j]] <- TPM_data[[j]] / as.numeric(geneLen[row.names(TPM_data), 2])# 除以基因长度
    # 对数据量标准化
    RPK_sum[[j]]<-sum(TPM_data[[j]])
    TPM_data[[j]]<-TPM_data[[j]]/RPK_sum[[j]]*1000000
  }
  data <- ceiling(TPM_data) # 转换为整数
  return(data)
}
##===============================GO and KEGG function===========================
library(clusterProfiler)
process.go <- function(entrez_id, data_base, pvalue = 0.05, qvalue = 0.05, title_name = "GO Term", ...){
  result <- list()
  if (!is.null(entrez_id) && length(entrez_id) > 0) {
    result$ALL <- enrichGO(gene = entrez_id, OrgDb = data_base, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = pvalue, qvalueCutoff = qvalue, readable = T, ...)
    # result$BP <- enrichGO(gene = entrez_id, OrgDb = data_base, keyType = "ENTREZID", ont = "BP", pvalueCutoff = pvalue, qvalueCutoff = qvalue, readable = T, ...)
    # result$MF <- enrichGO(gene = entrez_id, OrgDb = data_base, keyType = "ENTREZID", ont = "MF", pvalueCutoff = pvalue, qvalueCutoff = qvalue, readable = T, ...)
    # result$CC <- enrichGO(gene = entrez_id, OrgDb = data_base, keyType = "ENTREZID", ont = "CC", pvalueCutoff = pvalue, qvalueCutoff = qvalue, readable = T, ...)
    
  }
  return(result)
}

process.kegg <- function(entrez_id, species, pvalue = 0.05, qvalue = 0.05, title_name = "KEGG Term", ...){
  result <- list()
  if (!is.null(entrez_id) && length(entrez_id) > 0) {
    result$KEGG <- enrichKEGG(gene = entrez_id, organism = species, pvalueCutoff = pvalue, qvalueCutoff = qvalue, ...)
    result$KEGG.plot <- dotplot(result$KEGG, title = title_name, showCategory = 30)
  }
}
##==========================================================================================================================================================
parameter <- matrix(c(
  "Input", "i", 1, "character", "<file> \t input file, gene experiment count matrix",
  "output", "o", 2, "character", "<dir> \t\t output dir (default ./)",
  "Prefix", "p", 2, "character", "<string> \t output prefix",
  "Control", "C", 1, "character", "<string> \t control name (must exist in group name)",
  "Treat", "T", 2, "character", "<string> \t treat names (defalut all group names exclude control name)",
  "Group", "g", 2, "character", "<file> \t group setting (include header and sample names must same as input file header)",
  "Species", "s", 2, "character", "<string> \t species database, used to process GO rich and KEGG rich (\"xtr\" or \"xla)\"",
  "GO", "G", 0, "logical", "\t\t execute GO term",
  "KEGG", "K", 0, "logical", "\t\t execute KEGG term",
  "BatchEffects", "B", 0, "logical", "\t\t remove Batch Effects",
  "TPM", "M", 0, "logical", "\t\t use tpm to replace raw count",
  "GenesLen", "L", 2, "character", "<file>\t\t genes transcripts length file",
  "Genes", "E", 2, "character", "<file> \t related gene list",
  "Config", "c", 2, "character", "<file> \t configure file"
),
byrow = T, ncol = 5
)
default_value <- list(OutDir="./",Prefix="DES_OUT",Species="xtr")
## options----------------------------------------------------------------------------------------------------------------------------------------------------
## ---------------------------
max_pvalue <- 0.01 # 最大p值
min_log2_fold_change <- 1 # 最小倍数差异
SignIndicator <- "pvalue" #使用p值做差异分析
## ---------------------------
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
# 从命令行接收参数
opt <- getopt(spec = parameter)
#------------------测试用参数设置-----------
# opt=list()
# opt$Input="f:/Project/RRS-RNA/透明Xenopu-Bulk RNAsequencing/result/TransparentFeatureCount.add.matrix"
# opt$output="f:/Project/RRS-RNA/透明Xenopu-Bulk RNAsequencing/result"
# opt$Prefix="DE_Analysis"
# opt$Control="mitf_PN_DS,mitf_PN_VS,mitf_PN_VS,mitf_ko6_NN_VS,mitf_ko6_PN_Eye"
# opt$Treat="mitf_ko6_NN_DS,mitf_ko6_NN_VS,mitf_PN_DS,mitf_ko6_NN_DS,mitf_ko6_NN_Eye"
# opt$group="f:/Project/RRS-RNA/透明Xenopu-Bulk RNAsequencing/scripts/sample_name.txt"
# # opt$TPM=F
# # opt$GenesLen="f:/Project/Xenopus/XENTR_10.0_geneLen.txt"
# opt$Species="xtr"
# opt$GO=T
# opt$KEGG=T
# opt$Config <- "f:/Project/RRS-RNA/透明Xenopu-Bulk RNAsequencing/scripts/config.txt"
if (!is.null(opt$Config) && file.exists(opt$Config)) {
  opt <- get_configure(opt,config_file = opt$Config)
}else{
  warning("No or incorrect configure file: ", opt$Config)
}
show_opt(opt)
options(nwarnings = 1000)
#-------------------------------------------
if (is.null(opt$Input) || is.null(opt$Control)) {
  stop(getopt(spec = parameter, usage = T))
}
## ---------------------------
library("DESeq2", quietly = TRUE, warn.conflicts = F)
library("ggplot2", quietly = TRUE, warn.conflicts = F)
library("pheatmap", quietly = TRUE, warn.conflicts = F)
library("openxlsx", quietly = TRUE, warn.conflicts = F)
## ---------------------------
InFile <- opt$Input
OutPrefix <- if (!is.null(opt$Prefix)) opt$Prefix else "DES_OUT" #输出前缀
OutDir <- if (!is.null(opt$OutDir)) opt$OutDir else "./" #输出目录
Species <- if (!is.null(opt$Species)) opt$Species else "xtr" #物种名 xtr xla hsa
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
if (!is.null(opt$TPM) && !is.null(GeneLen)) {
  CountData <- CalculateTPM(data = CountData, geneLen = GeneLen)
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf(paste(OutDir, "/", OutPrefix, ".plot.pdf", sep = ""), height = 10, width = 16) # 暂时注释
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
    if (!is.null(data_base) && !is.null(opt$GO)) {
      cat("process GO enrich:\t", aSample, "\t")
      
      # go_term <- enrichGO(gene = na.omit(entrez_id[up_gene_index]), OrgDb = data_base, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.01, qvalueCutoff = 0.05, minGSSize = 5, readable = T,universe = na.omit(entrez_id))
      go_result <- process.go(entrez_id = na.omit(entrez_id[up_gene_index]), data_base = data_base, pvalue = 0.05, qvalue = 0.05, universe = na.omit(entrez_id))
      go_term_list[[paste(title_name,"_Up_gene",sep = "")]] <- go_result$ALL
      if (!is.null(go_result$ALL)) {
        print(dotplot(go_result$ALL, title = paste(aSample, " vs ", aControl, " Up genes\tGO term", sep = ""), showCategory = 30))
      }
      # go_term <- enrichGO(gene = na.omit(entrez_id[down_gene_index]), OrgDb = data_base, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.01, qvalueCutoff = 0.05, minGSSize = 5, readable = T,universe = na.omit(entrez_id))
      go_result <- process.go(entrez_id = na.omit(entrez_id[down_gene_index]), data_base = data_base, pvalue = 0.05, qvalue = 0.05, universe = na.omit(entrez_id))
      go_term_list[[paste(title_name,"_Down_gene",sep = "")]] <- go_result$ALL
      # go_term <- enrichGO(gene = na.omit(entrez_id[down_gene_index]), OrgDb = data_base, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T, universe = na.omit(entrez_id))
      if (!is.null(go_result$ALL)) {
        print(dotplot(go_result$ALL, title = paste(aSample, " vs ", aControl, " Down genes\tGO term", sep = ""), showCategory = 30))
      }
    }
    if (!is.null(opt$KEGG) && !is.null(Species)) {
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
write.xlsx(go_term_list, file = paste(OutDir, "/", OutPrefix, ".go.xlsx", sep = ""))
dev.off()
quit()

# changed_name <- Changed_gene$gene
# changed_entrez_id <- as.character(gene_id_trans_table$Human_Entrez_ID[match(changed_name, gene_id_trans_table$XB_GENEPAGE_NAME)])
# go_result <- process.go(entrez_id = na.omit(changed_entrez_id[Changed_gene[[2]]==1 & Changed_gene[[3]]==0 & Changed_gene[[4]]==1]), data_base = data_base, pvalue = 0.05, qvalue = 0.05, universe = na.omit(changed_entrez_id))
# print(dotplot(go_result$ALL, title = "1,0,0", showCategory = 30))


