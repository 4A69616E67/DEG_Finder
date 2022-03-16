###数据预处理
PreProcess <- function(data, condition) {
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
  if (!is.null(condition)) {
    if (dim(condition)[1] <= 0 || (!"sampleid" %in% names(condition))) {
      cat("can not find corrected group content, use file header as sample name\n")
      condition <- data.frame(sample = colnames(data), condition = gsub(x = colnames(data), pattern = "_[^_]*$", replacement = "", perl = T))
    }
    # 若无group，则依据sample值设置
    if (!"group" %in% names(condition)) {
      cat("set group depend on sample name\n")
      condition <- data.frame(condition, condition = gsub(x = condition[["sampleid"]], pattern = "_[^_]*$", replacement = "", perl = T))
    }
  } else {
    # 如果group为空则使用输入数据的header作为sample
    cat("Empty group, use file header as sample name\n")
    condition <- data.frame(sample = colnames(data), group = gsub(x = colnames(data), pattern = "_[^_]*$", replacement = "", perl = T))
  }
  condition$group <- factor(condition$group)
  # 判断sample名是否一致
  if (!all(condition[["sampleid"]] %in% names(data))) {
    # 不一致，退出
    cat("sample name include unkonow contnet\n")
    quit(status = 1)
  }
  data <- data[, match(condition$sampleid, names(data))] # 只保留sample中含有的样本,match返回的是前一个数组的值在后一个数组的索引值
  row.names(condition) <- as.character(condition$sampleid)
  if ("name" %in% names(condition)) {
    # 如果有name词条，则将样本名改成name
    names(data) <- as.character(condition$name)
    row.names(condition) <- as.character(condition$name)
  }
  data <- data[which(rowSums(data) > dim(data)[2]), ] # 去除count值较小的行
  cat("Sample name:\t", as.character(condition[["sampleid"]]), "\n")
  cat("Group name:\t", as.character(condition[["group"]]), "\n")
  res <- list(data = data, condition = condition, gene_length = GeneLength)
  return(res)
}

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
  length_set <- as.numeric(geneLen[row.names(TPM_data), 2])
  # 计算TPM
  for (j in names(TPM_data)) {
    TPM_data[[j]] <- TPM_data[[j]] / length_set * 10^3# 除以基因长度
    # 对数据量标准化
    RPK_sum[[j]]<- sum(TPM_data[[j]])
    TPM_data[[j]]<- TPM_data[[j]] / RPK_sum[[j]] * 10^6
  }
  # data <- ceiling(TPM_data) # 转换为整数
  return(TPM_data)
}

###GO富集
process.go <- function(entrez_id, data_base, pvalue = 0.05, qvalue = 0.05, ...){
  result <- list()
  if (!is.null(entrez_id) && length(entrez_id) > 0) {
    result$ALL <- enrichGO(gene = entrez_id, OrgDb = data_base, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = pvalue, qvalueCutoff = qvalue, readable = T, ...)
  }
  return(result)
}

plot.go <- function(go_result, ont=c("BP","MF","CC"), ...){
  plot_data <- go_result@result
    if (dim(plot_data)[1] <= 0) {
     return()
  }
  
  
}

###KEGG富集
process.kegg <- function(entrez_id, species, pvalue = 0.05, qvalue = 0.05, title= "",  show_category =20, ...){
  result <- list()
  if (!is.null(entrez_id) && length(entrez_id) > 0) {
    result$term <- enrichKEGG(gene = entrez_id, organism = species, pvalueCutoff = pvalue, qvalueCutoff = qvalue, ...)
    if (!is.null(result$term ) && dim(result$term )[1] > 0) {
      result$plot <- barplot(result$term, title = title, showCategory = show_category)
    }
  }
  return(result)
}

###获取配置文件的参数
get_configure <- function(opt, config_file) {
  config_data <- read.table(config_file, sep = "=", header = F, strip.white = T, encoding="UTF-8", colClasses = "character")
  for (i in 1:dim(config_data)[1]) {
    if(is.null(opt[[config_data[i,1]]]) && !config_data[i,2]==""){
      opt[[config_data[i,1]]] <- gsub("\\\\","/",config_data[i,2])
    }
  }
  return(opt)
}

###参数初始化
opt_init <- function(opt, default_value){
  for(name in names(default_value)){
    if (is.null(opt[[name]])) {
      opt[[name]] <- default_value[[name]]
    }
  }
  opt$GO <- if(is.na(as.logical(opt$GO)) || is.null(opt$GO)) default_value$GO else as.logical(opt$GO)
  opt$KEGG <- if(is.na(as.logical(opt$KEGG)) || is.null(opt$KEGG)) default_value$KEGG else as.logical(opt$KEGG)
  opt$TPM <- if(is.na(as.logical(opt$TPM)) || is.null(opt$TPM)) default_value$TPM else as.logical(opt$TPM)
  opt$MaxP <- if(is.na(as.numeric(opt$MaxP)) || is.null(opt$MaxP)) default_value$MaxP else as.numeric(opt$MaxP)
  opt$MinLog2FC <- if(is.na(as.numeric(opt$MinLog2FC)) || is.null(opt$MinLog2FC)) default_value$MinLog2FC else as.numeric(opt$MinLog2FC)
  return(opt)
}

###显示参数值
show_opt <- function(opt){
  for(name in names(opt)){
    cat(name," : ", opt[[name]], "\n")
  }
}
############趋势分析
trendAnalyse <- function(data, condition, k=10, algo="cm", sed="12580"){
  library(TCseq)
  res <- list()
  if ("timegroup" %in% names(condition)) {
    group_list <- levels(factor(condition[["timegroup"]]))
  }else{
    group_list <- c("1")
  }
  for(time_group in group_list){
    
  }
  # tca <- TCseq::TCA(design = condition, genomicFeature = gene_list, counts = input_data)
  # tca <- TCseq::DBanalysis(tca)
  # tca <- TCseq::timecourseTable(tca, value = "expression", norm.method = "none", filter = TRUE)
  set.seed(sed)
  tca <- TCseq::timeclust(as.matrix(data), algo = algo, k = k, standardize = TRUE)
  res$tca <- tca
  res$plot <- timeclustplot(tca, value = "z-score(RPKM)", cols = 4)
  # a <- as.data.frame(tca@clusterRes@cluster)
}
