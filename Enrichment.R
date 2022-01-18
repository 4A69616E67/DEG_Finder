

GOenrich <- function(){
  clusterProfiler::enrichGO()
}

library(openxlsx)
M6B_P53_changed_gene_list  <-openxlsx::read.xlsx("f:/Project/RRS-RNA/M6B_P53/res.3.6/DES_OUT.res.xlsx",sheet = "ChangedGene")
P53_CDKN2B_changed_gene_list <- openxlsx::read.xlsx("f:/Project/RRS-RNA/P53_cdkn2b/res.3.6/DES_OUT.res.xlsx",sheet = "ChangedGene")
merged_changed_gene_list <- merge(M6B_P53_changed_gene_list, P53_CDKN2B_changed_gene_list, by = "gene", all = TRUE)
merged_changed_gene_list[is.na(merged_changed_gene_list)] <- 0
col_names <- names(merged_changed_gene_list)
out_list<-list()
process_list <- merged_changed_gene_list[,c(1,grep("nevi",col_names))]
nevi_union_gene_list <- process_list[rowSums(abs(process_list[,2:dim(process_list)[2]]))>0,]
nevi_intersect_gene_list <- process_list[rowSums(abs(process_list[,2:dim(process_list)[2]]))>=(dim(process_list)[2]-1),]
out_list[["nevi_union"]] <- nevi_union_gene_list
out_list[["nevi_intersect"]] <- nevi_intersect_gene_list
process_list <- merged_changed_gene_list[,c(1,grep("M6B",col_names))]
M6B_union_gene_list <- process_list[rowSums(abs(process_list[,2:dim(process_list)[2]]))>0,]
M6B_intersect_gene_list <- process_list[rowSums(abs(process_list[,2:dim(process_list)[2]]))>=(dim(process_list)[2]-1),]
out_list[["M6B_union"]] <- M6B_union_gene_list
out_list[["M6B_intersect"]] <- M6B_intersect_gene_list
########################1
intersect_gene <- intersect(nevi_union_gene_list$gene,M6B_union_gene_list$gene)
out_list[["or_or_and"]] <- intersect_gene
result.go <- enrichGO(gene = mapIds(x = data_base, keys = intersect_gene, keytype = "SYMBOL", column = "ENTREZID"), keyType = "ENTREZID", ont = "ALL", OrgDb = data_base, readable = T)
out_list[["or_or_and_go"]] <- data.frame(result.go)
result.kegg <- enrichKEGG(gene = mapIds(x = data_base, keys = intersect_gene, keytype = "SYMBOL", column = "ENTREZID"),organism = "xtr")
out_list[["or_or_and_kegg"]] <- data.frame(result.kegg)
######################2
intersect_gene <- intersect(nevi_union_gene_list$gene,M6B_intersect_gene_list$gene)
out_list[["or_and_and_intersect"]] <- intersect_gene
result.go <- enrichGO(gene = mapIds(x = data_base, keys = intersect_gene, keytype = "SYMBOL", column = "ENTREZID"), keyType = "ENTREZID", ont = "ALL", OrgDb = data_base, readable = T)
out_list[["or_and_and__go"]] <- data.frame(result.go)
result.kegg <- enrichKEGG(gene = mapIds(x = data_base, keys = intersect_gene, keytype = "SYMBOL", column = "ENTREZID"),organism = "xtr")
out_list[["or_and_and__kegg"]] <- data.frame(result.kegg)
openxlsx::write.xlsx(out_list,file = "f:/Project/RRS-RNA/P53_cdkn2b/result.xlsx")


#########################################################趋势分析
library(TCseq)
gene_list <- read.table("f:/Project/Xenopus/XENTR_10.0_gene.bed",col.names = c("chr","start","end","id","score","strand"),colClasses = "character")
overlap_gene <- intersect(rownames(CountData),gene_list$id)
input_data <- as.matrix(na.omit(CountData[overlap_gene,]))
gene_list <- gene_list[match(overlap_gene,gene_list$id),]
gene_list$start <- as.numeric(gene_list$start)
gene_list$end <- as.numeric(gene_list$end)
input_Group <- data.frame(Group,timepoint=c("0d","0d","1d","1d","3d","3d","4d","4d","7d","7d","14d","14d","30d","30d"),group=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7))
names(Group)[1] <- "sampleid"
tca <- TCseq::TCA(design = input_Group,genomicFeature = gene_list,counts = input_data)
tca <- TCseq::DBanalysis(tca)
tca <- TCseq::timecourseTable(tca, value = "expression", norm.method = "rpkm", filter = TRUE)
tca <- TCseq::timeclust(tca, algo = "cm", k = 10, standardize = TRUE)
p <- timeclustplot(tca, value = "z-score(RPKM)", cols = 4)
a <- as.data.frame(tca@clusterRes@cluster)
