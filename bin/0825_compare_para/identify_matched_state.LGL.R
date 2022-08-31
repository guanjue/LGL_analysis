library(pheatmap)
library(dynamicTreeCut)

### read state para files
d0 = read.table('/Users/guanjuexiang/Documents/projects/Joint_Human_Mouse_IDEAS_State/VBSJ_052021_outputs_para_pdf/ES/06a_S3V2_IDEAS_hg38_r3_withHg38Mm10prior.para.modified.para', header=T, comment.char='~')
d1 = read.table('ideasLglMar2022.para', header=T, comment.char='~')

### get state signal
d1s = d1[,2:6]/d1[,1]
rownames(d1s) = paste0('LGL:',as.numeric(rownames(d1s))-1)
d0s = d0[,2:9]/d0[,1]
rownames(d0s) = paste0('J:',as.numeric(rownames(d0s))-1)

### get shared features
d0s1 = d0s[,is.element(colnames(d0s), colnames(d1s))]

### get pooled matrix
d01s = rbind(d0s1, d1s)
rownames_d01s_OD = rownames(d01s)

### hclust pooled matrix
d01s_hclust = hclust(dist(d01s), method='complete')

### cutreeDynamic
hclust_cCRE_DTC = cutreeDynamic(d01s_hclust, minClusterSize=1, deepSplit=4, dist= as.matrix(dist(d01s)), method='hybrid')

### modifications
hclust_cCRE_DTC[rownames(d01s)=='J:0'] = 0
hclust_cCRE_DTC[rownames(d01s)=='LGL:0'] = 0
hclust_cCRE_DTC[rownames(d01s)=='J:3'] = 17
hclust_cCRE_DTC[rownames(d01s)=='LGL:1'] = 17
hclust_cCRE_DTC[rownames(d01s)=='J:4'] = 18
hclust_cCRE_DTC[rownames(d01s)=='LGL:2'] = 18
hclust_cCRE_DTC[rownames(d01s)=='J:6'] = 19
hclust_cCRE_DTC[rownames(d01s)=='LGL:11'] = 19
hclust_cCRE_DTC[rownames(d01s)=='J:21'] = 20
hclust_cCRE_DTC[rownames(d01s)=='LGL:19'] = 20


table(hclust_cCRE_DTC)
cbind(rownames(d01s), hclust_cCRE_DTC)[order(hclust_cCRE_DTC),]
rownames(d01s) = paste(hclust_cCRE_DTC, rownames(d01s), sep='_')

### plot pooled state heatmap
breaksList = seq(0, 11, by = 0.001)
my_colorbar=colorRampPalette(c('white', 'red'))(n = length(breaksList))
pdf('test.LGL.VISION.state.DTC.pdf', height=18, width=6)
pheatmap(d01s[order(hclust_cCRE_DTC),], cex=1.5, color=my_colorbar, breaks = breaksList, cluster_rows=F)
dev.off()

### write output
output_mat = cbind(rownames_d01s_OD, hclust_cCRE_DTC)[order(hclust_cCRE_DTC),]
colnames(output_mat) = c('IDEAS_State_ID', 'merged_cluster_ID')
write.table(output_mat, 'LGL_state.VS.JointVISION_state.txt', quote=F, col.names=T, row.names=F, sep='\t')



