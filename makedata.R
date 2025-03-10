library(maftools)
meta <- read.csv('/home/yincheng23/MultOmic/LIHC/Mutation/gdc_sample_sheet.2023-04-12.tsv',sep='\t')
meta$Case.ID <- unlist(lapply(meta$Case.ID, function(x) unique(strsplit(x, ', ')[[1]])))

#https://portal.gdc.cancer.gov/projects/TCGA-LIHC
meta2 <- read.csv('/home/yincheng23/MultOmic/LIHC/Lnet/clinical.tsv',sep='\t')
sub <- meta2[meta2$case_submitter_id %in% meta$Case.ID ,]
sub <- sub[,c('case_submitter_id','age_at_index','vital_status','days_to_death','ajcc_pathologic_m','ajcc_pathologic_n','ajcc_pathologic_stage')]
sub <- sub[!duplicated(sub),]
sub$days_to_death[sub$vital_status == 'Alive'] <- 0
sub$days_to_death <- as.numeric(sub$days_to_death)

keepdf <- data.frame()
dddf <- data.frame()
for(i in 1:nrow(sub)){
  if(nrow(sub[sub$case_submitter_id==sub$case_submitter_id[i],]) == 1){
    keepdf <- rbind(keepdf, sub[i,])
  } else {
    sub_d <- sub[sub$case_submitter_id==sub$case_submitter_id[i],]
    sub_d$nm <- rowSums(sub_d[,5:7] != "'--")
    sub_d <- sub_d[sub_d$nm == max(sub_d$nm), 1:7]
    if(nrow(sub_d) > 1){sub_d <- sub_d[1,]}
    dddf <- rbind(dddf, sub_d)
  }
}
clinical <- rbind(keepdf, dddf)
clinical <- clinical[!duplicated(clinical),]
rownames(clinical) <- clinical$case_submitter_id
rm(dddf,keepdf,meta2,sub,sub_d)

myRowSums <- function(x) {
  if(is.null(ncol(x))){
    return(x)
  } else {
    return(rowSums(x))
  }
}

Frame_Shift <- data.frame()
In_Frame <- data.frame()
Other_Type <- data.frame()
Missense <- data.frame()
Frame_Shift_kw <- c("Frame_Shift_Del","Frame_Shift_Ins")
In_Frame_kw <- c("In_Frame_Del","In_Frame_Ins")
Other_Type_kw <- c("Splice_Site","Nonstop_Mutation")
root <- '~/MultOmic/LIHC/Mutation'
for(i in 1:nrow(meta)){
  #meta$Case.ID
  filepath <- paste(c(root,meta$File.ID[i],meta$File.Name[i]),collapse = '/')
  data <- read.maf(filepath)
  mutinfo <- data.frame(data@gene.summary)
  gid <- mutinfo$Hugo_Symbol
  mutinfo <- mutinfo[,2:ncol(mutinfo)]
  tmp <- myRowSums(mutinfo[,colnames(mutinfo) %in% Frame_Shift_kw])
  Frame_Shift <- rbind(Frame_Shift, data.frame('Case_ID'=meta$Case.ID[i],'Gene'=gid, 'Count'=tmp))
  
  tmp <- myRowSums(mutinfo[,colnames(mutinfo) %in% In_Frame_kw])
  In_Frame <- rbind(In_Frame, data.frame('Case_ID'=meta$Case.ID[i],'Gene'=gid, 'Count'=tmp))
  
  tmp <- myRowSums(mutinfo[,colnames(mutinfo) %in% Other_Type_kw])
  Other_Type <- rbind(Other_Type, data.frame('Case_ID'=meta$Case.ID[i],'Gene'=gid, 'Count'=tmp))
  
  tmp <- myRowSums(mutinfo[,colnames(mutinfo) %in% c("Missense_Mutation")])
  Missense <- rbind(Missense, data.frame('Case_ID'=meta$Case.ID[i],'Gene'=gid, 'Count'=tmp))
}

Frame_Shift <- as.data.frame.matrix(xtabs(Count ~ Case_ID + Gene, data = Frame_Shift))
In_Frame    <- as.data.frame.matrix(xtabs(Count ~ Case_ID + Gene, data = In_Frame))
Other_Type  <- as.data.frame.matrix(xtabs(Count ~ Case_ID + Gene, data = Other_Type))
Missense    <- as.data.frame.matrix(xtabs(Count ~ Case_ID + Gene, data = Missense))

sp <- intersect(rownames(clinical), rownames(Frame_Shift))
clinical <- clinical[sp,]
Frame_Shift <- Frame_Shift[sp,]
In_Frame <- In_Frame[sp,]
Other_Type <- Other_Type[sp,]
Missense <- Missense[sp,]

write.csv(clinical,'/home/yincheng23/MultOmic/LIHC/Lnet/my_clinical.csv')
write.csv(Frame_Shift,'/home/yincheng23/MultOmic/LIHC/Lnet/frame_shift.csv')
write.csv(In_Frame,'/home/yincheng23/MultOmic/LIHC/Lnet/in_frame.csv')
write.csv(Other_Type,'/home/yincheng23/MultOmic/LIHC/Lnet/other_type.csv')
write.csv(Missense,'/home/yincheng23/MultOmic/LIHC/Lnet/missense.csv')

rowSums(Frame_Shift[Frame_Shift])
rowSums(In_Frame)
rowSums(Other_Type)
rowSums(Missense)

plot(rowSums(Frame_Shift), rowSums(Missense))
plot(rowSums(In_Frame), rowSums(Missense))
plot(rowSums(Other_Type), rowSums(Missense))
plot(rowSums(Missense), rowSums(Missense))


meta <- read.csv('/home/yincheng23/MultOmic/LIHC/RNA/gdc_sample_sheet.2023-04-12.tsv',sep='\t')
meta$Case.ID <- unlist(lapply(meta$Case.ID, function(x) unique(strsplit(x, ', ')[[1]])))
meta <- meta[meta$Sample.Type == 'Primary Tumor',]
meta <- meta[meta$Case.ID %in% sp, ]
root <- '~/MultOmic/LIHC/RNA'
RNA <- data.frame()
for(i in 1:nrow(meta)){
  filepath <- paste(c(root, gsub("-", replacement="_", meta$Sample.ID[i]), meta$File.Name[i]),collapse = '/')
  data <- read.delim2(filepath, sep='\t', comment.char = '#')
  data <- data[data$gene_type != "",]
  RNA <- rbind(RNA, data.frame('Case_ID'=meta$Case.ID[i],'Gene'=data$gene_name, 'TPM'=as.numeric(data$tpm_unstranded)))
}
RNA <- as.data.frame.matrix(xtabs(TPM ~ Case_ID + Gene, data = RNA))
write.csv(RNA,'/home/yincheng23/MultOmic/LIHC/Lnet/mRNA.csv')
