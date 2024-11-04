####做acnes和相关散点图 r大于0.3的#####
CIBERSORTx_Job28_Results_DM22_ADPLUSPSO <- read.csv("~/KI-POSTDOC/MAARS data-microarry/cibersort/for MAARS/CLEAN DATA RERUN/RESULTS/CIBERSORTx_Job28_Results_DM22_ADPLUSPSO.csv")
CIBERSORTx_Job29_Results_LM22_ADPLUSPSO <- read.csv("~/KI-POSTDOC/MAARS data-microarry/cibersort/for MAARS/CLEAN DATA RERUN/RESULTS/CIBERSORTx_Job29_Results_LM22_ADPLUSPSO.csv")
#AD
load("step2_css_ad_tight_speciesl_funga_0505.Rdata")
funga_ad_tight_species_norm_log
otu_funga_logcss_species
df_input_data<- t(otu_funga_logcss_species)
df_input_data2<-data.frame(df_input_data)
df_input_data2$scilife_id<-rownames(df_input_data2)
#DM22
pheno_microarray_all <- read.csv("~/KI-POSTDOC/Funga_maars/pheno_microarray_all.csv")
sampledf<-merge(sampledf,pheno_microarray_all,by="scilife_id")

datTraits_PSO_LB<-merge(df_input_data2,sampledf,by="scilife_id")
datTraits_PSO_LB_select_LM22<-merge(datTraits_PSO_LB,CIBERSORTx_Job29_Results_LM22_ADPLUSPSO,by.x="sample_id",by.y="Mixture")
datTraits_PSO_LB_select_DM22<-merge(datTraits_PSO_LB,CIBERSORTx_Job28_Results_DM22_ADPLUSPSO,by.x="sample_id",by.y="Mixture")
rownames(datTraits_PSO_LB_select_DM22)<-datTraits_PSO_LB_select_DM22$scilife_id
rownames(datTraits_PSO_LB_select_LM22)<-datTraits_PSO_LB_select_LM22$scilife_id
if(T){
  #cortest_psy <- corr.test(CIBERSORTx_Job29_Results_LM22_PSO_LB[,c(12:285)], CIBERSORTx_Job29_Results_LM22_PSO_LB[,c(286:310)], method = "spearman")
  cortest_psy_adj <- corr.test(datTraits_PSO_LB_select_DM22[,c(3:98)], datTraits_PSO_LB_select_DM22[,c(120:141)], method = "spearman", adjust = "fdr")
  moduleTraitCor<-round(cortest_psy_adj$r,2)
  moduleTraitPvalue<-round(cortest_psy_adj$p,2) 
  moduleTraitPvalue_adjust<-round(cortest_psy_adj$p.adj,digits = 2)
  # moduleTraitCor<-moduleTraitCor[,-c(2,17,22)]
  # moduleTraitPvalue<-moduleTraitPvalue[,-c(2,17,22)]
  # moduleTraitPvalue_adjust<-moduleTraitPvalue_adjust[,-c(2,17,22)]
  write.csv(moduleTraitCor,file = "fungi_ciber_dm22_association_ad.csv")
  write.csv(moduleTraitPvalue,file = "fungi_ciber_dm22_ad_pvalue.csv")
  write.csv(moduleTraitPvalue_adjust,file = "fungi_ciber_dm22_ad_padjust.csv")
  
  moduleTraitCor_res<-filter_all(data.frame(moduleTraitCor),any_vars(abs(.)>0.3))
  moduleTraitPvalue_adjust_res<-filter_all(data.frame(moduleTraitPvalue_adjust),any_vars(.<0.05))
  res_id<-rownames(moduleTraitCor_res)
  moduleTraitPvalue_adjust_res<-moduleTraitPvalue_adjust_res[res_id,]
  max(moduleTraitCor_res)
  min(moduleTraitCor_res)
  moduleTraitPvalue_adjust_res<-moduleTraitPvalue_adjust_res[,-c(2,15)]
  moduleTraitCor_res<-moduleTraitCor_res[,-c(2,15)]
  write.csv(moduleTraitCor_res,file = "fungi_ad_dm22_asso.csv")
  write.csv(moduleTraitPvalue_adjust_res,file = "fungi_ad_dm22_padjust.csv")
 
 
  png(filename = "mpa_cibersort_dm22_correlation_les_ctrl.png",width = 6000,height = 6000,res = 300)
  col<-colorRampPalette(c("blue","white","red"))(10)
  corrplot(as.matrix(moduleTraitCor_res), method = "square",tl.cex = 0.4,col = col,tl.col = "black",
           p.mat=as.matrix(moduleTraitPvalue_adjust_res),#传入相关系数的P值
           sig.level = 0.05,#设置显著的水平
           insig = "blank",
  ) #0505
  dev.off()
}
save(pheno_microarray_all,sampledf,
     datTraits_PSO_LB_select_DM22,cortest_psy_adj,file="fungi_ad_DM22_asso.rdata")
# moduleTraitCor<-data.frame(t(moduleTraitCor))
# #做acnes和相关散点图 r大于0.3的
# colnames(datTraits_PSO_LB_select_DM22)[111]<-"C.simulans"
# colnames(datTraits_PSO_LB_select_DM22)[114]<-"C.acnes"
# mydf_simulans_keratino_dm22<-datTraits_PSO_LB_select_DM22[,c(111,145)]
# mydf_acnes_keratino_dm22<-datTraits_PSO_LB_select_DM22[,c(114,145)]
# p<-ggscatter(mydf_simulans_keratino_dm22, x = "C.simulans", y ="keratinocytes" ,  #keratinocytes myeloid.dendritic.cells sc.adipocytes
#              color = "black", size = 3, # 点的颜色与大小
#              add = "reg.line", # 添加回归线
#              add.params = list(color = "red", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
#              conf.int = TRUE, # 添加回归线的置信区间
#              cor.coef = TRUE, # 添加相关系数
#              cor.coef.size = 7,
#              cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")) + #选择Pearson相关
#   theme(axis.title.x=element_text(face="bold",size=25),
#         axis.title.y = element_text(face="bold",size=25))           
# #theme(text=element_text(face="bold",size=20))
#partial spearman
if(T){
  
  library(ppcor)
  datExpr_hubgene_les_ctrl$sex<-ifelse(datExpr_hubgene_les_ctrl$sex.x=="Male",0,1)
  datTraits_PSO_LB_partial<-cbind(MEs_rmgrey,datTraits_PSO_LB)
  colnames(datExpr_hubgene_les_ctrl)[17]<-"age"
  datExpr_hubgene_les_ctrl<-datExpr_hubgene_les_ctrl[,c(1:24,559,72:88,25:71,89:558)]
  #bac 43:89 , gene 90:559 47 470
  datExpr_hubgene_les_ctrl_female<-subset(datExpr_hubgene_les_ctrl,sex.x=="Female")
  moduleTraitCor_LB_partial_temp<-as.data.frame(matrix(nrow=47,ncol = 470)) #nroe是module数 ncol是细胞数
  moduleTraitPvalue_LB_partial_temp<-as.data.frame(matrix(nrow=47,ncol = 470))
  #rm(list = c("moduleTraitPvalue_LB_partial_temp"))
  #moduleTraitPvalue_LB_partial_temp_test<-as.data.frame(matrix(nrow=7,ncol = 22))
  for (i in c(43:89)) {
    
    for (j in c(90:559)) {
      
      moduleTraitCor_LB_partial_temp[i-42,j-89]<-pcor.test(datExpr_hubgene_les_ctrl_female[,i], datExpr_hubgene_les_ctrl_female[,j], datExpr_hubgene_les_ctrl_female[,c( "age")],method = "spearman")[1,1]
      moduleTraitPvalue_LB_partial_temp[i-42,j-89]<-pcor.test(datExpr_hubgene_les_ctrl_female[,i], datExpr_hubgene_les_ctrl_female[,j], datExpr_hubgene_les_ctrl_female[,c("age")],method = "spearman")[1,2]
    }
    
  }
  table(datExpr_hubgene_les_ctrl$sex)
  rownames(moduleTraitCor_LB_partial_temp)<-colnames(datExpr_hubgene_les_ctrl)[43:89]
  colnames(moduleTraitCor_LB_partial_temp)<-colnames(datExpr_hubgene_les_ctrl)[90:559]
  moduleTraitCor_LB_partial_temp<-data.frame(t(moduleTraitCor_LB_partial_temp))
  write.csv(moduleTraitCor_LB_partial_temp,file = "mpa_tur_hubgenes_partial_spearman_sexonly.csv")
  
  rownames(moduleTraitPvalue_LB_partial_temp)<-colnames(datExpr_hubgene_les_ctrl)[43:89]
  colnames(moduleTraitPvalue_LB_partial_temp)<-colnames(datExpr_hubgene_les_ctrl)[90:559]
  moduleTraitPvalue_LB_partial_temp<-data.frame(t(moduleTraitPvalue_LB_partial_temp))
  write.csv(moduleTraitPvalue_LB_partial_temp,file = "mpa_tur_hubgenes_partial_spearman_pvalue_sexonly.csv")
  
  
  
  
  
}
#LM22
if(T){
  #cortest_psy <- corr.test(CIBERSORTx_Job29_Results_LM22_PSO_LB[,c(12:285)], CIBERSORTx_Job29_Results_LM22_PSO_LB[,c(286:310)], method = "spearman")
  cortest_psy_adj <- corr.test(datTraits_PSO_LB_select_LM22[,c(3:98)], datTraits_PSO_LB_select_LM22[,c(120:141)], method = "spearman", adjust = "fdr")
  moduleTraitCor<-round(cortest_psy_adj$r,2)
  moduleTraitPvalue<-round(cortest_psy_adj$p,2) 
  moduleTraitPvalue_adjust<-round(cortest_psy_adj$p.adj,digits = 2)
  # moduleTraitCor<-moduleTraitCor[,-c(2,17,22)]
  # moduleTraitPvalue<-moduleTraitPvalue[,-c(2,17,22)]
  # moduleTraitPvalue_adjust<-moduleTraitPvalue_adjust[,-c(2,17,22)]
  write.csv(moduleTraitCor,file = "fungi_ciber_lm22_association_ad.csv")
  write.csv(moduleTraitPvalue,file = "fungi_ciber_lm22_ad_pvalue.csv")
  write.csv(moduleTraitPvalue_adjust,file = "fungi_ciber_lm22_ad_padjust.csv")
  
  moduleTraitCor_res<-filter_all(data.frame(moduleTraitCor),any_vars(abs(.)>0.3))
  moduleTraitPvalue_adjust_res<-filter_all(data.frame(moduleTraitPvalue_adjust),any_vars(.<0.05))
  res_id<-rownames(moduleTraitCor_res)
  moduleTraitPvalue_adjust_res<-moduleTraitPvalue_adjust_res[res_id,]
  max(moduleTraitCor_res)
  min(moduleTraitCor_res)

  write.csv(moduleTraitCor_res,file = "fungi_ad_lm22_asso.csv")
  write.csv(moduleTraitPvalue_adjust_res,file = "fungi_ad_lm22_padjust.csv")
  
  
  png(filename = "mpa_cibersort_lm22_correlation_les_ctrl.png",width = 6000,height = 6000,res = 300)
  col<-colorRampPalette(c("blue","white","red"))(10)
  corrplot(as.matrix(moduleTraitCor_res), method = "square",tl.cex = 0.4,col = col,tl.col = "black",
           p.mat=as.matrix(moduleTraitPvalue_adjust_res),#传入相关系数的P值
           sig.level = 0.05,#设置显著的水平
           insig = "blank",
  ) #0505
  dev.off()
}
save(pheno_microarray_all,sampledf,
     datTraits_PSO_LB_select_LM22,cortest_psy_adj,file="fungi_ad_LM22_asso.rdata")
#PSO
load("step2_css_pso_lb_speciesl_funga_0505.Rdata")
funga_lb_pso_species_norm_log
otu_funga_logcss_species<-data.frame(funga_lb_pso_species_norm_log@otu_table)
df_input_data<- t(otu_funga_logcss_species)
sampledf<-data.frame(sample_data(funga_lb_pso_species_norm_log))
sampledf[is.na(sampledf)]<-0
df_input_data<- t(otu_funga_logcss_species)
df_input_data2<-data.frame(df_input_data)
df_input_data2$scilife_id<-rownames(df_input_data2)
#DM22
pheno_microarray_all <- read.csv("~/KI-POSTDOC/Funga_maars/pheno_microarray_all.csv")
sampledf<-merge(sampledf,pheno_microarray_all,by="scilife_id")

datTraits_PSO_LB<-merge(df_input_data2,sampledf,by="scilife_id")
datTraits_PSO_LB_select_LM22<-merge(datTraits_PSO_LB,CIBERSORTx_Job29_Results_LM22_ADPLUSPSO,by.x="sample_id",by.y="Mixture")
datTraits_PSO_LB_select_DM22<-merge(datTraits_PSO_LB,CIBERSORTx_Job28_Results_DM22_ADPLUSPSO,by.x="sample_id",by.y="Mixture")
rownames(datTraits_PSO_LB_select_DM22)<-datTraits_PSO_LB_select_DM22$scilife_id
rownames(datTraits_PSO_LB_select_LM22)<-datTraits_PSO_LB_select_LM22$scilife_id
if(T){
  #cortest_psy <- corr.test(CIBERSORTx_Job29_Results_LM22_PSO_LB[,c(12:285)], CIBERSORTx_Job29_Results_LM22_PSO_LB[,c(286:310)], method = "spearman")
  cortest_psy_adj <- corr.test(datTraits_PSO_LB_select_DM22[,c(3:14)], datTraits_PSO_LB_select_DM22[,c(36:57)], method = "spearman", adjust = "fdr")
  moduleTraitCor<-round(cortest_psy_adj$r,2)
  moduleTraitPvalue<-round(cortest_psy_adj$p,2) 
  moduleTraitPvalue_adjust<-round(cortest_psy_adj$p.adj,digits = 2)

  write.csv(moduleTraitCor,file = "fungi_ciber_dm22_association_pso.csv")
  write.csv(moduleTraitPvalue,file = "fungi_ciber_dm22_pso_pvalue.csv")
  write.csv(moduleTraitPvalue_adjust,file = "fungi_ciber_dm22_pso_padjust.csv")
  
  moduleTraitCor_res<-filter_all(data.frame(moduleTraitCor),any_vars(abs(.)>0.3))
  moduleTraitPvalue_adjust_res<-filter_all(data.frame(moduleTraitPvalue_adjust),any_vars(.<0.05))
  res_id<-rownames(moduleTraitCor_res)
  moduleTraitPvalue_adjust_res<-moduleTraitPvalue_adjust_res[res_id,]
  max(moduleTraitCor_res)
  min(moduleTraitCor_res)
  moduleTraitPvalue_adjust_res<-moduleTraitPvalue_adjust_res[,-c(2,16,22)]
  moduleTraitCor_res<-moduleTraitCor_res[,-c(2,16,22)]
  write.csv(moduleTraitCor_res,file = "fungi_pso_dm22_asso.csv")
  write.csv(moduleTraitPvalue_adjust_res,file = "fungi_pso_dm22_padjust.csv")
  
  
  #png(filename = "mpa_cibersort_dm22_correlation_les_ctrl.png",width = 6000,height = 6000,res = 300)
  col<-colorRampPalette(c("blue","white","red"))(10)
  corrplot(as.matrix(moduleTraitCor_res), method = "square",tl.cex = 0.4,col = col,tl.col = "black",
           p.mat=as.matrix(moduleTraitPvalue_adjust_res),#传入相关系数的P值
           sig.level = 0.05,#设置显著的水平
           insig = "blank",
  ) #0505
  dev.off()
}
# save(pheno_microarray_all,sampledf,
#      datTraits_PSO_LB_select_DM22,cortest_psy_adj,file="fungi_pso_DM22_asso.rdata")

#partial spearman
if(T){
  
  library(ppcor)
  datExpr_hubgene_les_ctrl$sex<-ifelse(datExpr_hubgene_les_ctrl$sex.x=="Male",0,1)
  datTraits_PSO_LB_partial<-cbind(MEs_rmgrey,datTraits_PSO_LB)
  colnames(datExpr_hubgene_les_ctrl)[17]<-"age"
  datExpr_hubgene_les_ctrl<-datExpr_hubgene_les_ctrl[,c(1:24,559,72:88,25:71,89:558)]
  #bac 43:89 , gene 90:559 47 470
  datExpr_hubgene_les_ctrl_female<-subset(datExpr_hubgene_les_ctrl,sex.x=="Female")
  moduleTraitCor_LB_partial_temp<-as.data.frame(matrix(nrow=47,ncol = 470)) #nroe是module数 ncol是细胞数
  moduleTraitPvalue_LB_partial_temp<-as.data.frame(matrix(nrow=47,ncol = 470))
  #rm(list = c("moduleTraitPvalue_LB_partial_temp"))
  #moduleTraitPvalue_LB_partial_temp_test<-as.data.frame(matrix(nrow=7,ncol = 22))
  for (i in c(43:89)) {
    
    for (j in c(90:559)) {
      
      moduleTraitCor_LB_partial_temp[i-42,j-89]<-pcor.test(datExpr_hubgene_les_ctrl_female[,i], datExpr_hubgene_les_ctrl_female[,j], datExpr_hubgene_les_ctrl_female[,c( "age")],method = "spearman")[1,1]
      moduleTraitPvalue_LB_partial_temp[i-42,j-89]<-pcor.test(datExpr_hubgene_les_ctrl_female[,i], datExpr_hubgene_les_ctrl_female[,j], datExpr_hubgene_les_ctrl_female[,c("age")],method = "spearman")[1,2]
    }
    
  }
  table(datExpr_hubgene_les_ctrl$sex)
  rownames(moduleTraitCor_LB_partial_temp)<-colnames(datExpr_hubgene_les_ctrl)[43:89]
  colnames(moduleTraitCor_LB_partial_temp)<-colnames(datExpr_hubgene_les_ctrl)[90:559]
  moduleTraitCor_LB_partial_temp<-data.frame(t(moduleTraitCor_LB_partial_temp))
  write.csv(moduleTraitCor_LB_partial_temp,file = "mpa_tur_hubgenes_partial_spearman_sexonly.csv")
  
  rownames(moduleTraitPvalue_LB_partial_temp)<-colnames(datExpr_hubgene_les_ctrl)[43:89]
  colnames(moduleTraitPvalue_LB_partial_temp)<-colnames(datExpr_hubgene_les_ctrl)[90:559]
  moduleTraitPvalue_LB_partial_temp<-data.frame(t(moduleTraitPvalue_LB_partial_temp))
  write.csv(moduleTraitPvalue_LB_partial_temp,file = "mpa_tur_hubgenes_partial_spearman_pvalue_sexonly.csv")
  
  
  
  
  
}
#LM22
if(T){
  #cortest_psy <- corr.test(CIBERSORTx_Job29_Results_LM22_PSO_LB[,c(12:285)], CIBERSORTx_Job29_Results_LM22_PSO_LB[,c(286:310)], method = "spearman")
  cortest_psy_adj <- corr.test(datTraits_PSO_LB_select_LM22[,c(3:14)], datTraits_PSO_LB_select_LM22[,c(36:57)], method = "spearman", adjust = "fdr")
  moduleTraitCor<-round(cortest_psy_adj$r,2)
  moduleTraitPvalue<-round(cortest_psy_adj$p,2) 
  moduleTraitPvalue_adjust<-round(cortest_psy_adj$p.adj,digits = 2)
  write.csv(moduleTraitCor,file = "fungi_ciber_lm22_association_pso.csv")
  write.csv(moduleTraitPvalue,file = "fungi_ciber_lm22_pso_pvalue.csv")
  write.csv(moduleTraitPvalue_adjust,file = "fungi_ciber_lm22_pso_padjust.csv")
  
  moduleTraitCor_res<-filter_all(data.frame(moduleTraitCor),any_vars(abs(.)>0.3))
  moduleTraitPvalue_adjust_res<-filter_all(data.frame(moduleTraitPvalue_adjust),any_vars(.<0.05))
  res_id<-rownames(moduleTraitCor_res)
  moduleTraitPvalue_adjust_res<-moduleTraitPvalue_adjust_res[res_id,]
  max(moduleTraitCor_res)
  min(moduleTraitCor_res)
  
  write.csv(moduleTraitCor_res,file = "fungi_pso_lm22_asso.csv")
  write.csv(moduleTraitPvalue_adjust_res,file = "fungi_pso_lm22_padjust.csv")
  
  
  #png(filename = "mpa_cibersort_lm22_correlation_les_ctrl.png",width = 6000,height = 6000,res = 300)
  col<-colorRampPalette(c("blue","white","red"))(10)
  corrplot(as.matrix(moduleTraitCor_res), method = "square",tl.cex = 0.4,col = col,tl.col = "black",
           p.mat=as.matrix(moduleTraitPvalue_adjust_res),#传入相关系数的P值
           sig.level = 0.05,#设置显著的水平
           insig = "blank",
  ) #0505
  dev.off()
}
# save(pheno_microarray_all,sampledf,
#      datTraits_PSO_LB_select_LM22,cortest_psy_adj,file="fungi_ad_LM22_asso.rdata")
