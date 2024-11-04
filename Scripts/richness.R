rm(list = ls())  
options(stringsAsFactors = F)

library(vegan)
library(BiodiversityR)

#ps
daf<- data.frame(ps@otu_table)
sam<-data.frame(ps@sam_data)
# Shannon
daf<-data.frame(t(daf))
daf2<-round(daf*1000)

shannon <- diversity(daf, index = "shannon")

# Chao1
chao1 <- vegan::estimateR(daf2)

# Pielou
richness <- specnumber(daf)
pielou <- shannon / log(richness)

# dataframe
a<-ps1[["data"]]
df_aplha <- data.frame(
  Shannon = shannon,
  Chao1 = chao1[1,],
  Pielou = pielou
)
df_aplha$scilife_id<-rownames(df_aplha)
df_aplha<-merge(df_aplha,sam,by="scilife_id")


# boxplot

df_aplha$group <- factor(df_aplha$group, levels = c("a","b","CTRL"))
p_alpha<-ggboxplot(df_aplha, x = 'group', y = 'Shannon', color = 'group', palette = 'jco', add = 'jitter') +
  geom_pwc(
    aes(group = group), 
    tip.length = 0,
    method = "wilcox_test", 
    p.adjust.method="BH",
    label = "p.adj.format",
    hide.ns = FALSE,
    bracket.nudge.y = -0.08
  ) + theme(text = element_text(size = 15))
p_alpha
ggsave(p_alpha,file="Shannon.pdf")

p2<-ggboxplot(df_aplha, x = 'group', y = 'Chao1', color = 'group', palette = 'jco', add = 'jitter') +
  geom_pwc(
    aes(group = group), 
    tip.length = 0,
    method = "wilcox_test", 
    p.adjust.method="BH",
    label = "p.adj.format",
    hide.ns = FALSE,
    bracket.nudge.y = -0.08
  ) + theme(text = element_text(size = 15))
p2
ggsave(p2,file="p2.pdf")
graphics.off()

p3<-ggboxplot(df_aplha, x = 'group', y = 'Pielou', color = 'group', palette = 'jco', add = 'jitter') +
  geom_pwc(
    aes(group = group), 
    tip.length = 0,
    method = "wilcox_test", 
    p.adjust.method="BH",
    label = "p.adj.format",
    hide.ns = FALSE,
    bracket.nudge.y = -0.08
  ) + theme(text = element_text(size = 15))
p3
ggsave(p3,file="Pielou.pdf")
graphics.off()

