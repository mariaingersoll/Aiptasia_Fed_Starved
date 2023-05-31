nfkb_apo = bp_apo %>%
  filter(str_detect(term, "GO:1901222|GO:0043123|GO:0051092|GO:0043122|GO:0042345|GO:0007249")) %>%
  left_join(rlog_apo_10) %>%
  left_join(gene) %>%
  mutate(gene_symbol = make.names(gene_symbol, unique = TRUE)) %>%
  column_to_rownames(var = "gene_symbol") %>%
  dplyr::select(-name, -term, -lev, -value, -description, -gene_ID, -accession) %>%
  drop_na() %>%
  dplyr::select(sort(current_vars()))
head(nfkb_apo)
nrow(nfkb_apo)



test_df2 = bp_apo %>%
  filter(str_detect(term, "GO:1901222|GO:0043123|GO:0051092|GO:0043122|GO:0042345|GO:0007249")) %>% # select desired GO terms (nfkb)
  left_join(rlog_apo_10) %>% 
  drop_na() %>%
  distinct(gene_ID)#%>% # left_join the sig DEGs (contains rlog values)
  # left_join(gene) %>%
  # mutate(gene_symbol = make.names(gene_symbol, unique = TRUE)) %>%
  # column_to_rownames(var = "gene_symbol") %>%
  # dplyr::select(-name, -term, -lev, -value, -description, -gene_ID, -accession) %>%
  # drop_na() %>%
  # dplyr::select(sort(current_vars()))
  
  
View(test_df2)
nrow(test_df2)
#217

test_df3 = bp_apo %>%
  filter(str_detect(term, "GO:1901222|GO:0043123|GO:0051092|GO:0043122|GO:0042345|GO:0007249")) %>% # select desired GO terms (nfkb)
  left_join(rlog_apo_10) %>% 
  drop_na()
nrow(test_df3)
view(test_df3)

test4 = bp_apo %>%
  filter(str_detect(term, "GO:1901222|GO:0043123|GO:0051092|GO:0043122|GO:0042345|GO:0007249")) %>% # select desired GO terms (nfkb)
  left_join(rlog_apo_10) %>% 
  drop_na() %>%
  distinct(gene_ID, .keep_all = TRUE)
view(test4)

nrow(bp_apo %>%
       filter(str_detect(term, "GO:1901222|GO:0043123|GO:0051092|GO:0043122|GO:0042345|GO:0007249")))
#217
nrow(rlog_apo_10)


View(rlog_apo_10)



head(nfkb_apo)
nrow(nfkb_apo)

sym_dens <- read.delim("sym_dens.txt")
SF <- sym_dens[c(1,5,9,13,17,21),]
SF
SF2 <- sym_dens[c(2,6,10,14,18,22),]
SS <- sym_dens[c(3,7,11,15,19,23),]
SS
SS2 <- sym_dens[c(4,8,12,16,20,24),]

fed <- sym_dens[c(1,2,5,6,9,10,13,14,17,18,21,22),]
fed %>% ggplot(aes(Week, Sym_Density)) + 
  geom_boxplot(fill=c("lightsteelblue3", "gray50"),alpha=0.3, color=c("gray40","gray40")) + 
  geom_point() +
  geom_line(aes(group=Pair), linetype="longdash") +
  theme(legend.position = "none")+
  theme_bw() +
  coord_cartesian(ylim = c(30, 100))

purpsym <- c("violet", "purple4")
title="Starved Clones"

starved <- sym_dens[c(3,4,7,8,11,12,15,16,19,20,23,24),]
starved
starved %>% ggplot(aes(Week, Sym_Density)) + 
  geom_boxplot(fill=c("plum2", "purple1"),alpha=0.3, color=c("purple3","purple3")) + 
  geom_point() +
  geom_line(aes(group=Pair), linetype="longdash") +
  theme(legend.position = "none")+
  theme_bw() +
  coord_cartesian(ylim = c(30, 100))

oc_seq2iso <- read.delim("O_arbuscula_sequenceID_to_isogroup.tab", header=FALSE)
head(oc_seq2iso)
b <- c("seq", "genes")
colnames(oc_seq2iso)=paste(b)
head(oc_seq2iso)
oc_kog_2<-read.delim("Oculina_gene2kogClass_final.tab", header=TRUE)
head(oc_kog_2)
m <- c("seq", "term")
colnames(oc_kog_2)=paste(m)
head(oc_kog_2)

oc_kog_2 <- oc_kog_2 %>%
  left_join(oc_seq2iso)
head(oc_kog_2)
oc_kog_2 <- oc_kog_2[,2:3]
head(oc_kog_2)
oc_kog_2 <- oc_kog_2 %>% distinct()
head(oc_kog_2)
n <- c("term", "seq")
colnames(oc_kog_2)=paste(n)
head(oc_kog_2)
oc_kog_2 <- oc_kog_2[,c(2,1)]
head(oc_kog_2)
write.table(oc_kog_2, file="ocu_kog_mvi_final.txt", quote=F, sep="\t")

aip_kog_2<-read.delim("aiptasia_gene2kogClass_final.tab", header=TRUE)
head(aip_kog_2)
colnames(aip_kog_2)=paste(m)
head(aip_kog_2)
write.table(aip_kog_2, file="aip_kog_mvi_final.txt", quote=F, sep="\t")

nv_kog_2 <- read.delim("nematostella_gene2kogClass_final.tab", header=TRUE)
head(nv_kog_2)
colnames(nv_kog_2)=paste(m)
head(nv_kog_2)
write.table(nv_kog_2, file="nv_kog_mvi_final.txt", quote=F, sep="\t")


