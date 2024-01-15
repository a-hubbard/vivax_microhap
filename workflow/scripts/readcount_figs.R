# Original Author: Elizabeth Hemming-Schroeder
# Revised By: Alfred Hubbard

setwd("~/pv_microhaplotype_2023-06")


ella <- read.table("ella_readcount.txt")
vera <- read.table("vera_readcount.txt")

library(ggpubr)
library(tidyverse)


teal="#006E82"
purple="#8214A0"
blue="#005AC8"
azure="#00A0FA"
pink="#FA78FA"
aqua="#14D2DC"
raspberry="#AA0A3C"
green="#0A9B4B"
vermillion="#FF825F"
yellow="#EAD644"
lgreen="#A0FA82"
banana="#FAE6BE"

############################################
####### 1st sequencing run ################
####### May 2022 ###########################
############################################

# get meta data
m22_metad <- read.csv("Samples_GTseek_UCI_Feb2022_Pv.csv") %>%
  mutate(Source = paste0(Source,"+", Extraction)) %>%
  select(-Extraction)

# counts per primer
may2022 <- read.table("readcount.txt") %>%
  group_by(V1) %>%
  summarise(Primer = V1, Sample = V2, Reads = V3, mean_reads = median(V3)) %>%
  ungroup() %>%
  select(-V1) %>%
  mutate(treatment=ifelse(grepl("preAmped",Sample)==TRUE,"PREAMP","NONE")) %>%
  separate(Sample,into=c("junk","Wellx"),sep="_Feb2022_") %>%
  mutate(WellL=substr(Wellx, 1,1)) %>%
  mutate(WellN=substr(Wellx,2,3)) %>%
  mutate(WellN=as.numeric(WellN)) %>%
  mutate(Well=paste0(WellL,WellN)) %>%
  select(-Wellx,-WellL,-WellN,-junk) %>%
  left_join(m22_metad,.) %>%
  mutate(Treatment = paste0(Enrichment,"+",treatment)) %>%
  select(-treatment) %>%
  filter(Treatment != "sWGA+PREAMP") %>%
  mutate(Enrichmentx = ifelse(Treatment == "sWGA+NONE","sWGA",
                            ifelse(Treatment=="none+PREAMP", "Targ. Pre-amp.",
                                   ifelse(Treatment =="none+NONE", "None",NA)))) 
  
rc_prim <- ggplot(data=may2022) +
  geom_boxplot() +
  facet_wrap(~factor(Enrichmentx, levels=c('None', 'Targ. Pre-amp.','sWGA'))) +
  scale_y_continuous(trans="log10") +
  aes(x=reorder(Primer,-mean_reads),y=Reads) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  theme_bw() +
  xlab("Locus") +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=0.5,color="black", size =6),
        axis.text.y = element_text(hjust=0,color="black", size=8),
        axis.title.x = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=8),
        legend.title = element_text(color="black", size=8),
        legend.text = element_text(face="italic",color="black", size=8),
        legend.spacing.x= unit(.1,"cm"),
        legend.spacing.y= unit(.2,"cm"),
        legend.position = "right",
        strip.background = element_rect(color="black", fill="gray95", linewidth=0.4, linetype="solid"),
        strip.text.x = element_text(size = 8)) 

rc_prim
ggsave("figures/readcount_primer.png", rc_prim, device="png", width=8, height=5, units = c("in"), dpi=600)
# total on target reads
m22ot <- may2022 %>%
  group_by(Well,Treatment) %>%
  summarise(OT_reads = sum(Reads)) %>%
  distinct()  %>%
  left_join(m22_metad,.)  %>%
  ungroup()
  
# total reads
m22tr <- read.table("may2022_totalreadcount.txt") %>%
  mutate(treatment=ifelse(grepl("preAmped",V1)==TRUE,"PREAMP","NONE")) %>%
  separate(V1,into=c("junk","Wellx"),sep="_Feb2022_") %>%
  mutate(WellL=substr(Wellx, 1,1)) %>%
  mutate(WellN=substr(Wellx,2,3)) %>%
  mutate(WellN=as.numeric(WellN)) %>%
  mutate(Well=paste0(WellL,WellN)) %>%
  select(-Wellx,-WellL,-WellN,-junk) %>%
  dplyr::rename(Total_reads=V2) %>%
  left_join(m22_metad,.) %>%
  mutate(Treatment = paste0(Enrichment,"+",treatment)) %>%
  select(-Enrichment,-treatment)  %>%
  filter(Treatment != "sWGA+PREAMP") %>%
  mutate(Enrichmentx = ifelse(Treatment == "sWGA+NONE","sWGA",
                            ifelse(Treatment=="none+PREAMP", "Targ. Pre-amp.",
                                   ifelse(Treatment =="none+NONE", "None",NA)))) 
str(m22tr)
str(m22ot)
head(m22tr)
tail(m22tr)
head(m22ot)
tail(m22ot)
#compare OT reads to total reads
m22tr %>% filter(Well=="A1")
OT_TR <-left_join(m22ot,m22tr) %>%
  mutate(Perc_OT = (OT_reads/Total_reads)*100) %>%
  #rename(Enrichment = Treatment) %>%
  mutate(Enrichment = factor(Enrichmentx,levels=c("None","Targ. Pre-amp.","sWGA"))) %>%
  mutate(source= as.factor(Source)) %>%
  mutate(ParasiteDensity= ifelse(Source=="Whole blood+Kit",(10^((CT_Pv - 41.663)/-3.289))/5, (10^((CT_Pv - 41.663)/-3.289)))) %>%
  select(-Well,-Species,-CT,-CT_Pv) %>%
  mutate(Source2 = ifelse(str_detect(Source,"DBS"),"DBS","Whole blood"))

library(ggpubr)
library(rstatix)
library(tidyverse)


#hist(comp_ot_tr$Perc_OT)

Reads_summary <- OT_TR %>%
  group_by(Source,Enrichment) %>%
  dplyr::summarise(meanPE=mean(Perc_OT),medPE=median(Perc_OT),
            maxPE=max(Perc_OT), minPE=min(Perc_OT),
            meanOT=mean(OT_reads),medOT=median(OT_reads),
            maxOT=max(OT_reads), minOT=min(OT_reads),n=n(),
            meanTR=mean(Total_reads),medTR=median(Total_reads),
            maxTR=max(Total_reads), minTR=min(Total_reads))
Reads_summary
write.csv(Reads_summary,"May2022_reads_summary.csv")

# pairwise
my_comparisons <- list(c("None","Targ. Pre-amp."),c("None","sWGA"))

OT_TR2 <- OT_TR %>%
  ungroup() %>%
  #filter(Enrichment != "Targ. Pre-amp.") %>%
  #filter(Enrichment != "sWGA") %>%
  mutate(SID_Source=paste0(SID,"_",Source)) %>%
  #select(-SID) %>%
  #mutate(SID=as.factor(SID_Source)) %>%
  mutate(SID_Enrichment = paste0(SID,"_",Enrichment))
 # group_by(SID)
head(OT_TR2)
#comp_ot_tr$Enrichment <- as.factor(comp_ot_tr$Enrichment)
#comp_ot_tr$Source <- as.factor(comp_ot_tr$Source)
stat.testx <- OT_TR2 %>%
  ungroup() %>%
  #filter(Source=="DBS+Chelex") %>%
  arrange(Enrichment) %>%
  #group_by(Source) %>%
  wilcox_test(Perc_OT~Enrichment, paired=TRUE) %>%
  adjust_pvalue(method="holm") %>%
  add_significance() %>%
  add_xy_position(x = "Enrichment", dodge = 0.8) %>%
  mutate(y.position = ifelse(statistic==32,y.position+9,
                             ifelse(xmax==3 &xmin==1,y.position+3.5,y.position-2)))


########################################
### MS Figure 1 ########################
########################################

rc_KW_X <- ggplot() +
 # aes(x=Enrichment ,y=Perc_OT,fill=Enrichment ) +
  #facet_wrap(~Source) +
  # scale_y_continuous(trans="log10") +
  geom_boxplot(data=OT_TR2,aes(x=Enrichment ,y=Perc_OT,fill=Enrichment)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw() +
  scale_fill_manual(values=c(vermillion,azure,green)) +
  #xlab("Locus") +
  ylab("OT reads (%)") + 
  scale_y_continuous(limits=c(0,105)) +
  stat_pvalue_manual(stat.testx, label = "{p.adj}{p.adj.signif}", 
                     tip.length = 0,hide.ns = FALSE,size=2.5) +
  #stat_compare_means(comparisons = my_comparisons,label="p.format")+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50)  +   # Add global p-value
  theme(axis.text.x = element_text(angle=60,vjust=0.5,hjust=0.5,color="black", size =8),
        axis.text.y = element_text(hjust=0,color="black", size=8),
        axis.title.x = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=8),
        legend.title = element_text(color="black", size=8),
        legend.text = element_text(color="black", size=8),
        legend.spacing.x= unit(.1,"cm"),
        legend.spacing.y= unit(.2,"cm"),
        legend.position = "none",
        strip.background = element_rect(color="black", fill="gray95", size=0.4, linetype="solid"),
        strip.text.x = element_text(size=8)) 
rc_KW_X
#ggsave("figures/readcount_Wilcox_Enrichment.png", rc_KW_E, device="png", width=6, height=4, units = c("in"), dpi=600)
#hist(comp_ot_tr$Perc_OT)


#comp_ot_tr$Enrichment <- as.factor(comp_ot_tr$Enrichment)
#comp_ot_tr$Source <- as.factor(comp_ot_tr$Source)
my_comparisons2 <- list(c("DBS+Chelex","DBS_Kit"),c("DBS+Kit","Whole blood+Kit"))
stat.test1 <- OT_TR2 %>%
  ungroup() %>%
  select(-SID_Source) %>%
  filter(Source !="Whole blood+Kit") %>%
  #group_by(Enrichment) %>%
  wilcox_test(Perc_OT~Source,paired=TRUE) %>%
  adjust_pvalue(method="holm") %>%
  add_significance() %>%
  add_xy_position(x = "Source", dodge = 0.8)

stat.test2 <- OT_TR2 %>%
  ungroup() %>%
  select(-SID_Source) %>%
  filter(Source !="DBS+Chelex") %>%
  #group_by(Enrichment) %>%
  wilcox_test(Perc_OT~Source) %>%
  adjust_pvalue(method="holm") %>%
  add_significance() %>%
  add_xy_position(x = "Source", dodge = 0.8) %>%
  mutate(xmin = xmin +1) %>%
  mutate(xmax = xmax+1)

OT_TR3 <- OT_TR2 %>% filter(Source !="DBS+Chelex")
stat.test3 <- OT_TR3 %>%
  ungroup() %>%
  #group_by(Enrichment) %>%
  wilcox_test(ParasiteDensity~Source) %>%
  adjust_pvalue(method="holm") %>%
  add_significance() %>%
  add_xy_position(x = "Source", dodge = 0.8) 


stat.test <- rbind(stat.test1,stat.test2) %>%
  mutate(y.position= ifelse(statistic==56,y.position+12.5,y.position+4))


rc_KW_E <- ggplot() +
  #facet_wrap(~Enrichment) +
  geom_boxplot(data=OT_TR,aes(x=Source ,y=Perc_OT,fill=Source)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw() +
  scale_fill_manual(values=c(vermillion,azure,green)) +
  scale_y_continuous(limits=c(0,105)) +
  #xlab("Locus") +
  ylab("OT reads (%)") + 
  stat_pvalue_manual(stat.test, label = "{p.adj}{p.adj.signif}", 
                     tip.length = 0,hide.ns = FALSE,size=2.5) +
  theme(axis.text.x = element_text(angle=60,vjust=0.5,hjust=0.5,color="black", size =8),
        axis.text.y = element_text(hjust=0,color="black", size=8),
        axis.title.x = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=8),
        legend.title = element_text(color="black", size=8),
        legend.text = element_text(color="black", size=8),
        legend.spacing.x= unit(.1,"cm"),
        legend.spacing.y= unit(.2,"cm"),
        legend.position = "none",
        strip.background = element_rect(color="black", fill="gray95", size=0.4, linetype="solid"),
        strip.text.x = element_text(size = 8)) 
rc_KW_E
#ggsave("figures/readcount_Wilcox_Source.png", rc_KW_E, device="png", width=6, height=4, units = c("in"), dpi=600)
#hist(comp_ot_tr$Perc_OT)

rc_KW_PD <- ggplot() +
  #facet_wrap(~Enrichment) +
  geom_boxplot(data=OT_TR3,aes(x=Source ,y=ParasiteDensity,fill=Source)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw() +
  scale_fill_manual(values=c(banana,green)) +
  #xlab("Locus") +
  ylab("ParasiteDensity") + 
  scale_y_continuous(limits=c(0,19000)) +
  scale_x_discrete(labels=c("DBS","Whole blood")) +
  stat_pvalue_manual(stat.test3, label = "{p.adj}{p.adj.signif}", 
                     tip.length = 0,hide.ns = FALSE,size=2.5) +
  theme(axis.text.x = element_text(angle=60,vjust=0.5,hjust=0.5,color="black", size =8),
        axis.text.y = element_text(hjust=0,color="black", size=8),
        axis.title.x = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=8),
        legend.title = element_text(color="black", size=8),
        legend.text = element_text(color="black", size=8),
        legend.spacing.x= unit(.1,"cm"),
        legend.spacing.y= unit(.2,"cm"),
        legend.position = "none",
        strip.background = element_rect(color="black", fill="gray95", size=0.4, linetype="solid"),
        strip.text.x = element_text(size = 8)) 
rc_KW_PD

ggplot() +
  geom_histogram(data=OT_TR,aes(x=Perc_OT)) +
  facet_wrap(~Enrichment)

my.formula=y~x
library(ggpmisc)
ot_line <- ggplot(data=OT_TR,aes(x=ParasiteDensity,y=Perc_OT)) +
  #facet_grid(cols=vars(Source),rows=vars(Enrichment),scales="free_x") +
  facet_grid(cols=vars(Enrichment),rows=vars(Source2),scales="free_x") +
  geom_point(color="gray40") +
  geom_smooth(method="lm",se=TRUE,color="gray40",fill="gray60") +
  stat_poly_eq(formula = my.formula, rr.digits=2, label.x=0.51, label.y=0.96, hjust=0.5,coef.digits = 2,
               aes(label = paste(..eq.label.., sep = "~~~")), size=2.7, 
               parse = TRUE, color="black", rsquared.conf.level = 0.99) +    
  stat_poly_eq(formula = my.formula, rr.digits=2, label.x=0.51, label.y=0.84, hjust=0.5,coef.digits = 2,
               aes(label = paste( ..rr.label.., sep = "~~~")), size=2.7, 
               parse = TRUE, color="black", rsquared.conf.level = 0.99) +    
  theme_bw() +
  scale_y_continuous(limits=c(-15,130)) +
  #scale_color_manual(values=c(vermillion,azure,green)) +
  xlab("Parasite Density") +
  ylab("OT reads (%)") + 
  theme(axis.text.x = element_text(angle=60,vjust=0.5,hjust=0.5,color="black", size =8),
        axis.text.y = element_text(hjust=0,color="black", size=8),
        axis.title.x = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=8),
        legend.title = element_text(color="black", size=8),
        legend.text = element_text(color="black", size=8),
        legend.spacing.x= unit(.1,"cm"),
        legend.spacing.y= unit(.2,"cm"),
        legend.position = "none",
        strip.background = element_rect(color="black", fill="gray95", size=0.4, linetype="solid"),
        strip.text.x = element_text(size = 8)) 
ot_line
#ggsave("figures/readcount_line.png", ot_line, device="png", width=5.5, height=3, units = c("in"), dpi=600)


ggarrange(rc_KW_E,rc_KW_X,ot_line)


MS_READ_PLOT1 <- ggarrange(ggarrange( rc_KW_X, rc_KW_E,rc_KW_PD, ncol = 3, labels = c("A", "B","C"),vjust=1.5), # Second row with box and dot plots
          ot_line,
          nrow = 2, 
          labels = c("","D"),
          heights = c(1, 1.1)# Labels of the scatter plot
) 


MS_READ_PLOT1
ggsave("MS_READ_PLOT1.png",MS_READ_PLOT1,height=4.8, width=4.8,units="in",dpi=600)

############################################
####### 2nd sequencing run ################
####### ella ###########################
############################################

# get metadata
ella_samples <- read.csv("Pv_samples_Idaho.csv") %>%
  mutate(Well = substr(SeqPlate_Well,4,6)) %>%
  select(-SeqPlate_Well) %>%
  mutate(ParasiteDensity= as.numeric(ParasiteDensity))

# get read counts by primer
ella <- read.table("ella_readcount.txt") %>%
  mutate(treatment=ifelse(startsWith(V2, "PV_PREAMP"),"PREAMP","sWGA")) %>%
  separate(V2, into=c("j1","j2","Well","j3","j4"),sep="_") %>%
  select(-j1,-j2,-j3,-j4) %>%
  group_by(V1) %>%
  summarise(Primer = V1, Well=Well, treatment = treatment,Reads = V3, mean_reads = 2*median(V3)+mean(V3)) %>%
  ungroup() %>%
  select(-V1) %>%
  left_join(ella_samples,.)

# get on target read counts
ella2 <- ella %>%
  group_by(Well,treatment) %>%
  summarise(OT_reads = sum(Reads)) %>%
  ungroup() %>%
  distinct() %>%
  left_join(ella_samples,.)

# get total read counts
ella_tr <- read.table("ella_totalreadcount.txt") %>%
  mutate(treatment=ifelse(startsWith(V1, "PV_PREAMP"),"PREAMP","sWGA")) %>%
  separate(V1, into=c("j1","j2","Well","j3","j4"),sep="_") %>%
  select(-j1,-j2,-j3,-j4) %>%
  rename(Total_reads = V2)

# compare OT to total reads
ella_comp <- left_join(ella2,ella_tr) %>%
  mutate(Perc_OT = (OT_reads/Total_reads)*100) %>%
  arrange(desc(Perc_OT))

ggplot(data=ella_comp) +
  aes(x=CT_Pv,y=OT_reads, color=treatment,fill=treatment) +
  geom_point() +
  geom_smooth(method="lm") +
  theme(axis.text.x = element_text(angle = 90)) 
 # scale_y_continuous(trans="log10") 
  #scale_x_continuous(trans="log10")


ggplot(data=ella) +
  geom_boxplot() +
  facet_wrap(~treatment) +
  scale_y_continuous(trans="log10") +
  aes(x=reorder(Primer,-mean_reads),y=Reads) +
  theme(axis.text.x = element_text(angle = 90))

############################################
####### 3rd sequencing run ################
####### vera ###########################
############################################

# metadata
meta_data <- read.csv("Pv_samples_Idaho.csv") %>%
  mutate(Well = substr(SeqPlate_Well,4,6)) %>%
  select(-SeqPlate_Well) %>%
  mutate(ParasiteDensity= as.numeric(ParasiteDensity)) %>%
  mutate(Mixed=ifelse(!is.na(CT_Pf18s) | !is.na(CT_PfvarATS),"Yes","No"))

# number of reads by primer by sample
vera <- read.table("vera_readcount.txt") %>%
  dplyr::filter(V3>0) %>%
  mutate(treatment = ifelse(startsWith(V2,"PV_SWGA_OLD"),"OLD","NEW")) %>%
  group_by(V1) %>% 
  summarise(Primer = V1, Sample = V2, treatment = treatment, Reads = V3, mean_reads = (median(V3)+(mean(V3)*0.01))) %>%
  ungroup() %>%
  filter(grepl("R1|R2",Primer)==FALSE)

verax <- vera %>% 
  group_by(Primer) %>%
  summarise(Reads = quantile(Reads, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75)) %>%
  filter(q==0.5) %>%
  arrange(Reads)

  
write.csv(verax,"vera_primer_check.csv")
verax
ggplot(data=vera) +
  geom_boxplot() +
  facet_wrap(~treatment,scales="free_x") +
  scale_y_continuous(trans="log10") +
  aes(x=reorder(Primer,-mean_reads),y=Reads) +
  theme(axis.text.x = element_text(angle = 90))

ggplot(data=vera,aes(x = fct_infreq(Primer))) +
  geom_bar() +
  labs(x = "Primer") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~treatment,nrow=2)

# number of on target reads per sample
n_reads <- vera %>%
  separate(Sample,into=c("j1","j2","SeqID","j3","j4")) %>%
  select(-starts_with("j")) %>%
  group_by(SeqID,treatment) %>%
  summarise(SeqID = SeqID,treatment=treatment,OT_reads = sum(Reads)) %>%
  ungroup() %>%
  distinct()



# total number of reads per sample
vera_tr <- read.table("vera_totalreadcount.txt")  %>%
  mutate(treatment = ifelse(startsWith(V1,"PV_SWGA_OLD"),"OLD","NEW")) %>%
  filter(V1 != "Undetermined_S0_L001") %>%
  separate(V1,into=c("j1","j2","SeqID","j3","j4")) %>%
  select(-starts_with("j")) %>%
  rename(Total_reads = V2)

# compare on target reads to total reads
vera_comp <- left_join(n_reads,vera_tr) %>%
  mutate(Perc_OT = (OT_reads/Total_reads)*100) 

# connects barcode to our naming system
names_barcodes <- read.csv("vera_library_2023-03-18_HemmingSchroeder_UCDavis.csv") %>%
  mutate(treatment=ifelse(grepl("new",Library.Name)==TRUE,"NEW",ifelse(grepl("old",Library.Name)==TRUE,"OLD",NA))) %>%
  select(X,treatment,i7.index.Sequence,i5.index.Sequence) %>%
  rename(Well=X)

# connects barcode to sequence id
seqid_barcodes <- read.csv("vera_lane_summary.csv") %>%
  separate(Barcode.sequence, into=c("i7.index.Sequence","i5.index.Sequence")) %>%
  separate(Sample,into=c("j1","j2","SeqID")) %>%
  select(-starts_with("j"))


# connect meta data, to barcode, to sequence name, to read count
vera_combined <- left_join(meta_data,names_barcodes) %>%
  left_join(.,seqid_barcodes) %>%
  left_join(.,vera_comp) %>%
  filter(treatment=="OLD")

ggplot(data=vera_combined) +
  aes(x=ParasiteDensity,y=Perc_OT) +
  #facet_wrap(~treatment) +
  geom_point() +
  geom_smooth(method="lm") +
  #scale_y_continuous(trans="log10") +
  scale_x_continuous(trans="log10")
my.formula=y~x
ot_line_vera <- ggplot(data=vera_combined,aes(x=ParasiteDensity,y=Perc_OT)) +
  #facet_grid(cols=vars(Source),rows=vars(Enrichment),scales="free_x") +
 # facet_grid(cols=vars(Enrichment),rows=vars(Source2),scales="free_x") +
  geom_point(color="gray70") +
  geom_smooth(method="lm",se=TRUE,color="gray40",fill="gray80") +
  scale_x_continuous(trans="log10") +
  stat_poly_eq(formula = my.formula, rr.digits=2, label.x=0.05, label.y=0.96, hjust=0,coef.digits = 2,
               aes(label = paste(..eq.label.., sep = "~~~")), size=2.7, 
               parse = TRUE, color="black", rsquared.conf.level = 0.99) +    
  stat_poly_eq(formula = my.formula, rr.digits=2, label.x=0.05, label.y=0.84, hjust=0,coef.digits = 2,
               aes(label = paste( ..rr.label.., sep = "~~~")), size=2.7, 
               parse = TRUE, color="black", rsquared.conf.level = 0.99) +    
  theme_bw() +
  scale_y_continuous(limits=c(0,100)) +
  #scale_color_manual(values=c(vermillion,azure,green)) +
  xlab("Parasite Density") +
  ylab("OT reads (%)") + 
  theme(axis.text.x = element_text(angle=0,vjust=0.5,hjust=0.5,color="black", size =8),
        axis.text.y = element_text(hjust=0,color="black", size=8),
        axis.title.x = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=8),
        legend.title = element_text(color="black", size=8),
        legend.text = element_text(color="black", size=8),
        legend.spacing.x= unit(.1,"cm"),
        legend.spacing.y= unit(.2,"cm"),
        legend.position = "none",
        strip.background = element_rect(color="black", fill="gray95", size=0.4, linetype="solid"),
        strip.text.x = element_text(size = 8)) 
ot_line_vera


ggplot(data=vera_combined) +
  aes(x=ParasiteDensity,y=Perc_OT) +
  facet_wrap(~treatment) +
  geom_point() +
  geom_smooth(method="lm") +
 # scale_y_continuous(trans="log10") +
  scale_x_continuous(trans="log10")

ggplot(data=vera_combined) +
  aes(x=CT_Pv,y=Perc_OT,color=Mixed,fill=Mixed) +
  facet_wrap(~treatment) +
  geom_point() +
  geom_smooth(method="lm") 
  #scale_y_continuous(trans="log10") 



