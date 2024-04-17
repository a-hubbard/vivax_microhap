# Original Author: Elizabeth Hemming-Schroeder
# Revised By: Alfred Hubbard

library(ggpmisc)
library(ggpubr)
library(optparse)
library(rstatix)
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(make_option("--out_base", help = "Basename for output figure"))
arg <- parse_args(OptionParser(option_list = opts))

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
m22_metad <- read.csv("resources/Samples_GTseek_UCI_Feb2022_Pv.csv") %>%
  mutate(Source = paste0(Source,"+", Extraction)) %>%
  select(-Extraction)

# counts per primer
may2022 <- read.table("resources/readcount.txt") %>%
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

# total on target reads
m22ot <- may2022 %>%
  group_by(Well,Treatment) %>%
  summarise(OT_reads = sum(Reads)) %>%
  distinct()  %>%
  left_join(m22_metad,.)  %>%
  ungroup()
  
# total reads
m22tr <- read.table("resources/may2022_totalreadcount.txt") %>%
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

my.formula=y~x
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

MS_READ_PLOT1 <- ggarrange(ggarrange( rc_KW_X, rc_KW_E,rc_KW_PD, ncol = 3, labels = c("A", "B","C"),vjust=1.5), # Second row with box and dot plots
          ot_line,
          nrow = 2, 
          labels = c("","D"),
          heights = c(1, 1.1)# Labels of the scatter plot
) 


w <- 4.8
h <- 4.8
ggsave(
  str_c(arg$out_base, ".pdf"), 
  plot = MS_READ_PLOT1, 
  width = w, 
  height = h, 
  units = "in"
)
ggsave(
  str_c(arg$out_base, ".png"), 
  plot = MS_READ_PLOT1, 
  width = w, 
  height = h, 
  units = "in"
)
