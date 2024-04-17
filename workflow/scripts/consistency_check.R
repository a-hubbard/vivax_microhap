# Original Author: Elizabeth Hemming-Schroeder
# Revised By: Alfred Hubbard

library(tidyverse)
library(ggpubr)
library(rstatix)
library(optparse)

# Parse arguments ------------------------------------------------------
opts <- list(make_option("--out_base", help = "Basename for output figure"))
arg <- parse_args(OptionParser(option_list = opts))

# Create mode() function to calculate mode
mode <- function(x, na.rm = FALSE) {
  
  if(na.rm){ #if na.rm is TRUE, remove NA values from input x
    x = x[!is.na(x)]
  }
  
  val <- unique(x)
  tab.val <- tabulate(match(x,val))
  max.val <- which(tab.val==max(tab.val))
  out.val <- val[max.val]
  out.col <- paste(out.val, sep="", collapse=":") 
  return(out.col)
}


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
new_df <- read.csv("resources/Microhap_PV_combined-data_2023-06-21") %>%
  select(-X)

# reps are indicated by
# 1) extraction method
# 2) Dilution
# 3) Primer_set
# 4) Treatment

# split by SID then identify 
# may2022 first
splitx <- new_df %>%
  mutate(mh_allele_orig = mh_allele) %>%
  mutate(mh_allele = LETTERS[mh_allele]) %>%
  filter(!startsWith("PvSNP",Primer)) %>%
  filter(!is.na(mh_allele)) %>%
  filter(SeqRun=="may2022") %>%
  filter(Treatment!="sWGA+PREAMP") %>%
  select(-Dilution,-Primer_set) %>%
  group_by(SID, Primer,Source) %>%
  select(Source,Treatment,ParasiteDensity,
         SeqRun,
         Total_reads_sample,max_MOI,Reads,MOI,Prop_rank,Prop_value,
         mh_allele_snps,mh_allele,mh_allele_orig) %>%
  mutate(rep = paste0(Treatment,"+",Source)) %>%
  select(-Treatment,-Source) %>%
  group_split()

#splitx[[1300]]
#splitx[[1331]]
#i=1300
#i=1331
#j=2
new_listx <- c()
new_list <- c()
for (i in 1:length(splitx)) {
  
  # check if there are replicates
  check_reps <- length(unique(splitx[[i]]$rep))
  
  this_df <- splitx[[i]] %>%
    mutate(N_reps = check_reps)
  
  this_df_split <- this_df %>%
    group_by(rep) %>%
    group_split
  mhaps_list <- c()
  for (j in 1:length(this_df_split)) {
    this_df2 <- this_df_split[[j]]
    
    this_df2_allele_list <- this_df2 %>%
      group_by(mh_allele) %>%
      group_split
    
    # ignore length polymorphic variants
    # combine them into one haplotype if have same snps
    for (z in 1:length(this_df2_allele_list)) {
      if (nrow(this_df2_allele_list[[z]])>1) {
        new_prop_value <- sum(this_df2_allele_list[[z]]$Prop_value)
        new_df2 <- this_df2_allele_list[[z]][1,]
        new_df2$Prop_value <- new_prop_value
        this_df2_allele_list[[z]] <- NA
        this_df2_allele_list[[z]] <- new_df2
      }
    }
    this_df2 <- NA
    this_df2 <- do.call(rbind, this_df2_allele_list)
    
    mhaps <- toString(sort(this_df2$mh_allele))
    mhaps_list[[j]] <- mhaps
    this_df_split[[j]] <- NA
    this_df_split[[j]] <- this_df2
  }
  this_df <- NA
  this_df <- do.call(rbind,this_df_split)
  
  
  mhaps_df <- do.call(rbind,mhaps_list)
  mhap_mode <- mode(mhaps_df)
  
  mhaps_df_long <- as.data.frame(mhaps_df) %>%
    mutate(haps = gsub(" ","",V1)) %>%
    mutate(nhaps = ifelse(str_count(haps,",")>=1,str_count(haps,",")+1,1))
  
  max_mhaps <- max(mhaps_df_long$nhaps)
  
  filter_mhaps_long <- mhaps_df_long %>%
    filter(nhaps == max_mhaps)
  
  mhap_long_mode <- mode(filter_mhaps_long$haps)
    
  
  moi_mode_df <- this_df %>%
    select(rep,MOI) %>%
    distinct()

  
  moi_mode <- mode(moi_mode_df$MOI)
  moi_max <- max(moi_mode_df$MOI)
  #moi_2ndmax <- sort(moi_mode_df$MOI)[2]
  
  
  this_df_out <- this_df %>%
    mutate(mhaps_mode=mhap_mode[1]) %>%
    mutate(mhaps_long_mode=mhap_long_mode[1]) %>%
    mutate(mh_allele=as.character(mh_allele)) %>%
    mutate(mhap_max_ambig = str_count(mhaps_long_mode,":")+1) %>%
    mutate(mhap_mode_ambig = str_count(mhaps_mode,":")+1) %>%
    mutate(mhap_mode_check = ifelse(str_count(mhaps_mode,mh_allele)==mhap_mode_ambig,"Pass",
                                    ifelse(str_count(mhaps_mode,mh_allele)>=1,"Fail:ambiguous",
                                           ifelse(str_count(mhaps_mode)==0,"Fail:unambiguous",NA)))) %>%
    mutate(mhap_max_check = ifelse(str_count(mhaps_long_mode,mh_allele)==mhap_max_ambig,"Pass",
                                    ifelse(str_count(mhaps_long_mode,mh_allele)>=1,"Fail:ambiguous",
                                           ifelse(str_count(mhaps_long_mode)==0,"Fail:unambiguous",NA)))) %>%
    mutate(MOI_mode = moi_mode) %>%
    mutate(MOI_max = moi_max) %>%
    #mutate(MOI_2ndmax = moi_2ndmax) %>%
    mutate(MOI_mode_check = ifelse(MOI==MOI_mode,"Pass","Fail")) %>%
    mutate(MOI_max_check  = ifelse(MOI==MOI_max,"Pass","Fail"))  %>%
    mutate(rep= gsub("\\+",":",rep)) %>%
    tidyr::separate(col=rep,into=c("Enrichment","Source","Extraction"),sep=":")
    
  Minor_allele <- this_df_out %>%
    filter(Prop_rank>0) %>%
    filter(Enrichment=="None") %>%
    select(mh_allele)
  
  this_df_out_split <- this_df_out %>% 
    mutate(is_minor_allele = ifelse(is.na(Minor_allele[1,1]),NA,
                                           ifelse(mh_allele==Minor_allele[1,1],TRUE,FALSE))) %>%
    group_by(mh_allele) %>%
    group_split()
  
  lsplit_list <- c()
  for (k in 1:length(this_df_out_split)) {
    dfx_tmp <- this_df_out_split[[k]] %>%
      filter(MOI_max_check == "Pass") 
    n_treats <- nrow(dfx_tmp)
    dfx_split <- this_df_out_split[[k]] %>%
      mutate(n_treatments = n_treats) 
    if(nrow(dfx_split)<=1) {
      dfx_split_out <- dfx_split %>%
        mutate(Prop_sd = NA) %>%
        mutate(Prop_ref = NA) %>%
        mutate(Prop_diff = NA) }
    if (nrow(dfx_split)>1) {
    sd_prop <- sd(dfx_split$Prop_value)
    dfx_split_none <- dfx_split %>%
      filter(Enrichment=="None") %>%
      filter(Reads >10) %>%
      filter(MOI_max_check == "Pass")
    if (nrow(dfx_split_none) >= 1) {
    prop_ref <- mean(dfx_split_none$Prop_value,na.rm=TRUE)
    dfx_split_out <- dfx_split %>%
      mutate(Prop_sd = sd_prop) %>%
      mutate(Prop_ref = prop_ref) %>%
      mutate(Prop_diff = abs(Prop_ref - Prop_value)) }
    if (nrow(dfx_split_none) <1) {
      dfx_split_out <- dfx_split %>%
        mutate(Prop_sd = sd_prop) %>%
        mutate(Prop_ref = NA) %>%
        mutate(Prop_diff = NA) }
    }
      
    lsplit_list[[k]] <- dfx_split_out
  }
    OUT_DF <- do.call(rbind,lsplit_list)
    new_list[[i]] <- OUT_DF
}


match_df <- do.call(rbind,new_list) %>%
  mutate(Reads_bin = cut(Reads, breaks=c(0,10,100,1000,51000))) %>%
  mutate(ParasiteDensity_bin =cut(ParasiteDensity,breaks=c(0,100,1000,10000,20000))) 

match_df1 <- match_df %>% 
  filter(mhap_max_check=="Pass")  %>%
  filter(!is.na(Prop_ref)) %>%
  filter(n_treatments>=2) %>%
  filter(MOI_max >1)

match_df2 <- match_df1 %>%
  mutate(comp_id = paste0(SID,":",Extraction,":",Primer,":",mh_allele)) %>%
  filter(is_minor_allele == TRUE) %>%
  filter(Enrichment !="None") 


match_df3 <- as.data.frame(match_df2) %>%
  ungroup() %>%
  select(comp_id,Enrichment,Prop_diff) %>%
 # mutate(comp_id = as.factor(comp_id)) %>%
  mutate(Enrichment = as.factor(Enrichment)) %>%
  arrange(comp_id,Enrichment) 


####################################33
### split to do paired wilcox test all three groups at 

#########################################################
### microhaplotype consistency ##########################
########################################################

#MOI consistency
MOI_df <- match_df %>%
  mutate(comp_id = paste0(SID,":",Extraction,":",Primer)) %>%
  select(comp_id,Enrichment,ParasiteDensity,Reads,MOI,N_reps,Prop_diff,is_minor_allele, Prop_value, mh_allele) %>%
  filter(N_reps==3) %>%
  distinct() %>%
  #mutate(Enrichment = ifelse(Enrichment=="PREAMP","Targ. Pre-amp.", Enrichment)) %>%
  #mutate(Enrichment = factor(Enrichment,levels=c("None","Targ. Pre-amp.","sWGA"))) %>%
  arrange(comp_id,Enrichment) %>%
  dplyr::select(-Prop_value,-mh_allele,-is_minor_allele,-Prop_diff) %>%
  distinct() %>%
  mutate(MOI = as.numeric(MOI)) %>%
  mutate(Enrichment = factor(Enrichment,levels=c("None","PREAMP","sWGA"))) %>%
  mutate(Enrichment = as.numeric(Enrichment))

st_pd1 <- MOI_df %>%
  ungroup() %>%
  wilcox_test(MOI~Enrichment,paired=TRUE) %>%
  adjust_pvalue(method="holm") %>%
  add_significance() %>%
  add_xy_position(x = "Enrichment", dodge = 0.8) %>%
  mutate(y.position = ifelse(statistic==0,y.position+.2,
                             ifelse(statistic==96,y.position+.55,
                                    ifelse(statistic==976,y.position+0.9,NA))))

MOI_df_orig <- MOI_df

MOI_df$Enrichment <- jitter(MOI_df$Enrichment,amount=0.25)
MOI_df$MOI <- jitter(MOI_df$MOI,amount=0.25)
moi_consist <- ggplot(data=MOI_df,aes(x=Enrichment,y=MOI,group=comp_id)) +
  geom_line(colour="gray",alpha=0.4) +
  geom_point(shape = 16,alpha=0.4,aes(colour=Enrichment)) +
  stat_pvalue_manual(st_pd1, label = "{p.adj}{p.adj.signif}", 
                    tip.length = 0,hide.ns = FALSE,size=2.5) +

  theme_bw() +
  scale_colour_stepsn(colours=c(vermillion,vermillion,azure,azure,azure,green,green)) +
  #scale_fill_steps2(low=vermillion,mid=azure,high=green,midpoint=2) +
  scale_y_continuous(breaks=c(1,2,3),minor_breaks = c(1,2,3),limits=c(0.5,4.1)) +
  scale_x_continuous(breaks =c(1,2,3),minor_breaks=c(1,2,3), labels= c("None","Targ. Pre-amp.","sWGA")) +
  ylab("Alleles in\ninfection") + 
#  xlab("Haplotype proportion of reads") +
  theme(axis.text.x = element_text(angle=0,vjust=0.5,hjust=0.5,color="black", size =8),
        axis.text.y = element_text(hjust=0,color="black", size=8),
        axis.title.x = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=8),
        legend.title = element_text(color="black", size=8),
        legend.text = element_text(color="black", size=8),
        legend.spacing.x= unit(.1,"cm"),
        legend.spacing.y= unit(.2,"cm"),
        legend.position = "none",
        legend.background = element_rect(fill = "white",color="black"),
        strip.background = element_rect(color="black", fill="gray95", size=0.4, linetype="solid"),
        strip.text.x = element_text(size = 8)) 

################################
# allele consist #####################
#####################################

MOI_dfz <- match_df %>%
  mutate(comp_id = paste0(SID,":",Extraction,":",Primer)) %>%
  select(comp_id,Enrichment,ParasiteDensity,Reads,MOI,n_treatments,Prop_diff,is_minor_allele, Prop_value, mh_allele) %>%
  filter(n_treatments==3) %>%
  distinct() %>%
  # filter(Dilution != 1) %>%
  filter(is_minor_allele == TRUE) %>%
  #filter(Enrichment!="None") %>%
  filter(MOI>=2) %>%
  mutate(G2_id = paste0(comp_id,":",mh_allele)) %>%
  mutate(Enrichment = ifelse(Enrichment=="PREAMP","Targ. Pre-amp.", Enrichment)) %>%
  mutate(Enrichment = factor(Enrichment,levels=c("Targ. Pre-amp.","sWGA"))) 
write.csv(MOI_dfz, "results/MOI_dfz_check.csv")

st_pd2 <- MOI_dfz %>%
  ungroup() %>%
  select(-comp_id) %>%
  pairwise_wilcox_test(Prop_diff~Enrichment) %>%
  #friedman_test(MOI~Dilution| group_id) %>%
  adjust_pvalue(method="holm") %>%
  add_significance() %>%
  add_xy_position(x = "Enrichment", dodge = 1,group="Enrichment") 
st_pd2$y.position <- c(0.5)

minor_consist <- ggplot(data=MOI_dfz) +
  geom_jitter(aes(x=as.factor(Enrichment),y=Prop_diff, color = Enrichment),width=0.2,height=0.01) +
  #geom_point(aes(color=G2_id))
  #geom_line(aes(x=Enrichment,y=Prop_diff,group=G2_id, color=G2_id),alpha=0.5) +
  #geom_boxplot(width=0.7,aes(x=as.factor(Enrichment),y=Prop_diff,fill=Enrichment)) +
   stat_pvalue_manual(st_pd2, label = "{p.adj}{p.adj.signif}", 
                 tip.length = 0,hide.ns = FALSE,size=2.5,xmin="group1",xmax="group2",x=NULL,y.position="y.position") +
  #scale_y_continuous(trans="log10") +
  #stat_summary(aes(x=PD_log10_bin,y=MOI),fun.y=mean, geom="point", shape=17, size=5, color="black", fill="black") +
  #scale_x_continuous(trans="log10")+
  #stat_pvalue_manual(st_pd1, label = "{p.adj}{p.adj.signif}", 
  #                   tip.length = 0,hide.ns = FALSE,size=2.5) +
  theme_bw() +
  scale_fill_manual(values=c(azure,green)) +
  #scale_fill_steps2(low=vermillion,mid=azure,high=green,midpoint=2) +
   scale_y_continuous(limits=c(0,0.55)) +
  #scale_x_continuous(breaks =c(1,2,3),minor_breaks=c(1,2,3), labels= c("None","Targ. Pre-amp.","sWGA")) +
  ylab("Minor allele\nprevalence\ndifference") + 
  xlab("Enrichment") +
  theme(axis.text.x = element_text(angle=0,vjust=0.5,hjust=0.5,color="black", size =8),
        axis.text.y = element_text(hjust=0,color="black", size=8),
        axis.title.x = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=8),
        legend.title = element_text(color="black", size=8),
        legend.text = element_text(color="black", size=8),
        legend.spacing.x= unit(.1,"cm"),
        legend.spacing.y= unit(.2,"cm"),
        legend.position = "none",
        legend.background = element_rect(fill = "white",color="black"),
        strip.background = element_rect(color="black", fill="gray95", size=0.4, linetype="solid"),
        strip.text.x = element_text(size = 8)) 


# moi_none <- MOI_df %>%
#   filter(Enrichment=="None") %>%
#   select(MOI)
# 
# moi_tpa <- MOI_df %>%
#   filter(Enrichment=="PREAMP") %>%
#   select(MOI)
# 
# moi_swga <- MOI_df %>%
#   filter(Enrichment=="sWGA") %>%
#   select(MOI)
# d_moi <- data.frame(None=moi_none[,1],Preamp=moi_tpa[,1], swga=moi_swga[,1]) %>%
#   rename(None = MOI,`Targ. Pre-amp` = MOI.1, sWGA = MOI.2)
# 
# 
# ggpaired(d_moi,cond1 = "Targ. Pre-amp", cond2 = "sWGA", cond3
#          fill = "condition", line.color = "gray", line.size = 0.4,
#          xlab="Enrichment",ylab="MOI",
#          palette = c(azure,green,vermillion),ggtheme = theme_bw(),legend="none")+
#   stat_compare_means(paired = TRUE)
# 
# 
# Preamp <- match_df2 %>%
#   filter(Enrichment == "PREAMP") %>%
#   select(Prop_diff)
# 
# swga <- match_df2 %>%
#   filter(Enrichment == "sWGA") %>%
#   select(Prop_diff)
# d <- data.frame(Preamp=Preamp, swga=swga) %>%
#   rename(`Targ. Pre-amp` = Prop_diff, sWGA = Prop_diff.1)
# str(match_df2)
# check_match <- match_df2 %>%
#   select(comp_id,is_minor_allele,Reads,MOI,Prop_value,Prop_ref,Prop_diff,Enrichment,Source,Extraction,mh_allele_snps) 
# check_match
# st_pd1 <- match_df2 %>%
#   ungroup() %>%
#   wilcox_test(Prop_diff~Enrichment,paired=TRUE) %>%
#   add_significance() %>%
#   add_xy_position(x = "Enrichment", dodge = 0.8)
# st_pd1
# minor_consist <- ggpaired(d,cond1 = "Targ. Pre-amp", cond2 = "sWGA",
#          fill = "condition", line.color = "gray", line.size = 0.4,
#          xlab="Enrichment",ylab="Minor allele\nprevalence\ndifference (%)",
#          palette = c(azure,green),ggtheme = theme_bw(),legend="none", 
#          font.x = list(size = 8, face = "plain", color ="black"),
#          font.tickslab = list(size = 8, face = "plain", color ="black"),
#          font.y = list(size = 8, face = "plain", color ="black"),)+
#   stat_compare_means(paired = TRUE,  aes(label = paste0("p =", ..p.format..," ",p.signif)), 
#                      tip.length = 0,hide.ns = FALSE,size=2.5,label.x=1,label.y=90)
# 
# minor_consist
# mean(d$`Targ. Pre-amp`)
# mean(d$sWGA)
# nrow(d)
match_filtered <- match_df %>%
  filter(!is.na(mhap_max_check)) %>%
  mutate(Reads_hap = Reads*Prop_value) %>%
  mutate(Enrichment = ifelse(Enrichment == "PREAMP","Targ. Pre-amp.",Enrichment)) %>%
  mutate(`Haplotype match` = ifelse(mhap_max_check =="Fail:ambiguous","Fail",mhap_max_check))

match_pass <- match_filtered %>%
  filter(`Haplotype match` == "Pass")

mfchck <- match_df %>%
  #filter(!is.na(mhap_max_check)) %>%
  mutate(Reads_hap = Reads*Prop_value) %>%
  mutate(Enrichment = ifelse(Enrichment == "PREAMP","Targ. Pre-amp.",Enrichment)) %>%
  mutate(`Haplotype match` = ifelse(mhap_max_check =="Fail:ambiguous","Fail",mhap_max_check)) %>%
  dplyr::filter(SID == "GBPCD_P04_D04" & Primer == "PvP01_14_v1_2261601_2261800") %>%
  dplyr::filter(SID == "GBPCD_P04_D04") %>%
  #dplyr::filter(SID=="AH-46" & Primer =="PvP01_07_v1_1417201_1417400") %>%
  #dplyr::filter(SID=="AH-46") %>%
  #dplyr::filter(SID =="GBPCD_P08_D11" & Primer == "PvP01_09_v1_1814501_1814700") %>%
  #dplyr::filter(SID =="GBPCD_P08_D11") %>%
  #filter(Prop_value < 0.5) %>%
 # filter(max_MOI==3) %>%
  #filter(Prop_rank ==2) %>%
  select(SID,Primer,Prop_rank,max_MOI,Reads,MOI, Prop_value,mh_allele, mh_allele_snps,
         Enrichment, Source,Extraction,n_treatments,Reads_hap,`Haplotype match`)

write.csv(mfchck, "results/m_check.csv")

match_fail <- match_filtered %>%
  filter(`Haplotype match` == "Fail") 
write.csv(match_fail,"results/failed_matches.csv")
n_pass <- match_filtered %>%
  #filter(Prop_rank==0) %>%
  group_by(SID,Primer, Source,Extraction) %>%
  group_split()

comps_vect <- c()
pass_out <- c()
for (i in 1:length(n_pass)) {
  comps_per <- nrow(n_pass[[i]]) 
  comps_vect[i] <- comps_per
  pass_out[[i]] <- n_pass[[i]] %>%
    mutate(n_comps = comps_per)
}
# 900 comparisons source+SID+primer
pass_out_df <- do.call(rbind,pass_out) %>%
  filter(n_comps >1)

pass_out_df$Enrichment = factor(pass_out_df$Enrichment, levels=c("None","Targ. Pre-amp.","sWGA"))
match_fail$Enrichment = factor(match_fail$Enrichment, levels=c("None","Targ. Pre-amp.","sWGA"))

 hap_consist <- ggplot(data=pass_out_df,aes(x=Prop_value,y=Reads_hap),color="white",alpha=0) +
  geom_point(data=pass_out_df,aes(color=`Haplotype match`,shape=`Haplotype match`),alpha=0.4,size=1.2) +
  #geom_point(data=pass_out_df, alpha=0.3, color=azure,size=1.2) +
   scale_shape_manual(values=c(17,16)) +
  geom_point(data=match_fail,alpha=1,color=vermillion,size=1.4,shape=17) +
  #geom_point(aes(alpha=`Haplotype match`),size=2) +
  #scale_alpha_discrete(range=c(1,0.2)) +
  facet_grid(cols=vars(Enrichment)) +
  scale_y_continuous(trans="log2") +
  #scale_x_continuous(trans="log10") +
  theme_bw() +
  scale_color_manual(values=c(vermillion,azure)) +
  ylab("Haplotype\nreads") + 
  xlab("Haplotype proportion of reads") +
  theme(axis.text.x = element_text(angle=60,vjust=0.5,hjust=0.5,color="black", size =8),
        axis.text.y = element_text(hjust=0,color="black", size=8),
        axis.title.x = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=8),
        legend.title = element_text(color="black", size=8),
        legend.text = element_text(color="black", size=8),
        legend.spacing.x= unit(.1,"cm"),
        legend.spacing.y= unit(.2,"cm"),
        legend.position = "right",
        legend.background = element_rect(fill = "gray95",color="black",size=0.4),
        strip.background = element_rect(color="black", fill="gray95", size=0.4, linetype="solid"),
        strip.text.x = element_text(size = 8)) 

 MS_READ_PLOT2 <- ggarrange(ggarrange(moi_consist, minor_consist, ncol = 2, labels = c("A", "B"),vjust=1.5), # Second row with box and dot plots
                            hap_consist,
                            nrow = 2, 
                            labels = c("","C"),
                            heights = c(1.2, 1)# Labels of the scatter plot
 ) 

w <- 5
h <- 3.5
ggsave(
  str_c(arg$out_base, ".pdf"), 
  plot = MS_READ_PLOT2, 
  width = w, 
  height = h, 
  units = "in"
)
ggsave(
  str_c(arg$out_base, ".png"), 
  plot = MS_READ_PLOT2, 
  width = w, 
  height = h, 
  units = "in"
)

#chck_matches <- match_filtered %>%
#  filter(mhap_max_check!="Pass") %>%
#  select(SID,Primer,Prop_rank,Prop_value,mh_allele,mhaps_long_mode)
#chck_matches

#chck_match_out <- c()
#for (i in 1:length(chck_matches)) {
#  match2chck <- chck_matches[i,]
#  chck_matches2 <- match_filtered %>%
#    filter(SID %in% match2chck$SID) %>%
#    filter(Primer %in% match2chck$Primer) %>%
#    select(SID,Primer,Enrichment,Source,Extraction,MOI,Reads,Prop_rank,Prop_value,mh_allele,mhaps_long_mode,mhap_max_check,mhap_max_ambig) %>%
#    mutate(Reads_hap = Reads*Prop_value) %>%
#  chck_match_out[[i]] <- chck_matches2
#}
##chck_out_df <- do.call(rbind,chck_match_out)
#length(chck_match_out)
#chck_match_out[[6]]

#chk_df <- do.call(rbind,chck_match_out) %>%
#  filter(mhap_max_check!="Pass")


# 
# match_df %>% filter(is.na(mhap_max_check))
# 
# ggplot(data=match_df, aes(x=MOI,fill=mhap_max_check)) +
#   geom_histogram()
# ggplot(data=match_df, aes(x=ParasiteDensity,fill=mhap_max_check)) +
#   geom_histogram()
# 
# ggplot(data=match_df, aes(x=Reads,fill=mhap_max_check)) +
#   geom_histogram()  +
#   scale_x_continuous(trans="log10")
# 
# ggplot(data=match_df,aes(x=Primer,fill=mhap_max_check)) +
#   geom_bar()
# 
# #########################################################
# ####### SD of prop reads ################################
# ########################################################
# 
# propsd <- match_df %>%
#   select(SID,Primer,ParasiteDensity,Reads,MOI,Prop_sd) %>%
#   distinct() %>%
#   filter(MOI>1) %>%
#   mutate(sd_bin = cut(Prop_sd, breaks=c(-0.01,0.05,0.15,0.25,1))) 
#   
# 
# ggplot(data=propsd,aes(x=ParasiteDensity,y=Prop_sd)) +
#   #facet_wrap(~MOI) +
#   geom_point()
# 
# ggplot(data=propsd,aes(x=Primer,y=Prop_sd)) +
#   geom_boxplot()
# 
# ggplot(data=propsd, aes(x=sd_bin,fill=as.factor(MOI))) +
#   geom_histogram(stat="count") 
# 
# ggplot(data=propsd, aes(x=as.factor(sd_bin), y=Reads)) +
#   geom_boxplot() +
#   scale_y_continuous(trans="log10")
# 
# ggplot(data=propsd, aes(x=as.factor(MOI), y=Reads)) +
#   geom_boxplot() +
#   scale_y_continuous(trans="log10")
# 
# ggplot(data=propsd, aes(x=Prop_sd)) +
#   geom_histogram() 
# 
# quantile(propsd$Prop_sd,seq(0,1,0.05))
# 
# ggplot(data=propsd,aes(x=ParasiteDensity,y=Prop_sd)) +
#   geom_point() +
#   scale_y_continuous(trans="log10") +
#   geom_smooth(method="lm")
# 
# ggplot(data=propsd,aes(x=as.factor(MOI),y=Prop_sd)) +
#   geom_boxplot() 
#   
# 
# 
# ggplot(data=propsd,aes(x=Reads,y=Prop_sd)) +
#   geom_point() +
#   scale_y_continuous(trans="log10") +
#   scale_x_continuous(trans="log10") +
#   geom_smooth(method="lm")
# 
# ch1 <- match_df %>% filter(mhap_max_check == "Fail" & Prop_value ==1) %>% arrange(SID, Primer) %>%
#   select(-Total_reads_sample,-ParasiteDensity) %>% filter(SID=="GBPCD_P01_A01")
# 
# ch2 <- match_df %>% filter(SID=="GBPCD_P01_A01" & Primer == "PvP01_13_v1_258501_258700") 
# 
# 
# ggplot(data=match_df, aes(x=mhap_mode_check,y=Prop_value)) +
#   geom_boxplot()
# 
# ggplot(data=match_df, aes(x=mhap_max_check,y=Prop_value)) +
#   geom_boxplot()
# 
# ggplot(data=match_df, aes(x=mhap_mode_check,y=Reads)) +
#   geom_boxplot()
# 
# ggplot(data=match_df, aes(x=mhap_max_check,y=Reads)) +
#   geom_boxplot()
# 
# ggplot(data=match_df, aes(x=mhap_mode_check,y=Prop_sd)) +
#   geom_boxplot()
# 
# ggplot(data=match_df, aes(x=mhap_max_check,y=Prop_sd)) +
#   geom_boxplot()
# 
# ggplot(data=match_df, aes(x=mhap_mode_check,y=MOI)) +
#   geom_boxplot()
# 
# ggplot(data=match_df, aes(x=mhap_max_check,y=MOI)) +
#   geom_boxplot()
# 
# ggplot(data=match_df, aes(x=MOI_mode_check,y=Reads)) +
#   geom_boxplot()
# 
# ggplot(data=match_df, aes(x=MOI_max_check,y=Reads)) +
#   geom_boxplot()
# 
# #########################################################
# ### MOI consistency #######################################
# ##########################################################
# 
# match_df_unique_moi <- match_df %>%
#   select(SID,Primer,rep,MOI,Reads,ParasiteDensity,MOI_mode_check,MOI_max_check,N_reps) %>%
#   distinct()
# 
# ggplot(data=match_df_unique_moi,aes(x=Reads,fill=MOI_max_check)) +
#   geom_histogram() +
#   scale_x_continuous(trans="log10")
# 
# ggplot(data=match_df_unique_moi,aes(x=Reads,fill=MOI_mode_check)) +
#   geom_histogram() +
#   scale_x_continuous(trans="log10")
# 
# ggplot(data=match_df_unique_moi,aes(x=ParasiteDensity,fill=MOI_max_check)) +
#   geom_histogram() +
#   scale_x_continuous(trans="log10")
# ggplot(data=match_df,aes(x=MOI_max_check,fill=rep)) +
#   geom_bar()
# 
# # get distribution of replicates
# 
# reps_df <- match_df_unique_moi %>%
#   select(SID,Primer,N_reps) %>%
#   distinct()
#   
# ggplot(data=reps_df,aes(x=N_reps)) +
#   geom_histogram() +
#   scale_x_continuous(trans="log10")
# 
# 
# 
# ggplot
# 
# ggplot(data=match_df, aes(x=Reads,y=Prop_sd)) +
#   geom_point()
