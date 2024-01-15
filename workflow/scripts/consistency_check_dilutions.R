# Original Author: Elizabeth Hemming-Schroeder
# Revised By: Alfred Hubbard

library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggpmisc)
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

new_df %>% filter(str_detect(Primer, "09_v1_1564401")) %>% filter(str_detect(SID,"A03")) %>% filter(SeqRun=="vera") %>% filter(ParasiteDensity >27175) %>% select(-Sequence)

# reps are indicated by
# 1) extraction method
# 2) Dilution
# 3) Primer_set
# 4) Treatment

# split by SID then identify 
# may2022 first
splitx <- new_df %>%
  mutate(mh_allele = LETTERS[mh_allele]) %>%
  filter(Primer_set == "OLD") %>%
  filter(!startsWith("PvSNP",Primer)) %>%
  filter(SeqRun=="vera") %>%
  select(-Source,-Treatment) %>%
  group_by(SID, Primer) %>%
  select(Dilution,Primer_set,ParasiteDensity,
         Total_reads_sample,max_MOI,Reads,MOI,Prop_rank,Prop_value,
         mh_allele_snps,mh_allele) %>%
  mutate(rep = paste0(Primer_set,"+",Dilution)) %>%
  #select(-Treatment,-Source) %>%
  group_split()
# 
# splitx[[1300]]
# 
# getidx <- c()
# 
# for (i in 1:length(splitx)) {
#   chck1 <- splitx[[i]] %>% filter(str_detect(SID,"A03")) %>% filter(str_detect(Primer,"09_v1_1564401"))
#   chck2 <- nrow(chck1)
#   if(chck2 >= 1) {
#     getidx[i] <- i
#   }
# }
# 
 i=130
j=1

new_listx <- c()
new_list <- c()
for (i in 1:length(splitx)) {
  
  # check if there are replicates
  check_reps <- length(unique(splitx[[i]]$rep))
  
 MAX_MOI_allREPS <- max(splitx[[i]]$max_MOI)
  
  this_df <- splitx[[i]] %>%
    mutate(N_reps = check_reps) %>%
    mutate(MAX_MOI_allREPS)
  
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
  
  
  mhaps_twicex <- paste(mhaps_df[,1],collapse=",")
  mhaps_twicexx <- strsplit(mhaps_twicex,",")[[1]]
  #mhaps_twicexx <- c("B","B","D","D")
  mhaps_twice_idx <- which(mhaps_twicexx>1)
  mhaps_twice <- unique(mhaps_twicexx[mhaps_twice_idx])
  
  # get reference ("Neat") haplotype
  ref_hapsx <- this_df %>%
    dplyr::filter(Dilution =="Neat") 
  ref_haps <- paste(ref_hapsx$mh_allele, collapse = ":")
 # ref_haps <- ref_hapsx[,"mh_allele"]
  
  ref_haps1 <- str_split(ref_haps,":")[[1]]
  twice_or_refx <- c(mhaps_twice,ref_haps1)
  twice_or_refxx <- unique(twice_or_refx)
  twice_or_ref <- paste(twice_or_refxx, collapse = ":")
  
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
    mutate(mhap_ref_pres = ifelse(length(ref_haps)==0,FALSE,TRUE)) %>%
    mutate(mhap_ref_check = ifelse(mhap_ref_pres==TRUE & str_count(twice_or_ref,mh_allele) >=1,"Pass","Fail")) %>%
    mutate(mhap_max_ambig = str_count(mhaps_long_mode,":")+1) %>%
    mutate(mhap_mode_ambig = str_count(mhaps_mode,":")+1) %>%
    mutate(mhap_mode_check = ifelse(str_count(mhaps_mode,mh_allele)==mhap_mode_ambig,"Pass",
                                    ifelse(str_count(mhaps_mode,mh_allele)>=1,"Fail:ambiguous",
                                           ifelse(str_count(mhaps_mode,mh_allele)==0,"Fail:unambiguous",NA)))) %>%
    mutate(mhap_max_check = ifelse(str_count(mhaps_long_mode,mh_allele)==mhap_max_ambig,"Pass",
                                   ifelse(str_count(mhaps_long_mode,mh_allele)>=1,"Fail:ambiguous",
                                          ifelse(str_count(mhaps_long_mode,mh_allele)==0,"Fail:unambiguous",NA)))) %>%
    mutate(MOI_mode = moi_mode) %>%
    mutate(MOI_max = moi_max) %>%
    #mutate(MOI_2ndmax = moi_2ndmax) %>%
    mutate(MOI_mode_check = ifelse(MOI==MOI_mode,"Pass","Fail")) %>%
    mutate(MOI_max_check  = ifelse(MOI==MOI_max,"Pass","Fail"))  %>%
    mutate(rep= gsub("\\+","@",rep)) %>%
    tidyr::separate(col=rep,into=c("Primer_set","Dilution"),sep="@")
  
  Minor_allele <- this_df_out %>%
    filter(Prop_rank>0) %>%
    filter(Dilution=="Neat") %>%
    filter(Primer_set=="OLD") %>%
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
        mutate(Prop_ref = NA) }
    if (nrow(dfx_split)>1) {
      sd_prop <- sd(dfx_split$Prop_value)
      dfx_split_none <- dfx_split %>%
        filter(Dilution=="Neat") %>%
        filter(Primer_set=="OLD") %>%
        filter(Reads >10) %>%
        filter(MOI_max_check == "Pass")
      if (nrow(dfx_split_none) >= 1) {
        prop_ref <- mean(dfx_split_none$Prop_value,na.rm=TRUE)
        dfx_split_out <- dfx_split %>%
          mutate(Prop_sd = sd_prop) %>%
          mutate(Prop_ref = prop_ref) }
      if (nrow(dfx_split_none) <1) {
        dfx_split_out <- dfx_split %>%
          mutate(Prop_sd = sd_prop) %>%
          mutate(Prop_ref = NA) }
    }
    
    lsplit_list[[k]] <- dfx_split_out
  }
  OUT_DF <- do.call(rbind,lsplit_list)
  new_list[[i]] <- OUT_DF
}


match_df <- do.call(rbind,new_list) %>%
  mutate(Prop_diff = abs(Prop_ref - Prop_value)) %>%
  mutate(Reads_bin = cut(Reads, breaks=c(0,10,100,1000,51000))) %>%
  mutate(ParasiteDensity_bin =cut(ParasiteDensity,breaks=c(0,100,1000,10000,20000))) 

match_df1 <- match_df %>% 
  filter(mhap_max_check=="Pass")  %>%
  filter(!is.na(Prop_ref)) %>%
  filter(n_treatments>=2) %>%
  filter(MOI_max >1)

match_df2 <- match_df1 %>%
  mutate(comp_id = paste0(SID,":",Primer_set,":",Dilution,":",Primer,":",mh_allele)) %>%
  filter(is_minor_allele == TRUE) 


# match_df3 <- as.data.frame(match_df2) %>%
#   ungroup() %>%
#   select(comp_id,Enrichment,Prop_diff) %>%
#   # mutate(comp_id = as.factor(comp_id)) %>%
#   mutate(Enrichment = as.factor(Enrichment)) %>%
#   arrange(comp_id,Enrichment) 


####################################33
### split to do paired wilcox test all three groups at 

#########################################################
### microhaplotype consistency ##########################
########################################################

#MOI consistency
MOI_dfa <- match_df %>%
  mutate(comp_id = paste0(SID,":",Primer_set,":",Dilution)) %>%
  mutate(group_id = paste0(SID,":",Primer_set,":",Primer)) %>%
  select(comp_id,group_id,Primer,ParasiteDensity,Reads,MOI,N_reps,Dilution,Primer_set) %>%
  filter(Primer_set=="OLD") %>%
  filter(N_reps>=3) %>%
  filter(ParasiteDensity >=10) %>%
  mutate(Dilution = factor(Dilution,levels = c("Neat", "1:10", "1:100", "1:1000"))) %>%
  mutate(Dilution = as.numeric(Dilution)) %>%
  mutate(PD_log10 = log10(ParasiteDensity)) %>%
  mutate(PD_log10_bin = cut(PD_log10, breaks=c(1,2,3,4,5,6,7))) %>%
  mutate(PD_bin = cut(ParasiteDensity,breaks=c(10,100,1000,10000,31000))) %>%
  mutate(MOI_bin = ifelse(MOI <=2,MOI,ifelse(MOI>=3,3,NA))) %>%
  mutate(group_id = as.factor(group_id)) %>%
  distinct()
 # mutate(Enrichment = ifelse(Enrichment=="PREAMP","Targ. Pre-amp.", Enrichment)) %>%
  #mutate(Enrichment = factor(Enrichment,levels=c("None","Targ. Pre-amp.","sWGA"))) %>%
 # arrange(comp_id,Enrichment) %>%
  # mutate(Enrichment = factor(Enrichment,levels=c("None","PREAMP","sWGA"))) %>%
#  mutate(Enrichment = as.numeric(Enrichment))
#library(rstatix)
#library(ggpubr)

st_pd1 <- MOI_dfa %>%
  ungroup() %>%
  select(-comp_id) %>%
  pairwise_wilcox_test(MOI~PD_bin) %>%
  #friedman_test(MOI~Dilution| group_id) %>%
  adjust_pvalue(method="holm") %>%
  add_significance() %>%
  add_xy_position(x = "PD_bin", dodge = 1,step.increase=1,fun="max",group="PD_bin") 
st_pd1$y.position <- c(5,5.7,6.1,5.3,6.5,5)
MOI_dfb <- MOI_dfa
MOI_dfb$PD_bin <- as.numeric(MOI_dfb$PD_bin)
MOI_dfb$MOI <- jitter(MOI_dfb$MOI,amount=0.15)
MOI_dfb$PD_bin <- jitter(MOI_dfb$PD_bin,amount=0.25)
#MOI_df$Enrichment <- jitter(MOI_df$Enrichment,amount=0.25)
#MOI_df$MOI <- jitter(MOI_df$MOI,amount=0.25)
moi_consist <- ggplot(data=MOI_dfb,aes(x=PD_bin,y=MOI)) +
  #geom_point(aes(group=group_id),alpha=0.1, color="black") +
  geom_point(shape = 16,alpha=0.4,aes(colour=PD_bin)) +
  geom_line(aes(group=group_id),alpha=0.03, color="black") +
  #scale_colour_stepsn(colours=c(vermillion,azure,green,purple),n.breaks=4) +
  binned_scale(aesthetics = "colour", scale_name = "stepsn",
                palette = function(x) c(vermillion,azure,green,purple),
                breaks = c(0.5,1.5,2.5,3.5),
                limits = c(0.5,4.5),
                show.limits = TRUE,
                guide = "colorsteps") +
  #geom_boxplot(data=MOI_dfa, width=0.1) +
  stat_pvalue_manual(st_pd1, label = "{p.adj}{p.adj.signif}", 
                     tip.length = 0,hide.ns = FALSE,size=2.5) +
  #stat_summary(aes(x=PD_log10_bin,y=MOI),fun.y=mean, geom="point", shape=17, size=5, color="black", fill="black") +
  #scale_x_continuous(trans="log10")+
  #stat_pvalue_manual(st_pd1, label = "{p.adj}{p.adj.signif}", 
  #                   tip.length = 0,hide.ns = FALSE,size=2.5) +
  theme_bw() +
  #scale_colour_manual(values=c(vermillion,vermillion,azure,azure,azure,green,green)) +
  #scale_fill_steps2(low=vermillion,mid=azure,high=green,midpoint=2) +
  scale_y_continuous(limits=c(0.5,6.9), minor_breaks = c(1,2,3,4,5), breaks=c(1,2,3,4,5)) +
  scale_x_continuous(breaks =c(1,2,3,4),minor_breaks=c(1,2,3,4), labels= c("(10,100]","(100,1000]","(1000,10000]","(10000,31000]")) +
  ylab("Alleles in\ninfection") + 
  xlab("Parasite density bin") +
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
MOI_dfz <- match_df %>%
  mutate(comp_id = paste0(SID,":",Primer_set,":",Dilution)) %>%
  mutate(group_id = paste0(SID,":",Primer_set,":",Primer)) %>%
  select(comp_id,group_id,Primer,ParasiteDensity,Reads,MOI,N_reps,Dilution,Primer_set,is_minor_allele,Prop_diff,Prop_value,mh_allele,mhap_max_check) %>%
  filter(Primer_set=="OLD") %>%
  filter(mhap_max_check=="Pass") %>%
  filter(N_reps>=2) %>%
  filter(ParasiteDensity >=10) %>%
  filter(Dilution != "Neat") %>%
  mutate(Dilution = factor(Dilution,levels = c("1:10", "1:100", "1:1000"))) %>%
 # mutate(Dilution = as.numeric(Dilution)) %>%
  mutate(PD_log10 = log10(ParasiteDensity)) %>%
  mutate(PD_log10_bin = cut(PD_log10, breaks=c(1,2,3,4,5,6,7))) %>%
  mutate(PD_bin = cut(ParasiteDensity,breaks=c(10,100,1000,10000,50000))) %>%
  mutate(MOI_bin = ifelse(MOI <=2,MOI,ifelse(MOI>=3,3,NA))) %>%
  mutate(group_id = as.factor(group_id)) %>%
  distinct() %>%
  filter(is_minor_allele == TRUE) %>%
  filter(MOI>=2) %>%
  mutate(G2_id = paste0(group_id,":",mh_allele))
#MOI_dfz %>% arrange(desc(Prop_diff))  %>% select(group_id,mh_allele,Dilution,Prop_value,Prop_diff) %>% filter(group_id =="GBPCD_P01_B07:OLD:PvP01_12_v1_2032801_2033000")
#MOI_dfz %>% arrange(group_id,mh_allele,Dilution,Prop_value,Prop_diff) %>% select(group_id,mh_allele,Dilution,Prop_value,Prop_diff) %>% filter( Dilution == 1)
st_pd2 <- MOI_dfz %>%
  ungroup() %>%
  select(-comp_id) %>%
  pairwise_wilcox_test(Prop_diff~Dilution) %>%
  #friedman_test(MOI~Dilution| group_id) %>%
  adjust_pvalue(method="holm") %>%
  add_significance() %>%
  add_xy_position(x = "Dilution", dodge = 1,step.increase=1,fun="max",group="Dilution") 
st_pd2$y.position = c(0.78,0.85,0.93)
#st_pd2$y.position <- c(4700,4800,4900)

prop_consist <- ggplot(data=MOI_dfz,aes(x=as.factor(PD_bin),y=Prop_diff)) +
 # geom_line(aes(group=G2_id,x=as.factor(PD_bin)),alpha=0.1, color="black") +
 geom_boxplot( width=0.5,aes(fill = as.factor(PD_bin))) +
 # stat_pvalue_manual(st_pd2, label = "{p.adj}{p.adj.signif}", 
  #                   tip.length = 0,hide.ns = FALSE,size=2.5) +
  #scale_y_continuous(trans="log10") +
  #stat_summary(aes(x=PD_log10_bin,y=MOI),fun.y=mean, geom="point", shape=17, size=5, color="black", fill="black") +
  #scale_x_continuous(trans="log10")+
  stat_pvalue_manual(st_pd2, label = "{p.adj}{p.adj.signif}", 
                    tip.length = 0,hide.ns = FALSE,size=2.5) +
  theme_bw() +
  scale_fill_manual(values=c(azure,green,purple)) +
  scale_y_continuous(limits=c(0,1)) +
  #scale_fill_steps2(low=vermillion,mid=azure,high=green,midpoint=2) +
  # scale_y_continuous(breaks=c(1,2,3),minor_breaks = c(1,2,3),limits=c(0.5,4.1)) +
  scale_x_discrete(labels=  c("(10,100]","(100,1000]","(1000,10000]","(10000,31000]")) +
  ylab("Minor allele\nprevalence\ndifference") + 
  xlab("Parasite density bin") +
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

# MOI_dfzz <- MOI_dfz %>%
#   filter(Dilution != 1) 
#   
# 
# st_pd2 <- MOI_dfzz %>%
#   ungroup() %>%
#   select(-comp_id) %>%
#   pairwise_wilcox_test(Prop_diff~Dilution) %>%
#   #friedman_test(MOI~Dilution| group_id) %>%
#   adjust_pvalue(method="holm") %>%
#   add_significance() %>%
#   add_xy_position(x = "Dilution", dodge = 1,step.increase=1,fun="max",group="Dilution") 
# st_pd2
# st_pd2$y.position <- c(0.8,0.95,1.05)

# prop_consist <- ggplot(data=MOI_dfzz,aes(x=as.factor(Dilution),y=Prop_diff)) +
#   geom_line(aes(group=G2_id,x=as.factor(Dilution)),alpha=0.1, color="black") +
#   geom_point(aes(group=G2_id,x=as.factor(Dilution)),alpha=0.1, color="black") +
#   geom_boxplot( width=0.1) +
#   stat_pvalue_manual(st_pd2, label = "{p.adj}{p.adj.signif}", 
#                     tip.length = 0,hide.ns = FALSE,size=2.5) +
#   #scale_y_continuous(trans="log10") +
#   #stat_summary(aes(x=PD_log10_bin,y=MOI),fun.y=mean, geom="point", shape=17, size=5, color="black", fill="black") +
#   #scale_x_continuous(trans="log10")+
#   #stat_pvalue_manual(st_pd1, label = "{p.adj}{p.adj.signif}", 
#   #                   tip.length = 0,hide.ns = FALSE,size=2.5) +
#   theme_bw() +
#   #scale_colour_manual(values=c(vermillion,vermillion,azure,azure,azure,green,green)) +
#   #scale_fill_steps2(low=vermillion,mid=azure,high=green,midpoint=2) +
#   # scale_y_continuous(breaks=c(1,2,3),minor_breaks = c(1,2,3),limits=c(0.5,4.1)) +
#   scale_x_discrete(labels= c("1:10","1:100","1:1000")) +
#   ylab("Difference from undiluted\nsample in proportion\nof reads for minor alleles\n") + 
#   xlab("Dilution") +
#   theme(axis.text.x = element_text(angle=0,vjust=0.5,hjust=0.5,color="black", size =8),
#         axis.text.y = element_text(hjust=0,color="black", size=8),
#         axis.title.x = element_text(color="black", size=8),
#         axis.title.y = element_text(color="black", size=8),
#         legend.title = element_text(color="black", size=8),
#         legend.text = element_text(color="black", size=8),
#         legend.spacing.x= unit(.1,"cm"),
#         legend.spacing.y= unit(.2,"cm"),
#         legend.position = "none",
#         legend.background = element_rect(fill = "white",color="black"),
#         strip.background = element_rect(color="black", fill="gray95", size=0.4, linetype="solid"),
#         strip.text.x = element_text(size = 8)) 
# prop_consist

#LEFT HERE

match_filtered <- match_df %>%
  #filter(mhap_ref_pres == TRUE) %>%
 # filter(!is.na(mhap_max_check)) %>%
  mutate(Reads_hap = Reads*Prop_value) %>%
  #mutate(MOI_diff = max_MOI-MOI) %>%
  mutate(MOI_diff = MAX_MOI_allREPS-MOI) %>%
  #mutate(MOI_diff = MAX_MOI_allREPS-max_MOI) %>%
  mutate(MOI_diff2 = ifelse(MOI_diff >= 2,"ICD >1",
                            ifelse(MOI_diff == 0, " ICD 0",
                                   ifelse(MOI_diff == 1," ICD 1",MOI_diff)))) %>%
  mutate(Dilution2 = paste0("DF ",Dilution)) %>%
  filter(N_reps>=4) %>%
 # mutate(Enrichment = ifelse(Enrichment == "PREAMP","Targ. Pre-amp.",Enrichment)) %>%
 filter(Dilution != "Neat") %>%
  mutate(`Haplotype match` = ifelse(mhap_ref_check =="Fail:ambiguous" |mhap_ref_check =="Fail:unambiguous" ,"Fail",mhap_ref_check)) %>%
  filter(!is.na(`Haplotype match`))

match_pass <- match_filtered %>%
  #filter(mhap_max_check== "Pass")
  filter(`Haplotype match` == "Pass")
match_fail <- match_filtered %>%
 # filter(mhap_max_check== "Fail")
 filter(`Haplotype match` == "Fail") 

write.csv(match_fail,"results/match_fail.csv")

n_pass <- match_filtered %>%
  #filter(Prop_rank==0) %>%
  group_by(SID,Primer) %>%
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


hap_consist <- ggplot(data=pass_out_df,aes(x=Prop_value,y=Reads),color="white",alpha=0) +
  geom_point(data=pass_out_df,aes(color=`Haplotype match`,shape=`Haplotype match`),alpha=0.4,size=1.2) +
  #geom_point(data=pass_out_df, alpha=0.3, color=azure,size=1.2) +
 # scale_shape_manual(values=c(17,16)) +
  geom_point(data=match_fail,alpha=1,color=vermillion,size=1.4,shape=17) +
  #geom_point(aes(alpha=`Haplotype match`),size=2) +
  #scale_alpha_discrete(range=c(1,0.2)) +
  facet_grid(cols=vars(Dilution2),rows = vars(MAX_MOI_allREPS),scales="free_y") +
  scale_y_continuous(trans="log2") +
  #scale_x_continuous(trans="log10") +
  theme_bw() +
  #scale_color_manual(values=c(vermillion,azure)) +
  ylab("Haplotype\nreads") + 
  xlab("Haplotype proportion of reads") +
  theme(axis.text.x = element_text(angle=0,vjust=0.5,hjust=0.5,color="black", size =8),
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


hap_consist <- ggplot(data=pass_out_df,aes(x=Prop_value,y=Reads),color="white",alpha=0) +
  geom_point(data=pass_out_df,aes(color=`Haplotype match`,shape=`Haplotype match`),alpha=0.4,size=1.2) +
  #geom_point(data=pass_out_df, alpha=0.3, color=azure,size=1.2) +
   scale_shape_manual(values=c(17,16)) +
  geom_point(data=match_fail,alpha=1,color=vermillion,size=1.4,shape=17) +
  #geom_point(aes(alpha=`Haplotype match`),size=2) +
  #scale_alpha_discrete(range=c(1,0.2)) +
  facet_grid(cols=vars(Dilution2),rows = vars(MOI_diff2)) +
  scale_y_continuous(trans="log2") +
  #scale_x_continuous(trans="log10") +
  theme_bw() +
  scale_color_manual(values=c(vermillion,azure)) +
  ylab("Haplotype\nreads") + 
  xlab("Haplotype proportion of reads") +
  theme(axis.text.x = element_text(angle=0,vjust=0.5,hjust=0.5,color="black", size =8),
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

# number of reads by primer by sample
vera <- read.table("resources/vera_readcount.txt") %>%
  dplyr::filter(V3>0) %>%
  mutate(treatment = ifelse(startsWith(V2,"PV_SWGA_OLD"),"OLD","NEW")) %>%
  group_by(V1) %>% 
  summarise(Primer = V1, Sample = V2, treatment = treatment, Reads = V3, mean_reads = (median(V3)+(mean(V3)*0.01))) %>%
  ungroup() %>%
  filter(grepl("R1|R2",Primer)==FALSE)

# number of on target reads per sample
n_reads <- vera %>%
  separate(Sample,into=c("j1","j2","SeqID","j3","j4")) %>%
  select(-starts_with("j")) %>%
  group_by(SeqID,treatment) %>%
  summarise(SeqID = SeqID,treatment=treatment,OT_reads = sum(Reads)) %>%
  ungroup() %>%
  distinct()

# total number of reads per sample
vera_tr <- read.table("resources/vera_totalreadcount.txt")  %>%
  mutate(treatment = ifelse(startsWith(V1,"PV_SWGA_OLD"),"OLD","NEW")) %>%
  filter(V1 != "Undetermined_S0_L001") %>%
  separate(V1,into=c("j1","j2","SeqID","j3","j4")) %>%
  select(-starts_with("j")) %>%
  rename(Total_reads = V2)

# compare on target reads to total reads
vera_comp <- left_join(n_reads,vera_tr) %>%
  mutate(Perc_OT = (OT_reads/Total_reads)*100) 

# metadata
meta_data <- read.csv("resources/Pv_samples_Idaho.csv") %>%
  mutate(Well = substr(SeqPlate_Well,4,6)) %>%
  select(-SeqPlate_Well) %>%
  mutate(ParasiteDensity= as.numeric(ParasiteDensity)) %>%
  mutate(Mixed=ifelse(!is.na(CT_Pf18s) | !is.na(CT_PfvarATS),"Yes","No"))

# connects barcode to our naming system
names_barcodes <- read.csv("resources/vera_library_2023-03-18_HemmingSchroeder_UCDavis.csv") %>%
  mutate(treatment=ifelse(grepl("new",Library.Name)==TRUE,"NEW",ifelse(grepl("old",Library.Name)==TRUE,"OLD",NA))) %>%
  select(X,treatment,i7.index.Sequence,i5.index.Sequence) %>%
  rename(Well=X)

# connects barcode to sequence id
seqid_barcodes <- read.csv("resources/vera_lane_summary.csv") %>%
  separate(Barcode.sequence, into=c("i7.index.Sequence","i5.index.Sequence")) %>%
  separate(Sample,into=c("j1","j2","SeqID")) %>%
  select(-starts_with("j"))


# connect meta data, to barcode, to sequence name, to read count
vera_combined <- left_join(meta_data,names_barcodes) %>%
  left_join(.,seqid_barcodes) %>%
  left_join(.,vera_comp) %>%
  filter(treatment=="OLD")

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

MS_READ_PLOT2 <- ggarrange(ggarrange(ot_line_vera,moi_consist,prop_consist, ncol = 3, labels = c("A", "B","C"),vjust=1.5,widths=c(0.63,1.1,0.98)), # Second row with box and dot plots
                           hap_consist,
                           nrow = 2, 
                           labels = c("","D"),
                           heights = c(0.8, 1)# Labels of the scatter plot
) 

w <- 8
h <- 5
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
#    chck_match_out[[i]] <- chck_matches2
#}
##chck_out_df <- do.call(rbind,chck_match_out)
#length(chck_match_out)
#chck_match_out[[6]]

#chk_df <- do.call(rbind,chck_match_out) %>%
#  filter(mhap_max_check!="Pass")

