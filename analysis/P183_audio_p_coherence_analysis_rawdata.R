library(tidyverse)
library(glue)
library(readr)
library(magrittr)
library(lmPerm)
library(broom)
library(ggforce)
library(scales)


theme_set(theme_bw())

#PATH_ANALYSIS = "/Volumes/Nexus/Users/busha/Analysis/2020-08-05-audio_p-coherence"
SUBJECT = "DBS3001"
PATH_ANALYSIS = paste0("Z:/Commits/Vibration_artifacts/denoising_method/coherence-analysis/AllSubjects","/",SUBJECT)
PATH_DB = paste0(PATH_ANALYSIS,"/data")
PATH_FIG = paste0(PATH_ANALYSIS,"/fig")


setwd(PATH_ANALYSIS)
coherence_db <- read_tsv('data/coherence_audio_p.txt',
                         col_types=cols(
                           subject = col_character(),
                           starts = col_double(),
                           duration = col_double(),
                           session_id = col_double(),
                           trial_id = col_double(),
                           type = col_character(),
                           syl_id = col_double(),
                           stim = col_character(),
                           stim_volume = col_double(),
                           rms_audio_p = col_double(),
                           channel1 = col_character(),
                           coherence = col_character(),
                           norm1 = col_double(),
                           norm2 = col_double()
                         )) %>%
  rename(electrode = channel1) %>%
  mutate(electrode_type = tolower(str_extract(electrode, "^[a-zA-Z0-9]+"))) %>%
  mutate(coherence = str_replace(coherence,fixed('+-'),'-')) %>%
  mutate(coherence = as.complex(coherence)) %>%
  filter(!is.na(norm1)) %>%
  filter(!is.na(coherence))

# 
# coherence_db %>%
#   group_by(electrode) %>%
#   summarise(N_subjects=length(unique(subject)),electrode_type=unique(electrode_type)[1]) %>%
#   write_tsv('data/electrodes_rename.tsv')
# 
# coherence_db %>%
#   filter(electrode=='dbs_Lb') %>%
#   with(unique(subject))



coherence_db_norm_db <- coherence_db %>%
  mutate(norm_coherence = coherence / (norm1 * norm2))

# 
# agg_xcorr_db <- xcorr_norm_db %>%
#   group_by(subject,session_id,type,electrode_type,electrode) %>%
#   summarise(median_lag = median(lag),
#             wmean_lag = sum(lag*norm_maxcorr)/sum(norm_maxcorr),
#             sd_lag = sd(lag),
#             median_revlag = median(revlag),
#             wmean_revlag = sum(revlag*norm_maxrevcorr)/sum(norm_maxrevcorr),
#             sd_revlag = sd(revlag),
#             norm_maxcorr_q90 = quantile(norm_maxcorr,0.9,na.rm=TRUE),
#             norm_maxrevcorr_q90 = quantile(norm_maxrevcorr,0.9,na.rm=TRUE),
#             norm_maxcorr_median = median(norm_maxcorr,na.rm=TRUE),
#             norm_maxrevcorr_median = median(norm_maxrevcorr,na.rm=TRUE)
#             )
# 


# === Fig 01c corr vs lag for each channel, epoch, session, subject ========
se <- function(x) sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x)))

plot_coherence <- function(data){
  mdata <- data %>%
    group_by(type,electrode) %>%
    summarise(mRNC = mean(Re(norm_coherence),na.rm=TRUE), 
              mINC = mean(Im(norm_coherence),na.rm=TRUE),
              seRNC = se(Re(norm_coherence)), 
              seINC = se(Im(norm_coherence)))
  
  p<-ggplot(data) +
    aes(x=Re(norm_coherence),y=Im(norm_coherence),color=type,shape=type)+
    geom_point(size=0.5,alpha=0.5) +
    geom_vline(xintercept=0,color='grey',size=0.2) +
    geom_hline(yintercept=0,color='grey',size=0.2) +
    # geom_point(aes(x=mRNC*10,y=mINC*10,shape=type),data=mdata,size=3,color='white') +   
    # geom_point(aes(x=mRNC*10,y=mINC*10,color=type,shape=type),data=mdata,size=1.5) +
    # geom_errorbar(aes(x=mRNC*10,y=NULL,color=type,ymax=(mINC+seINC)*10,ymin=(mINC-seINC)*10),data=mdata) +
    # geom_errorbarh(aes(x=NULL,y=mINC*10,color=type,xmax=(mRNC+seRNC)*10,xmin=(mRNC-seRNC)*10),data=mdata) +
    facet_wrap(~electrode) +
    coord_fixed(ratio = 1, xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
  
  return(p)
}


coherence_db_norm_db %>%
  #filter(subject=='DBS3004') %>%
  filter(electrode_type %in% c("ecog","dbs")) %>%
  group_by(subject, session_id) %>%
  group_walk(~ggsave(paste0(PATH_FIG,'/01c_coherence_',.y[[1]],'_S',.y[[2]],'.png'),
                     plot=plot_coherence(.),
                     width=12,height=12,units='in')
  )
# ========================================================================
# 
# coherence_db_norm_db %>%
#   filter(subject=='DBS3004' & session_id==1) %>%
#   filter(electrode_type %in% c("ecog","dbs")) %>%
#   plot_coherence()
# 



# === Fig 01d corr vs lag for each channel, epoch, session, subject ========
se <- function(x) sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x)))

plot_coherence2 <- function(data){
  mdata <- data %>%
    group_by(type,electrode) %>%
    summarise(mRNC = mean(Re(norm_coherence),na.rm=TRUE), 
              mINC = mean(Im(norm_coherence),na.rm=TRUE),
              seRNC = se(Re(norm_coherence)), 
              seINC = se(Im(norm_coherence)),
              seNC = (seRNC + seINC)/2)
  
  p<-ggplot(mdata) +
    aes(x=mRNC/seNC,y=mINC/seNC,color=type,shape=type) +
    geom_point(size=1.5) +
    geom_errorbar(aes(x=mRNC/seNC,y=NULL,color=type,ymax=(mINC+seINC)/seNC,ymin=(mINC-seINC)/seNC)) +
    geom_errorbarh(aes(x=NULL,y=mINC/seNC,color=type,xmax=(mRNC+seRNC)/seNC,xmin=(mRNC-seRNC)/seNC)) +
    geom_point(size=0.2,alpha=0.5) +
    geom_vline(xintercept=0,color='grey',size=0.2) +
    geom_hline(yintercept=0,color='grey',size=0.2) +
    geom_circle(aes(x0=x0,y0=y0,r=r,x=NULL,y=NULL,color=NULL,shape=NULL),data=tibble(x0=0,y0=0,r=1),color='grey',size=0.4)+
    geom_circle(aes(x0=x0,y0=y0,r=r,x=NULL,y=NULL,color=NULL,shape=NULL),data=tibble(x0=0,y0=0,r=3),color='grey',size=0.4)+
    facet_wrap(~electrode) +
    coord_fixed(ratio = 1)
  
  return(p)
}


coherence_db_norm_db %>%
  #filter(subject=='DBS3004') %>%
  filter(electrode_type %in% c("ecog","dbs")) %>%
  group_by(subject, session_id) %>%
  group_walk(~ggsave(paste0(PATH_FIG,'/01d_coherence_',.y[[1]],'_S',.y[[2]],'.png'),
                     plot=plot_coherence2(.),
                     width=12,height=12,units='in')
  )
# ========================================================================


#20201019

se <- function(x) sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x)))

coherence_agg <- coherence_db_norm_db %>%
  filter(electrode_type %in% c('ecog','dbs','macro','micro','v0')) %>%
  group_by(subject,session_id,type,electrode_type,electrode) %>%
  summarise(mRNC = mean(Re(norm_coherence),na.rm=TRUE), 
            mINC = mean(Im(norm_coherence),na.rm=TRUE),
            seRNC = se(Re(norm_coherence)), 
            seINC = se(Im(norm_coherence)),
            seNC = sqrt(seRNC^2 + seINC^2)) %>%
  mutate(rNC = sqrt(mRNC^2 + mINC^2)/seNC)


coherence_agg %>%
  filter(rNC < 20) %>%
  ggplot() +
  aes(x=rNC) + 
  geom_histogram(binwidth = 0.2) + 
  facet_grid(electrode_type~type)


coherence_agg %>%
  filter(rNC < 20) %>%
  ggplot() +
  aes(x=rNC,y=..density..) + 
  geom_histogram(binwidth = 0.2) + 
  geom_vline(xintercept = 3,color='red') +
  facet_grid(electrode_type~type)
ggsave('fig/02_norm_coherence_histogram_by_type.png',width=5,height=6)

coherence_agg %>%
  filter(electrode_type == 'ecog') %>%
  filter(type == 'stim') %>%
  filter(rNC > 3) %>%
  View()


coherence_agg %>%
  group_by(electrode_type,type) %>%
  summarise(signif_coherence=sum(rNC>3,na.rm=TRUE),tot=n()) %>%
  mutate(f = 100*signif_coherence/tot)


# electrode_type type  signif_coherence   tot      f
# <chr>          <chr>            <int> <int>  <dbl>
#   1 dbs            ITI                 10   532  1.88 
# 2 dbs            prod               295   540 54.6  
# 3 dbs            stim                66   532 12.4  
# 4 ecog           ITI                 55 11660  0.472
# 5 ecog           prod              5776 11791 49.0  
# 6 ecog           stim               472 11705  4.03 
# 7 macro          ITI                  0   244  0    
# 8 macro          prod                99   244 40.6  
# 9 macro          stim                12   244  4.92 
# 10 micro          ITI                  0   242  0    
# 11 micro          prod               136   242 56.2  
# 12 micro          stim                17   242  7.02 
# 13 v0             ITI                 24   247  9.72 
# 14 v0             prod               166   251 66.1  
# 15 v0             stim                34   248 13.7  



coherence_agg %>%
  filter(type=='prod')%>%
  filter(str_sub(subject,4,4)=='3') %>%
  ggplot() +
  aes(x=session_id,y=electrode,fill=rNC)+
  geom_tile()+
  #scale_fill_gradient(limits=c(3,10),oob=squish)+
  scale_fill_gradient(limits=c(3,15),oob=squish,low='white',high=muted('red'))+
  facet_grid(electrode_type~subject,scales='free',space='free')
ggsave('fig/03_norm_coherence_raster_3000_series.png',width=30,height=18)


coherence_agg %>%
  filter(as.numeric(str_sub(subject,6,7))<=14) %>%
  filter(electrode_type %in% c('ecog')) %>%  
  filter(str_sub(electrode,6,6)=='1') %>%
  filter(type=='prod')%>%
  filter(str_sub(subject,4,4)=='3') %>%
  ggplot() +
  aes(x=session_id,y=electrode,fill=rNC)+
  geom_tile()+
  #scale_fill_gradient(limits=c(3,10),oob=squish)+
  scale_fill_gradient(limits=c(3,15),oob=squish,low='white',high=muted('red'))+
  facet_grid(electrode_type~subject,scales='free',space='free')
ggsave('fig/04_norm_coherence_raster_upto3014_ecog.png',width=12,height=10)


coherence_agg %>%
  filter(type=='prod')%>%
  filter(str_sub(subject,4,4)=='4') %>%
  ggplot() +
  aes(x=session_id,y=electrode,fill=rNC)+
  geom_tile()+
  #scale_fill_gradient(limits=c(3,10),oob=squish)+
  scale_fill_gradient(limits=c(3,15),oob=squish,low='white',high=muted('red'))+
  facet_grid(electrode_type~subject,scales='free',space='free')
ggsave('fig/05_norm_coherence_raster_4000_series.png',width=30,height=18)




#2020.10.21
se <- function(x) sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x)))

coherence_agg2 <- coherence_db_norm_db %>%
  filter(electrode_type %in% c('ecog','dbs','macro','micro','v0')) %>%
  filter(type=='prod') %>%
  group_by(subject,session_id,electrode_type,electrode) %>%
  summarise(mRNC = mean(Re(norm_coherence),na.rm=TRUE), 
            mINC = mean(Im(norm_coherence),na.rm=TRUE),
            seRNC = se(Re(norm_coherence)), 
            seINC = se(Im(norm_coherence)),
            seNC = sqrt(seRNC^2 + seINC^2)) %>%
  mutate(rNC = sqrt(mRNC^2 + mINC^2)/seNC) %>%
  select(-mRNC,-mINC,-seRNC,-seINC,-seNC)

write_tsv(coherence_agg2,"data/coherence_audio_p_by_electrode.tsv")

coherence_agg3 <- coherence_agg2 %>%
  group_by(subject,session_id) %>%
  summarise(median_rNC = median(rNC,na.rm=TRUE))
write_tsv(coherence_agg3,"data/coherence_audio_p_by_session.tsv")


coherence_agg4 <- coherence_db_norm_db %>%
  filter(electrode_type %in% c('ecog','dbs','macro','micro','v0')) %>%
  filter(type=='prod') %>%
  group_by(subject,session_id,electrode_type,electrode,stim_volume) %>%
  summarise(mRNC = mean(Re(norm_coherence),na.rm=TRUE), 
            mINC = mean(Im(norm_coherence),na.rm=TRUE),
            seRNC = se(Re(norm_coherence)), 
            seINC = se(Im(norm_coherence)),
            seNC = sqrt(seRNC^2 + seINC^2)) %>%
  mutate(rNC = sqrt(mRNC^2 + mINC^2)/seNC) %>%
  select(-mRNC,-mINC,-seRNC,-seINC,-seNC)  %>%
  group_by(subject,session_id,stim_volume) %>%
  summarise(median_rNC = median(rNC,na.rm=TRUE))

write_tsv(coherence_agg4,"data/coherence_audio_p_by_session_stim_volume.tsv")




