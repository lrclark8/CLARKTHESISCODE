
library(googledrive)
library(data.table)
library(dplyr)
library(MARSS)
library(tidyverse)
library(ggplot2)

drive_download('surv7.fin.csv', overwrite = T)
surv7 <- read.csv("surv7.fin.csv")

drive_download('DMRfin.csv', overwrite = T)
DMRfin <- read.csv ("DMRfin.csv")

drive_download('ecomon.fin.csv', overwrite = T)
ecomon_fin <- read.csv("ecomon.fin.csv")


drive_download('stratum_area.csv', overwrite = T)
stratum_area <- read.csv('stratum_area.csv')

#shrimp abun/tow from shrimp survey
n6=read.csv('c:/FWDP2020/nshrimp/data/nshrimp.survey/n6.survdat.shrimp.csv', h=T, sep=',', na.strings = c(NA))


#subset so that date only goes from 2000-2021 for all datasets 

surv7_sub <- surv7 %>% filter(year >= 2000 & year <= 2021)

DMR_sub <- DMRfin %>% filter(year >= 2000 & year <= 2021)

ecomon_sub <- ecomon_fin %>% filter (year >= 2000 & year <= 2021)

n6_sub= n6 %>% filter(year>=2000 &year <=2021)
n6_sub2=merge(data.frame('year'=2000:2021), n6_sub, by='year', all.x=T)


master=data.frame('year'=2000:2021, 'btsabun.fa'=surv7_sub$abundance[surv7_sub$season=='FALL'], 'btsabun.sp'=surv7_sub$abundance[surv7_sub$season=='SPRING'],'btsbio.fa'= surv7_sub$biomass[surv7_sub$season == "FALL"], 
                  'btsbio.sp'= surv7_sub$biomass[surv7_sub$season == 'SPRING'], eco.fa = ecomon_sub$m.caridea[ecomon_sub$season=='FALL'], eco.sp = ecomon_sub$m.caridea[ecomon_sub$season=='SPRING'],
                  'dmrweight.fa'= DMR_sub$stratmean_weight[DMR_sub$season == 'FALL'], 'dmrweight.sp' = DMR_sub$stratmean_weight[DMR_sub$season == 'SPRING'], 'dmrcatch.fa'= DMR_sub$stratmean_catch[DMR_sub$season== 'FALL'],
                  'dmrcatch.sp'= DMR_sub$stratmean_catch[DMR_sub$season == 'SPRING'], 'nshr.sum'=n6_sub2$abundance) 

#zscore()

master$z.btsabun.fa=zscore(master$btsabun.fa)
master$z.btsabun.sp = zscore(master$btsabun.sp)
master$z.btsbio.fa=zscore(master$btsbio.fa)
master$z.btsbio.sp=zscore(master$btsbio.sp)
master$z.eco.fa = zscore (master$eco.fa)
master$z.eco.sp = zscore (master$eco.sp)
master$z.dmrweight.fa = zscore (master$dmrweight.fa)
master$z.dmrweight.sp = zscore (master$dmrweight.sp)
master$z.dmrcatch.fa = zscore (master$dmrcatch.fa)
master$z.dmrcatch.sp = zscore (master$dmrcatch.sp)
master$z.nshr.sum= zscore(master$nshr.sum)

#graph of master ggplot.
library(reshape)

m.=melt(master[,c(1,13,14,21:23)], value.name='z', id='year' )
m.$s1=substr(m.$variable, nchar(as.character(m.$variable))-1, nchar(as.character(m.$variable)))
m.$Season=ifelse(m.$s1=='sp', 'SPRING', ifelse(m.$s1=='fa', 'FALL', 'SUMMER')) 
  
m.$var=ifelse(m.$variable=='z.btsabun.fa'|m.$variable=='z.btsabun.sp', 'BTS', ifelse(m.$variable=='z.dmrcatch.fa'|m.$variable=='z.dmrcatch.sp', 'DMR', ifelse(m.$variable=='z.nshr.sum', 'N.SHRIMP','X'))) 
  
png('c:/FWDP2020/nshrimp/data/nshrimp.survey/master.vars.png', height = 4, width = 4, res=600, units = 'in')
masvar=ggplot()+ylab('Abundance series')+xlab('')+
  geom_line(m., mapping=aes(x=year, y=value, color=Season))+
  geom_point(m.,mapping=aes(x=year, y=value, color=Season), size=0.5)+
  geom_hline(yintercept=0, linetype="dashed")+facet_wrap(~var, ncol=1)+
  scale_color_manual(values = c('#F8766D', '#00BFCF', 'black'))+
  theme_bw()
print(masvar)

dev.off()


#graph of master zscores?
zvars <- master %>%
  select(year, starts_with("z.")) %>%
  pivot_longer(cols = -year, names_to = "variable", values_to = "zscore")

ggplot(zvars, aes(x = year, y = zscore)) +
  geom_line(na.rm = TRUE, linewidth = 1, color = "steelblue") +
  geom_point(na.rm = TRUE, color = "black", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  theme_minimal(base_size = 12) +
  labs(title = " Master Z-scored Time Series ", x = "Year", y = "Z-score") +
  theme(strip.text = element_text(face = "bold")) 


##limiting zscore graph to bts abun, dmrcatch and eco and nshr

zcols <- master %>%
  select(year,
         z.btsabun.fa, z.btsabun.sp,
         z.dmrcatch.fa, z.dmrcatch.sp,
         z.eco.fa, z.eco.sp, z.nshr.sum)

plot_data <- zcols %>%
  pivot_longer(
    cols = -year,
    names_to = "variable",
    values_to = "zscore"
  )

plot_data <- plot_data %>%
  mutate(variable = recode(variable,
                           "z.btsabun.fa" = "BTS abundance (Fall)",
                           "z.btsabun.sp"= "BTS abundance (Spring)",
                           "z.dmrcatch.fa"= "DMR catch (Fall)",
                           "z.dmrcatch.sp"= "DMR catch (Spring)",
                           "z.eco.fa" = "EcoMon abundance (Fall)",
                           "z.eco.sp" = "Ecomon abundance (Spring)",
                           "z.nshr.sum" = "Northern shrimp abundance (Summer)"
                           ))
                           


ggplot(plot_data, aes(x = year, y = zscore)) +
  geom_line(na.rm = TRUE, linewidth = 1, color = "steelblue") +
  geom_point(na.rm = TRUE, color = "black", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  theme_minimal(base_size = 12) +
  labs(title = "Master Z-scored Time Series Plot 2", x = "Year", y = "Z-score") +
  theme(strip.text = element_text(face = "bold"))


#final transposed data for modeling
master.abun=t(master[,c(13,14,17,18,21,22,23)])
#master.biom=t(master[,c()])

#master.abun[is.na(master.abun)] <- 0


#running MARSS model
 survey.data=data.frame()
 
 
 levels.R=c('diagonal and equal', 'diagonal and unequal', 'unconstrained')
 cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
 for (R in levels.R)  {
    for (m in 1:3) {
      survey.dfa.model= list(A='zero', R=R, m=m)
      kemz= MARSS(master.abun[-c(3,4),], model=survey.dfa.model, form='dfa', control=cntl.list)#, covariates=t(master_exp[,c(70)]))
      survey.data=rbind(survey.data,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))
      assign(paste('kemz', m, R, sep='.'), kemz)
    } # end of m loop
  } #end end of R loop
 
  survey.model.tbl=survey.data[order(survey.data$AICc),]
  survey.model.tbl

 write.csv(survey.model.tbl, 'All7.minusecomon.surveys.model.tbl.csv')


#explanatory variables#
drive_download('panbor.csv', overwrite = T)
panbor= read.csv('panbor.csv')

drive_download('panbor.fin.csv', overwrite = T)
panborfin <- read.csv("panbor.fin.csv")

drive_download('allwt.fin.csv', overwrite = T)
allwtfin <- read.csv("allwt.fin.csv")

drive_download('fo.wtmfin.csv', overwrite = T)
fo_fin <- read.csv("fo.wtmfin.csv")

drive_download('nao_means.csv', overwrite = T)
nao_means<- read.csv('nao_means.csv')

drive_download('dmr.env.vars.csv', overwrite = T)
dmr_env <- read.csv("dmr.env.vars.csv")

drive_download('amo_means.csv', overwrite = T)
amo_means <- read.csv ('amo_means.csv')


library(dplyr)
library(tidyr)

#alternative way to set up for zscores?

dmr_env <- mutate(dmr_env, season = toupper(season))


surv7_wide <- surv7_sub %>%
  pivot_wider(
    id_cols = year,
    names_from = season,
    values_from = c(depth, surftemp, bottemp, surfsalin, botsalin, abundance, biomass),
    names_sep = "."
  )

fo_wide <- fo_fin %>%
  pivot_wider(
    id_cols = year,
    names_from = season,
    values_from = perfrq,
    names_prefix = "fo.perfrq."
  )

allwt_wide <- allwtfin %>%
  pivot_wider(
    id_cols = year,
    names_from = season,
    values_from = relmsw_s,
    names_prefix = "allwt.relmsw_s."
  )

dmr_env_wide <- dmr_env %>%
  pivot_wider(
    id_cols = year,
    names_from = season,
    values_from = c(dmr.bt, dmr.bsal, dmr.sst, dmr.ssal),
    names_sep = "."
  )


master_exp <- surv7_wide %>%
  left_join(fo_wide, by = "year") %>%
  left_join(allwt_wide, by = "year") %>%
  left_join(dmr_env_wide, by = "year") %>%
  left_join(nao_means %>% select(year, mean_value), by = "year") %>%
  left_join(amo_means %>% select(year, mean_value), by = "year")

colnames(master_exp)[c(28, 29)] <- c("nao", "amo")

#zscores


master_exp$z.depth.fa        <- zscore(master_exp$depth.FALL)
master_exp$z.depth.sp        <- zscore(master_exp$depth.SPRING)
master_exp$z.surftemp.fa     <- zscore(master_exp$surftemp.FALL)
master_exp$z.surftemp.sp     <- zscore(master_exp$surftemp.SPRING)
master_exp$z.bottemp.fa      <- zscore(master_exp$bottemp.FALL)
master_exp$z.bottemp.sp      <- zscore(master_exp$bottemp.SPRING)
master_exp$z.surfsalin.fa    <- zscore(master_exp$surfsalin.FALL)
master_exp$z.surfsalin.sp    <- zscore(master_exp$surfsalin.SPRING)
master_exp$z.botsalin.fa     <- zscore(master_exp$botsalin.FALL)
master_exp$z.botsalin.sp     <- zscore(master_exp$botsalin.SPRING)

master_exp$z.fo.perfrq.fa    <- zscore(master_exp$fo.perfrq.FALL)
master_exp$z.fo.perfrq.sp    <- zscore(master_exp$fo.perfrq.SPRING)

master_exp$z.allwt.relmsw_s.fa <- zscore(master_exp$allwt.relmsw_s.FALL)
master_exp$z.allwt.relmsw_s.sp <- zscore(master_exp$allwt.relmsw_s.SPRING)

master_exp$z.dmr.env.bt.fa     <- zscore(master_exp$dmr.bt.FALL)
master_exp$z.dmr.env.bt.sp     <- zscore(master_exp$dmr.bt.SPRING)
master_exp$z.dmr.env.bsal.fa   <- zscore(master_exp$dmr.bsal.FALL)
master_exp$z.dmr.env.bsal.sp   <- zscore(master_exp$dmr.bsal.SPRING)
master_exp$z.dmr.env.sst.fa    <- zscore(master_exp$dmr.sst.FALL)
master_exp$z.dmr.env.sst.sp    <- zscore(master_exp$dmr.sst.SPRING)
master_exp$z.dmr.env.ssal.fa   <- zscore(master_exp$dmr.ssal.FALL)
master_exp$z.dmr.env.ssal.sp   <- zscore(master_exp$dmr.ssal.SPRING)




#explanatory marss model will be a one trend unconstrained model




#create lead/lag covariates
master_exp$lead1.fo.sp=lead(master_exp$z.fo.perfrq.sp, 1)
master_exp$lead2.fo.sp=lead(master_exp$z.fo.perfrq.sp, 2)
master_exp$lead3.fo.sp=lead(master_exp$z.fo.perfrq.sp, 3)
master_exp$lead4.fo.sp=lead(master_exp$z.fo.perfrq.sp, 4)
master_exp$lead5.fo.sp=lead(master_exp$z.fo.perfrq.sp, 5)
master_exp$lead1.fo.fa=lead(master_exp$z.fo.perfrq.fa, 1)
master_exp$lead2.fo.fa=lead(master_exp$z.fo.perfrq.fa, 2)
master_exp$lead3.fo.fa=lead(master_exp$z.fo.perfrq.fa, 3)
master_exp$lead4.fo.fa=lead(master_exp$z.fo.perfrq.fa, 4)
master_exp$lead5.fo.fa=lead(master_exp$z.fo.perfrq.fa, 5)

master_exp$lead1.allwt.sp=lead(master_exp$z.allwt.relmsw_s.sp, 1)
master_exp$lead2.allwt.sp=lead(master_exp$z.allwt.relmsw_s.sp, 2)
master_exp$lead3.allwt.sp=lead(master_exp$z.allwt.relmsw_s.sp, 3)
master_exp$lead4.allwt.sp=lead(master_exp$z.allwt.relmsw_s.sp, 4)
master_exp$lead5.allwt.sp=lead(master_exp$z.allwt.relmsw_s.sp, 5)
master_exp$lead1.allwt.fa=lead(master_exp$z.allwt.relmsw_s.fa, 1)
master_exp$lead2.allwt.fa=lead(master_exp$z.allwt.relmsw_s.fa, 2)
master_exp$lead3.allwt.fa=lead(master_exp$z.allwt.relmsw_s.fa, 3)
master_exp$lead4.allwt.fa=lead(master_exp$z.allwt.relmsw_s.fa, 4)
master_exp$lead5.allwt.fa=lead(master_exp$z.allwt.relmsw_s.fa, 5)



#lag
master_exp$lag1.fo.sp=lag(master_exp$z.fo.perfrq.sp,1)
master_exp$lag2.fo.sp=lag(master_exp$z.fo.perfrq.sp,2)
master_exp$lag3.fo.sp=lag(master_exp$z.fo.perfrq.sp,3)
master_exp$lag1.fo.fa=lag(master_exp$z.fo.perfrq.fa,1)
master_exp$lag2.fo.fa=lag(master_exp$z.fo.perfrq.fa,2)
master_exp$lag3.fo.fa=lag(master_exp$z.fo.perfrq.fa,3)

master_exp$lag1.allwt.sp=lag(master_exp$z.allwt.relmsw_s.sp,1)
master_exp$lag2.allwt.sp=lag(master_exp$z.allwt.relmsw_s.sp,2)
master_exp$lag3.allwt.sp=lag(master_exp$z.allwt.relmsw_s.sp,3)
master_exp$lag1.allwt.fa=lag(master_exp$z.allwt.relmsw_s.fa,1)
master_exp$lag2.allwt.fa=lag(master_exp$z.allwt.relmsw_s.fa,2)
master_exp$lag3.allwt.fa=lag(master_exp$z.allwt.relmsw_s.fa,3)




#changing NA values to 0

master_exp[is.na(master_exp)] <- 0

#add 3x4 plot of expl vars sp/fal shown in 1 fig nao through dmr surf sal (24 total); nao and amo annual.
master_exp2=data.frame(master_exp[,c(1,28:51)]) 

library(reshape)
m3=melt(master_exp2, value.name='z', id='year' )
m3$s1=substr(m3$variable, nchar(as.character(m3$variable))-1, nchar(as.character(m3$variable)))

m3$Season=ifelse(m3$s1== 'sp', 'SPRING', ifelse(m3$s1== 'fa', 'FALL', ifelse(m3$s1=='ao', 'WINTER NAO', 'ANNUAL AMO')))
m3$var=ifelse(m3$variable=='z.depth.fa'|m3$variable=='z.depth.sp', 'Depth', ifelse(m3$variable=='z.surftemp.fa'|m3$variable=='z.surftemp.sp', 'SST.BTS', ifelse(m3$variable=='z.bottemp.fa'|m3$variable=='z.bottemp.sp', 'BT.BTS', ifelse(m3$variable=='z.surfsalin.fa'|m3$variable=='z.surfsalin.sp', 'SSal.BTS', ifelse(m3$variable=='z.botsalin.fa'|m3$variable=='z.botsalin.sp', 'BSal.BTS', ifelse(m3$variable=='z.fo.perfrq.fa'|m3$variable=='z.fo.perfrq.sp', 'Diet.FO', ifelse(m3$variable=='z.allwt.relmsw_s.fa'|m3$variable=='z.allwt.relmsw_s.sp', 'Diet.MSW', 
ifelse(m3$variable=='z.dmr.env.bt.fa'|m3$variable=='z.dmr.env.bt.sp', 'BT.DMR', ifelse(m3$variable=='z.dmr.env.bsal.fa'|m3$variable=='z.dmr.env.bsal.sp', 'BSal.DMR', ifelse(m3$variable=='z.dmr.env.sst.fa'|m3$variable=='z.dmr.env.sst.sp', 'SST.DMR', ifelse(m3$variable=='z.dmr.env.ssal.fa'|m3$variable=='z.dmr.env.ssal.sp', 'SSal.DMR', toupper(as.character(m3$variable)))))))))))))

#plot expl vars

m3$Season=factor(m3$Season, levels=c('FALL', 'SPRING','WINTER NAO', 'ANNUAL AMO'))

png('c:/FWDP2020/nshrimp/data/nshrimp.survey/expl.vars2.png', height = 6, width = 8, res=600, units = 'in')
expvar=ggplot()+ylab('Environmental series')+xlab('')+
  geom_line(m3, mapping=aes(x=year, y=value, color=Season))+
  geom_point(m3, mapping=aes(x=year, y=value, color=Season), size=0.5)+
  #geom_smooth(m3, mapping=aes(x=year, y=value, color=Season), alpha=0.15 )+
  geom_hline(yintercept=0, linetype="dashed")+facet_wrap(~var, ncol=4)+
  scale_color_manual(values = c('#F8766D', '#00BFCF', 'black', 'gray55'))+
  theme_bw()
print(expvar)

dev.off()


png(paste('c:/FWDP2020/WGECO/2025/tord/figs/size/', 'temp.seas.png', sep=''), height=6, width = 8, res=600, units='in')
temp=ggplot()+ylab('Temperature (C)') + xlab('')+
  geom_line(temp6822, mapping=aes(x=year, y=bt), color='black')+
  geom_line(temp6822, mapping=aes(x=year, y=sst), color='gray')+facet_wrap(~season, ncol=2)#+geom_smooth(aes(year, med.lon, color=pred.size),method='loess', span=0.8)

print(temp)
dev.off()

#explanatory MARSS models

#11/24  best model-- depth fall, surftemp fall, diet msw fall
levels.R.expl=c('diagonal and equal')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame()

  survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
  kemz.expl= MARSS(master.abun[-c(3:4),], model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(30,32, 42)]))
  
  survey.data.expl=rbind(survey.data.expl,
                         data.frame(cov='3.best', R='diagonal and equal', m=1, logLik=kemz.expl$logLik, K=kemz.expl$num.params, AICc=kemz.expl$AICc, stringsAsFactors=F))
  assign(paste('kemz.expl', names(master_exp[,c(30,32,42)]), sep='.'), kemz.expl)



survey.model.expl.tbl=survey.data.expl[order(survey.data.expl$AICc),]
survey.model.expl.tbl
write.csv(survey.model.expl.tbl, 'c:/FWDP2020/nshrimp/data/nshrimp.survey/BEST.3.tbl.csv')

#fitted data
nshrimp.best.model=survey.model.expl.tbl[1,]  #diagonal and equal
nshrimp.best.model.fitname= paste('kemz.expl')
nshrimp.best.fit=get(nshrimp.best.model.fitname)

nshrimp.best.fit2=fitted(nshrimp.best.fit)
nshrimp.best.fit2$year=c(2000:2021)

#nshrimp.best.fit2$Season=ifelse(nshrimp.best.fit2$.rownames=='z.btsabun.fa'|nshrimp.best.fit2$.rownames=='z.btsabun.', 'Fall', ifelse(nshrimp.best.fit2$.rownames=='z.btsabun.sp', 'Spring', ifelse(nshrimp.best.fit2$.rownames=='x'))
nshrimp.best.fit2$Season=ifelse(substr(nshrimp.best.fit2$.rownames, nchar(nshrimp.best.fit2$.rownames)-1, nchar(nshrimp.best.fit2$.rownames))=='fa', 'FALL', ifelse(substr(nshrimp.best.fit2$.rownames, nchar(nshrimp.best.fit2$.rownames)-1, nchar(nshrimp.best.fit2$.rownames))=='sp', 'SPRING', ifelse(substr(nshrimp.best.fit2$.rownames, nchar(nshrimp.best.fit2$.rownames)-1, nchar(nshrimp.best.fit2$.rownames))=='um', 'SUMMER', 'x')))

nshrimp.best.fit2$var=ifelse(nshrimp.best.fit2$.rownames=='z.btsabun.fa'|nshrimp.best.fit2$.rownames=='z.btsabun.sp', 'BTS', ifelse(nshrimp.best.fit2$.rownames=='z.dmrcatch.fa'|nshrimp.best.fit2$.rownames=='z.dmrcatch.sp', 'DMR', ifelse(nshrimp.best.fit2$.rownames=='z.nshr.sum', 'N.SHRIMP','X'))) 


#plot
library(ggplot2)

png('c:/FWDP2020/nshrimp/data/nshrimp.survey/REDO.best.11-24-25.fitted.png', height = 4, width = 4, res=600, units = 'in')
masvar=ggplot()+ylab('Abundance series')+xlab('')+
  geom_point(nshrimp.best.fit2, mapping=aes(x=year, y=y, color=Season), size=0.5)+
  geom_line(nshrimp.best.fit2, mapping=aes(x=year, y=.fitted, color=Season))+
  geom_hline(yintercept=0, linetype="dashed")+facet_wrap(~var, ncol=1)+
  scale_color_manual(values = c('#F8766D', '#00BFCF', 'black'))+
  theme_bw()
print(masvar)

dev.off()


png(paste('c:/FWDP2020/nshrimp/data/', 'REDO.best.11-24-25.fitted.png', sep=''), height=4, width = 4, res=600, units='in')
bts.fitted2=ggplot()+ylab('Shrimp abundance series') + xlab('')+theme_bw()+
  geom_point(bts.best.fit2, mapping=aes(x=year, y=y), shape=21, color='black', fill='white')+geom_line(bts.best.fit2, mapping=aes(x=year, y=.fitted))+geom_hline(yintercept=0, linetype="dashed")+facet_wrap(~Season, ncol=2)#+geom_smooth(aes(year, med.lon, color=pred.size),method='loess', span=0.8)

print(bts.fitted2)
dev.off()

write.csv(bts.best.fit2, 'bts.best.fit2.csv')




levels.R.expl=c('diagonal and equal')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame()

for (i in c(28:51)){  #all expl vars plus lead/lag for diet

survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun[-c(3:4),], model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(i)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(cov=names(master_exp)[i], R='diagonal and equal', m=1, logLik=kemz.expl$logLik, K=kemz.expl$num.params, AICc=kemz.expl$AICc, stringsAsFactors=F))
assign(paste('kemz.expl', names(master_exp)[i], sep='.'), kemz.expl)
}


survey.model.expl.tbl=survey.data.expl[order(survey.data.expl$AICc),]
survey.model.expl.tbl
write.csv(survey.model.expl.tbl, 'All7.minusecomon.minusleadlag.survey.model.single.expl.tbl.csv')


#test nmfs (2) abundance with combos of expl vars: lead1.fo.fa, amo, allwt.fa, allwt.sp.

levels.R.expl=c('diagonal and equal')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame()
survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)  #A= 'zero' 1. because needs to be fixed when estimating Z and also because data are already demeaned (zscored)
kemz.expl= MARSS(master.abun[-c(3:4),], model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(32, 42)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(model='kemz.expl', R='diagonal and equal', m=1, logLik=kemz.expl$logLik, K=kemz.expl$num.params, AICc=kemz.expl$AICc, stringsAsFactors=F))
assign(paste('kemz.expl', 'combo.vars', sep='.'), kemz.expl)

survey.model.expl.tbl=survey.data.expl[order(survey.data.expl$AICc),]
survey.model.expl.tbl
write.csv(survey.model.expl.tbl, 'minusecomon.combotesttop2.survey.model.expl.tbl.csv') 
#so far fall surftemp and fall allwt_relmsw as combo is best.
#test other combos

#fitted data
bts.best.model=survey.model.expl.tbl[1,]  #unconstrained did not converge
bts.best.model.fitname= paste('kemz.expl')#, bts.best.model$m,bts.best.model$R, sep='.')
bts.best.fit=get(bts.best.model.fitname)

bts.best.fit2=fitted(bts.best.fit)
bts.best.fit2$year=c(2000:2021)

bts.best.fit2$Season=ifelse(bts.best.fit2$.rownames=='z.btsabun.fa', 'Fall', ifelse(bts.best.fit2$.rownames=='z.btsabun.sp', 'Spring', 'x'))
#plot
library(ggplot2)
png(paste('c:/FWDP2020/nshrimp/data/', 'bts.fitted.png', sep=''), height=4, width = 8, res=600, units='in')
bts.fitted2=ggplot()+ylab('Shrimp abundance series') + xlab('')+theme_bw()+
  geom_point(bts.best.fit2, mapping=aes(x=year, y=y), shape=21, color='black', fill='white')+geom_line(bts.best.fit2, mapping=aes(x=year, y=.fitted))+geom_hline(yintercept=0, linetype="dashed")+facet_wrap(~Season, ncol=2)#+geom_smooth(aes(year, med.lon, color=pred.size),method='loess', span=0.8)

print(bts.fitted2)
dev.off()

write.csv(bts.best.fit2, 'bts.best.fit2.csv')

#need to also do fitted plot for dmr


#lead1.fo.fa+allwt.fa+allwt.sp+amo AICc= 108.6786
#lead1.fo.fa+amo = 91.19585 !!!
#allwt.sp+allwt.fa = 100.7707
#amo test = 94.63812
#lead1.fo.fa test = 95.26368
#allwt.fa+amo= 98.63864


#####is this useful?


levels.R.expl <- c('diagonal and equal')
cntl.list <- list(minit = 1000, maxit = 100000, allow.degen = FALSE, safe = TRUE)
survey.data.expl <- data.frame()

covariate_indices <- 28:51  

combo_2vars <- combn(covariate_indices, 2, simplify = FALSE)

for (combo in combo_2vars) {
  
  cov_names <- names(master_exp)[combo]
  
  survey.dfa.model.expl <- list(A = 'zero', R = levels.R.expl, m = 1)
  kemz.expl <- MARSS(master.abun[-c(3:4), ],
                     model = survey.dfa.model.expl,
                     form = 'dfa',
                     control = cntl.list,
                     covariates = t(master_exp[, combo]))
  
  survey.data.expl <- rbind(survey.data.expl,
                            data.frame(cov = paste(cov_names, collapse = "+"),
                                       R = 'diagonal and equal',
                                       m = 1,
                                       logLik = kemz.expl$logLik,
                                       K = kemz.expl$num.params,
                                       AICc = kemz.expl$AICc,
                                       stringsAsFactors = FALSE))
  
  assign(paste('kemz.expl', paste(cov_names, collapse = "_"), sep = "."), kemz.expl)
}

survey.model.expl.tbl2 <- survey.data.expl[order(survey.data.expl$AICc), ]
survey.model.expl.tbl2

write.csv(survey.model.expl.tbl2, "nshrimp.survey.minusecomon.combo2.loop1.csv", row.names = FALSE)


#try three combo-- depth fa, allwt fa, surftemp fa

levels.R.expl=c('diagonal and equal')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame()
survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)  #A= 'zero' 1. because needs to be fixed when estimating Z and also because data are already demeaned (zscored)
kemz.expl= MARSS(master.abun[-c(3:4),], model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(28,30, 32, 40,42)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(model='kemz.expl', R='diagonal and equal', m=1, logLik=kemz.expl$logLik, K=kemz.expl$num.params, AICc=kemz.expl$AICc, stringsAsFactors=F))
assign(paste('kemz.expl', 'combo.vars', sep='.'), kemz.expl)

survey.model.expl.tbl=survey.data.expl[order(survey.data.expl$AICc),]
survey.model.expl.tbl
write.csv(survey.model.expl.tbl, 'minusecomon.combotesttop4test.survey.model.expl.tbl.csv') 


#try all 3 combos...
combo_3vars <- combn(covariate_indices, 3, simplify = FALSE)

for (combo3 in combo_3vars) {
  
  cov_names <- names(master_exp)[combo3]
  
  survey.dfa.model.expl <- list(A = 'zero', R = levels.R.expl, m = 1)
  kemz.expl <- MARSS(master.abun[-c(3:4), ],
                     model = survey.dfa.model.expl,
                     form = 'dfa',
                     control = cntl.list,
                     covariates = t(master_exp[, combo3]))
  
  survey.data.expl <- rbind(survey.data.expl,
                            data.frame(cov = paste(cov_names, collapse = "+"),
                                       R = 'diagonal and equal',
                                       m = 1,
                                       logLik = kemz.expl$logLik,
                                       K = kemz.expl$num.params,
                                       AICc = kemz.expl$AICc,
                                       stringsAsFactors = FALSE))
  
  assign(paste('kemz.expl', paste(cov_names, collapse = "_"), sep = "."), kemz.expl)
}

survey.model.expl.tbl3 <- survey.data.expl[order(survey.data.expl$AICc), ]
survey.model.expl.tbl3

write.csv(survey.model.expl.tbl3, "nshrimp.survey.minusecomon.combo3all.loop2.csv", row.names = FALSE)




##gave a AICc of  311.76 (fa surf salin?)

levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame()


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(40)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=levels.R.expl, m=1, logLik=kemz.expl$logLik, K=kemz.expl, AICc=kemz$AICc, stringsAsFactors=F))


##AICc 349.52 (fo perfrq fall)

levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(41)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)

##AICc 351.11 ( fo perfrq sp)


levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(42)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)


##allwt fa AICc 346.35 

levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(43)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)

##allwt sp AICc 353.16

levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(34)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)


##fa bottom temp AICc 357.09


levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(35)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)

##sp bottom temp AICc 357.80


levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(48)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)


##dmr fa sst AICc 353.42


levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(49)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)

##dmr spring sst AICc 355.29


levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(44)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)


##dmr fa bt AICc 345.53



levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(45)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)


##dmr sp bt AICc 357.00

levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(37)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)


##sp suf salin AICc 342.14


levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(28)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)

##nao AICc 351.53


levels.R.expl=c('unconstrained')
cntl.list=list(minit=1000, maxit=100000, allow.degen=F, safe=TRUE)
survey.data.expl <- data.frame(R = character(), m = numeric(), logLik = numeric(), K = numeric(), AICc = numeric(), stringsAsFactors = FALSE)


survey.dfa.model.expl= list(A='zero', R=levels.R.expl, m=1)
kemz.expl= MARSS(master.abun, model=survey.dfa.model.expl, form='dfa', control=cntl.list, covariates=t(master_exp[,c(29)]))

survey.data.expl=rbind(survey.data.expl,
                       data.frame(R=R, m=m, logLik=kemz$logLik, K=kemz$num.params, AICc=kemz$AICc, stringsAsFactors=F))

assign(paste('kemz', m, R, sep='.'), kemz)

##amo AICc 353.71





#graph of explanatory zscores 

zscore2 <- names(master_exp)[grepl("^z\\.", names(master_exp)) & !grepl("panbor", names(master_exp))]

plot_data2 <- master_exp %>%
  select(year, all_of(zscore2), ) %>%
  pivot_longer(
    cols = -year,
    names_to = "variable",
    values_to = "zscore"
  )

plot_data2 <- plot_data2 %>%
  mutate(variable = recode(variable,
                           "z.depth.fa" = "Surv7 Depth (Fall)",
                           "z.depth.sp" = "Surv7 Depth (Spring)",
                           "z.surftemp.fa" = "Surv7 Surface Temp (Fall)",
                           "z.surftemp.sp" = "Surv7 Surface Temp (Spring)",
                           "z.bottemp.fa" = "Surv7 Bottom Temp (Fall)",
                           "z.bottemp.sp" = "Surv7 Bottom Temp (Spring)",
                           "z.surfsalin.fa" = "Surv7 Surface Salinity (Fall)",
                           "z.surfsalin.sp" = "Surv7 Surface Salinity (Spring)",
                           "z.botsalin.fa" = "Surv7 Bottom Salinity (Fall)",
                           "z.botsalin.sp" = "Surv7 Bottom Salinity (Spring)",
                           "z.fo.perfrq.fa" = "FO perfrq (Fall)",
                           "z.fo.perfrq.sp" = "FO perfrq (Spring)",
                           "z.allwt.relmsw_s.fa" = "AllWT relMSW (Fall)",
                           "z.allwt.relmsw_s.sp" = "AllWT relMSW (Spring)",
                           "z.dmr.env.bt.fa" = "DMR Bottom Temp (Fall)",
                           "z.dmr.env.bt.sp" = "DMR Bottom Temp (Spring)",
                           "z.dmr.env.bsal.fa" = "DMR Bottom Salinity (Fall)",
                           "z.dmr.env.bsal.sp" = "DMR Bottom Salinity (Spring)",
                           "z.dmr.env.sst.fa" = "DMR SST (Fall)",
                           "z.dmr.env.sst.sp" = "DMR SST (Spring)",
                           "z.dmr.env.ssal.fa" = "DMR Surface Salinity (Fall)",
                           "z.dmr.env.ssal.sp" = "DMR Surface Salinity (Spring)"
                           
                           
  ))

amo_plot <- master_exp %>%
  select(year, amo) %>%
  mutate(variable = "AMO (Unscaled)") %>%
  rename(zscore = amo)  

nao_plot <- master_exp %>%
  select(year, nao) %>%
  mutate(variable = "Nao winterized") %>%
  rename(zscore = nao)  

plot_data2_combined <- bind_rows(plot_data2, nao_plot, amo_plot)

ggplot(plot_data2_combined, aes(x = year, y = zscore)) +
  geom_line(na.rm = TRUE, linewidth = 1, color = "lightblue") +
  geom_point(na.rm = TRUE, color = "black", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ variable, scales = "free_y") +
  labs(
    title = "Z-Scores of Environmental Variables",
    x = "Year",
    y = "Z-Score / Value"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 10, face = "bold"))

# in master_exp, you can get rid of:
## abun.fall.abun.spring, biofall, biospring

#mapping survey areas#


#create a map of the gulf of maine

library(marmap)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ggspatial)

##GOM coordinates?
lon1 <- -72
lon2 <- -64
lat1 <- 38
lat2 <- 45

##testing a gom bathymetry map (useful?)

gom_bathy <- getNOAA.bathy(lon1, lon2, lat1, lat2, resolution = 1)
gom_df <- fortify (gom_bathy)

#m100=subset(gom_df, gom_df$z>=-100.4&gom_df$z<=-99.5)

gommap=ggplot() +
  #geom_raster(data = gom_df, aes(x = x, y = y, fill = z)) +
  #scale_fill_viridis(option = "D", name = "Depth (m)", direction = -1) +
  geom_contour(data = gom_df, aes(x = x, y = y, z = z),
               breaks = c(-100), color = "black", size = 0.3) +
  borders("world", fill = "gray80", colour = "black") +
  coord_fixed(xlim = c(-72, -64), ylim = c(38, 45), expand = FALSE) +  
  labs(title = "", x = "", y = "") +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  annotation_scale(location = "bl",  plot_unit = 'km')+
  annotate("text", x = -70.5, y = 44.5, label = "ME", size = 4, fontface = "bold") +
  annotate("text", x = -70.3, y = 42, label = "MA", size = 4, fontface = "bold") +
  annotate("text", x = -66.5, y = 45.5, label = "New Brunswick", size = 4, fontface = "bold") +
  annotate("text", x = -65, y = 44.5, label = "NS", size = 4, fontface = "bold")

#create 100m isobath from bathymetry data used in gomap (ggplot2)
gg_build=ggplot_build(gommap)
m100=data.frame(gg_build$data[1])
m1002=split(m100, m100$piece)

add_na_row <- function(m100) {
  # Create a new row of NAs matching the structure/column names of the input data frame
  na_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(m100)))
  colnames(na_row) <- colnames(m100)
  
  combined_df <- rbind(m100, na_row)
  return(combined_df)
}

dfs_with_nas <- lapply(m1002, add_na_row)

# 4. Unsplit (combine) the list of data frames into a single data frame
m1003 <- do.call(rbind, dfs_with_nas)

#csv with declon declat for 100 m isobath

#write.csv(m1003, 'c:/FWDP2012/basemaps/bathy/m1003.csv')


##gom base map

gom_map <- ggplot() +
  borders("world", xlim = c(lon1, lon2), ylim = c(lat1, lat2), fill = "gray80")+#, colour = "black") +
  #borders("state", xlim = c(lon1, lon2), ylim = c(lat2, lat2)) +
  coord_fixed(xlim = c(lon1, lon2), ylim = c(lat1, lat2), expand = FALSE) +
  labs(title = "Gulf of Maine", x = "Degrees West", y = "Degrees North") +
  theme_minimal()+

annotate("text", x = -70.5, y = 44.5, label = "ME", size = 4, fontface = "bold") +
  annotate("text", x = -70.8, y = 42, label = "MA", size = 4, fontface = "bold") +
  #annotate("text", x = -66.5, y = 45.5, label = "New Brunswick", size = 4, fontface = "bold") +
  annotate("text", x = -65, y = 44.5, label = "NS", size = 4, fontface = "bold")



gom_map

#mapping panbor sample sites

convert_to_decimal_degrees <- function(dm) {
  deg <- floor(dm / 100)
  min <- dm %% 100
  return(deg + (min / 60))
}

panbor_sub$lat_dd <- convert_to_decimal_degrees(panbor_sub$beglat)
panbor_sub$lon_dd <- -convert_to_decimal_degrees(panbor_sub$beglon)  

#new map 12/19/25
predkp=c(155, 102, 73, 101, 22, 156, 104, 197, 74, 26, 163, 75, 77, 164, 72, 27, 15, 28, 76, 108, 23)

pb2=subset(panbor, panbor$svspp%in%predkp& panbor$geoarea=='GoM')

library(maps)
library(mapdata)
#library(GISTools)
library(mapplots)
library(shape)

m200=read.csv('c:/FWDP2012/basemaps/bathy/m200.csv', header=T, sep=',', na.strings=c(NA))
m100=read.csv('c:/FWDP2012/basemaps/bathy/m1003.csv', header=T, sep=',', na.strings=c(NA))
#names(us100)[1:2]=c('declon', 'declat')

gom=read.csv('c:/FWDP2012/basemaps/gom.csv', header=T, sep=',', na.strings=c(NA))
gom

library(plotrix)


png(file=paste('c:/FWDP2020/nshrimp/data/nshrimp.survey/map.gom4.png', sep=""), width=6, height=8, units='in',res=800)
#setEPS()
#postscript('fig.1.eps')
#par(bg=NA)#,
map('worldHires', interior=F,  xlim=c(-72,-64), ylim=c(38,45.0), lwd=1, myborder=0, mar=c(2,4,2,2))
map('state', boundary=T,  col='grey', add=T, lwd=1, fill=T)
map('worldHires', interior=F, xlim=c(-72,-64), ylim=c(38, 45), lwd=1, col='grey', fill=T, add=T)
#lines(strat2$declon, strat2$declat, col='gray65', lty=1, lwd=1.5)


#map('worldHires', 'Canada', xlim=c(-78,-64), ylim=c(33,45.0), lwd=1, col='grey', add=T, fill=T)
#map('rivers', add=T, col='white')
#map('rivers', boundary=F, col='white', add=T, lwd=1, fill=T)
#map.cities(x=world.cities, country="usa", capitals=2, add=T)
#mtext('Degrees West', side=1, line=2.5)
#mtext('Degrees North', side=2, line=2.5)
#ax1=seq(-78, -66, 2)
#mtext(paste('Cruise', i, sep=" "), line=1.5, font=2)
axis(1, seq(-72, -64,2), labels=c(72, 70,68, 66,"64°W"))
axis(2, seq(37,45,2), labels=c(37, 39, 41, 43, "45°N"), las=1  )
#title(main=i)
maps::map.scale(-71,38.6, relwidth = .25, cex=0.75, ratio=F)

#strat2=strat
#strat2
#strat2[1:800, 2:3]=NA
#strat2

#lines(strat2$declon, strat2$declat, col='gray', lty=1, lwd=1)
#lines(iso100$lon, iso100$lat, col='black', lty=1, lwd=1.5)

#lines(sne$lon, sne$lat, col='black', lty=2, lwd=2)
#lines(mab2$lon, mab2$lat, col='black', lty=1, lwd=2)
#lines(mab3$lon, mab3$lat, col='black', lty=2, lwd=2)
#lines(gb$lon, gb$lat, col='black', lty=2, lwd=2)
#lines(hague$declon, hague$declat, col='black', lty=2, lwd=1)
#lines(scs$lon, scs$lat, col='black', lty=1, lwd=2)
#polygon(mab2$lon, mab2$lat, density=-100, col='#0099FF80',lty=1, border=NA)

#polygon(gb$lon, gb$lat, density=-100, col='#FF000080',lty=1,border=NA)
#polygon(sne$lon, sne$lat, density=-100, col='#0099FF80',lty=1,border=NA)
#polygon(gom$lon, gom$lat, density=-100, col='#0099FF80',lty=1,border=NA)
#polygon(scs$lon, scs$lat, density=-100, col='#0099FF80',lty=1,border=NA)
#lines(sab$lon, sab$lat, col='black', lty=1,  lwd=3)
#lines(sab$lon, sab$lat, col='black', lty=1, lwd=3)
#north.arrow(-65, 33.7, len=0.20,  cex.lab=0.75, col='black')
#text(-70,38, 'Atlantic Ocean', cex=1.25)

#recolor islands
map('state', boundary=F,  col='grey', add=T, lwd=1, fill=T)
#map('worldHires', 'Canada', xlim=c(-78,-62), ylim=c(33,45.0), lwd=1, col='grey', add=T, fill=T)
#axis(1, seq(-78, -64,2), labels=c(78,76, 74, 72,  70, 68, 66,"64°W"))
#axis(2, seq(33,45,2), labels=c(33, 35, 37, 39, 41, 43, "45°N"), las=1  )

lines(m100$x, m100$y, lwd=1, col='black')
#lines(us100$declon, us100$declat, lwd=1, col='black')
text(-70.5,44.5, 'ME',  cex=1)
text(-71.6, 43.4, 'NH', cex=1)
text(-71, 42, 'MA', cex=1)
text(-65, 44.5, 'NS',  cex=1)

points(pb2$declon*-1, pb2$declat, pch=1, col='red') 
lines(gom$lon, gom$lat, col='black', lty=2, lwd=1)
#draw.ellipse(x = -69.2, y = 43.7, a = -1.5, b = 2, angle = 45, border = "black", lwd = 2)


#arrows(-75.55, 34.5, -75.55, 35.1, length=0.1, angle=30, lwd=2)
#text(-74.35, 34.8, 'Cape Hatteras', cex=1)
#text(-65.4,44.5, 'Nova', cex=1, srt=35)
#text(-65.1, 44.3, 'Scotia', cex=1, srt=35)
#Arrowhead( -75.55, 35.1,  arr.length=0.5, arr.type='ellipse', angle=90)
#arrows( 35.2, -73, 35.2, -74, length=0.25, angle=30, code=5)




par(plt=c(.6, .9, .3, .45), new=T, xaxs='i', yaxs='i') # bg=NA, 0.6, 0.85, 0.3, 0.5
plot.window(xlim=c(-82,-50),  ylim=c(28,50))
axis(1, seq(-80, -50,10),labels=c(80,70, 60, "50°W"), padj=-.8, cex.axis=1,tck=-0.02)  #join axes!! start
axis(1, seq(-82, -50,10), labels=c(NA, NA, NA, NA), tck=0)
axis(2, seq(30,50,10), las=1, labels=c(30, 40, "50°N"), hadj=0.8, cex.axis=1,tck=-0.02)
axis(2, seq(28,50, 10), labels=c(NA, NA, NA), tck=0)
axis(3, seq(-80,-50,10), labels=c(NA,NA,NA,NA), tck=-0.02)
axis(3, seq(-82,-50,10), labels=c(NA,NA,NA,NA), tck=0)
axis(4, seq(30, 50,10), labels=c(NA, NA, NA, NA), las=1, tck=-0.02)
axis(4, seq(28, 50,10), labels=c(NA, NA, NA, NA), las=1, tck=0)
map('worldHires', interior=T,  xlim=c(-82,-50),  ylim=c(28,50), lwd=0.5, myborder=0, mar=c(4,7,0,6), add=T)
map('usa', boundary=F, add=T, col='grey', lwd=0.5, fill=T)
map('state', boundary=F, xlim=c(-100,-40), ylim=c(20,60), lwd=0.5, lty=1, add=T)
map('worldHires', 'canada', interior=T, xlim=c(-100,-40), ylim=c(20, 60), lwd=0.5, col='grey', fill=T, add=T)


#map('worldHires', interior=F,  xlim=c(-82,-50),  ylim=c(28,50), lwd=0.5, myborder=0, mar=c(4,7,0,6), add=T)
#map('worldHires', interior=F,  xlim=c(-78,-64), ylim=c(33,45.0), lwd=1, myborder=0)
#map('worldHires', interior=F, xlim=c(-82,-50), ylim=c(28, 50), lwd=0.5, col='grey', fill=T, add=T)
#plot.window(xlim=c(-82,-50),  ylim=c(28,50))
#map('usa', boundary=F, xlim=c(-82,-50),  ylim=c(28,50), col='grey', lwd=0.5, fill=T, add=T)
#map('state', boundary=F, xlim=c(-82,-50), ylim=c(28,50), col='grey', fill=T, lwd=0.5, lty=1, add=T)
#map('worldHires', 'canada', interior=T, xlim=c(-82,-40), ylim=c(28, 50), lwd=0.5, col='grey', fill=T, add=T)
axis(1, seq(-80, -50,10),labels=c(NA, NA, NA, NA), cex.axis=1,tck=-0.02)  # padj=.3,
axis(1, seq(-82, -50,10), labels=c(NA, NA, NA, NA), tck=0)
axis(2, seq(30,50,10), las=1, labels=c(NA, NA, NA), hadj=1, cex.axis=1,tck=-0.02)
axis(2, seq(28,50, 10), labels=c(NA, NA, NA), tck=0)
axis(3, seq(-80,-50,10), labels=c(NA,NA,NA,NA), tck=-0.02)
axis(3, seq(-82,-50,10), labels=c(NA,NA,NA,NA), tck=0)
axis(4, seq(30, 50,10), labels=c(NA, NA, NA, NA), las=1, tck=-0.02)
axis(4, seq(28, 50,10), labels=c(NA, NA, NA, NA), las=1, tck=0)
box=data.frame('declat'=c(37,45,45,37, 37), 'declon'=c(-72,-72,-64,-64, -72) )
lines(box$declon, box$declat, col='black', lty=1, lwd=2)
text(-58, 35.0, 'Atlantic',  cex=0.75)
text(-58, 32.7, 'Ocean', cex=0.75)



dev.off()



#mapping stratum area

stratum_area$mid_lat_dd <- convert_to_decimal_degrees(stratum_area$MIDLAT)
stratum_area$mid_lon_dd <- -convert_to_decimal_degrees(stratum_area$MIDLON) 


gom_map +
  geom_point(data = stratum_area, 
             aes(x = mid_lon_dd, y = mid_lat_dd, color = STRATUM), 
             size = 2) +
  labs(title = "stratum area map", color = "stratum number") +
  theme_minimal()+
  annotate("text", x = -70.5, y = 44.5, label = "ME", size = 4, fontface = "bold") +
  annotate("text", x = -70.3, y = 42, label = "MA", size = 4, fontface = "bold") +
  annotate("text", x = -66.5, y = 45.5, label = "New Brunswick", size = 4, fontface = "bold") +
  annotate("text", x = -65, y = 44.5, label = "NS", size = 4, fontface = "bold")

  

#plot ecomon_pandalidea_count

ecomon_p <-filter(ecomon_p, REGION == "GOM")

gom_map +
  geom_point(data = ecomon_p, aes( x = LONGITUDE, y = LATITUDE))+
  labs(title = "ecomon survey sites")+
  theme_minimal()+
  annotate("text", x = -70.5, y = 44.5, label = "ME", size = 4, fontface = "bold") +
  annotate("text", x = -70.2, y = 41.5, label = "MA", size = 4, fontface = "bold") +
  annotate("text", x = -66.5, y = 45.5, label = "New Brunswick", size = 4, fontface = "bold") +
  annotate("text", x = -65, y = 44.5, label = "NS", size = 4, fontface = "bold")



