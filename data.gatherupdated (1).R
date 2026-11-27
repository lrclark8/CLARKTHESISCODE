#Data for Lydia/NERTO 2025 northern shrimp

#### stratified mean abundance and biomass/tow of n shrimp, spring and fall time series.
library(googledrive)
library(data.table); #library(rgdal)

drive_download('survdat.RData', overwrite = T)
load('survdat.RData')


setwd("C:/Users/lydia.clark/Downloads/data")


names(survdat)=tolower(names(survdat))

survdat$cruise6=as.numeric(survdat$cruise6)
survdat$station=as.numeric(survdat$station)
survdat$stratum=as.numeric(survdat$stratum)
survdat$tow=as.numeric(survdat$tow)
survdat$svspp=as.numeric(survdat$svspp)
survdat$year=as.numeric(survdat$year)


#assign geoareas
survdat$geoarea=ifelse((survdat$stratum>=1010 & survdat$stratum <=1120)|(survdat$stratum >=3010 & survdat$stratum<=3140)|(survdat$stratum>=3450&survdat$stratum<=3550),'SNE',#|(survdat$stratum >=3470 & survdat$stratum <=3510)|(survdat$stratum ==3910 | survdat$stratum == 3990)|(survdat$stratum ==3540 | survdat$stratum==3530),'SNE',
                       ifelse((survdat$stratum >=1610 & survdat$stratum <=1760)| (survdat$stratum >=3150 & survdat$stratum<=3440),'MAB',
                              ifelse((survdat$stratum >=1130 & survdat$stratum <=1230),'GB',#| (survdat$stratum ==1250|survdat$stratum==3460 | survdat$stratum==3520 | survdat$stratum==3550 | survdat$stratum==3560|survdat$stratum== 2160 | survdat$stratum == 2170 | survdat$stratum == 2210), 'GB',
                                     ifelse((survdat$stratum == 1240)| (survdat$stratum >=1260 & survdat$stratum <=1300)| (survdat$stratum >=1360 & survdat$stratum <=1400)|(survdat$stratum >=3560 & survdat$stratum<=3900), 'GoM',
                                            ifelse((survdat$stratum >=1310 & survdat$stratum <=1340)|(survdat$stratum==1350), 'ScS',#(survdat$stratum >=1410 & survdat$stratum <=1600)|(survdat$stratum ==3920|survdat$stratum == 5470 | survdat$stratum == 5480)|(survdat$stratum >= 5410 & survdat$stratum <=5440), 'ScS',
                                                   ifelse((survdat$stratum >=7000 & survdat$stratum <=8990),'SAB','X')))))) 


#bring in stratum area from google drive
drive_download('stratum_area.csv', overwrite = T)

strarea1=read.csv('stratum_area.csv', sep=',', header=T, na.strings=c(NA)) 
names(strarea1)=tolower(names(strarea1)) 

strarea=strarea1[,c(2,4)] 

survdat2=merge(survdat, strarea, by=c('stratum'), all.x=T) 


#agg data by stratum season geoarea year remove length by row; each row is station n/tow and svspp etc.
survdat3=aggregate(survdat2$numlen, by=list('svspp'=survdat2$svspp, 'cruise6'=survdat2$cruise6, 'station'=survdat2$station, 'stratum'=survdat2$stratum, 'geoarea'=survdat2$geoarea, 'season'=survdat2$season, 'stratum_area'=survdat2$stratum_area, 'year'=survdat2$year,'depth'=survdat2$depth, 'surftemp'=survdat2$surftemp, 'bottemp'=survdat2$bottemp, 'surfsalin'=survdat2$surfsalin, 'botsalin'=survdat2$botsalin, 'abundance'=survdat2$abundance, 'biomass'=survdat2$biomass ), length )

#remove 0 abundances ~600 AND LIMIT TO FALL AND SVSPP OF INTEREST
SVSPP2=306
survdat4=subset(survdat3, survdat3$abundance>0&(survdat3$season=='FALL'|survdat3$season=='SPRING')&survdat3$svspp%in%SVSPP2&survdat3$geoarea=='GoM')


survdat4a=data.table(survdat4)
#find out how many data pts per block here BY YEAR!
towdata=aggregate(survdat4a$abundance, by=list('svspp'=survdat4a$svspp, 'geoarea'=survdat4a$geoarea, 'season'=survdat4a$season, 'year'=survdat4a$year), length)
towdata2=towdata[order(towdata$svspp, towdata$geoarea, towdata$season, towdata$year),]
names(towdata2)[ncol(towdata2)]='ntows'

head(survdat4)
head(survdat4a)

#weighted average n/tow and bio/tow


survdat5=survdat4a[,lapply(.SD, weighted.mean, stratum_area, na.rm=TRUE), by=list(svspp, geoarea, season, year)] 
library(Hmisc)
survdat5var=survdat4a[, lapply(.SD, wtd.var, stratum_area, na.rm=TRUE), by=list(svspp, geoarea, season, year)]

#create df of unique strata, strata area and sum per svspp, geoarea, year
strat=aggregate(survdat4a$abundance, by=list('svspp'=survdat4a$svspp, 'stratum'=survdat4a$stratum, 'geoarea'=survdat4a$geoarea, 'season'=survdat4a$season, 'year'=survdat4a$year, 'stratum_area'=survdat4a$stratum_area),sum)
strat2=aggregate(strat$stratum_area, by=list('svspp'=strat$svspp, 'geoarea'=strat$geoarea,'season'=strat$season, 'year'=strat$year),sum)
names(strat2)[max(ncol(strat2))]='sweptarea'

#need to remove unnecessary columns
survdat5[,':='(cruise6=NULL,stratum=NULL, x=NULL,station=NULL,stratum_area=NULL)]

survdat5var[,':='(cruise6=NULL,stratum=NULL, x=NULL,station=NULL,stratum_area=NULL)]
names(survdat5var)[5:ncol(survdat5var)]=c('depth.var', 'surftemp.var', 'bottemp.var', 'surfsalin.var', 'botsalin.var', 'abundance.var', 'biomass.var')

survdat5a=merge(merge(merge(survdat5, survdat5var,by=c('svspp', 'geoarea','season', 'year'),all.x=T), towdata2, by=c('svspp', 'geoarea' ,'season', 'year'), all.x=T), strat2, by=c('svspp', 'geoarea','season', 'year'), all.x=T)

#final 
survdat6=survdat5a[order(survdat5a$svspp, survdat5a$geoarea, survdat5a$season, survdat5a$year),] 
#spring and fall;  need to add NA for year 2020 and possibly others.


#survdat6$abunpop=(survdat6$abundance*survdat6$sweptarea)/0.0112
#survdat6$abunpop.sd=(survdat6$abundance.sd*survdat6$sweptarea)/0.0112

write.csv(survdat6, 'survdat.shrimp.csv') 



#%FO and % mass northern shrimp 
#first determine predators

#select prey occurrences (aka preyinpred.sas)

#setwd("c:/FWDP2020/nshrimp/data")

#read in allfh
drive_download('lc.allfh.RData', overwrite = T)
load('lc.allfh.RData')

#allfh

#column names to lowercase 
names(allfh)=tolower(names(allfh))

#preyinpred  Limit data here and select predators (svspp). 

allfh2=subset(allfh, allfh$pynam!='BLOWN'& allfh$pynam != 'PRESERVED' & allfh$pynam != ' ' & allfh$purcode == 10&allfh$year!=''&(allfh$season=='SPRING'|allfh$season=='FALL')&(allfh$geoarea=='GoM'))

prey=subset(allfh2, allfh2$pyabbr =='PANBOR' ) #select PREY 

#LIMIT COLUMNS and reorder
keep=c('pynam',	'cruise',	'station',	'surftemp',	'bottemp',	'beglat',	'beglon',	'declat',	'declon',	'month',	'day',	'year',	'season',	'geoarea',	'pyamtv',	'perpyv',	'pdcomnam',	'svspp',	'pdlen')

prey2=prey[keep]

prey3=prey2[order(prey2$svspp, prey2$year),]

write.csv(prey3, 'panbor.csv')

shrimppd=unique(prey3$svspp)



#%FO  grab from dietdash
drive_download('allfreqp2025-06-05.csv', overwrite = T)  #this automatically goes to c:Users/your.name/Documents
fo=read.csv('allfreqp2025-06-05.csv', h=T, sep=',', na.strings = c(NA))

fo2=subset(fo, fo$prey=='pandalid'&fo$geoarea=='GoM'&(fo$season=='SPRING'|fo$season=='FALL')&fo$nstom>=40)

#weighted mean fo across predators (weighted by n stomachs).
#first simple mean
fo.m=aggregate(list('m.pctfo'=fo2$perfrq), by=list('season'=fo2$season, 'year'=fo2$year), mean, na.rm=T)
fo.sd=aggregate(list('sd.pctfo'=fo2$perfrq), by=list('season'=fo2$season, 'year'=fo2$year), sd, na.rm=T)
fo.n=aggregate(list('n'=fo2$perfrq), by=list('season'=fo2$season, 'year'=fo2$year), length)

fo.m2=merge(merge(fo.m, fo.sd, by=c('season', 'year'), all.x=T), fo.n, by=c('season', 'year'), all.x=T)
fo.m2$ci=2*(fo.m2$sd.pctfo/sqrt(fo.m2$n))

#wt mean--reference code starts here 
library(data.table)
fo3=data.table(fo2[,c(3:8)]) #remove pred and prey columns to allow for wt.mean

fo.wtm=fo3[,lapply(.SD, weighted.mean, nstom, na.rm=TRUE), by=list(geoarea, season, year)] 
#rm cols
fo.wtm[,':='(freq=NULL,nstom=NULL)]



library(Hmisc)
fo.var=fo3[, lapply(.SD, wtd.var, nstom, na.rm=TRUE), by=list(geoarea, season, year)]
#rm cols
fo.var[,':='(freq=NULL, nstom=NULL)]

names(fo.var)[ncol(fo.var)]='perfrq.var'



fo.wtm2=merge(merge(fo.wtm, fo.var, by=c('geoarea', 'season', 'year'), all.x=T), fo.n, by=c('season', 'year'), all.x=T)
fo.wtm2$ci=2*(sqrt(fo.wtm2$perfrq.var/fo.wtm2$n))
#need to add 2020=NA for spring and fall

write.csv(fo.wtm2, file = "/Users/lydia.clark/Downloads/data/fo.wtm2.csv", row.names = FALSE)
#######

#%mass

#calc weighted diet composition by weight (aka allwt.sas) USE 64-BIT R

#setwd("c:/FWDP2020/nshrimp/data")

#read in allfh
#allfh=read.csv("all_raw.csv", header=T, sep=',', na.strings=c(NA))
fo.wtm2 <- read.csv("/Users/lydia.clark/Downloads/data/fo.wtm2.csv")
#allfh

names(allfh)=tolower(names(allfh))


#CHOOSE byvars and byvar names (e.g. geoarea, sizecat, season, year, etc.).  Order according to output and fill as necessary or leave 'x'. Allows for max of 4.
#1
allfh$byvar1=allfh$season
name_byvar1=c('season')

#2
allfh$byvar2=allfh$year
name_byvar2=c('year')

#3
allfh$byvar3=allfh$geoarea
name_byvar3=c('geoarea')

#4
allfh$byvar4='x'
name_byvar4=c('x')


#Define pyvar; pyamt covers all years; global var
#pyvar=('pyamt') #for now use pyamt throughout pr pyvar

########Define taxnam; analcat, collcat, gencat, or create your own with ifelse etc.

allfh$pycat=ifelse(allfh$collcat=='PANFAM', 'PANFAM', 'OTHER')

allfh$taxnam=allfh$pycat


########Define svspp;

pred=allfh$svspp%in%shrimppd#|allfh$svspp==23
pred

#Define byvars;
#byvars=c("year")#(season, geoarea, year...)
#byvars

# DON'T CHANGE. This allows for the appropriate number of columns when setting up data.
nbyvar=4

allfh1=subset(allfh, pred & (allfh$pynam!='BLOWN'& allfh$pynam != 'PRESERVED' & allfh$pynam != ' ' & allfh$purcode == 10&(allfh$season=='SPRING'|allfh$season=='FALL')&allfh$geoarea=='GoM')) #add other data limitations here


allfh1$tax=allfh1$taxnam
allfh1$perpy=allfh1$perpyw
allfh1$pyamt=allfh1$pyamtw

#pyamt is the variable for prey weight

keep1=c('cruise6', 'stratum', 'station', 'byvar1','byvar2', 'byvar3', 'byvar4', 'svspp', 'pdid', 'pdsex', 'pdlen', 'tax', 'pyamt', 'catnum', 'numlen', 'tot_catnum_stratum', 'tot_catwgt_stratum', 'tot_tows_spp_stratum', 'stratum_area')  #add byvars before svspp
allfh2=allfh1[keep1]
#head(allfh2,10)

allwt1=allfh2[order(allfh2$cruise6, allfh2$station, allfh2$byvar1, allfh2$byvar2, allfh2$byvar3, allfh2$byvar4, allfh2$svspp, allfh2$pdid, allfh2$pdsex, allfh2$pdlen),] #add byvars before svspp
#allwt1

new_numlen1=aggregate(allwt1$numlen, list('cruise6'=allwt1$cruise6, 'station'=allwt1$station, 'byvar1'=allwt1$byvar1, 'byvar2'=allwt1$byvar2, 'byvar3'=allwt1$byvar3, 'byvar4'=allwt1$byvar4, 'svspp'=allwt1$svspp, 'pdid'=allwt1$pdid, 'pdsex'=allwt1$pdsex, 'pdlen'=allwt1$pdlen), length) #add byvars before svspp
#head(new_numlen1,10)
colnames(new_numlen1)[7+nbyvar]='xxxx' #change column number if byvars added

new_numlen1a=new_numlen1[order(new_numlen1$cruise6, new_numlen1$station, new_numlen1$byvar1, new_numlen1$byvar2, new_numlen1$byvar3, new_numlen1$byvar4, new_numlen1$svspp, new_numlen1$pdsex, new_numlen1$pdlen),] #add byvars before svspp
#new_numlen1a

#create numlen2 (ie count of true numbers at pdlen); eg in some cases (2) 40cm svspp individuals sampled per station
new_numlen2=aggregate(new_numlen1a$pdlen, list('cruise6'=new_numlen1a$cruise6, 'station'=new_numlen1a$station, 'byvar1'=new_numlen1a$byvar1, 'byvar2'=new_numlen1a$byvar2, 'byvar3'=new_numlen1a$byvar3, 'byvar4'=new_numlen1a$byvar4, 'svspp'=new_numlen1a$svspp, 'pdsex'=new_numlen1a$pdsex, 'pdlen'=new_numlen1a$pdlen), length) #add byvars before svspp
#head(new_numlen2,10)
colnames(new_numlen2)[6+nbyvar]='numlen2' #change column number if byvars added

#write.csv(new_numlen2,'new_numlen2.csv')

allfh_prep=merge(allwt1, new_numlen2, by=c('cruise6', 'station','byvar1', 'byvar2', 'byvar3', 'byvar4',  'svspp', 'pdsex', 'pdlen'), all.x=T) #add byvars before svspp
#head(allfh_prep,10)

allfh=allfh_prep
allfh$numlen_fin=ifelse(allfh$numlen< allfh$numlen2, allfh$numlen2, ifelse(is.na(allfh$numlen), allfh$numlen, allfh$numlen))
#head(allfh,10)



keep2=c('cruise6', 'stratum', 'station', 'byvar1', 'byvar2', 'byvar3', 'byvar4', 'svspp', 'pdid', 'pdsex', 'pdlen', 'tax', 'pyamt', 'catnum', 'tot_catnum_stratum', 'tot_catwgt_stratum', 'tot_tows_spp_stratum', 'stratum_area', 'numlen', 'numlen_fin')  #add byvars before svspp
#edit
allfh=allfh[keep2]

#numlen_fin created but only for non spiny dog and smoothdog; create rnumlen_fin: use instead of numlen; numlen_fin includes null numlen
#fix numlens for 15 and 13
allfh=allfh[order(allfh$cruise6, allfh$station, allfh$byvar1, allfh$byvar2, allfh$byvar3, allfh$byvar4, allfh$svspp, allfh$pdsex, allfh$pdlen, allfh$numlen_fin),] #add byvars before svspp
#head(allfh,10)


allfh_numlen_fin=aggregate(allfh$numlen, list('cruise6'=allfh$cruise6, 'station'=allfh$station, 'byvar1'=allfh$byvar1, 'byvar2'=allfh$byvar2,'byvar3'=allfh$byvar3,'byvar4'=allfh$byvar4, 'svspp'=allfh$svspp, 'pdsex'=allfh$pdsex, 'pdlen'=allfh$pdlen, 'numlen_fin'=allfh$numlen_fin), sum, na.rm=T) #add byvars before svspp
#head(allfh_numlen_fin,10)
colnames(allfh_numlen_fin)[7+nbyvar]='dumby_var' #change column number if byvars added

allfh_numlen_fin=allfh_numlen_fin[order(allfh_numlen_fin$cruise6, allfh_numlen_fin$station, allfh_numlen_fin$byvar1, allfh_numlen_fin$byvar2, allfh_numlen_fin$byvar3, allfh_numlen_fin$byvar4, allfh_numlen_fin$svspp, allfh_numlen_fin$pdlen),] #add byvars before svspp
#head(allfh_numlen_fin,10)

allfh_rnumlen_fin=aggregate(allfh_numlen_fin$numlen_fin, list('cruise6'=allfh_numlen_fin$cruise6, 'station'=allfh_numlen_fin$station, 'byvar1'=allfh_numlen_fin$byvar1, 'byvar2'=allfh_numlen_fin$byvar2, 'byvar3'=allfh_numlen_fin$byvar3, 'byvar4'=allfh_numlen_fin$byvar4, 'svspp'=allfh_numlen_fin$svspp, 'pdlen'=allfh_numlen_fin$pdlen),sum, na.rm=T)  #add byvars before svspp
#head(allfh_rnumlen_fin,10)
colnames(allfh_rnumlen_fin)[5+nbyvar]='rnumlen_fin' #change column number if byvars added


allfh=allfh[order(allfh$cruise6, allfh$station, allfh$byvar1, allfh$byvar2, allfh$byvar3, allfh$byvar4, allfh$svspp, allfh$pdlen),] #add byvars before svspp

allfh=merge(allfh, allfh_rnumlen_fin, by=c('cruise6', 'station', 'byvar1', 'byvar2', 'byvar3', 'byvar4', 'svspp', 'pdlen'), all.x=T) #add byvars before svspp
#allfh_rnumlen_fin missing 2014

allfh$rnumlen_fin=ifelse(allfh$svspp!=13 & allfh$svspp!=15, allfh$numlen_fin, allfh$rnumlen_fin)
#head(allfh,10)

#numlen_fin=NA??

#end fixes for numlen 13 and 15

#sort and sum stomach data (volume) over prey taxa - (e.g. over analcat)
allfh=allfh[order(allfh$cruise6, allfh$station, allfh$byvar1, allfh$byvar2, allfh$byvar3, allfh$byvar4, allfh$svspp, allfh$pdsex, allfh$pdid, allfh$pdlen, allfh$tax),] #add byvars before svspp

allfh=allfh[order(allfh$cruise6, allfh$station, allfh$byvar1, allfh$byvar2, allfh$byvar3, allfh$byvar4, allfh$svspp, allfh$pdsex, allfh$catnum),] #add byvars before svspp

allfh_catnum=aggregate(allfh$numlen, list('cruise6'=allfh$cruise6, 'station'=allfh$station, 'byvar1'=allfh$byvar1, 'byvar2'=allfh$byvar2, 'byvar3'=allfh$byvar3, 'byvar4'=allfh$byvar4, 'svspp'=allfh$svspp, 'pdsex'=allfh$pdsex, 'catnum'=allfh$catnum), sum, na.rm=T) #add byvars before svspp
#head(allfh_catnum,10)
colnames(allfh_catnum)[6+nbyvar]='dumby_var' #change column number if byvars added

allfh_catnum=allfh_catnum[order(allfh_catnum$cruise6, allfh_catnum$station, allfh_catnum$byvar1, allfh_catnum$byvar2, allfh_catnum$byvar3, allfh_catnum$byvar4, allfh_catnum$svspp),] #add byvars before svspp

allfh_rcatnum=aggregate(allfh_catnum$catnum, list('cruise6'=allfh_catnum$cruise6, 'station'=allfh_catnum$station, 'byvar1'=allfh_catnum$byvar1, 'byvar2'=allfh_catnum$byvar2, 'byvar3'=allfh_catnum$byvar3, 'byvar4'=allfh_catnum$byvar4, 'svspp'=allfh_catnum$svspp), sum,na.rm=T) #add byvars before svspp
#head(allfh_rcatnum,10)
colnames(allfh_rcatnum)[4+nbyvar]='rcatnum' #change column number if byvars added

allfh=merge(allfh, allfh_rcatnum, by=c('cruise6', 'station', 'byvar1', 'byvar2', 'byvar3', 'byvar4', 'svspp'), all.x=T) #add byvars before svspp
#head(allfh,10)

#fix rcatnum for svspp ne 015 and ne 013
allfh$rcatnum=ifelse(allfh$svspp!=13 & allfh$svspp!=15, allfh$catnum, allfh$rcatnum)

#015 013 tot_catnum_stratum tot_catwgt_stratum fix, catch and tow numbers by stratum separated by sex
allfh=allfh[order(allfh$cruise6, allfh$stratum, allfh$byvar1, allfh$byvar2, allfh$byvar3, allfh$byvar4, allfh$svspp, allfh$pdsex, allfh$tot_tows_spp_stratum, allfh$tot_catnum_stratum, allfh$tot_catwgt_stratum),]  #add byvars before svspp


allfh_catstratum=aggregate(allfh$rcatnum, list('cruise6'=allfh$cruise6, 'stratum'=allfh$stratum, 'byvar1'=allfh$byvar1, 'byvar2'=allfh$byvar2, 'byvar3'=allfh$byvar3, 'byvar4'=allfh$byvar4, 'svspp'=allfh$svspp, 'pdsex'=allfh$pdsex, 'tot_tows_spp_stratum'=allfh$tot_tows_spp_stratum, 'tot_catnum_stratum'=allfh$tot_catnum_stratum, 'tot_catwgt_stratum'=allfh$tot_catwgt_stratum),sum,na.rm=T) #add byvars before svspp
#head(allfh_catstratum,10)
colnames(allfh_catstratum)[8+nbyvar]='dum_var'

allfh_catstratum=allfh_catstratum[order(allfh_catstratum$cruise6, allfh_catstratum$stratum, allfh_catstratum$byvar1, allfh_catstratum$byvar2, allfh_catstratum$byvar3, allfh_catstratum$byvar4, allfh_catstratum$svspp),] #add byvars before svspp

#keeps tot_tows_spp_stratum ala id statement in proc means
max_tot_tows_spp_stratum=aggregate(allfh_catstratum$tot_tows_spp_stratum, list('cruise6'=allfh_catstratum$cruise6, 'stratum'=allfh_catstratum$stratum, 'byvar1'=allfh_catstratum$byvar1, 'byvar2'=allfh_catstratum$byvar2, 'byvar3'=allfh_catstratum$byvar3,  'byvar4'=allfh_catstratum$byvar4, 'svspp'=allfh_catstratum$svspp),max) #add byvars before svspp
#head(max_tot_tows_spp_stratum,10)
colnames(max_tot_tows_spp_stratum)[4+nbyvar]='tot_tows_spp_stratum'

allfh_rcatstratum=aggregate(cbind(allfh_catstratum$tot_catnum_stratum, allfh_catstratum$tot_catwgt_stratum), list('cruise6'=allfh_catstratum$cruise6, 'stratum'=allfh_catstratum$stratum, 'byvar1'=allfh_catstratum$byvar1, 'byvar2'=allfh_catstratum$byvar2, 'byvar3'=allfh_catstratum$byvar3, 'byvar4'=allfh_catstratum$byvar4, 'svspp'=allfh_catstratum$svspp),sum,na.rm=T) #add byvars before svspp
#head(allfh_rcatstratum,10)
colnames(allfh_rcatstratum)[4+nbyvar]='rtot_catnum_stratum' #change column numbers if adding byvars
colnames(allfh_rcatstratum)[5+nbyvar]='rtot_catwgt_stratum' #change column numbers if adding byvars


allfh_rcatstratum=merge(allfh_rcatstratum, max_tot_tows_spp_stratum, by=c('cruise6', 'stratum', 'byvar1', 'byvar2', 'byvar3', 'byvar4', 'svspp'), all.x=T) #add byvars before svspp

allfh_rcatstratum=allfh_rcatstratum[order(allfh_rcatstratum$cruise6, allfh_rcatstratum$stratum, allfh_rcatstratum$byvar1, allfh_rcatstratum$byvar2, allfh_rcatstratum$byvar3, allfh_rcatstratum$byvar4, allfh_rcatstratum$svspp),] #add byvars before svspp
allfh_rcatstratum$rtot_tows_spp_stratum=allfh_rcatstratum$tot_tows_spp_stratum


allfh_rcatstratum2=subset(allfh_rcatstratum, select=-c(tot_tows_spp_stratum))

#merge in final rcatstratum numbers
allfh=merge(allfh, allfh_rcatstratum2, by=c('cruise6', 'stratum', 'byvar1', 'byvar2', 'byvar3', 'byvar4', 'svspp'), all.x=T) #add byvars before svspp

#fix rcatstratum for svspp ne 015 and ne 013, also created rtot_tows_spp_stratum need to keep until for diet calc
allfh$rtot_catnum_stratum=ifelse(allfh$svspp!=13&allfh$svspp!=15, allfh$tot_catnum_stratum, allfh$rtot_catnum_stratum)
allfh$rtot_catwgt_stratum=ifelse(allfh$svspp!=13&allfh$svspp!=15, allfh$tot_catwgt_stratum, allfh$rtot_catwgt_stratum)
allfh$rtot_tows_spp_stratum=ifelse(allfh$svspp!=13&allfh$svspp!=15, allfh$tot_tows_spp_stratum, allfh$rtot_tows_spp_stratum)

#end of fixes for 13 and 15 catnums and catwgts and numbers per stratum

allfh=allfh[order(allfh$cruise6, allfh$station, allfh$byvar1, allfh$byvar2, allfh$byvar3, allfh$byvar4, allfh$svspp, allfh$pdid, allfh$pdlen, allfh$tax),]

fish=aggregate(allfh$pyamt, list('cruise6'=allfh$cruise6, 'station'=allfh$station, 'byvar1'=allfh$byvar1, 'byvar2'=allfh$byvar2, 'byvar3'=allfh$byvar3, 'byvar4'=allfh$byvar4, 'svspp'=allfh$svspp, 'pdid'=allfh$pdid, 'pdlen'=allfh$pdlen, 'tax'=allfh$tax),sum,na.rm=T) #add byvars before svspp
#head(fish,10)
colnames(fish)[7+nbyvar]='pysum'

#count the number of prey taxa

fish=fish[order(fish$tax),]

taxcount=data.frame('tax'=unique(fish$tax))
taxcount

numpy=data.frame(cbind(taxcount,'num_prey'=nrow(taxcount)))
numpy

keep3=c('cruise6', 'station', 'pdid', 'pysum', 'byvar1',  'byvar2', 'byvar3', 'byvar4', 'svspp') #add byvars before svspp
pred= fish[keep3]
pred=pred[order(pred$cruise6, pred$station, pred$byvar1, pred$byvar2, pred$byvar3, pred$byvar4, pred$svspp, pred$pdid),] # add byvars before svspp
pred
one=aggregate(pred$pysum, list('cruise6'=pred$cruise6, 'station'=pred$station, 'byvar1'=pred$byvar1, 'byvar2'=pred$byvar2, 'byvar3'=pred$byvar3, 'byvar4'=pred$byvar4, 'svspp'=pred$svspp, 'pdid'=pred$pdid),head,1) #add byvars before svspp
#head(one,10)
colnames(one)[5+nbyvar]='pysum' #change column number if byvars added

one=one[order(one$cruise6, one$station, one$byvar1, one$byvar2, one$byvar3, one$byvar4, one$svspp),] #add byvars before svspp
#head(one,15)


new=aggregate(one$pysum, list('cruise6'=one$cruise6, 'station'=one$station, 'byvar1'=one$byvar1, 'byvar2'=one$byvar2, 'byvar3'=one$byvar3, 'byvar4'=one$byvar4, 'svspp'=one$svspp),length)  #add byvars before svspp
#head(new,10)
colnames(new)[4+nbyvar]='nostom' #change column number if byvars added

new=new[order(new$cruise6, new$station, new$byvar1, new$byvar2, new$byvar3, new$byvar4, new$svspp),] # add byvars before svspp

fish=fish[order(fish$cruise6, fish$station, fish$byvar1, fish$byvar2, fish$byvar3, fish$byvar4, fish$svspp),] # add byvars before svspp

fish=merge(fish, new, by=c('cruise6', 'station', 'byvar1', 'byvar2', 'byvar3', 'byvar4', 'svspp')) #add byvars before svspp


fish=fish[order(fish$cruise6, fish$station, fish$byvar1, fish$byvar2, fish$byvar3, fish$byvar4, fish$svspp),] # add byvars before svspp


#Dataset fish now has the foodhabits data, plus the number
#of stomachs per cruise6 station combination and the number of prey taxa
#in data set

#Create prey codes to set up for the macro loops

pylist=cbind(taxcount,'pycode'=paste('prey',1:nrow(taxcount),sep=''))
pylist


fish=merge(pylist, fish, by=c('tax'), all.x=T)
#head(fish,10)

fish=fish[order(fish$cruise6, fish$station, fish$byvar1, fish$byvar2, fish$byvar3, fish$byvar4, fish$svspp, fish$pdid, fish$pdlen),] #add byvars before svspp


#Merge in the catch at length data
length1=allfh[order(allfh$cruise6, allfh$station, allfh$byvar1, allfh$byvar2, allfh$byvar3, allfh$byvar4, allfh$svspp, allfh$pdid, allfh$pdlen),] #add byvars before svspp

length2=aggregate(length1$stratum_area, list('cruise6'=length1$cruise6, 'station'=length1$station, 'byvar1'=length1$byvar1, 'byvar2'=length1$byvar2, 'byvar3'=length1$byvar3, 'byvar4'=length1$byvar4, 'svspp'=length1$svspp, 'pdid'=length1$pdid, 'pdlen'=length1$pdlen),head,1) #add byvars before svspp
#head(length2,10)
colnames(length2)[6+nbyvar]='stratum_area' #change column number if adding byvars

length2=merge(length2, length1, by=c('cruise6', 'station', 'byvar1', 'byvar2', 'byvar3', 'byvar4', 'svspp', 'pdid', 'pdlen', 'stratum_area'), all.x=T) #byvars before svspp
#head(length2,10)

keep4=c('cruise6', 'station', 'stratum', 'byvar1', 'byvar2', 'byvar3', 'byvar4','svspp',  'pdid', 'pdlen', 'catnum', 'rcatnum', 'tot_catnum_stratum', 'rtot_catnum_stratum', 'tot_catwgt_stratum', 'rtot_catwgt_stratum', 'tot_tows_spp_stratum', 'rtot_tows_spp_stratum', 'stratum_area', 'numlen', 'numlen_fin', 'rnumlen_fin') #byvars before svspp
keep4

length2=length2[keep4]
#head(length2,10)

length2=length2[order(length2$cruise6, length2$station, length2$byvar1, length2$byvar2, length2$byvar3, length2$byvar4, length2$svspp, length2$pdid, length2$pdlen),] #add byvars before svspp

fish=subset(fish, select=-c(tax))

fish=merge(fish, length2, by=c('cruise6','station',  'byvar1', 'byvar2', 'byvar3', 'byvar4', 'svspp', 'pdid', 'pdlen'), all.x=T) #byvars before svspp
#head(fish,10)

fish=fish[order(fish$cruise6, fish$stratum, fish$station, fish$byvar1, fish$byvar2, fish$byvar3, fish$byvar4, fish$svspp, fish$pdid,    fish$pdlen, fish$catnum, fish$rcatnum, fish$tot_catnum_stratum, fish$rtot_catnum_stratum, fish$tot_catwgt_stratum, fish$rtot_catwgt_stratum, fish$tot_tows_spp_stratum, fish$rtot_tows_spp_stratum, fish$stratum_area, fish$numlen, fish$numlen_fin, fish$rnumlen_fin, fish$nostom),]
#head(fish,10)

#transpose cast pycode data and merge table to fish
library(reshape2)
library(data.table)

fish=data.table(fish)

tran1=data.frame(dcast.data.table(fish,cruise6+stratum+station+byvar1+byvar2+byvar3+byvar4+svspp+pdid+pdlen+catnum+rcatnum+tot_catnum_stratum+rtot_catnum_stratum+tot_catwgt_stratum+rtot_catwgt_stratum+tot_tows_spp_stratum+rtot_tows_spp_stratum+stratum_area+numlen+rnumlen_fin+nostom~pycode, fun.aggregate=median, value.var='pysum') )   #add byvars before svspp

#head(tran1,10)
#set prey NA to zero in fish

#preycol=dput(as.character(pylist$pycode))
preycol=dput(as.character(pylist$pycode))

tran1[,preycol][is.na(tran1[,preycol])]=0




tran2=data.table(tran1)

new1=tran2[,lapply(.SD,mean, na.rm=T),by=list(cruise6, stratum, station, byvar1, byvar2, byvar3, byvar4, svspp, pdid, pdlen, catnum, rcatnum, tot_catnum_stratum, rtot_catnum_stratum, tot_catwgt_stratum,  #add byvars before svspp
                                              rtot_catwgt_stratum, tot_tows_spp_stratum, rtot_tows_spp_stratum, stratum_area, numlen, rnumlen_fin, nostom)]

#NA for many rnumlen_fin records??

#remove nas  here
new1rm=new1[!is.na(new1$rnumlen_fin)]   
new1rm
new11=new1rm[,lapply(.SD,  weighted.mean,rnumlen_fin,na.rm=TRUE), by=list(cruise6,stratum, station, byvar1, byvar2, byvar3, byvar4, svspp,rcatnum, rtot_catnum_stratum, rtot_tows_spp_stratum, stratum_area)]  #add byvars before svspp

#need to remove unnecessary columns
new11[,':='(pdid=NULL, pdlen=NULL ,catnum=NULL, tot_catnum_stratum=NULL, tot_catwgt_stratum=NULL, rtot_catwgt_stratum=NULL, tot_tows_spp_stratum=NULL, numlen=NULL, rnumlen_fin=NULL, nostom=NULL)]


#new11 already sorted by cruise6, stratum, station byvars svspp
#and py amts already means by station so no need for proc means here

#addition of new weighted estimator

new2=new11[,lapply(.SD,mean, na.rm=T), by=list(cruise6, stratum, station, byvar1, byvar2, byvar3, byvar4, svspp, rcatnum, rtot_catnum_stratum, rtot_tows_spp_stratum, stratum_area)]


musw2=new11[,lapply(.SD,function(x) x*rcatnum),by=list(cruise6, stratum, station, byvar1, byvar2, byvar3, byvar4, svspp,rcatnum, rtot_catnum_stratum, rtot_tows_spp_stratum, stratum_area)] #add byvars before svspp


new2s=musw2[,lapply(.SD,sum ), by=list(cruise6, stratum, byvar1, byvar2, byvar3, byvar4, svspp,rtot_catnum_stratum, rtot_tows_spp_stratum, stratum_area)] #add byvars before svspp


musw_strat=new2s[,lapply(.SD, function(x) x/rtot_tows_spp_stratum),by=list(cruise6, stratum, byvar1, byvar2, byvar3, byvar4, svspp,rcatnum, rtot_catnum_stratum, rtot_tows_spp_stratum, stratum_area)] #add byvars before svspp
musw_strat$munfish_strat=musw_strat$rtot_catnum_stratum/musw_strat$rtot_tows_spp_stratum
#retain munfish_strat
munfish_stratdat=musw_strat[, c('cruise6', 'stratum',  'byvar1', 'byvar2', 'byvar3', 'byvar4', 'svspp', 'munfish_strat'), with=F] #add byvars before 


musw_strat=musw_strat[order(musw_strat$svspp, musw_strat$byvar1, musw_strat$byvar2,musw_strat$byvar3,musw_strat$byvar4),] #byvars AFTER svspp


#weight by stratum area
#first remove NA stratum areas
musw_stratum=musw_strat[!is.na(musw_strat$stratum_area)]

new3s=musw_stratum[,lapply(.SD, weighted.mean,stratum_area, na.rm=T), by=list(svspp,byvar1,byvar2,byvar3,byvar4)] #add byvars AFTER svspp   #msw_strat

new3s[,':='(cruise6=NULL, stratum=NULL ,station=NULL, rcatnum=NULL, rtot_catnum_stratum=NULL, rtot_tows_spp_stratum=NULL, stratum_area=NULL)]

#new3s$num_stra=NROW(musw_stratum) #fix for multiple svspp
num_stra=musw_stratum[, .N,by=list(svspp,byvar1,byvar2,byvar3,byvar4)] #add byvars AFTER svspp
names(num_stra)=gsub('N','num_stra',names(num_stra))

new3s=merge(new3s, num_stra, by=c('svspp', 'byvar1', 'byvar2', 'byvar3', 'byvar4'), all.x=T)  #add byvars AFTER svspp




meansw_s=new3s[,lapply(.SD,function(x) x/munfish_strat),by=list(svspp,byvar1,byvar2,byvar3,byvar4)] #add byvars AFTER svspp 

meansw_s[,':='(munfish_strat=NULL, num_stra=NULL)]

#CHANGE labels to meansw_s by prey number

three=new2[order(new2$svspp,new2$byvar1,new2$byvar2,new2$byvar3,new2$byvar4),] #add byvars AFTER svspp

meansw=three[,lapply(.SD,weighted.mean, rcatnum, na.rm=T), by=list(svspp,byvar1,byvar2,byvar3,byvar4)] #add byvars AFTER svspp
meansw[,':='(cruise6=NULL, stratum=NULL ,station=NULL, rcatnum=NULL, rtot_catnum_stratum=NULL, rtot_tows_spp_stratum=NULL, stratum_area=NULL)]


#meansw$num_tows=NROW(three) fix for multiple svspp
num_tows=three[, .N,by=list(svspp,byvar1,byvar2,byvar3,byvar4)] #add byvars AFTER svspp
names(num_tows)=gsub('N','num_tows',names(num_tows))

meansw=merge(meansw, num_tows, by=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4'), all.x=T)  #add byvars AFTER svspp

#melt three to start prod
master_musw=melt(three, id.vars=c("cruise6", "stratum", "station", "svspp",  "byvar1", "byvar2", "byvar3", "byvar4", "rcatnum", "rtot_catnum_stratum", "rtot_tows_spp_stratum", "stratum_area"),
                 measure.vars=preycol, variable.name='pycode', value.name='musw')

master_tmsw_strat=melt(new2s, id.vars=c("cruise6", "stratum", "svspp", "byvar1", "byvar2", "byvar3", "byvar4"),
                       measure.vars=preycol, variable.name='pycode', value.name='tmsw_strat')

master_musw2=melt(musw2, id.vars=c("cruise6", "stratum", "station", "svspp", "byvar1", "byvar2", "byvar3", "byvar4"),
                  measure.vars=preycol, variable.name='pycode', value.name='musw2')

master_musw_strat=melt(musw_strat, id.vars=c("cruise6", "stratum",  "svspp", "byvar1", "byvar2", "byvar3", "byvar4" ),
                       measure.vars=preycol, variable.name='pycode', value.name='musw_strat')


master_meansw=melt(meansw, id.vars=c("svspp", "byvar1", "byvar2", "byvar3", "byvar4"),     #add byvars AFTER svspp
                   measure.vars=preycol, variable.name='pycode', value.name='meansw')

master_meansw_s=melt(meansw_s, id.vars=c("svspp", "byvar1", "byvar2", "byvar3", "byvar4"),   #add byvars AFTER svspp
                     measure.vars=preycol, variable.name='pycode', value.name='meansw_s')

master_msw_strat=melt(new3s, id.vars=c("svspp", "byvar1", "byvar2", "byvar3", "byvar4"),   #add byvars AFTER svspp
                      measure.vars=preycol, variable.name='pycode', value.name='msw_strat')

keep7=c('cruise6', 'stratum',  'svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4', 'munfish_strat') #add byvars AFTER svspp
master_munfish_strat=musw_strat[,keep7, with=F]

num_stra_fish=new3s[,c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4','munfish_strat', 'num_stra'),with=F] #add byvars AFTER svspp
names(num_stra_fish)=gsub('munfish_strat','m_nfish_strat',names(num_stra_fish))

master_num_tows=meansw[,c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4','num_tows'), with=F] #add byvars AFTER svspp

#merge masters
#memory.size(4095)
merge1=merge(master_musw, master_musw_strat, by=c('cruise6', 'stratum',  'svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4', 'pycode'), all.x=T)  #add byvars AFTER svspp
merge2=merge(merge1, master_tmsw_strat, by=c('cruise6', 'stratum',  'svspp',   'byvar1', 'byvar2', 'byvar3', 'byvar4', 'pycode'), all.x=T)   #add byvars AFTER svspp
merge3=merge(merge2, master_munfish_strat, by=c('cruise6', 'stratum',  'svspp', 'byvar1', 'byvar2', 'byvar3', 'byvar4'), all.x=T)      #add byvars AFTER svspp
merge4=merge(merge3, master_musw2, by=c('cruise6', 'stratum', 'station', 'svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4','pycode'), all.x=T)   #add byvars AFTER svspp
merge5=merge(merge4, master_meansw, by=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4','pycode'), all.x=T)   #add byvars AFTER svspp
merge6=merge(merge5, master_num_tows, by=c('svspp', 'byvar1', 'byvar2', 'byvar3', 'byvar4'), all.x=T)          #add byvars AFTER svspp
merge7=merge(merge6, master_meansw_s, by=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4','pycode'), all.x=T)  #add byvars AFTER svspp
merge8=merge(merge7, num_stra_fish, by=c('svspp', 'byvar1', 'byvar2', 'byvar3', 'byvar4'), all.x=T)           #add byvars AFTER svspp
master=merge(merge8, master_msw_strat, by=c('svspp', 'byvar1', 'byvar2', 'byvar3', 'byvar4', 'pycode'), all.x=T)  #add byvars AFTER svspp


#prod    
master$prod=master$rcatnum^2*((master$musw-master$meansw))^2
master$prodf= (master$rcatnum-master$munfish_strat)^2
master$prodd= (master$musw2-master$musw_strat)^2
master$prod_cov=((master$rcatnum-master$munfish_strat)*(master$musw2-master$musw_strat))

#master=master[order(master$svspp),] #add byvars AFTER svspp ##causes crash?


#sprod   better way?  remove cols except svspp and pycode?
#fix?  master=data.table(master); master needs to be a data.table
master=data.table(master)
sprod=master[,lapply(.SD, sum, na.rm=T), by=list(svspp,byvar1,byvar2,byvar3,byvar4, pycode)]#, num_tows, meansw, meansw_s, msw_strat, m_nfish_strat, num_stra)] #add byvars AFTER svspp
keep8=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4','pycode', 'prod')#'num_tows', 'meansw', 'meansw_s', 'msw_strat', 'm_nfish_strat', 'num_stra', 'prod') #add byvars AFTER svspp
sprod2=sprod[,keep8, with=F]
names(sprod2)=gsub('prod','sprod',names(sprod2))

#mprod and mnumfish
mprod_mnumfish=master[,lapply(.SD, mean, na.rm=T), by=list(svspp, byvar1,byvar2,byvar3,byvar4, pycode, pycode)]#, num_tows, meansw, meansw_s, msw_strat, m_nfish_strat, num_stra)] #add byvars AFTER svspp
keep9=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4', 'pycode', 'prod', 'rcatnum')#'num_tows', 'meansw', 'meansw_s', 'msw_strat', 'm_nfish_strat', 'num_stra', 'prod', 'rcatnum')        #add byvars AFTER svspp
mprod_mnumfish2=mprod_mnumfish[,keep9,with=F]
names(mprod_mnumfish2)=gsub('prod', 'mprod', names(mprod_mnumfish2))
names(mprod_mnumfish2)=gsub('rcatnum', 'mnumfish', names(mprod_mnumfish2))

new4=merge(sprod2, mprod_mnumfish2, by=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4','pycode'), all.x=T)#, 'num_tows', 'meansw', 'meansw_s', 'msw_strat', 'm_nfish_strat', 'num_stra'), all.x=T)
#add byvars?

sprodf_d_cov=master[,lapply(.SD, sum, na.rm=T), by=list(svspp, byvar1,byvar2,byvar3,byvar4, pycode, cruise6, stratum, stratum_area, rtot_tows_spp_stratum, m_nfish_strat)]      #add byvars after svspp
keep10=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4', 'pycode', 'cruise6', 'stratum', 'stratum_area', 'rtot_tows_spp_stratum', 'm_nfish_strat', 'prodf', 'prodd', 'prod_cov')  #add byvars after svspp
sprodf_d_cov2=sprodf_d_cov[,keep10, with=F]
names(sprodf_d_cov2)=gsub('prodf', 'sprodf', names(sprodf_d_cov2))
names(sprodf_d_cov2)=gsub('prodd', 'sprodd', names(sprodf_d_cov2))
names(sprodf_d_cov2)=gsub('prod_cov', 'sprod_cov', names(sprodf_d_cov2))

sprodf_d_cov2$dfntows_strat= sprodf_d_cov2$rtot_tows_spp_stratum-1
sprodf_d_cov2$varprodf=ifelse(sprodf_d_cov2$rtot_tows_spp_stratum>1, (sprodf_d_cov2$stratum_area^2)*((sprodf_d_cov2$sprodf/sprodf_d_cov2$dfntows_strat)/sprodf_d_cov2$rtot_tows_spp_stratum), 0)
sprodf_d_cov2$varprodd= ifelse(sprodf_d_cov2$rtot_tows_spp_stratum>1, (sprodf_d_cov2$stratum_area^2)*((sprodf_d_cov2$sprodd/sprodf_d_cov2$dfntows_strat)/sprodf_d_cov2$rtot_tows_spp_stratum), 0)
sprodf_d_cov2$varprod_cov=ifelse(sprodf_d_cov2$rtot_tows_spp_stratum>1, (sprodf_d_cov2$stratum_area^2)*((sprodf_d_cov2$sprod_cov/sprodf_d_cov2$dfntows_strat)/sprodf_d_cov2$rtot_tows_spp_stratum), 0)

new6=sprodf_d_cov2[order(sprodf_d_cov2$svspp,sprodf_d_cov2$byvar1,sprodf_d_cov2$byvar2,sprodf_d_cov2$byvar3,sprodf_d_cov2$byvar4),] #add byvars AFTER svspp

#sumvarprod     some vars still zero?
sumvarprod1=new6[,lapply(.SD,sum, na.rm=T), by=list(svspp, byvar1,byvar2,byvar3,byvar4,pycode)] #add byvars AFTER svspp
keep11=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4', 'pycode', 'varprodf', 'varprodd', 'varprod_cov', 'stratum_area')   #add byvars AFTER svspp
sumvarprod=sumvarprod1[,keep11, with=F]
names(sumvarprod)=gsub('varprodf', 'svarprodf', names(sumvarprod))
names(sumvarprod)=gsub('varprodd', 'svarprodd', names(sumvarprod))
names(sumvarprod)=gsub('varprod_cov', 'svarprod_cov', names(sumvarprod))
names(sumvarprod)=gsub('stratum_area', 'sstratum_area', names(sumvarprod))

sumvarprod$varf=sumvarprod$svarprodf/sumvarprod$sstratum_area^2
sumvarprod$vard=sumvarprod$svarprodd/sumvarprod$sstratum_area^2
sumvarprod$var_cov=sumvarprod$svarprod_cov/sumvarprod$sstratum_area^2

#new4_strat
new4_strat=merge(new4, sumvarprod, by=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4', 'pycode'), all.x=T) #add byvars AFTER svspp

#merge with  select master columns


six1=merge(new4_strat, master_meansw, by=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4','pycode'), all.x=T)  #add byvars AFTER svspp
six2= merge(six1, master_meansw_s, by=c('svspp', 'byvar1', 'byvar2', 'byvar3', 'byvar4', 'pycode'), all.x=T)      #add byvars AFTER svspp
six3=merge(six2, master_msw_strat, by=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4','pycode'), all.x=T)     #add byvars AFTER svspp
six4=merge(six3, master_num_tows, by=c('svspp', 'byvar1', 'byvar2', 'byvar3', 'byvar4'), all.x=T) #add byvars AFTER svspp
six=merge(six4, num_stra_fish, by=c('svspp', 'byvar1', 'byvar2', 'byvar3', 'byvar4'), all.x=T) #add byvars AFTER svspp 


six$variance=ifelse(six$num_tows>1&six$mnumfish!=0&six$meansw!=0&six$m_nfish_strat!=0&six$meansw_s!=0, 1/(six$num_tows*(six$mnumfish)^2)*(six$sprod/(six$num_tows-1)),0)
six$cv=ifelse(six$num_tows>1&six$mnumfish!=0&six$meansw!=0&six$m_nfish_strat!=0&six$meansw_s!=0,((six$variance)^0.5)/six$meansw,0)
six$var_s=ifelse(six$num_tows>1&six$mnumfish!=0&six$meansw!=0&six$m_nfish_strat!=0&six$meansw_s!=0, (six$meansw_s^2)*((six$varf/(six$m_nfish_strat^2))+(six$vard/(six$msw_strat^2))-(2*six$var_cov/six$m_nfish_strat/six$msw_strat)),0)  
six$cv_s=ifelse(six$num_tows>1&six$mnumfish!=0&six$meansw!=0&six$m_nfish_strat!=0&six$meansw_s!=0,((six$var_s)^0.5)/six$meansw_s,0)

#renamed meansw since used earlier
meanswf=six[,':='(sprod=NULL, mprod=NULL, mnumfish=NULL, msw_strat=NULL, svarprodf=NULL, svarprodd=NULL, svarprod_cov=NULL, sstratum_area=NULL, varf=NULL, vard=NULL, var_cov=NULL)]

#merge tax names back in
output=merge(meanswf, pylist, by=c('pycode'), all.x=T)


names(output)=gsub('tax', 'prey', names(output))

totwt1=output[,lapply(.SD,sum, na.rm=T), by=list(svspp,byvar1,byvar2,byvar3,byvar4), .SDcols=c('meansw', 'meansw_s')] #add byvars AFTER svspp
keep12=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4', 'meansw', 'meansw_s')  #add byvars AFTER svspp
totwt=totwt1[,keep12,with=F]
names(totwt)=gsub('meansw', 'totwt', names(totwt))
names(totwt)=gsub('meansw_s', 'totwt_s', names(totwt))

total=merge(output, totwt, by=c('svspp',  'byvar1', 'byvar2', 'byvar3', 'byvar4'), all.x=T)  #add byvars AFTER svspp

final=total[,':='(pycode=NULL, m_nfish_strat=NULL)]

final$relmsw=100*(final$meansw/final$totwt)
final$relmsw_s=100*(final$meansw_s/final$totwt_s)
final$ci=sqrt(final$variance/final$num_tows)*2
final$relci=(final$ci/final$totwt)*100
final$ci_s=sqrt(final$var_s/final$num_tows)*2
final$relci_s=(final$ci_s/final$totwt_s)*100

#remove strata-based metrics for now; only minor differences from original weighted metrics by tow  
#final[,':=' (meansw_s=NULL, num_stra=NULL, var_s=NULL, cv_s=NULL, totwt_s=NULL, relmsw_s=NULL, ci_s=NULL, relci_s=NULL )]


if(unique(final$byvar1=='x')){final$byvar1=NULL} 
if(unique(final$byvar2=='x')){final$byvar2=NULL}
if(unique(final$byvar3=='x')){final$byvar3=NULL}
if(unique(final$byvar4=='x')){final$byvar4=NULL}

names(final)[names(final)=='byvar1']=name_byvar1                                  
names(final)[names(final)=='byvar2']=name_byvar2
names(final)[names(final)=='byvar3']=name_byvar3                                  
names(final)[names(final)=='byvar4']=name_byvar4                                  


final

#export csv
write.csv(final, 'allwt_panfamGOM.csv')

#similar to fo, come up with wtd mean % over preds, wtd by nstrat
#save for lydia
#column of interest: relmsw_s which is the stratified mean % mass of nshrimp by svspp (predator), season, and year.

library(data.table)
drive_download('allwt_panfamGOM.csv', overwrite = T) 

allwt=read.csv('allwt_panfamGOM.csv', h=T, sep = ',', na.strings = c(NA))

allwt2=subset(allwt, allwt$prey=='PANFAM' )

allwt3=data.table(allwt2[,c(3:5, 9,18 )]) #remove pred and prey columns to allow for wt.mean

m.wtm=allwt3[,lapply(.SD, weighted.mean, num_stra, na.rm=TRUE), by=list(geoarea, season, year)] 
#rm cols
m.wtm[,':='(num_stra=NULL)]
m.wtm2=m.wtm[order(m.wtm$geoarea, m.wtm$season, m.wtm$year,)]
m.wtm2

#start here
library(Hmisc)
wtmean.var=allwt3[, lapply(.SD, wtd.var, num_stra, na.rm=TRUE), by=list(geoarea, season, year)]
#rm cols
wtmean.var[,':='(num_stra=NULL)]

names(wtmean.var)[ncol(wtmean.var)]='wtmn.var'

head(wtmean.var)

#merge wtmean.var and m.wtm datasets 
m.n=aggregate(list('n'=allwt2$relmsw_s), by=list('geoarea'=allwt2$geoarea, 'season'=allwt2$season, 'year'=allwt2$year), length)


allwt.wtm=merge(merge(m.wtm2, wtmean.var, by=c('geoarea', 'season', 'year'), all.x=T), m.n, by=c('season', 'year'), all.x=T)
allwt.wtm$ci = 2 * sqrt(allwt.wtm$wtmn.var / allwt.wtm$n)
#need to add 2020=NA for spring and fall


head(allwt.wtm)

write.csv(allwt.wtm, file = "/Users/lydia.clark/Downloads/data/allwt.wtm.csv", row.names = FALSE)

library(data.table)
library(Hmisc)
library(ggplot2)
library(dplyr)


#ts of stratified number of shrimp per tows in both seasons

shrimp <- filter(allwt, prey == "PANFAM")
shrimp <- as.data.table(shrimp) 

ts_d <- shrimp[, .(meanshrimptow = mean(meansw, na.rm = TRUE)), by = .(season, year)]

ggplot( ts_d, aes ( x = year, y = meanshrimptow, color = season))+
  geom_line(linewidth = 1)+
  geom_point()+
  labs(
    title = "time series of stratified mean number of shrimp per tow", 
    x = "year", y = "mean number of shrimp per tow", 
    color = "season")+
  theme_minimal()


#DMR
library(readr)
DMRindices <- read_csv("MaineNShrimp_Indices_thru2024.csv")
View(DMRindices)

setDT(DMRindices)

##stratified catch (#/tow) 


##plot kg per tow

p1 <- ggplot(DMRindices, aes (x = Year, y = StratMean_Weight, color = Season))+
  geom_line()+
  geom_point()+
  labs(title = "stratified mean kg of shrimp per tow", x = "year", y = "kg per tow")+
  theme_minimal()

p2<- ggplot(DMRindices, aes (x = Year, y = StratMean_Catch, color = Season))+
  geom_line()+
  geom_point()+
  labs(title = "stratified mean number of shrimp per tow", x = "year", y = "numbers per tow")+
  theme_minimal()

library(gridExtra)
grid.arrange(p1, p2, ncol = 2)


#weighted mean over time
##use dataset allwt_wtm (relmsw_s)

ggplot(allwt.wtm, aes ( x = year, y = relmsw_s, color = season))+
  geom_line(linewidth = 1)+
  geom_point()+
  labs(
    title = "weighted mean over time", 
    x = "year", y = "weighted mean", 
    color = "season")+
  theme_minimal()


#percent frequency over time
##use fo.m2 (m.pctfo)

ggplot(fo.m2, aes ( x = year, y = m.pctfo, color = season))+
  geom_line(linewidth = 1)+
  geom_point()+
  labs(
    title= "mean %FO over time", 
    x = "year", y = "%FO", 
    color = "season")+
  theme_minimal()


#amo and nao datasets

library(readr)
library(dplyr)
amo <- read_csv("amo.csv")

nao <- read.csv("nao.csv")

##now to "winterize" nao dataset

nao_w <- filter(nao, month %in% c(1, 2, 3, 12))
nao_w

library(tibble)

##nao-making the date a column so month and year aren't separate 

##note- first amo and nao graphs don't use mean values

nao_date <- nao_w %>% 
  mutate(
    Month = base::match(month, base::month.name),
    Date = as.Date(paste(year, '-', month, '-01', sep=''))
  ) %>% 
  select(-c(month, year))

nao_date <- subset( nao_date, select = -Month )

head(nao_date) #worked

ggplot(nao_date, aes( x = Date, y = nao))+
  geom_line(color = "steelblue")+
  labs( title = "nao values per month per year")+
  theme_minimal()

#amo- reformat the data so it's date and value vs how it is now?

library (tidyr)

##amo graph 1

amo_long <- amo %>%
  pivot_longer(cols = jan:dec, names_to = "month", values_to = "value") %>%
  mutate(month_num = match(tolower(month), tolower(month.abb)),
         date = as.Date(paste(year, month_num, "01", sep = "-"))) %>%
  arrange(date)

amo_clean <- amo_long[, c(1, 2, 3, 5)]
colnames(amo_clean) <- c("year", "month", "value", "date")

amo_clean <- amo_clean %>%
  filter(value > -99)


ggplot(amo_clean, aes(x = date, y = value, group = 1)) +
  geom_line(color = "steelblue") +
  labs(title = "AMO Monthly Data Over Time",
       x = "Date",
       y = "Value") +
  theme_minimal()

##positive amo - higher SST; negative AMO-lower SST

##amo means?

amo_means <- amo_clean %>%
  group_by( year) %>%
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%
  arrange(year)

ggplot(amo_means, aes ( x = year, y = mean_value))+
  geom_line(linewidth = 1)+
  geom_point()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  labs(
    title= "AMO mean values", 
    x = "year", y = " Means")+
  theme_minimal()

write.csv(amo_means,'amo_means.csv')

##nao means

nao_means <- nao_w %>%
  group_by( year) %>%
  summarise(mean_value = mean(nao, na.rm = TRUE)) %>%
  arrange(year)


ggplot(nao_means, aes(x = year, y = mean_value)) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  
  labs(
    title = "NAO mean values", 
    x = "Year", 
    y = "Means"
  ) +
  theme_minimal()

##strong positive phases of nao associated with warmer weather in the US



###ECOMON DATA####
drive_download('EcoMon_Pandalidae_Conc_Count.csv', overwrite = T)
em=read.csv('EcoMon_Pandalidae_Conc_Count.csv', head=T, sep=',', na.strings=c(NA))

names(em)=tolower(names(em))

#add season
em$season=ifelse((em$month>=3&em$month<=5), 'SPRING', ifelse((em$month>=9&em$month<=11), 'FALL', 'OTHER'))

#use caridea_conc for now as standardized n true shrimp/100m3 water

#limit to gom
em2=subset(em, em$region=='GOM'&(em$season=='FALL'|em$season=='SPRING'))

em3=aggregate(list('m.caridea'=em2$caridea_conc), by=list('season'=em2$season, 'year'=em2$year), mean)
em4=em3[order(em3$season, em3$year),]

write.csv(em4, 'em4.csv')

#graphing
ggplot(em4, aes ( x = year, y = m.caridea, color = season))+
  geom_line(linewidth = 1)+
  geom_point()+
  labs(
    title= "ECOMON graph", 
    x = "year", y = "M.caridea per 100m3 water", 
    color = "season")+
  theme_minimal()


wilcox.test(m.caridea ~ season, data = em4,
            exact = FALSE)

###wilcox test indicates that the difference between weights for each season is not significant

kruskal.test(m.caridea~season, data = em4)
###kruskal test indicates the same as wilcox (pvalues >0.05)

###looks like more shrimp are caught in the spring

mean(em4$m.caridea) #18.3
median(em4$m.caridea) #7.7
max(em4$m.caridea) #146.4

t.test(data = em4, m.caridea ~ season, alternative = "greater" )
###ttest shows that spring mean is higher than fall mean
###pvalue futher indicates difference isn't significant


#TS WORKING#

##want to use em4, panbor, survdat6, DMRindices, allwt_wtm, fo.m2
##should split dataframes by season

#panbor fall

library(readr)
panbor <- read_csv("panbor.csv")
View(panbor)

library(dplyr)
library(lubridate)
library(zoo)

panbor_fall <- panbor %>%
  filter(season == "FALL") %>%
  mutate(
    month = as.numeric(month),
    Date = as.Date(paste(year, month, "01", sep = "-"))
  )

panbor_subset <- panbor_fall %>%
  select(Date, season, pyamtv, svspp) %>%
  arrange(Date)

library(dplyr)
library(ggplot2)
library(forecast)

panbor_yearly <- panbor_subset %>%
  mutate(year = as.numeric(format(Date, "%Y"))) %>%
  group_by(year) %>%
  summarise(mean_pyamtv = mean(pyamtv, na.rm = TRUE)) %>%
  arrange(year)


panborfall_ts <- ts(panbor_yearly$mean_pyamtv, start = min(panbor_yearly$year), frequency = 1)

autoplot(panborfall_ts) +
  labs(title = "Fall-Yearly Mean  pyamtv Over Time", x = "Year", y = "Mean pyamtv") +
  theme_minimal()

#panbor spring
library(dplyr)

panbor_spring <- panbor %>%
  filter(season == "SPRING") %>%
  mutate(
    month = as.numeric(month),
    Date = as.Date(paste(year, month, "01", sep = "-"))
  )

panbor_subset2 <- panbor_spring %>%
  select(Date, season, pyamtv, svspp) %>%
  arrange(Date)

panbor_yearly2 <- panbor_subset2 %>%
  mutate(year = as.numeric(format(Date, "%Y"))) %>%
  group_by(year) %>%
  summarise(mean_pyamtv = mean(pyamtv, na.rm = TRUE)) %>%
  arrange(year)


panborspring_ts <- ts(panbor_yearly2$mean_pyamtv, start = min(panbor_yearly2$year), frequency = 1)

autoplot(panborspring_ts) +
  labs(title = "Spring-Yearly Mean  pyamtv Over Time", x = "Year", y = "Mean pyamtv") +
  theme_minimal()


#em4 fall
em4_fall <- em4 %>%
  filter(season == "FALL") %>%
  select(c(-season, -year))

em4_fallts <- ts(em4_fall, start = 1977)
autoplot(em4_fallts)


#em4 spring

em4_spring <-em4 %>%
  filter (season == "SPRING") %>%
  select(c(-season, -year))

em4_springts <- ts (em4_spring, start = 1977)

autoplot(em4_springts)


#survdat6 fall

survdat6fall <- filter (survdat6, season == "FALL")

survdat6_fall <- select(survdat6fall, year, abundance, biomass)

survdat6ts_fall <- ts(data = survdat6_fall, start = c(1992), frequency = 1)

autoplot(survdat6ts_fall, facets = TRUE)

#survdat6 spring
survdat6spring <- filter (survdat6, season == "SPRING")

survdat6_spring <- select(survdat6spring, year, abundance, biomass)
survdat6ts_spring <- ts(data = survdat6_spring, start = c(1992))
autoplot(survdat6ts_spring, facets = TRUE)
##^for some reason when graphing this, the year only goes to 2015

library(gridExtra)
##solution?
p3 <- ggplot(survdat6_spring, aes(  x = year,   y = abundance))+
  geom_line()+
  labs(title = "survdat6 spring abundance and biomass over time")

p4 <- ggplot(survdat6_spring, aes ( x = year, y = biomass))+
  geom_line()+
  labs()

grid.arrange(p3, p4)


#DMR ts

DMRfiltered <- select(DMRindices, StratMean_Catch, StratMean_Weight)

DMRfall <- DMRfiltered[1:25, ] #only fall

DMRts_fall <- ts (data = DMRfall, start = c(2001), frequency = 1)
autoplot(DMRts_fall, facets = TRUE)


DMRspring <- DMRfiltered [26:48, ] #only spring

DMRts_spring <- ts(data = DMRspring, start = c(2001), frequency = 1)

autoplot(DMRts_spring, facets = TRUE)
##both time series include year and survey as a graph, unsure how to fix


#allwt_wtm
##relmsw_s is weighted mean value 

allwt_wtm_filtered <- select(allwt.wtm, season, year, relmsw_s)
allwt_fall <- filter (allwt_wtm_filtered, season == "FALL") #only fall 
allwt_spring <- filter (allwt_wtm_filtered, season == "SPRING") #only spring 

allwt_fall <- select(allwt_fall, relmsw_s)
allwt_spring <- select (allwt_spring, relmsw_s)

allwtts_fall <- ts(data = allwt_fall, start = c(1973))
autoplot(allwtts_fall)

allwtts_spring <- ts (data = allwt_spring, start = c(1973))
autoplot(allwtts_spring)


##graph fix

allwt_filtered <- allwt.wtm %>%
  select(season, year, relmsw_s) %>%
  mutate(season = toupper(season)) %>%
  filter(year >= 1973 & year <= 2021)

allwt_filled <- allwt_filtered %>%
  complete(year = 1973:2021, season = c("FALL", "SPRING"))

allwt_fall <- allwt_filled %>%
  filter(season == "FALL") %>%
  arrange(year)

allwt_spring <- allwt_filled %>%
  filter(season == "SPRING") %>%
  arrange(year)

allwt_ts_fall <- ts(allwt_fall$relmsw_s, start = 1973)
allwt_ts_spring <- ts(allwt_spring$relmsw_s, start = 1973)

autoplot(allwt_ts_fall) + ggtitle("Allwt Fall Predation Index")
autoplot(allwt_ts_spring) + ggtitle("Allwt Spring Predation Index")



#fo.m2

fo_sub <- fo.m2 %>%
  mutate(season = toupper(season)) %>%
  filter(year >= 1973 & year <= 2021) %>%
  complete(year = 1973:2021, season = c("FALL", "SPRING"))  

fo_fall <- fo_sub %>%
  filter(season == "FALL") %>%
  arrange(year)

fo_spring <- fo_sub %>%
  filter(season == "SPRING") %>%
  arrange(year)

fo_fall_ts <- ts(fo_fall$m.pctfo, start = 1973)
fo_spring_ts <- ts(fo_spring$m.pctfo, start = 1973)

autoplot(fo_fall_ts) + ggtitle("FO Fall - m.pctfo") + theme_minimal()
autoplot(fo_spring_ts) + ggtitle("FO Spring - m.pctfo") + theme_minimal()


write.csv(nao_means, 'nao_means.csv')
