##############################################################################
# FileName:		  R_Epilepsy_Wordcloud.R
# Description:	This file has the R code that can be used to get abstracts from PUBMED
#               and analyze them with text analytics 
# Author(s):    Jai Singh
# Version:      $Id$
# History:
# 2017/11/20  First Version
###############################################################################

#get current working dir
getwd()
options(max.print=1000000)

#Set current working directory to the clustering folder
setwd("//home/a159899/01-Projects/Pubmed_Analysis")
getwd()


##Install packages
# install.packages("/home/a159899/01-Projects/Pubmed_Analysis/felixfan-PubMedWordcloud-v0.3-14-g00674ad.tar.gz", repos = NULL, type="source")
# install.packages("/home/a159899/01-Projects/Pubmed_Analysis/slam_0.1-37.tar.gz", repos = NULL, type="source")
# install.packages("devtools")
# install.packages("twitteR")
# install.packages("reshape2")
# install.packages("sentimentr")
# install.packages("plyr")
# install.packages("ggplot2")
# install.packages("lazyeval")
# install.packages("wordcloud")
# install.packages("RColorBrewer")
# install.packages("ggplot2")
# install.packages("SnowballC")
# install.packages("devtools")
# install.packages("tm")
# install.packages("shiny")
# install.packages("shinythemes")
# install.packages("rsconnect")
# install.packages("NLP")
# install.packages("openNLP")
# install.packages("reshape2")
# install.packages("RColorBrewer")
# install.packages("plotly")
# install.packages("topicmodels")
# install.packages("tidytext")
# install.packages("DT")
# install.packages("sentimentr")
# install.packages("dplyr")
# install.packages("RWeka")
# install.packages("shiny")
# install.packages("shinythemes")
# install.packages("stringr")
# install.packages("RedditExtractoR")
# install.packages("scales")
# install.packages("qdap")
# install.packages("plotly")
# install.packages("magrittr")
# install.packages("shinydashboard")
# install.packages("textstem")
# install.packages("stringr")
# install.packages("splitstackshape")
# install.packages("tidyr")


#call libraries
library(PubMedWordcloud)
library(textstem)
library(RColorBrewer)
library(devtools)
library(twitteR)
library(reshape2)
library(sentimentr)
library(plyr)
library(ggplot2)
library(lazyeval)
library(wordcloud)
library(RColorBrewer)
library(ggplot2)
library(SnowballC)
library(devtools)
library(tm)
library(shiny)
library(dplyr)
library(shinythemes)
library(rsconnect)
library(NLP)
library(openNLP)
library(reshape2)
library(RColorBrewer)
library(plotly)
library(topicmodels)
library(tidytext)
library(DT)
library(sentimentr)
library(dplyr)
library(RWeka)
library(shiny)
library(shinythemes)
library(stringr)
library(RedditExtractoR)
library(scales)
library(qdap)
library(plotly)
library(magrittr)
library(shinydashboard)
library(stringr)
library(splitstackshape)
library(tidyr)

#search by keyword
epilepsy_drug<-getPMIDsByKeyWords(keys = "Antiepileptic Drug", journal = NULL, dFrom = 2007,
                                  dTo = 2017,n=10000, https = TRUE)
#Save epilepsy abstracts as rds
epilepsy_abs<-getAbstracts(epilepsy, https = TRUE, s = 100)
#saveRDS(epilepsy_abs,"epilepsy_abstracts.rds")

#Get epilepsy drug abstracts and save them
epilepsy_drug_abs<-getAbstracts(epilepsy_drug, https = TRUE, s = 100)
#saveRDS(epilepsy_drug_abs,"epilepsy_drug_abstracts.rds")

#Read previously saved rds data
#epilepsy_drug_abs<-readRDS("epilepsy_drug_abstracts.rds")
#epilepsy_abs<-readRDS("epilepsy_abstracts.rds")

#clean abstracts of stopwords, line nums, etc
clean_epilepsy=cleanAbstracts(epilepsy_abs)
clean_epilepsy_drug=cleanAbstracts(epilepsy_drug_abs)
dim(clean_epilepsy_drug)
head(clean_epilepsy_drug)
write.csv(clean_epilepsy_drug, file = "clean_epilepsy_drug_word_freq.csv")

#Plot wordclouds
plotWordCloud(clean_epilepsy_drug, scale = c(3, 0.3), min.freq = 1, max.words = 100,
random.order = FALSE, rot.per = 0.35, use.r.layout = FALSE,
colors = brewer.pal(8, "Dark2"))

#Read in drug names from https://catalog.data.gov/dataset/drugsfda-database
drugs <- read_delim("~/01-Projects/Pubmed_Analysis/drug_name_full.txt", 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
head(drugs)
names(drugs)
lapply(drugs,typeof)
drugs$stmarkdate<-as.Date(drugs$stmarkdate)


#convert to data.frame for inner join
drugs<-data.frame(drugs)
clean_epilepsy_drug<-data.frame(clean_epilepsy_drug)
epilepsy_drug<-merge(clean_epilepsy_drug, drugs, by.x="word",by.y="propname")
head(ftable(epilepsy_drug$pharmclas),2)
#Export this file for tableau
write.csv(epilepsy_drug, file = "EpilepsyDrugs.csv", row.names=FALSE)

#Split 'Pharmclas' into multiple columns
epilepsy_drug2<-cSplit(epilepsy_drug, "pharmclas", sep=",")
epilepsy_drug3 <- melt(epilepsy_drug2, id = c("word","freq","ndc","prodtype","npropname","dosename","routename","stmarkdatestr","stmarkdate","markname","appnum","labelname","subname","actnumstr","actingunit"))



pharmaclass<-data.frame(table(epilepsy_drug$pharmclas))
pharmaclass<-pharmaclass[order(pharmaclass[,2],decreasing=TRUE),]
barplot(table(epilepsy_drug$pharmclas),horiz = TRUE)

#remove duplicated rows
epilepsy_drug = epilepsy_drug[!duplicated(epilepsy_drug$word),]
epilepsy_drug$word

#Order by frequency
epilepsy_drug = epilepsy_drug[order(-epilepsy_drug[,'freq']),]
dim(epilepsy_drug)
head(epilepsy_drug,20)

#Convert to df
epilepsy_drug_abs<-data.frame(epilepsy_drug_abs)

# Keep only drugs 
epilepsy_drug_abs2<-dplyr::filter(epilepsy_drug_abs, grepl(
  'gabapentin|
levetiracetam|
  topiramate|
  lamotrigine|
  acetazolamide|
  analgesic|
  riluzole|
  zonisamide|
  propofol|
  ethosuximide|
  lidocaine|
  fluoxetine|
  duloxetine|
  haloperidol|
  tizanidine|
  acetaminophen|
  venlafaxine|
  adenosine|
  ibuprofen|
  risperidone|
  phosphate|
  olanzapine|
  aripiprazole|
  naproxen|
  antibacterial|
  celecoxib|
  nicotine|
  indomethacin|
  temozolomide|
  furosemide|
  sumatriptan|
  testosterone|
  bupropion|
  etomidate|
  nifedipine|
  modafinil|
  temazepam|
  ceftriaxone|
  diphenhydramine|
  pioglitazone|
  escitalopram|
  ondansetron|
  simvastatin|
  ketoconazole|
  paroxetine|
  cimetidine|
  mirtazapine|
  paclitaxel|
  amantadine|
  fenofibrate|
  zaleplon|
  ciprofloxacin|
  erythromycin|
  meloxicam|
  ropinirole|
  hydrochlorothiazide|
  oxaliplatin|
  aspirin|
  hydrocortisone|
  isoniazid|
  azithromycin|
  cefepime|
  cytarabine|
  decitabine|
  sirolimus|
  doxycycline|
  epinephrine|
  eszopiclone|
  letrozole|
  carboplatin|
  carisoprodol|
  clarithromycin|
  clopidogrel|
  menthol|
  methazolamide|
  ofloxacin|
  rifampin|
  telmisartan|
  atomoxetine|
  gemcitabine|
  metronidazole|
  acyclovir|
  fluorouracil|
  itraconazole|
  tetrabenazine|
  ampicillin|
  disposable|
  levofloxacin|
  valsartan|
  voriconazole|
  omeprazole|
  calcitriol|
  cefazolin|
  cortisone|
  fluconazole|
  linezolid|
  lovastatin|
  nevirapine|
  oxybutynin|
  bexarotene|
  cilostazol|
  clotrimazole|
  emilia|
  expectorant|
  guanfacine|
  hemorrhoidal|
  laxative|
  leflunomide|
  lisinopril|
  loratadine|
  misoprostol|
  piroxicam|
  ribavirin|
  almotriptan|
  bacitracin|
  budesonide|
  chlorzoxazone|
  cyanocobalamine|
  ganciclovir|
  gatifloxacin|
  lansoprazole|
  mercaptopurine|
  minoxidil|
  mupirocin|
  oxandrolone|
  propylthiouracil|
  rifabutin|
  suboxone|
  sulfacetamide|
  topotecan|'
  , ignore.case=TRUE,epilepsy_drug_abs))

#check dimensions
dim(epilepsy_drug_abs2)

#create drug_indicator variables
epilepsy_drug_abs2$gabapentin <- as.factor(ifelse(grepl('gabapentin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$levetiracetam <- as.factor(ifelse(grepl('levetiracetam',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$topiramate <- as.factor(ifelse(grepl('topiramate',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$lamotrigine <- as.factor(ifelse(grepl('lamotrigine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$acetazolamide <- as.factor(ifelse(grepl('acetazolamide',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$analgesic <- as.factor(ifelse(grepl('analgesic',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$riluzole <- as.factor(ifelse(grepl('riluzole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$zonisamide <- as.factor(ifelse(grepl('zonisamide',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$propofol <- as.factor(ifelse(grepl('propofol',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$ethosuximide <- as.factor(ifelse(grepl('ethosuximide',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$lidocaine <- as.factor(ifelse(grepl('lidocaine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$fluoxetine <- as.factor(ifelse(grepl('fluoxetine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$duloxetine <- as.factor(ifelse(grepl('duloxetine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$haloperidol <- as.factor(ifelse(grepl('haloperidol',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$tizanidine <- as.factor(ifelse(grepl('tizanidine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$acetaminophen <- as.factor(ifelse(grepl('acetaminophen',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$venlafaxine <- as.factor(ifelse(grepl('venlafaxine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$adenosine <- as.factor(ifelse(grepl('adenosine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$ibuprofen <- as.factor(ifelse(grepl('ibuprofen',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$risperidone <- as.factor(ifelse(grepl('risperidone',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$phosphate <- as.factor(ifelse(grepl('phosphate',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$olanzapine <- as.factor(ifelse(grepl('olanzapine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$aripiprazole <- as.factor(ifelse(grepl('aripiprazole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$naproxen <- as.factor(ifelse(grepl('naproxen',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$antibacterial <- as.factor(ifelse(grepl('antibacterial',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$celecoxib <- as.factor(ifelse(grepl('celecoxib',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$nicotine <- as.factor(ifelse(grepl('nicotine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$indomethacin <- as.factor(ifelse(grepl('indomethacin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$temozolomide <- as.factor(ifelse(grepl('temozolomide',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$furosemide <- as.factor(ifelse(grepl('furosemide',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$sumatriptan <- as.factor(ifelse(grepl('sumatriptan',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$testosterone <- as.factor(ifelse(grepl('testosterone',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$bupropion <- as.factor(ifelse(grepl('bupropion',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$etomidate <- as.factor(ifelse(grepl('etomidate',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$nifedipine <- as.factor(ifelse(grepl('nifedipine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$modafinil <- as.factor(ifelse(grepl('modafinil',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$temazepam <- as.factor(ifelse(grepl('temazepam',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$ceftriaxone <- as.factor(ifelse(grepl('ceftriaxone',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$diphenhydramine <- as.factor(ifelse(grepl('diphenhydramine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$pioglitazone <- as.factor(ifelse(grepl('pioglitazone',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$escitalopram <- as.factor(ifelse(grepl('escitalopram',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$ondansetron <- as.factor(ifelse(grepl('ondansetron',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$simvastatin <- as.factor(ifelse(grepl('simvastatin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$ketoconazole <- as.factor(ifelse(grepl('ketoconazole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$paroxetine <- as.factor(ifelse(grepl('paroxetine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$cimetidine <- as.factor(ifelse(grepl('cimetidine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$mirtazapine <- as.factor(ifelse(grepl('mirtazapine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$paclitaxel <- as.factor(ifelse(grepl('paclitaxel',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$amantadine <- as.factor(ifelse(grepl('amantadine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$fenofibrate <- as.factor(ifelse(grepl('fenofibrate',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$zaleplon <- as.factor(ifelse(grepl('zaleplon',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$ciprofloxacin <- as.factor(ifelse(grepl('ciprofloxacin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$erythromycin <- as.factor(ifelse(grepl('erythromycin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$meloxicam <- as.factor(ifelse(grepl('meloxicam',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$ropinirole <- as.factor(ifelse(grepl('ropinirole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$hydrochlorothiazide <- as.factor(ifelse(grepl('hydrochlorothiazide',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$oxaliplatin <- as.factor(ifelse(grepl('oxaliplatin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$aspirin <- as.factor(ifelse(grepl('aspirin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$hydrocortisone <- as.factor(ifelse(grepl('hydrocortisone',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$isoniazid <- as.factor(ifelse(grepl('isoniazid',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$azithromycin <- as.factor(ifelse(grepl('azithromycin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$cefepime <- as.factor(ifelse(grepl('cefepime',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$cytarabine <- as.factor(ifelse(grepl('cytarabine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$decitabine <- as.factor(ifelse(grepl('decitabine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$sirolimus <- as.factor(ifelse(grepl('sirolimus',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$doxycycline <- as.factor(ifelse(grepl('doxycycline',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$epinephrine <- as.factor(ifelse(grepl('epinephrine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$eszopiclone <- as.factor(ifelse(grepl('eszopiclone',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$letrozole <- as.factor(ifelse(grepl('letrozole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$carboplatin <- as.factor(ifelse(grepl('carboplatin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$carisoprodol <- as.factor(ifelse(grepl('carisoprodol',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$clarithromycin <- as.factor(ifelse(grepl('clarithromycin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$clopidogrel <- as.factor(ifelse(grepl('clopidogrel',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$menthol <- as.factor(ifelse(grepl('menthol',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$methazolamide <- as.factor(ifelse(grepl('methazolamide',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$ofloxacin <- as.factor(ifelse(grepl('ofloxacin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$rifampin <- as.factor(ifelse(grepl('rifampin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$telmisartan <- as.factor(ifelse(grepl('telmisartan',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$atomoxetine <- as.factor(ifelse(grepl('atomoxetine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$gemcitabine <- as.factor(ifelse(grepl('gemcitabine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$metronidazole <- as.factor(ifelse(grepl('metronidazole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$acyclovir <- as.factor(ifelse(grepl('acyclovir',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$fluorouracil <- as.factor(ifelse(grepl('fluorouracil',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$itraconazole <- as.factor(ifelse(grepl('itraconazole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$tetrabenazine <- as.factor(ifelse(grepl('tetrabenazine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$ampicillin <- as.factor(ifelse(grepl('ampicillin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$disposable <- as.factor(ifelse(grepl('disposable',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$levofloxacin <- as.factor(ifelse(grepl('levofloxacin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$valsartan <- as.factor(ifelse(grepl('valsartan',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$voriconazole <- as.factor(ifelse(grepl('voriconazole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$omeprazole <- as.factor(ifelse(grepl('omeprazole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$calcitriol <- as.factor(ifelse(grepl('calcitriol',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$cefazolin <- as.factor(ifelse(grepl('cefazolin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$cortisone <- as.factor(ifelse(grepl('cortisone',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$fluconazole <- as.factor(ifelse(grepl('fluconazole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$linezolid <- as.factor(ifelse(grepl('linezolid',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$lovastatin <- as.factor(ifelse(grepl('lovastatin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$nevirapine <- as.factor(ifelse(grepl('nevirapine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$oxybutynin <- as.factor(ifelse(grepl('oxybutynin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$bexarotene <- as.factor(ifelse(grepl('bexarotene',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$cilostazol <- as.factor(ifelse(grepl('cilostazol',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$clotrimazole <- as.factor(ifelse(grepl('clotrimazole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$emilia <- as.factor(ifelse(grepl('emilia',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$expectorant <- as.factor(ifelse(grepl('expectorant',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$guanfacine <- as.factor(ifelse(grepl('guanfacine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$hemorrhoidal <- as.factor(ifelse(grepl('hemorrhoidal',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$laxative <- as.factor(ifelse(grepl('laxative',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$leflunomide <- as.factor(ifelse(grepl('leflunomide',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$lisinopril <- as.factor(ifelse(grepl('lisinopril',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$loratadine <- as.factor(ifelse(grepl('loratadine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$misoprostol <- as.factor(ifelse(grepl('misoprostol',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$piroxicam <- as.factor(ifelse(grepl('piroxicam',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$ribavirin <- as.factor(ifelse(grepl('ribavirin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$almotriptan <- as.factor(ifelse(grepl('almotriptan',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$bacitracin <- as.factor(ifelse(grepl('bacitracin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$budesonide <- as.factor(ifelse(grepl('budesonide',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$chlorzoxazone <- as.factor(ifelse(grepl('chlorzoxazone',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$cyanocobalamine <- as.factor(ifelse(grepl('cyanocobalamine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$ganciclovir <- as.factor(ifelse(grepl('ganciclovir',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$gatifloxacin <- as.factor(ifelse(grepl('gatifloxacin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$lansoprazole <- as.factor(ifelse(grepl('lansoprazole',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$mercaptopurine <- as.factor(ifelse(grepl('mercaptopurine',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$minoxidil <- as.factor(ifelse(grepl('minoxidil',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$mupirocin <- as.factor(ifelse(grepl('mupirocin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$oxandrolone <- as.factor(ifelse(grepl('oxandrolone',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$propylthiouracil <- as.factor(ifelse(grepl('propylthiouracil',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$rifabutin <- as.factor(ifelse(grepl('rifabutin',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$suboxone <- as.factor(ifelse(grepl('suboxone',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$sulfacetamide <- as.factor(ifelse(grepl('sulfacetamide',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))
epilepsy_drug_abs2$topotecan <- as.factor(ifelse(grepl('topotecan',epilepsy_drug_abs2$epilepsy_drug_abs),'TRUE','FALSE'))


#check dimensions
dim(epilepsy_drug_abs2)
head(epilepsy_drug_abs2,1)
names(epilepsy_drug_abs2)


#Add average sentiment analysis
#Add a column with clean text
epilepsy_drug_abs2$sentiment<- sentiment_by(epilepsy_drug_abs2$epilepsy_drug_abs)$ave_sentiment 
drugname<-c("gabapentin",
            "levetiracetam",
            "topiramate",
            "lamotrigine",
            "acetazolamide",
            "analgesic",
            "riluzole",
            "zonisamide",
            "propofol",
            "ethosuximide",
            "lidocaine",
            "fluoxetine",
            "duloxetine",
            "haloperidol",
            "tizanidine",
            "acetaminophen",
            "venlafaxine",
            "adenosine",
            "ibuprofen",
            "risperidone",
            "phosphate",
            "olanzapine",
            "aripiprazole",
            "naproxen",
            "antibacterial",
            "celecoxib",
            "nicotine",
            "indomethacin",
            "temozolomide",
            "furosemide",
            "sumatriptan",
            "testosterone",
            "bupropion",
            "etomidate",
            "nifedipine",
            "modafinil",
            "temazepam",
            "ceftriaxone",
            "diphenhydramine",
            "pioglitazone",
            "escitalopram",
            "ondansetron",
            "simvastatin",
            "ketoconazole",
            "paroxetine",
            "cimetidine",
            "mirtazapine",
            "paclitaxel",
            "amantadine",
            "fenofibrate",
            "zaleplon",
            "ciprofloxacin",
            "erythromycin",
            "meloxicam",
            "ropinirole",
            "hydrochlorothiazide",
            "oxaliplatin",
            "aspirin",
            "hydrocortisone",
            "isoniazid",
            "azithromycin",
            "cefepime",
            "cytarabine",
            "decitabine",
            "sirolimus",
            "doxycycline",
            "epinephrine",
            "eszopiclone",
            "letrozole",
            "carboplatin",
            "carisoprodol",
            "clarithromycin",
            "clopidogrel",
            "menthol",
            "methazolamide",
            "ofloxacin",
            "rifampin",
            "telmisartan",
            "atomoxetine",
            "gemcitabine",
            "metronidazole",
            "acyclovir",
            "fluorouracil",
            "itraconazole",
            "tetrabenazine",
            "ampicillin",
            "disposable",
            "levofloxacin",
            "valsartan",
            "voriconazole",
            "omeprazole",
            "calcitriol",
            "cefazolin",
            "cortisone",
            "fluconazole",
            "linezolid",
            "lovastatin",
            "nevirapine",
            "oxybutynin",
            "bexarotene",
            "cilostazol",
            "clotrimazole",
            "emilia",
            "expectorant",
            "guanfacine",
            "hemorrhoidal",
            "laxative",
            "leflunomide",
            "lisinopril",
            "loratadine",
            "misoprostol",
            "piroxicam",
            "ribavirin",
            "almotriptan",
            "bacitracin",
            "budesonide",
            "chlorzoxazone",
            "cyanocobalamine",
            "ganciclovir",
            "gatifloxacin",
            "lansoprazole",
            "mercaptopurine",
            "minoxidil",
            "mupirocin",
            "oxandrolone",
            "propylthiouracil",
            "rifabutin",
            "suboxone",
            "sulfacetamide",
            "topotecan")
avg_sentiment<-c(
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$gabapentin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$levetiracetam=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$topiramate=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$lamotrigine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$acetazolamide=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$analgesic=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$riluzole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$zonisamide=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$propofol=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$ethosuximide=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$lidocaine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$fluoxetine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$duloxetine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$haloperidol=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$tizanidine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$acetaminophen=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$venlafaxine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$adenosine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$ibuprofen=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$risperidone=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$phosphate=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$olanzapine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$aripiprazole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$naproxen=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$antibacterial=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$celecoxib=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$nicotine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$indomethacin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$temozolomide=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$furosemide=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$sumatriptan=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$testosterone=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$bupropion=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$etomidate=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$nifedipine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$modafinil=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$temazepam=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$ceftriaxone=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$diphenhydramine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$pioglitazone=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$escitalopram=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$ondansetron=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$simvastatin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$ketoconazole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$paroxetine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$cimetidine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$mirtazapine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$paclitaxel=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$amantadine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$fenofibrate=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$zaleplon=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$ciprofloxacin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$erythromycin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$meloxicam=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$ropinirole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$hydrochlorothiazide=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$oxaliplatin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$aspirin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$hydrocortisone=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$isoniazid=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$azithromycin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$cefepime=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$cytarabine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$decitabine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$sirolimus=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$doxycycline=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$epinephrine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$eszopiclone=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$letrozole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$carboplatin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$carisoprodol=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$clarithromycin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$clopidogrel=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$menthol=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$methazolamide=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$ofloxacin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$rifampin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$telmisartan=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$atomoxetine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$gemcitabine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$metronidazole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$acyclovir=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$fluorouracil=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$itraconazole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$tetrabenazine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$ampicillin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$disposable=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$levofloxacin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$valsartan=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$voriconazole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$omeprazole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$calcitriol=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$cefazolin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$cortisone=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$fluconazole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$linezolid=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$lovastatin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$nevirapine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$oxybutynin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$bexarotene=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$cilostazol=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$clotrimazole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$emilia=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$expectorant=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$guanfacine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$hemorrhoidal=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$laxative=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$leflunomide=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$lisinopril=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$loratadine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$misoprostol=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$piroxicam=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$ribavirin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$almotriptan=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$bacitracin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$budesonide=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$chlorzoxazone=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$cyanocobalamine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$ganciclovir=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$gatifloxacin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$lansoprazole=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$mercaptopurine=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$minoxidil=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$mupirocin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$oxandrolone=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$propylthiouracil=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$rifabutin=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$suboxone=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$sulfacetamide=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE)),
  epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$topotecan=="TRUE") %>% summarize(Mean = mean(sentiment, na.rm=TRUE))
 )
drug_avg_sent <- data.frame(cbind(drugname, avg_sentiment))

#####topic modelling by drug
#gabapentin
sub<-epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$gabapentin=="TRUE")
sen_epi<-CleanDataframe(sub[1])
epi_drug_topics<-show_em_topics(sen_epi$clean.text,2)
epi_drug_topics

#levetiracetam
sub<-epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$levetiracetam=="TRUE")
sen_epi<-CleanDataframe(sub[1])
epi_drug_topics<-show_em_topics(sen_epi$clean.text,2)
epi_drug_topics

#topiramate
sub<-epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$topiramate=="TRUE")
sen_epi<-CleanDataframe(sub[1])
epi_drug_topics<-show_em_topics(sen_epi$clean.text,4)
epi_drug_topics

#lamotrigine
sub<-epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$lamotrigine=="TRUE")
sen_epi<-CleanDataframe(sub[1])
epi_drug_topics<-show_em_topics(sen_epi$clean.text,4)
epi_drug_topics

#acetazolamide
sub<-epilepsy_drug_abs2 %>% filter(epilepsy_drug_abs2$acetazolamide=="TRUE")
sen_epi<-CleanDataframe(sub[1])
epi_drug_topics<-show_em_topics(sen_epi$clean.text,3)
epi_drug_topics

####################################functions#####################################
# Topic models fucntion
show_em_topics  <- function(final_txt, count.topic) 
{
  # Create DTM
  myCorpus <- Corpus(VectorSource(final_txt))
  DTM <- DocumentTermMatrix(myCorpus)
  rowTotals <- apply(DTM , 1, sum) #Find the sum of words in each Document
  DTM   <- DTM[rowTotals > 0, ] 
  #Compute word frequencies
  freq <- sort(colSums(as.matrix(DTM)), decreasing = TRUE) 
  wf <- data.frame(word = names(freq), 
                   freq = freq)   
  # Set parameters for Gibbs sampling
  burnin <- 4000
  iter <- 2000
  thin <- 500
  seed <- list(2003, 5, 63, 100001, 765)
  nstart <- 5
  best <- TRUE
  #Number of topics
  #Run LDA using Gibbs sampling
  ldaOut <- LDA(DTM,as.numeric(count.topic))
  #Get top 5 words for each topic
  top_terms <- as.data.frame(terms(ldaOut,15))
  #Return the top topics 
  return(top_terms)
}

#Sentiment Polarity Function
SentimentPolarity <- function(sentiment) {    # Required for CleanDataframe & ConvertDataframe
  if (sentiment > 0) {
    return("positive")
  } else if (sentiment == 0) {
    return("neutral")
  } else {
    return("negative")
  }
}

#clean the text fucntion
CleanText <- function(some_txt) {      # Clean text function
  # remove html links
  some_txt <- gsub("http\\S+\\s*",
                   "",
                   some_txt)
  # remove retweet entities
  some_txt <- gsub("(RT|via)((?:\\b\\W*@\\w+)+)",
                   "",
                   some_txt)
  # remove at people
  some_txt <- gsub("@\\w+", 
                   "", 
                   some_txt)
  
  try.error <- function(x) {                    # "tolower error handling" function
    # create missing value
    y <- NA
    # tryCatch error
    try_error <- tryCatch(tolower(x), 
                          error = function(e) e)
    # if not an error
    if ( !inherits(try_error, "error") ) {
      y <- tolower(x)
    }
    # result
    return(y)
  }
  # lower case using try.error with sapply 
  some_txt <- sapply(some_txt, try.error)
  # remove NAs in some_txt
  some_txt <- some_txt[ !is.na(some_txt) ]
  names(some_txt) <- NULL
  myCorpus <- Corpus(VectorSource(some_txt))
  myCorpus <- tm_map(myCorpus, content_transformer(tolower))
  myCorpus <- tm_map(myCorpus, removePunctuation)
  myCorpus <- tm_map(myCorpus, content_transformer(strip), char.keep = ".")    # Keep period
  myCorpus <- tm_map(myCorpus, removeNumbers)
  #Add words to be excluded from the list of stop words here
  exceptions <- c("not","nor","neither","never")
  my_stopwords <- setdiff(stopwords("en"), exceptions)
  myCorpus <- tm_map(myCorpus, removeWords, my_stopwords)
  #myCorpus <- tm_map(myCorpus, stemDocument)
  some_txt_clean <- as.character(unlist(sapply(myCorpus, `[`, "content")))
  # remove trailing/leading spaces
  some_txt_clean <- str_trim(some_txt_clean)
  return(some_txt_clean)
}

#attach sentiment score to dataframe
CleanDataframe <- function(text.dataframe) {    # Clean text.dataframe and attach sentiment value
  text.dataframe$clean.text <- CleanText(text.dataframe[1,])
  text.dataframe$sentiment <- sentiment_by(text.dataframe$clean.text)$ave_sentiment
  text.dataframe$sentiment <- text.dataframe$sentiment
  text.dataframe$sentiment.pol <- lapply(text.dataframe$sentiment, SentimentPolarity)
  return(text.dataframe)
}

#################################################################################