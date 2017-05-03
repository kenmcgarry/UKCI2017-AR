# adr_comorbid.r     27/04/2017

library(igraph)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)
library(xtable)
library(arulesViz)
library(arules)
library(DOSE)

setwd("C:/R-files/FDA-ADR")  # for writing and reading files

INDQ1 <- file.path('C://R-files//FDA-ADR','INDI16Q1.txt') %>% read.delim(na.strings='',sep = "$",header=T)
INDQ2 <- file.path('C://R-files//FDA-ADR','INDI16Q2.txt') %>% read.delim(na.strings='',sep = "$",header=T)
INDQ3 <- file.path('C://R-files//FDA-ADR','INDI16Q3.txt') %>% read.delim(na.strings='',sep = "$",header=T)
INDQ4 <- file.path('C://R-files//FDA-ADR','INDI16Q4.txt') %>% read.delim(na.strings='',sep = "$",header=T)

DRUG1 <- file.path('C://R-files//FDA-ADR','DRUG16Q1.txt') %>% read.delim(na.strings='',sep = "$",header=T)
DRUG2 <- file.path('C://R-files//FDA-ADR','DRUG16Q2.txt') %>% read.delim(na.strings='',sep = "$",header=T)
DRUG3 <- file.path('C://R-files//FDA-ADR','DRUG16Q3.txt') %>% read.delim(na.strings='',sep = "$",header=T)
DRUG4 <- file.path('C://R-files//FDA-ADR','DRUG16Q4.txt') %>% read.delim(na.strings='',sep = "$",header=T)

DEMO1 <- file.path('C://R-files//FDA-ADR','DEMO16Q1.txt') %>% read.delim(na.strings='',sep = "$",header=T)
DEMO2 <- file.path('C://R-files//FDA-ADR','DEMO16Q2.txt') %>% read.delim(na.strings='',sep = "$",header=T)
DEMO3 <- file.path('C://R-files//FDA-ADR','DEMO16Q3.txt') %>% read.delim(na.strings='',sep = "$",header=T)
DEMO4 <- file.path('C://R-files//FDA-ADR','DEMO16Q4.txt') %>% read.delim(na.strings='',sep = "$",header=T)

REACT1 <- file.path('C://R-files//FDA-ADR','REAC16Q1.txt') %>% read.delim(na.strings='',sep = "$",header=T)
REACT2 <- file.path('C://R-files//FDA-ADR','REAC16Q2.txt') %>% read.delim(na.strings='',sep = "$",header=T)
REACT3 <- file.path('C://R-files//FDA-ADR','REAC16Q3.txt') %>% read.delim(na.strings='',sep = "$",header=T)
REACT4 <- file.path('C://R-files//FDA-ADR','REAC16Q4.txt') %>% read.delim(na.strings='',sep = "$",header=T)

IND2016 <- rbind(INDQ1,INDQ2,INDQ3,INDQ4)
DRUG2016 <- rbind(DRUG1,DRUG2,DRUG3,DRUG4)
DEMO2016 <- rbind(DEMO1,DEMO2,DEMO3,DEMO4)
REAC2016 <- rbind(REACT1,REACT2,REACT3,REACT4)

rm(REACT1,REACT2,REACT3,REACT4); rm(DEMO1,DEMO2,DEMO3,DEMO4);rm(DRUG1,DRUG2,DRUG3,DRUG4);rm(INDQ1,INDQ2,INDQ3,INDQ4)
#ontargets <- read.csv(file='C://R-files//sider//drugbank-proteins.tsv', header=TRUE, sep="\t")

# order dataframe
df_IND2016 <- IND2016[order(IND2016$primaryid),]
#convert member number to numeric
df_IND2016$primaryid <- as.numeric(df_IND2016$primaryid)

#convert item description to categorical format
df_IND2016$indi_pt <- as.factor(df_IND2016$indi_pt)

# group indications by same patient
df_itemlist <- ddply(df_IND2016, "primaryid", function(df1)paste(df1$indi_pt,collapse = ","))
#remove member number and date
df_itemlist$primaryid <- NULL
colnames(df_itemlist) <- c("itemlist")

# write to csv format
write.csv(df_itemlist,"indication-items.csv", quote = FALSE, row.names = TRUE)

# convert data in csv file to the basket format expected by arules: took 25 minutes for 3.1 million transactions
trans = read.transactions(file="indication-items.csv", rm.duplicates= FALSE, format="basket",sep=",",cols=1);

#remove quotes from transactions
trans@itemInfo$labels <- gsub("\"","",trans@itemInfo$labels)

#run apriori algorithm
comorbid_rules <- apriori(trans,parameter= list(minlen=2,sup = 0.0001, conf = 0.01, target="rules"))
comorbid_rules <- apriori(trans,parameter= list(minlen=2,sup = 0.001, conf = 0.01, target="rules"),appearance = list(lhs = "Atrial fibrillation"))

#view rules
inspect(comorbid_rules)

#convert to datframe and view; optional
df_comorbid <- as(comorbid_rules,"data.frame")
df_comorbid$confidence <- df_comorbid$confidence * 100
df_comorbid$support <- df_comorbid$support * nrow(df_comorbid)

write.csv(df_comorbid,"comorbid-rules.csv",row.names = FALSE)

#plot the rules
plot(comorbid_rules)
set.seed(8000)
plot(comorbid_rules[1:120], method = "grouped", control = list(k = 5))
plot(comorbid_rules[1:100,], method="graph", control=list(type="items"))
plot(comorbid_rules[1:50],method="graph",interactive=TRUE,shading=NA)
plot(comorbid_rules, method="paracoord",  control=list(alpha=.5, reorder=TRUE))
itemFrequencyPlot(trans, topN = 10)

trans.unknown <- trans[trans %in%  "Atrial fibrillation"]
trans.unknown <- head(sort(trans.unknown, decreasing = TRUE, na.last = NA,  by = "support"))

itemFrequencyPlot(trans.unknown[, 1:110],  population = trans[ ,1:110])
itemFrequencyPlot(trans.unknown[, 1:80],  population = trans[ ,1:80], confidence = 0.1, lift = TRUE, horiz = FALSE)

plot(trans.unknown[1:50],method="graph",interactive=TRUE,shading=NA)

itemsets <- apriori(trans.unknown, parameter = list(target = "rules", supp=0.001, minlen = 3, maxlen=4))
itemrule <- inspectDT(head(sort(itemsets), n=20))
plot(itemsets[1:50],method="graph",interactive=TRUE,shading=NA)

tli.table <- xtable(itemrule)

quality(itemsets)$lift <- interestMeasure(itemsets, measure="lift", trans = trans)
inspect(head(sort(itemsets, by = "lift"), n=50))

# Plot itemsets as a graph. Different subgroups with items that are related to each other can be identified.
plot(head(sort(itemsets, by = "lift"), n=20), method = "graph", control=list(cex=.8))

# xtable for LaTex file
tli.table <- xtable(itemrule[1:20,])
digits(tli.table)[c(3,4,5)] <- print(tli.table,floating=TRUE)

plot(comorbid_rules, shading="order", control = list(main = "Two-key plot",  col=rainbow(max(size(comorbid_rules))-1L)))

quality(itemsets) <- interestMeasure(itemsets, trans=trans)
head(quality(itemsets))
plot(itemsets, measure=c("support", "allConfidence"), shading="lift")

# starting to panic here....
newrules<-apriori(data=trans, parameter=list(supp=0.001,conf = 0.05), 
               appearance = list(default="lhs",rhs="Atrial fibrillation"),
               control = list(verbose=TRUE))
newrules<-sort(newrules, decreasing=TRUE,by="confidence")
inspect(newrules)

ruledf = data.frame(
  lhs = labels(lhs(newrules)),
  rhs = labels(rhs(newrules)), 
  newrules@quality)
tli.table <- xtable(ruledf)
digits(tli.table)[c(3,4,5)] <- print(tli.table,floating=TRUE)

# do it the other way
newrules<-apriori(data=trans, parameter=list(supp=0.001,conf = 0.05,minlen=2), 
               appearance = list(default="rhs",lhs="Atrial fibrillation"),
               control = list(verbose=TRUE))
newrules<-sort(newrules, decreasing=TRUE,by="confidence")
inspect(newrules)

ruledf = data.frame(
  lhs = labels(lhs(newrules)),
  rhs = labels(rhs(newrules)), 
  newrules@quality)
tli.table <- xtable(ruledf)
digits(tli.table)[c(3,4,5)] <- print(tli.table,floating=TRUE)

# now  do the semantic similarity
a <- c("DOID:6713","DOID:10763","DOID:0060224","DOID:11695")
names(a) <- c("cerebrovascular disease","hypertension","atrial fibrillation","portal vein thrombosis")

b <- a #c("DOID:9409", "DOID:2491", "DOID:4467", "DOID:3498", "DOID:11256")

# "DOID:6713" 	cerebrovascular disease
# "DOID:10763" 	hypertension
# "DOID:0060224" 	atrial fibrillation
# "DOID:11695" 	portal vein thrombosis
# 

s <- doSim(a, b, measure="Wang")
colnames(s)<-names(a)
rownames(s)<-names(a)
simplot(s,
        color.low="white", color.high="red", ylab=a, xlab=b,
        labs=TRUE, digits=2, labs.size=5,
        font.size=15)



