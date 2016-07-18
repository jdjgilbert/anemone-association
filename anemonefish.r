#Anemone associations

#   1) Did anemone associations evolve numerous times (if so, how many)?
#   2) Does phylogenetic relatedness, or adult max body size best explain the life history stage at which anemone-associated fishes are associated with anemones?
#   3) Is there any real difference in these analysis when comparing facultative and obligate anemone associated fishes?
#   
#   3 is less important, as all we're trying to formally show is whether there's a good and broad case for size-selective predation being a generalisable selection force for anemone-associations in fishes.


library(ape)
library(phytools)
library(caper)

nrow(dat<-read.csv('data/Anemonefish Comparative Data 22.06.2016.csv'))
#[1] 87
names(dat)
# [1] "Common.name"                     "Family"                         
# [3] "Scientific.name"                 "Location"                       
# [5] "Association"                     "Ontogentic.Stage.of.Association"
# [7] "Adult.Max.Length..mm."           "COI"  

#Ntip(phy<-read.tree('COI\ genbank\ alignment/tree.phy'))
#Ntip(phy<-read.tree('phylogenies/20160711_averaged_jModeltest.phy'))
Ntip(phy<-read.tree('phylogenies/20160714_averaged_jModeltest_withnames.phy'))

# [1] 64

### phy has a basal polytomy, so try creating an outgroup (worked)
plot(outg <- read.tree(text='(outgroup:0.05,TREE:0.00);'))
phy1 <- bind.tree(outg, phy, where=grep('TREE', outg$tip.label))
phy2 <- root(phy1, 'outgroup')

### Need to match genbank IDs in phy tip labels to actual species names
#   head(trans <- read.csv('COI\ genbank\ alignment/genbank_spp_translate.txt'))
#   ##  genbank                       sp
#   ##  1 gi|584295226|gb|KF929610.1|      Apogon aunerolineatus
#   ##  2 gi|614136154|gb|KJ466130.1|   Ostorhinchus cyanosoma
#   ##  3 gi|356996146|gb|JN827969.1|         Apogon maculatus
#   ##  4 gi|356996320|gb|JN828056.1|   Apogon quadrisquamatus
#   ##  5 gi|584296689|gb|KF930338.1|      Pterapogon kauderni
#   ##  6 gi|816375469|gb|KP194051.1| Plagiotremus tapeinosoma
#   
#   trans$gb <- substr(trans$genbank, 1, 10)
#   
#   head(intersect(trans$genbank <- substr(trans$genbank, 1, 10), phy2$tip.label))
#     
#   head(phy2$tip.label <- as.character(trans$sp[match(phy2$tip.label, trans$genbank)]))
#   
#   phy2$tip.label[1]<-'outgroup' 
  
dd <- subset(dat, select=c('Scientific.name','Family','Ontogentic.Stage.of.Association','Adult.Max.Length..mm.'))
names(dd)<-c('sp','family','stage','length')

dd$sp <- as.character(dd$sp)
dd$sp.trans <- paste(substr(dd$sp, 1, 2), '.', sub('^.* ','',dd$sp), sep='')

levels(dd$stage)[1] <- 'Blank (presume no data?)'

phy2$tip.label <- sub('Os.cyanosoma','Os.nanus', phy2$tip.label)

phy2$tip.label <- dd$sp[match(phy2$tip.label, dd$sp.trans)]
phy2$tip.label <- as.character(phy2$tip.label)
phy2$tip.label[1]<-'outgroup'
treedat <- comparative.data(phy2, dd, sp)


summary(pgls1 <- pgls(length ~ stage, data=treedat, lambda='ML'))
##  Call:
##    pgls(formula = length ~ stage, data = treedat, lambda = "ML")
##  
##  Residuals:
##    Min       1Q   Median       3Q      Max 
##  -108.265  -23.871    1.031   19.396  120.312 
##  
##  Branch length transformations:
##    
##    kappa  [Fix]  : 1.000
##  lambda [ ML]  : 0.637
##  lower bound : 0.000, p = 0.00048113
##  upper bound : 1.000, p = < 2.22e-16
##  95.0% CI   : (0.288, 0.849)
##  delta  [Fix]  : 1.000
##  
##  Coefficients:
##    Estimate Std. Error t value  Pr(>|t|)    
##  (Intercept)      151.8956     3.7715 40.2751 < 2.2e-16 ***
##    stageJuvAndAdult -23.1969    28.0817 -0.8260     0.412    
##  stageJuvenile    181.4165    34.1788  5.3079 1.696e-06 ***
##    ---
##    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##  
##  Residual standard error: 40.79 on 60 degrees of freedom
##  Multiple R-squared: 0.4018,	Adjusted R-squared: 0.3819 
##  F-statistic: 20.15 on 2 and 60 DF,  p-value: 2.02e-07  
  
pdf(paste('plots/','Length vs Stage', '-', Sys.Date(),'.pdf',sep=''), height=8, width=8)
with(treedat$data, boxplot(length~stage, las=1, xlab='Stage of association', ylab='Adult max length'))
plot(jitter(resid(pgls1))~jitter(fitted(pgls1)), las=1, main='PGLS residuals versus fitted values')
dev.off()


