# Analysis of genetic differentiation for megafauna-dependent species.
# Combining Duminil et al. 2007 and our own compiliation for brazilian species.
# Duminil, J., S. Fineschi, A. Hampe, P. Jordano, D. Salvini, G. G. Vendramin, and R. J. Petit. 2007. Can population genetic structure be predicted from life-history traits? The American Naturalist 169:662–672.
#
library(dplyr)
#
# Data input. --------------------------------------------------------------
# POOLED ANALYSIS WITH ALL DATA.
# Combining Duminil et al. 2007 and our own compiliation for brazilian species. Added Rosanne's dataset.
# NOTE: genetdiff as for 3Mar mixes Gst and Fst data and different nuclear markers (SSR, isozymes, etc.).
gst_full <-read.table("./datasets/gst_full.txt", header=TRUE, sep="\t", dec=".", na.strings="NA")
str(gst_full)
'data.frame':	615 obs. of  24 variables:
    $ id             : Factor w/ 481 levels "10","100","105",..: 127...
$ idpaper        : int  248 248 248 248 248 248 2 23 23 15 ...
$ data           : Factor w/ 3 levels "duminil","galetti",..: 3 3 ...
$ fam            : Factor w/ 82 levels "Acanthaceae",..: 1 1 1 1 ...
$ sp             : Factor w/ 369 levels "Abies alba","Abies fraseri",..: 32
$ spcode         : Factor w/ 225 levels "ACACUL","ACCAVE",..: 18 18 ...
$ marker         : Factor w/ 14 levels "AFLP","cAFLP",..: 9 3 14 9 3 14 ...
$ lifeform       : Factor w/ 8 levels "Cact","Creeper",..: 8 8 8 8 8 ...
$ dispersal      : Factor w/ 8 levels "ants","autochory",..: 5 5 5 5 ...
$ dispcode       : Factor w/ 8 levels "ant","au","bi",..: 5 5 5 ...
$ anacron        : Factor w/ 5 levels "ab","b","md",..: 1 1 1 1 ...
$ pollinator     : Factor w/ 12 levels "bat","beetle",..: 10 10 10 10 ...
$ pollicode      : Factor w/ 11 levels "bat","beet","bi",..: 9 9 9 ...
$ distgeo        : Factor w/ 2 levels "restric","wide": 2 2 2 2 2 2...
$ habitat        : Factor w/ 23 levels "A","AE","Andean",..: 10 10 ...
$ diverhaplot    : Factor w/ 59 levels "0.00000","0.01680",..:  NA NA ...
$ divernucleot   : num  NA NA NA NA NA NA ...
$ genetdiff      : num  0.57 NA 0.49 0.57 NA 0.49 NA 0.16 0.97 0.04 ...
$ AllelesLocus   : num  NA NA NA NA NA ...
$ Porcpoliloci   : num  NA NA NA NA NA NA NA NA NA 61.5 ...
$ AllelicRichness: num  NA NA 1.68 NA NA 3.77 1.45 6.46 NA NA ...
$ FIS            : num  NA NA 0.22 NA NA 0.01 0.317 0.077 NA -0.156 ...
$ phi_FT         : num  NA NA NA NA NA NA NA NA NA NA ...
$ seedwt         : num  NA NA NA NA NA NA NA NA NA NA ...

# Differences among all syndromes ------------------------------------------
m1<- lm(asin(gst_full$genetdiff)~ gst_full$dispersal)
summary(m1)

# Box plots
p <- ggplot(gst_full, aes(factor(dispersal), genetdiff, na.rm = TRUE))
p + geom_boxplot() +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

ggplot(gst_full, aes(factor(dispersal), genetdiff, fill=factor(dispersal))) +
    #    stat_boxplot(geom ='errorbar')+
    geom_boxplot()+
    scale_fill_manual(name= "Syndrome", 
        values=c("burlywood","lightblue",
            "lightblue","burlywood",
            "blue","lightblue",
            "burlywood"), 
        labels = c("autochory"= "Autochory", "birds"= "Birds",
            "epizoochory"="Epizoochory", "hidrochory"="Hidrochory",
            "mammals"= "Mammals", "mixed"= "Birds and Mammals", 
            "wind"= "Anemochory")) +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

# Megafauna dependent, megafauna independent, extant frugivores ------------
ppz<- gst_full %>% 
               dplyr::filter(anacron=="md" |
                      anacron=="mi"|
                      anacron=="b") %>% 
               dplyr::select(anacron, genetdiff)

m3<- lm(asin(ppz$genetdiff)~ ppz$anacron)
summary(m3)
anova(m3)

# Box plots. Just the zoochoric taxa.
p <- ggplot(ppz, aes(factor(anacron), genetdiff))
p + geom_boxplot()

p <- ggplot(ppz[ppz$anacron!="NA",], 
    aes(factor(anacron), genetdiff, fill=factor(anacron)))
p + geom_boxplot() +
    scale_fill_manual(name= "Syndrome", 
        values=c("darkblue","blue","lightblue")) +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

# Tables -------------------------------------------------------------------
table(gst_full$dispersal)
table(gst_full$dispcode)
table(gst_full$anacron)

# Differences among all syndromes ------------------------------------------
library(contrast)
m1<- lm(asin(gst_full$genetdiff)~ gst_full$dispersal)
summary(m1)

m2<- lm(asin(gst_full$genetdiff)~ gst_full$dispcode)
summary(m2)

m3<- lm(asin(gst_full$genetdiff)~ gst_full$anacron)
summary(m3)

# Box plots ----------------------------------------------------------------
p <- ggplot(gst_full, aes(factor(dispersal), genetdiff, na.rm = TRUE))
p + geom_boxplot() +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

p <- ggplot(gst_full, aes(factor(anacron), genetdiff, na.rm = TRUE))
p + geom_boxplot() +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

ggplot(gst_full, aes(factor(anacron), genetdiff, fill=factor(anacron))) +
    #    stat_boxplot(geom ='errorbar')+
    geom_boxplot()+
    scale_fill_manual(name= "Syndrome", 
        values=c("burlywood","lightblue",
            "lightblue","burlywood",
            "blue","burlywood",
            "burlywood"), 
        labels = c("ab"= "Abiotic", "b"= "Biotic",
            "md"="Megafauna dep.", mi="Megafauna indep.")) +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

# SUBSETS ------------------------------------------------------------------
# Only Trees.
#
str(gst_full)
test1 <- gst_full %>%
    dplyr::filter(gst_full$lifeform=="Tree")
m1<- lm(asin(test1$genetdiff)~ test1$dispersal)
summary(m1)

m2<- lm(asin(test1$genetdiff)~ test1$dispcode)
summary(m2)

m3<- lm(asin(test1$genetdiff)~ test1$anacron)
summary(m3)

ggplot(test1, aes(factor(anacron), genetdiff, fill=factor(anacron))) +
    #    stat_boxplot(geom ='errorbar')+
    geom_boxplot()+
    scale_fill_manual(name= "Syndrome", 
        values=c("burlywood","lightblue",
            "lightblue","burlywood",
            "blue","burlywood",
            "burlywood"), 
        labels = c("ab"= "Abiotic", "b"= "Biotic",
            "md"="Megafauna dep.", mi="Megafauna indep.")) +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

test1 %>% 
    group_by(anacron) %>%
    summarise(avg= mean(genetdiff, na.rm= T), 
        SD= sd(genetdiff, na.rm= T),
        min= min(genetdiff, na.rm= T), 
        max= max(genetdiff, na.rm= T),
        N= n())

# Megafauna dependent vs megafauna independent plus extant frugivores.------
# Create new column grouping extant frugivores (Birds + Mammals + Others).
gst_full$frugiv[gst_full$anacron=="b" | 
                gst_full$anacron=="mi"] <- "bmi"
gst_full$frugiv[gst_full$anacron=="md"] <- "md"
ppz<- gst_full %>%
    dplyr::filter(anacron!="ab")  %>%
    dplyr::select(lifeform, frugiv, genetdiff) %>% 
    dplyr::filter(lifeform=="Tree")
m3<- lm(asin(ppz$genetdiff)~ ppz$frugiv)
summary(m3)
anova(m3)

# Box plots
p <- ggplot(ppz, aes(factor(frugiv), genetdiff))
p + geom_boxplot()

ggplot(ppz[ppz$frugiv!="NA",], 
    aes(factor(frugiv), genetdiff, fill=factor(frugiv))) + 
    geom_boxplot() +
    scale_fill_manual(name= "Syndrome", 
        values=c("burlywood","lightblue"), 
        labels = c("md"="Megafauna dep.", bmi="Birds & Mammals")) +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

# Just Rosanne's dataset. --------------------------------------------------
#
# Data input.
# NOTE: gst as for 3Mar mixes Gst and Fst data and different nuclear markers (SSR, isozymes, etc.).
library(dplyr)
rosanne <-read.table("./datasets/megafauna-genetics datasets_Rosane/rosanne.txt", header=TRUE, sep="\t", dec=".", na.strings="NA")
str(rosanne)

# Differences among all syndromes
m1<- lm(asin(rosanne$FST)~ rosanne$disp)
summary(m1)

m2<- lm(asin(rosanne$FST)~ rosanne$anacron)
summary(m2)


# Box plots
p <- ggplot(rosanne, aes(factor(disp), FST, na.rm = TRUE))
p + geom_boxplot() +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

ggplot(rosanne, aes(factor(anacron), FST, fill=factor(anacron))) +
    #    stat_boxplot(geom ='errorbar')+
    geom_boxplot()+
    scale_fill_manual(name= "Syndrome", 
        values=c("burlywood","lightblue",
            "lightblue","burlywood",
            "blue","burlywood",
            "burlywood"), 
        labels = c("AA"= "Epizoochory", "AC"= "Synzoochory",
            "AI"="Endozoochory", G="Barochory",
            "M"= "Megafauna", "W"= "Anemochory")) +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

#---------------------------------------------------------------------------

# dplyr NOTES ==============================================================
# 




