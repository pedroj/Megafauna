# Analysis of genetic differentiation for megafauna-dependent species.
# Combining Duminil et al. 2007 and our own compiliation for brazilian species.
# Duminil, J., S. Fineschi, A. Hampe, P. Jordano, D. Salvini, G. G. Vendramin, and R. J. Petit. 2007. Can population genetic structure be predicted from life-history traits? The American Naturalist 169:662–672.
#
# Data input.
# NOTE: gst as for 3Mar mixes Gst and Fst data and different nuclear markers (SSR, isozymes, etc.).
library(dplyr)
gst <-read.table("./datasets/gst.txt", header=TRUE, sep="\t", dec=".", na.strings="NA")
str(gst)

# Differences among all syndromes
m1<- lm(asin(gst$gst_nr)~ gst$syndr)
summary(m1)

# Box plots
p <- ggplot(gst, aes(factor(syndr), gst_nr, na.rm = TRUE))
p + geom_boxplot() +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

ggplot(gst, aes(factor(syndr), gst_nr, fill=factor(syndr))) +
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

# Megafauna dependent, megafauna independent, extant frugivores
pp1<- gst %.%
    filter(contr2=="Megafauna dependent" |
           contr2=="Megafauna independent" |
           contr2== "Zoochoric") %.%
    select(contr2, gst_nr)
m3<- lm(asin(pp1$gst_nr)~ pp1$contr2)
summary(m3)
anova(m3)

# Box plots
p <- ggplot(pp1, aes(factor(contr2), gst_nr))
p + geom_boxplot()

p <- ggplot(pp1[pp1$contr2!="NA",], 
    aes(factor(contr2), gst_nr, fill=factor(contr2)))
p + geom_boxplot() +
    scale_fill_manual(name= "Syndrome", 
        values=c("darkblue","blue","lightblue")) +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

# Just the zoochoric taxa
pp<- gst %.%
     filter(contr1=="Megafauna" |
            contr1== "Zoochoric") %.%
     select(contr1, gst_nr)

m4<- lm(asin(pp$gst_nr)~ pp$contr1)
summary(m4)
anova(m4)

p <- ggplot(pp[pp$contr1!="NA",], 
            aes(factor(contr1), gst_nr, fill=factor(contr1)))
p + geom_boxplot() +
    scale_fill_manual(name= "Syndrome", 
        values=c("lightblue","burlywood")) +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")






