# Analysis of genetic differentiation for megafauna-dependent species.
# Combining Duminil et al. 2007 and our own compiliation for brazilian species.
# Duminil, J., S. Fineschi, A. Hampe, P. Jordano, D. Salvini, G. G. Vendramin, and R. J. Petit. 2007. Can population genetic structure be predicted from life-history traits? The American Naturalist 169:662–672.
#
# We now use Rosanne Collevatti's dataset, updated August 2014, with many more
# brazilian species.
# Data input.
# NOTE: gst as for 3Mar mixes Gst and Fst data and different nuclear markers (SSR, isozymes, etc.).
require(dplyr)
gst <-read.table("./datasets/gst.txt", header=TRUE, sep="\t", dec=".", na.strings="NA")
str(gst)
rosanne <-read.table("./datasets/rosanne.txt", header=TRUE, sep="\t", dec=".", na.strings="NA")
str(rosanne)
# Differences among all syndromes
m1<- lm(asin(gst$gst_nr)~ gst$syndr)
summary(m1)

# Zoochoric comparison.
table(rosanne$anacron)
md  mi   n 
62  23 350 

table(rosanne$dispersal)
autochory         bat       birds hidrochory      mammals       mixed 
23           2          71          17         163           3 
wind 
156 

table(rosanne$dispcode)
au bat  bi  hi  ma mix  wi 
23   2  71  17 163   3 156


m1<- lm(asin(rosanne$FST)~ rosanne$anacron)
summary(m1)

m2<- lm(asin(rosanne$mantel)~ rosanne$anacron)
summary(m2)


# Box plots - Rosanne Dataset
p <- ggplot(rosanne, aes(factor(anacron), FST, na.rm = TRUE))
p + geom_boxplot() +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation")

ggplot(rosanne, aes(factor(anacron), FST, fill=factor(anacron))) +
    #    stat_boxplot(geom ='errorbar')+
    geom_boxplot()+
    scale_fill_manual(name= "Syndrome", 
        values=c("burlywood","lightblue","blue"), 
        labels = c("md"= "Megafauna dependent", "mi"= "Megafauna independent",
            "n"="Other")) +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation, Fst")


# Megafauna dependent, megafauna independent, extant frugivores
ppz<- rosanne %.% 
      filter(anacron=="md" |anacron=="mi") %.%
      select(anacron, FST)
m3<- lm(asin(ppz$FST)~ ppz$anacron)
summary(m3)
anova(m3)

# Box plots
p <- ggplot(ppz, aes(factor(anacron), FST))
p + geom_boxplot()

p <- ggplot(ppz[ppz$anacron!="NA",], 
    aes(factor(anacron), FST, fill=factor(anacron)))
p + geom_boxplot() +
    scale_fill_manual(name= "Syndrome", 
        values=c("blue","lightblue")) +
    xlab("Seed dispersal syndrome") +
    ylab("Population genetic differentiation, Fst")

# Just the zoochoric taxa
He<- rosanne %.%
     filter(anacron=="md" | anacron=="mi") %.%
     select(anacron, He)

m4<- lm(asin(He$He)~ He$anacron)
summary(m4)
anova(m4)

# Zoochoric comparison.





