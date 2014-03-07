# Regression models. Body mass vs digestive variables. Megafauna.
# Pedro. Sevilla, 18 Feb 2014.
# I used the body mass (g) dataset from:
# Smith, F., S. Lyons, S. Ernest, K. Jones, D. Kaufman, T. Dayan, P. Marquet, J. Brown, and J. Haskell. 2003. Body mass of late quaternary mammals. Ecology 84:3403–3403.
# and from:
# Fariña, R. A., S. F. Vizcaíno, and M. S. Bargo. 2004. Body mass estimations in Lujanian (late Pleistocene-early Holocene of South America) mammal megafauna. Mastozoología Neotropical 5:87–108.
#
# Our objective is to fit gut capacity (wet weight of content, in g) and gut retention time (in h) to body mass across the full range of body mass for terrestrial mammals. Then we use both the confidence interval for slope and the prediction interval (95%) to estimate the value and range of expected gut capacity and gut retention time of extinct megafauna. 
# For these I used the regressions in van Soest (1996).
# Van Soest, P. J. 1982. Nutritional Ecology of the Ruminant. Cornell University Press, NY, USA.
# Van Soest, P. J. 1996. Allometry and ecology of feeding behavior and digestive capacity in herbivores: A review. Zoo Biology 15:455–479.
# Owen-Smith, R. N. 1988. Megaherbivores. Cambridge University Press, Cambridge, UK.

library(ggplot2)
library(Hmisc)
library(scales)
library(dplyr)
#
#-------------------------------------------------------------------------
# Body mass data for extinct taxa.
# Body mass (g) of Late Quaternary mammals MOM database (Smith et al 2003).
# Smith, F., S. Lyons, S. Ernest, K. Jones, D. Kaufman, T. Dayan, P. Marquet, J. Brown, and J. Haskell. 2003. Body mass of late quaternary mammals. Ecology 84:3403–3403.
# Body mass in g.
bmass_MOM <-read.table("./datasets/MOMv3.3.txt", header=TRUE, sep="\t",
    dec=".",na.strings="-999")

# Body mass range up to 15000 kg
extinct<- bmass_MOM %.%     # Selecting just the extinct taxa.
    filter(status=="extinct") %.%
    select(family, genus, species, meanbodymass)

# Selecting just the extinct frugivores.
extinct.frug.megafauna<- bmass_MOM %.% 
    filter(status=="extinct", 
           family== "Elephantidae" |
           family== "Equidae" |
           family== "Gomphotheriidae" |
           family== "Macraucheniidae" |
           family== "Mammutidae" | 
           family== "Megalonychidae" |
           family== "Megatheriidae" |
           family== "Mylodontidae") %.%
    select(family, genus, species, meanbodymass)
colnames(extinct.frug.megafauna)<- c("family", "genus", "species", "bmass")

# Selecting different genera with dplyr. *** Body masses in g ***.
eremoth<- bmass_MOM %.%
    filter(genus=="Eremotherium") %.%
    select(species, meanbodymass) %.%
    summarise(mean(meanbodymass))
toxodontids<- bmass_MOM %.%
    filter(family=="Toxodontidae") %.%
    select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass))
equids<- bmass_MOM %.%
    filter(family=="Equidae", status=="extinct") %.%
    select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass))
megatherids<- bmass_MOM %.%
    filter(family=="Megatheriidae", status=="extinct") %.%
    select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass))
gomphotherids<- bmass_MOM %.%
    filter(family=="Gomphotheriidae", status=="extinct") %.%
    select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass))
elephantids<- bmass_MOM %.%
    filter(family=="Elephantidae", status=="extinct") %.%
    select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass)) 
mylo<- bmass_MOM %.%
    filter(family=="Mylodontidae", status=="extinct") %.%
    select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass))
tags<- c("Eremotherium", "Toxodontidae", "Equidae", "Megatheriidae",
         "Gomphotheriidae", "Elephantidae", "Mylodontidae")
topredict<- data.frame(taxa=tags, 
                 bmass= rbind(eremoth, toxodontids, equids, 
                              megatherids, gomphotherids, 
                              elephantids, mylo))
colnames(topredict)<- c("taxa","bmass") # Body mass in g!

#-------------------------------------------------------------------------
# Body mass (kg) vs. digestive organ capacity (wet g).
# Dataset obtained from digitization of figure in van Soest (1982).
# Body mass in kg; digestive capacity in g wet content.
dig.organ <-read.table("./datasets/Body mass (kg) vs Gut capc (g wet).txt", 
                       header=TRUE,sep="\t",dec=".",na.strings=".")
str(dig.organ)
# dig.organ <- data.frame(cbind(logbmass=dig.organ$bmass,
#                               loggutcapac=dig.organ$gutcapac,
#                               bmass=10^dig.organ$bmass,
#                               gutcapac=10^dig.organ$gutcapac))
    
# Linear model with log-transformed values.
lm_fit0<- lm(log10(gutcapac) ~ log10(bmass), data= dig.organ)
summary(lm_fit0)

gut_with_pred<- data.frame(dig.organ,   # Predicted values for plotting.
    predict(lm_fit0, interval= 'prediction'))

# Numerical predictions of functional traits from body mass.
# We want predictions for the range of body masses between 1000 kg
# and 15000 kg for extinct megafauna; below this range we can use the 
# regression fitted to the extant species.
# I used the body mass data in Smith et al (2003) and Fariña et al. (2004).
# *** The data of body mass are in g ***.

# Datatset with body mass (kg) and estimated gut capacity (wet kg).
# Just the EXTINCT MEGAFAUNA FRUGIVOROUS species.
prediction<- data.frame(extinct.frug.megafauna, 
             predgutcap= 10^(2.0316 + 
                     1.0020*log10(extinct.frug.megafauna$bmass/1000)))
head(prediction)

# Dataset with predicted mean values by family.
predfammean<- data.frame(topredict, 
    predgutcap= 10^(2.0316 + 
            1.0020*log10(topredict$bmass/1000)))

#-------------------------------------------------------------------------
# PLOTS
# Without log scales
ggplot(dig.organ, aes(bmass, gutcapac)) + 
    stat_summary(fun.data= mean_cl_normal) + 
    scale_x_log10() + scale_y_log10() +
    geom_smooth(method= 'lm')

# Log10-scaled axes. Confidence interval for the regression slope.
p0<- ggplot(gut_with_pred, aes(x= log10(bmass), 
    y= log10(gutcapac))) + 
    geom_point(color="black", size=2) +
    geom_smooth(method= 'lm', aes(fill= 'confidence'), alpha= 0.5) +
    geom_ribbon(aes(y= fit, ymin= lwr, ymax= upr, fill= 'prediction'),
        alpha= 0.2) +
    scale_fill_manual('Interval', values= c('green', 'blue')) +
    xlab("log10(Body mass) (kg)") +
    ylab("log10(Gut capacity) (wet g)") +
    ggtitle("Gut capacity estimation for Megafauna frugivores") +
    coord_cartesian(xlim= c(-4.5,5), ylim= c(-3,7.5)) +
    theme(legend.position= c(0.20, 0.85))
p0 # The plot

p1<-p0 + geom_point(data= prediction, aes(x= log10(prediction$bmass/1000), 
                    y= log10(prediction$predgutcap/1000), 
                    color="white", size=2)) +
                    theme(legend.position = "none")

# Plots with averaged data for Extinct megafauna families
p0 + geom_point(data= predfammean, aes(x= log10(predfammean$bmass/1000), 
    y= log10(predfammean$predgutcap/1000))) + 
      geom_point(color="pink", size=1.5) +
    theme(legend.position = "none")
#-------------------------------------------------------------------------
# Body mass (kg) vs. mean retention time (h).

mean.retention <-read.table("./datasets/Body mass (kg)-Mean gut retention (h).txt",header=TRUE,sep="\t",dec=".",na.strings=".")
str(mean.retention)

# Linear model
#### Refitting with newly digitized data (7Mar2014).
# Regression equation for time and digestion-limited herbivores.
# Body mass (kg) vs. Mean gut retention time (h).
# Log10-scaled axes. Confidence interval for the regression slope.
# Body mass range up to 15000 kg
lm_fit1<- lm(log10(retention) ~ log10(bmass), data= mean.retention)
summary(lm_fit1)

# Datatset with body mass (kg) and estimated retention time (h).
# Using van Soest regression equation after I refitted it (lm_fit1).
# Just the EXTINCT MEGAFAUNA FRUGIVOROUS species.
predretfit1=  1.0193 + 
              0.2356 * log10(extinct.frug.megafauna$bmass)
prediction<- data.frame(prediction, 10^predretfit1)
head(prediction) # Now the dataset holds all variables.

ret_with_pred<- data.frame(mean.retention, 
    predict(lm_fit1, interval= 'prediction'))
p1<- ggplot(ret_with_pred, aes(x= log10(bmass), 
    y= log10(retention))) + 
    geom_point() +
    geom_smooth(method= 'lm', aes(fill= 'confidence'), alpha= 0.5) +
    geom_ribbon(aes(y= fit, ymin= lwr, ymax= upr, fill= 'prediction'),
        alpha= 0.2) +
    scale_fill_manual('Interval', values= c('green', 'blue')) +
    coord_cartesian(xlim= c(-4,5), ylim= c(-0.5,2.5)) +
    theme(legend.position= c(0.20, 0.85))
p1 # The plot

# New estimations
# Log10-scaled axes. Confidence interval for the regression slope.
p0<- ggplot(ret_with_pred, aes(x= log10(bmass), 
    y= log10(retention))) + 
    geom_smooth(method= 'lm', aes(fill= 'confidence'), alpha= 0.5) +
    geom_ribbon(aes(y= fit, ymin= lwr, ymax= upr, fill= 'prediction'),
        alpha= 0.2) +
    geom_point(color="white", size=4) +
    geom_point(color="black", size=2) +
    scale_fill_manual('Interval', values= c('green', 'blue')) +
    xlab("log10(Body mass) (g)") +
    ylab("log10(Gut retention time) (h)") +
    ggtitle("Gut retention time for Megafauna frugivores") +
    coord_cartesian(xlim= c(-4,5), ylim= c(-0.5,2.5)) +
    theme(legend.position= c(0.20, 0.85))
p0 # The plot

p1<-p0 + geom_point(data= prediction, aes(x= log10(prediction$bmass), 
    y= predretfit1-2.65)) 
p1

# Plots with averaged data for Extinct megafauna families
p1 + geom_point(data= predfammean, aes(x= log10(predfammean$bmass/1000), 
    y= log10(predfammean$predretfam)), size=4, color= "pink")  
  #  theme(legend.position = "none")

#--------------------------------------------------------------------------
# Datatset with body mass (g) and estimated gut capacity (wet g).
# Using van Soest regression equation.
# Regression equation for time and digestion-limited herbivores.
# gutret.h= 0.57 + 6.95*bmass.kg^0.33 + (0.33 + 0.43 * bmass^0.33)^0.5
# Predict specific values of gut retention from body mass data.
# Just the EXTINCT MEGAFAUNA FRUGIVOROUS species.
predret= 0.57 + 6.95 * ((extinct.frug.megafauna$bmass/1000)^0.33) + 
    (0.33 + 0.43 * (extinct.frug.megafauna$bmass/1000^0.33))^0.5
prediction<- data.frame(prediction, predret)
head(prediction) # Now the dataset holds all variables.

R> prediction
family            genus          species    bmass predgutcap predret
1          Equidae            Equus         capensis   350000   38085.19 172.702
2          Equidae        Hipparion          libycum   150000   16294.59 118.130
3     Elephantidae          Elephas         iolensis  6500000  711441.48 661.334
4             <NA>       Diprotodon            minor   900000   98118.52 265.169
5             <NA>       Diprotodon          optatum  1500000  163698.03 335.122
6             <NA>         Euowenia            grata   750000   81735.62 244.000
7             <NA>       Euryzygoma           dunese   500000   54446.25 202.929
8             <NA>      Nototherium        mitchelli   500000   54446.25 202.929
9             <NA>     Palorchestes            azeal   500000   54446.25 202.929
10            <NA>     Palorchestes           parvus   100000   10854.25  98.674
11            <NA>      Zygomaturus         trilobus   750000   81735.62 244.000
12            <NA> Plesiorycteropus germainepetterae     7810     843.41  32.812
13            <NA> Plesiorycteropus madagascariensis    13220    1429.14  40.988
14         Equidae            Equus          alaskae   372000   40484.06 177.519
15         Equidae            Equus       laurentius   648000   70598.94 228.287
16         Equidae            Equus         caballus   250000   27185.41 148.438
17         Equidae            Equus        fraternus   259000   28166.08 150.814
18         Equidae            Equus     niobrarensis   334000   36340.76 169.096
19         Equidae            Equus      complicatus   400000   43537.56 183.434
20         Equidae            Equus        giganteus   400000   43537.56 183.434
21         Equidae            Equus         hemionus   250000   27185.41 148.438
22         Equidae            Equus     conversidens   306000   33288.40 162.557
23         Equidae            Equus     occidentalis   574000   62521.55 216.044
24         Equidae            Equus           scotti   555000   60447.95 212.767
25    Elephantidae        Mammuthus        imperator 10000000 1095468.77 809.113
26    Elephantidae        Mammuthus      primigenius  5500000  601787.85 611.722
27    Elephantidae        Mammuthus          columbi  8000000  875983.99 728.775
28 Gomphotheriidae      Cuvieronius             spp.  5000000  546975.59 585.137
29      Mammutidae           Mammut       americanum  4523800  494782.58 558.489
30  Megalonychidae        Megalonyx      jeffersonii   600000   65359.32 220.437
31   Megatheriidae     Eremotherium         rusconii  3500000  382609.88 495.696
32   Megatheriidae   Nothrotheriops        shastense   300000   32634.39 161.114
33    Mylodontidae    Glossotherium          harlani  1587000  173212.04 343.925
34 Macraucheniidae     Macrauchenia      patachonica   988000  107732.43 276.722
35 Macraucheniidae     Windhausenia             spp.   700000   76276.06 236.450
36         Equidae            Equus     santaeelenae   350000   38085.19 172.702
37         Equidae            Equus         lasallei   350000   38085.19 172.702
38         Equidae            Equus          neogeus   378000   41138.34 178.806
39         Equidae            Equus           andium   220000   23917.05 140.169
40         Equidae            Equus        insulatus   351000   38194.23 172.924
41         Equidae        Hippidion         saldiasi   265000   28819.89 152.373
42         Equidae        Hippidion       principale   511000   55646.49 204.941
43         Equidae     Onohippidium             spp.   310700   33800.72 163.676
44 Gomphotheriidae      Cuvieronius             spp.  5000000  546975.59 585.137
45 Gomphotheriidae    Haplomastodon       chimborazi  6000000  656610.09 637.069
46 Gomphotheriidae    Notiomastodon             spp.  6193000  677773.96 646.556
47 Gomphotheriidae    Stegomastodon         superbus  7580000  829905.31 710.611
48  Megalonychidae        Nothropus             spp.   100000   10854.25  98.674
49  Megalonychidae    Nothrotherium             spp.   150000   16294.59 118.130
50  Megalonychidae          Ocnopus             spp.   300000   32634.39 161.114
51  Megalonychidae         Valgipes             spp.   200000   21738.62 134.314

# Dataset with predicted mean values by family.
predretfam= 0.57 + 6.95 * ((topredict$bmass/1000)^0.33) + 
    (0.33 + 0.43 * (topredict$bmass/1000^0.33))^0.5
predfammean<- data.frame(predfammean, predretfam)

# Predicted mean values for extinct megafauna.
R> predfammean
taxa   bmass predgutcap predretfam
1    Eremotherium 2600000     309754     431.90
2    Toxodontidae  830000      99007     255.54
3         Equidae  440001      52523     191.51
4   Megatheriidae  775000      92454     247.67
5 Gomphotheriidae  830000      99007     255.54
6    Elephantidae 2412500     287440     417.21
7    Mylodontidae  678571      80962     233.13


