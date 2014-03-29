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

require(ggplot2)
require(Hmisc)
require(scales)
require(grDevices)
require(dplyr)
#
#-------------------------------------------------------------------------
# BODY MASS SUMMARY DATASETS. EXTANT AND EXTINCT SPECIES.
#-------------------------------------------------------------------------
# Body mass data for extinct taxa.
# Body mass (g) of Late Quaternary mammals MOM database.
# Smith, F., S. Lyons, S. Ernest, K. Jones, D. Kaufman, T. Dayan, P. Marquet, J. Brown, and J. Haskell. 2003. Body mass of late quaternary mammals. Ecology 84:3403–3403.
# Body mass in g.
bmass_MOM <-read.table("./datasets/MOMv3.3.txt", header=TRUE, sep="\t",
    dec=".",na.strings="-999")

# Body mass range up to 15000 kg
extinct<- bmass_MOM %.%     # Selecting just the extinct taxa.
    dplyr::filter(status=="extinct") %.%
    dplyr::select(family, genus, species, meanbodymass)

# Selecting just the extinct frugivores. Body masses in g.
extinct.frug.megafauna<- bmass_MOM %.% 
    dplyr::filter(status=="extinct", 
           family== "Elephantidae" |
           family== "Equidae" |
           family== "Gomphotheriidae" |
           family== "Macraucheniidae" |
           family== "Mammutidae" | 
           family== "Megalonychidae" |
           family== "Megatheriidae" |
           family== "Mylodontidae") %.%
    dplyr::select(family, genus, species, meanbodymass)
colnames(extinct.frug.megafauna)<- c("family", "genus", 
                                     "species", "bmass")

# Selecting different genera with dplyr. *** Body masses in g ***.
eremoth<- bmass_MOM %.%
    dplyr::filter(genus=="Eremotherium") %.%
    dplyr::select(species, meanbodymass) %.%
    summarise(mean(meanbodymass))
toxodontids<- bmass_MOM %.%
    dplyr::filter(family=="Toxodontidae") %.%
    dplyr::select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass))
equids<- bmass_MOM %.%
    dplyr::filter(family=="Equidae", status=="extinct") %.%
    dplyr::select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass))
megatherids<- bmass_MOM %.%
    dplyr::filter(family=="Megatheriidae", status=="extinct") %.%
    dplyr::select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass))
gomphotherids<- bmass_MOM %.%
    dplyr::filter(family=="Gomphotheriidae", status=="extinct") %.%
    dplyr::select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass))
elephantids<- bmass_MOM %.%
    dplyr::filter(family=="Elephantidae", status=="extinct") %.%
    dplyr::select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass)) 
mylo<- bmass_MOM %.%
    dplyr::filter(family=="Mylodontidae", status=="extinct") %.%
    dplyr::select(genus, species, meanbodymass) %.%
    summarise(mean(meanbodymass))
tags<- c("Eremotherium", "Toxodontidae", "Equidae", "Megatheriidae",
         "Gomphotheriidae", "Elephantidae", "Mylodontidae")
topredict<- data.frame(taxa=tags, 
                 bmass= rbind(eremoth, toxodontids, equids, 
                              megatherids, gomphotherids, 
                              elephantids, mylo))
colnames(topredict)<- c("taxa","bmass") # Body mass in g!

#-------------------------------------------------------------------------
# DIGESTIVE ORGAN CAPACITY.
#-------------------------------------------------------------------------
# Body mass (kg) vs. digestive organ capacity (wet g).
# Dataset obtained from digitization of figure in van Soest (1982).
# Body mass in kg; digestive capacity in g wet content.
dig.organ <-read.table("./datasets/Body mass (kg) vs Gut capc (g wet).txt", 
                       header=TRUE,sep="\t",dec=".",na.strings=".")
# NOTE: dig.organ has Body mass (kg); Gut capc (g wet)
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
 #   stat_summary(fun.data= mean_cl_normal) + 
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
                    y= log10(prediction$predgutcap), 
                    color="white", size=2)) +
                    theme(legend.position = "none")

# Plots with averaged data for Extinct megafauna families
p0 + geom_point(data= predfammean, aes(x= log10(predfammean$bmass/1000), 
    y= log10(predfammean$predgutcap))) + 
      geom_point(color="pink", size=1.5) +
    theme(legend.position = "none")

#-------------------------------------------------------------------------
# RETENTION TIME IN THE GUT.
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
colnames(prediction)<- c("family","genus","species","bmass",
                         "predgutcap","predret") # Body mass in g!

# New estimations
# Log10-scaled axes. Confidence interval for both the regression slope
# and the prediction interval.

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
    coord_cartesian(xlim= c(-4,8), ylim= c(-0.5,3)) +
    theme(legend.position= c(0.20, 0.85))
p0 # The plot

p1<-p0 + geom_point(data= prediction, aes(x= log10(prediction$bmass), 
    y= predretfit1)) 
p1

# Plots with averaged data for Extinct megafauna families
p1 + geom_point(data= predfammean, aes(x= log10(predfammean$bmass), 
    y= log10(predfammean$predretfam)+0.45), size=4, color= "pink") +
     geom_vline(xintercept= log10(predfammean$bmass), 
                            linetype="dotted")
#-------------------------------------------------------------------------
# Datatset with body mass (g) and estimated gut capacity (wet g).
# Using van Soest regression equation.
# Regression equation for time and digestion-limited herbivores.
# gutret.h= 0.57 + 6.95*bmass.kg^0.33 + (0.33 + 0.43 * bmass^0.33)^0.5
# Predict specific values of gut retention from body mass data.
# Just the EXTINCT MEGAFAUNA FRUGIVOROUS species.
predret= 0.57 + 6.95 * ((extinct.frug.megafauna$bmass/1000)^0.33) + 
    (0.33 + 0.43 * (extinct.frug.megafauna$bmass/1000^0.33))^0.5
prediction<- data.frame(prediction, predret)
colnames(prediction)<-  c("family","genus","species","bmass",
    "predgutcap","predret.vanSoest","predret.vanSoest.timelim") 
head(prediction) # Now the dataset holds all variables.
# Body mass in g!
# predgutcap, in g
# predret.vanSoest: raw equation bmass vs retention time.
# predret.vanSoest.timelim: regr. equation van Soest, expanded with parameters
#                           for time and digestion-limited herbivores.

R> prediction
family           genus       species    bmass predgutcap
1          Equidae           Equus      capensis   350000      38085
2          Equidae       Hipparion       libycum   150000      16295
3     Elephantidae         Elephas      iolensis  6500000     711441
4          Equidae           Equus       alaskae   372000      40484
5          Equidae           Equus    laurentius   648000      70599
6          Equidae           Equus      caballus   250000      27185
7          Equidae           Equus     fraternus   259000      28166
8          Equidae           Equus  niobrarensis   334000      36341
9          Equidae           Equus   complicatus   400000      43538
10         Equidae           Equus     giganteus   400000      43538
11         Equidae           Equus      hemionus   250000      27185
12         Equidae           Equus  conversidens   306000      33288
13         Equidae           Equus  occidentalis   574000      62522
14         Equidae           Equus        scotti   555000      60448
15    Elephantidae       Mammuthus     imperator 10000000    1095469
16    Elephantidae       Mammuthus   primigenius  5500000     601788
17    Elephantidae       Mammuthus       columbi  8000000     875984
18 Gomphotheriidae     Cuvieronius          spp.  5000000     546976
19      Mammutidae          Mammut    americanum  4523800     494783
20  Megalonychidae       Megalonyx   jeffersonii   600000      65359
21   Megatheriidae    Eremotherium      rusconii  3500000     382610
22   Megatheriidae  Nothrotheriops     shastense   300000      32634
23    Mylodontidae   Glossotherium       harlani  1587000     173212
24 Macraucheniidae    Macrauchenia   patachonica   988000     107732
25 Macraucheniidae    Windhausenia          spp.   700000      76276
26         Equidae           Equus  santaeelenae   350000      38085
27         Equidae           Equus      lasallei   350000      38085
28         Equidae           Equus       neogeus   378000      41138
29         Equidae           Equus        andium   220000      23917
30         Equidae           Equus     insulatus   351000      38194
31         Equidae       Hippidion      saldiasi   265000      28820
32         Equidae       Hippidion    principale   511000      55646
33         Equidae    Onohippidium          spp.   310700      33801
34 Gomphotheriidae     Cuvieronius          spp.  5000000     546976
35 Gomphotheriidae   Haplomastodon    chimborazi  6000000     656610
36 Gomphotheriidae   Notiomastodon          spp.  6193000     677774
37 Gomphotheriidae   Stegomastodon      superbus  7580000     829905
38  Megalonychidae       Nothropus          spp.   100000      10854
39  Megalonychidae   Nothrotherium          spp.   150000      16295
40  Megalonychidae         Ocnopus          spp.   300000      32634
41  Megalonychidae        Valgipes          spp.   200000      21739
42   Megatheriidae    Eremotherium   laurillardi   800000      87196
43   Megatheriidae    Eremotherium      rusconii  3500000     382610
44   Megatheriidae     Megatherium    americanum  6265000     685670
45   Megatheriidae Paramegatherium          spp.  3500000     382610
46    Mylodontidae   Glossotherium      myloides  1200000     130900
47    Mylodontidae   Glossotherium      robustum  1713000     186993
48    Mylodontidae        Lestodon       armatus  3397000     371328
49    Mylodontidae         Mylodon        listai  1000000     109044
50    Mylodontidae      Scelidodon          spp.  1000000     109044
51    Mylodontidae  Scelidotherium leptocephalum  1119000     122047
predret.vanSoest predret.vanSoest.timelim
1            211.58                  172.702
2            173.30                  118.130
3            421.14                  661.334
4            214.64                  177.519
5            244.63                  228.287
6            195.46                  148.438
7            197.09                  150.814
8            209.26                  169.096
9            218.35                  183.434
10           218.35                  183.434
11           195.46                  148.438
12           204.99                  162.557
13           237.74                  216.044
14           235.86                  212.767
15           466.12                  809.113
16           404.88                  611.722
17           442.25                  728.775
18           395.89                  585.137
19           386.67                  558.489
20           240.23                  220.437
21           363.98                  495.696
22           204.04                  161.114
23           302.10                  343.925
24           270.19                  276.722
25           249.12                  236.450
26           211.58                  172.702
27           211.58                  172.702
28           215.46                  178.806
29           189.66                  140.169
30           211.73                  172.924
31           198.16                  152.373
32           231.32                  204.941
33           205.73                  163.676
34           395.89                  585.137
35           413.27                  637.069
36           416.36                  646.556
37           436.67                  710.611
38           157.51                   98.674
39           173.30                  118.130
40           204.04                  161.114
41           185.45                  134.314
42           257.08                  251.287
43           363.98                  495.696
44           417.50                  650.056
45           363.98                  495.696
46           282.85                  302.487
47           307.59                  356.235
48           361.43                  488.871
49           270.96                  278.254
50           270.96                  278.254
51           278.23                  292.952
R> 
#
# AVERAGED VALUES FOR FAMILIES. -----------------------------------------
# Dataset with predicted mean values by family.
str(prediction)
predfammean<- prediction %.% group_by(family) %.% summarise(mean(bmass), 
    mean(predgutcap), mean(predret.vanSoest), 
    mean(predret.vanSoest.timelim))
colnames(predfammean)<- c("family","bmass", "gutcapacity",
                      "retention","ret.timelim") 
head(predfammean) # Now the dataset holds all variables.

# Plots with averaged data for Extinct megafauna families
p1 + geom_point(data= predfammean, aes(x= log10(predfammean$bmass), 
    y= log10(predfammean$retention)), size=4, color= "pink", alpha=0.5) +
    geom_vline(xintercept= log10(predfammean$bmass), 
        linetype="dotted")

# Based on regressions. -------------------------------------------------
# Gut capacity (g).
predfammean<- data.frame(topredict, 
        predgutcap= 10^(2.0316 + 
                1.0020*log10(topredict$bmass/1000)))

# Gut retention time. Raw model, van Soest.
predfamret= 10^(1.0193 + 
    0.2356 * log10(topredict$bmass/1000))

# Gut retention time. Model for time-limited foragers, van Soest.
predfamrettimelim= 0.57 + 6.95 * ((topredict$bmass/1000)^0.33) + 
    (0.33 + 0.43 * (topredict$bmass/1000^0.33))^0.5

predfammean<- data.frame(predfammean, predfamret, predfamrettimelim)
colnames(predfammean)<-  c("family","genus","species","bmass",
    "predgutcap","predret.vanSoest","predret.vanSoest.timelim") 
head(predfammean) # Now the dataset holds all variables.

# -----------------------------------------------------------------------
# Body mass vs home range size (PanTHERIA dataset).
#-------------------------------------------------------------------------
# Body mass (g) vs. home-range area (km2).
# I have already trimmed the species with missing data.
pantheria <-read.table("./datasets/mypantheria.txt", 
                    header=TRUE, sep="\t", dec=".", na.strings="-999")

# Regression equation.
# Linear model with log-transformed values.
lm_fit5<- lm(log10(homerangekm2) ~ 
             log10(bmass/1000), 
             data= pantheria)
summary(lm_fit5)

hr_with_pred<- data.frame(pantheria,   # Predicted values for plotting.
                          predict(lm_fit5, interval= 'prediction'))
head(hr_with_pred)

prediction<- data.frame(prediction, 
                        predhomerange= 10^(-1.2070 + 
                                1.0435*log10(prediction$bmass/1000)))
head(prediction)

plot(pantheria$bmass/1000, pantheria$homerangekm2, log = "xy")

# Plots
# Log10-scaled axes. Confidence interval for the regression slope.
p0<- ggplot(hr_with_pred, aes(x= log10(bmass/1000), 
    y= log10(homerangekm2))) + 
    geom_point(color="black", size=2) +
    geom_smooth(method= 'lm', aes(fill= 'confidence'), alpha= 0.5) +
    geom_ribbon(aes(y= fit, ymin= lwr, ymax= upr, fill= 'prediction'),
        alpha= 0.2) +
    scale_fill_manual('Interval', values= c('green', 'blue')) +
    xlab("log10(Body mass) (kg)") +
    ylab("log10(Home range area) (km2)") +
    ggtitle("Home range area estimation for Megafauna frugivores") +
    coord_cartesian(xlim= c(-2.5,5), ylim= c(-6,6)) +
    theme(legend.position= c(0.20, 0.85))
p0 # The plot

p1<-p0 + geom_point(data= prediction, aes(x= log10(prediction$bmass/1000), 
    y= log10(prediction$predhomerange), 
    color="white", size=4)) +
    theme(legend.position = "none")
p1

prediction

# -----------------------------------------------------------------------
# Body mass vs population density (PanTHERIA dataset).
#-------------------------------------------------------------------------
# Body mass (g) vs. population density (no. ind/km2).
# I have already trimmed the species with missing data.
pantheria.popdens <-read.table("./datasets/mypantheria_popdens.txt", 
    header=TRUE, sep="\t", dec=".", na.strings="-999")

# Regression equation.
# Linear model with log-transformed values.
lm_fit6<- lm(log10(popdens.n.km2) ~ 
             log10(bmass/1000), 
             data= pantheria.popdens)
summary(lm_fit6)

popdens_with_pred<- data.frame(pantheria.popdens,   
                                    # Predicted values for plotting.
                                predict(lm_fit6, interval= 'prediction'))
head(popdens_with_pred)

prediction<- data.frame(prediction, 
                        predpopdens= 10^(1.6767 + 
                            -0.7976*log10(prediction$bmass/1000)))
head(prediction)

plot(pantheria.popdens$bmass/1000, pantheria.popdens$popdens.n.km2, log = "xy")

# Plots
# Log10-scaled axes. Confidence interval for the regression slope.
p0<- ggplot(popdens_with_pred, aes(x= log10(bmass/1000), 
    y= log10(popdens.n.km2))) + 
    geom_point(color="black", size=2) +
    geom_smooth(method= 'lm', aes(fill= 'confidence'), alpha= 0.5) +
    geom_ribbon(aes(y= fit, ymin= lwr, ymax= upr, fill= 'prediction'),
        alpha= 0.2) +
    scale_fill_manual('Interval', values= c('green', 'blue')) +
    xlab("log10(Body mass) (kg)") +
    ylab("log10(Population density) (n. ind km^-2)") +
    ggtitle("Population density estimation for Megafauna frugivores") +
    coord_cartesian(xlim= c(-2.5,5), ylim= c(-6,6)) +
    theme(legend.position= c(0.15, 0.15))
p0 # The plot

p1<-p0 + geom_point(data= prediction, aes(x= log10(prediction$bmass/1000), 
    y= log10(prediction$predpopdens), 
    color="white", size=4)) +
    theme(legend.position = "none")
p1

prediction

