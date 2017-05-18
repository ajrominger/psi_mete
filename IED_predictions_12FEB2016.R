##########
##########
# Clean up

par(mfrow=c(1,1))
rm(list=ls(all=TRUE))
##########
##########

##LOADING DATA##


BCI<- read.csv("BCI_PSI_test_csv.csv", head=TRUE, sep=",")
RMBL<- read.csv("RMBL_PSI_test.csv", head=TRUE, sep=",")

names(BCI)
names(RMBL)


#######################################################
##GRAPHING##

require(ggplot2)
library(ggplot2)
require(Rmisc)
library(Rmisc)
####################
# graphing setup

#par(mfrow=c(2,2))
#par(mar= c(4,4,2,2))
####################

x <- RMBL$log_rank
y1 <- RMBL$log_obs
y2 <- RMBL$log_predicted

df.rmbl <- data.frame(x, y1, y2) 
m1 <- ggplot(df.rmbl,aes(x, y = value)) +
      geom_point(aes(y = y1)) +
      geom_line(aes(y = y2))

m1 <- m1 + xlab("Ln(rank)")+ 
  ylab("Ln(metabolic rate)")

m1 <- m1 + theme_bw() 
m1 <- m1 + theme(axis.text = element_text(size = 12),
                 axis.title = element_text(size = 14))
m1 <- m1 + labs(title = "B                                                       ")
m1 <- m1 + theme(legend.position = c(0.2, 0.2))


##BCI DATA

x <- BCI$ln.rank.
y1 <- BCI$ln.dbh2.
y2 <- BCI$ln.PRED_METE.

df.bci <- data.frame(x, y1, y2) 

m2 <- ggplot(df.bci,aes(x, y = value)) +
  geom_point(aes(y = y1)) +
  geom_line(aes(y = y2))

m2 <- m2 + xlab("Ln(rank)")+ 
  ylab("Ln(metabolic rate)")

m2 <- m2 + theme_bw() 
m2 <- m2 + theme(axis.text = element_text(size = 12),
                 axis.title = element_text(size = 14))
m2 <- m2 + labs(title = "A                                                       ")
m2 <- m2 + theme(legend.position = c(0.2, 0.2), 
                 legend = c("data", "predicted"))


#Plot multiple graphs on one page
multiplot(m2, m1, cols=2)


##########
##########
# Clean up

par(mfrow=c(1,1))
rm(list=ls(all=TRUE))
##########
##########