library(meteR)
devtools::install_github('ajrominger/socorro')
library(socorro)

setwd('~/Dropbox/Research/psi_mete')

## load data

rmbl <- read.csv('RMBL_PSI_test.csv', as.is = TRUE)

bci <- read.csv('BCIS.csv', as.is = TRUE)
bci <- bci[bci$year == 1995, ]
bcitest <- read.csv('BCI_PSI_test.csv', as.is = TRUE)

arth <- read.csv('gruner_kohala.csv', as.is = TRUE)

## make METE objects

rmblIPD <- ipd(meteESF(S0 = 31, N0 = 877, E0 = sum(rmbl$observed)))
## need to custom add data to IPD for rmbl
rmblIPD$data <- sort(rmbl$observed, decreasing = TRUE)

bciIPD <- ipd(meteESF(spp = bci$spp, abund = bci$count, power = bci$dbh^2))
## metabolic scaling from Mori et al. PNAS paper
bciIPDmori <- ipd(meteESF(spp = bci$spp, abund = bci$count, power = (bci$dbh^(8/3))^0.805))

pdf('fig_bciPsiMori.pdf', width = 5, height = 5)
plot(bciIPD, log = 'xy', axes = FALSE)
logAxis(1)
logAxis(2)
box()
dev.off()

arthIPD <- ipd(meteESF(spp = arth$SpeciesCode, abund = arth$Abundance, 
                       power = arth$IND_BIOM^0.75))


## plotting theory for BCI
e <- exp(seq(log(1), log(500000), length = 1000))

jpeg('ms/fig_PsiThr.jpg', width = 3, height = 3, units = 'in', res = 380)

par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(e, bciIPD$d(e), log = 'xy', type = 'l', col = 'red',
     xaxt = 'n', yaxt = 'n', ylim = 10^c(-12, 0),
     xlab = 'Metabolic rate', ylab = 'Probability density',
     panel.first = {
         rect(xleft = 10, xright = 10000, ybottom = 10^par('usr')[3], ytop = 10^par('usr')[4],
              col = 'gray60', border = NA)
         rect(xleft = 10000, xright = 10^par('usr')[2], 
              ybottom = 10^par('usr')[3], ytop = 10^par('usr')[4],
              col = 'gray80', border = NA)
     })
logAxis(1, expLab = TRUE)
axis(2, at = 10^seq(-12, 0, by = 3), 
     labels = sapply(seq(-12, 0, by = 3), 
                     function(p) eval(substitute(expression(10^p), list(p = p)))))

dev.off()


## plotting data and theory

jpeg('ms/fig_PsiData.jpg', width = 8, height = 3, units = 'in', res = 380)

par(mfrow = c(1, 3), oma = c(3, 2, 0, 0) + 0.5, mar = c(0, 2, 1, 0) + 0.2, 
    cex = 1, mgp = c(2, 0.75, 0))

plot(exp(bcitest$ln.rank.), exp(bcitest$ln.dbh2.), log = 'xy',
     ylim = c(1, 300000), xaxt = 'n', yaxt = 'n',
     xlab = '', ylab = '')
points(exp(bcitest$ln.rank.), exp(bcitest$ln.PRED_METE.), type = 'l', col = 'red')
logAxis(1, expLab = TRUE)
logAxis(2, expLab = TRUE)
legend('topright', legend = c('Data', 'METE'), pch = c(1, NA), lwd = c(NA, 1), 
       col = c('black', 'red'), bty = 'n', cex = 0.9)
mtext('A', side = 3, at = 10^(par('usr')[1] + 0.05 * diff(par('usr')[1:2])), line = 0.2)

plot(rmblIPD, ptype = 'rad', log = 'xy', add.legend = FALSE, 
     ylim = c(1, 100000), xaxt = 'n', yaxt = 'n',
     xlab = '', ylab = '')
logAxis(1, expLab = TRUE)
logAxis(2, expLab = TRUE)
mtext('B', side = 3, at = 10^(par('usr')[1] + 0.05 * diff(par('usr')[1:2])), line = 0.2)

plot(arthIPD, ptype = 'rad', log = 'xy', add.legend = FALSE,
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
logAxis(1, expLab = TRUE)
logAxis(2, expLab = TRUE)
mtext('C', side = 3, at = 10^(par('usr')[1] + 0.05 * diff(par('usr')[1:2])), line = 0.2)

mtext('Metabolic rate', side = 2, line = 0.5, outer = TRUE)
mtext('Rank', side = 1, line = 2, outer = TRUE)

dev.off()
