dt1 = c(.056, .211, .258, .272, .334, .266,.264)
dt1 = rbind(dt1, c(-.004, -.009, -.005, -.008, -.007, -.008, -.005))
dt1 = rbind(dt1, c(.038, .147, .151, .154, .155, .162, .162))
dt1 = rbind(dt1, c(.005, .090, .102, .104, .109, .126, .136))
dt1 = rbind(dt1, c(.048, .137, .182, .212, .211, .210, .202))
dt1 = rbind(dt1, c(.103, .163, .169, .173, .175, .180, .178))
dt1 = rbind(dt1, c(.018, .090, .102, .109, .118, .122, .122))

colnames(dt1) = c('0 min', '5 min','10 min','15 min','20 min','25 min','30 min') 
rownames(dt1) = seq(1,7)

dt2 = c(NA, .004, .003, .005)
dt2 = rbind(dt2, c(.025, .001, .000, .005))
dt2 = rbind(dt2, c(.015, .007, .003, .003))
dt2 = rbind(dt2, c(.013, .003, .008, .006))
colnames(dt2) = c('0 min', '10 min','20 min','30 min')
rownames(dt2) = seq(3,7)[-3]


dt3 = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)
dt3 = cbind(dt3, c(.044, .045, .048, .049, .051, .053, .055, .057, .059, .061))
colnames(dt3) = c("Time (s)", "Abs 340 nm")


nn = ncol(t(dt1))
matplot(t(dt1), type = 'l')
legend("bottomright", colnames(t(dt1)),col=seq_len(nn),cex=0.8,fill=seq_len(nn))
