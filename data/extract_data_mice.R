load(file = "impc_quantitative_phen_matrix.RData")
cen.r = "10"  #Centre to extract
proc.r = "Hematology" #Procedure to extract

idv = meta[meta$cen == cen.r, "id"]
phenv = pmeta[pmeta$procnam == proc.r, "ph"]
dat = d[idv, phenv]
met = meta[match(idv, meta$id), ]
dat = cbind(met, dat)
pmet = pmeta[match(phenv, pmeta$ph), ]

# save(dat, met, pmet, file = paste(dat.dir, "/impc_quantitative_phen_submatrix.RData", sep = ""))
dim(dat)

#Note there are some missing values!!
#image(is.na(dat))
#image((dat))

# ignore any variables with > 50% missingness
pres <- apply(dat, 2, function(x) mean(is.na(x)) < 0.5)
final <- dat[complete.cases(dat[,pres]),pres]

#Note there are no more missing values!!
#image(is.na(final))
#image((final))

#geno <- as.integer(factor(final$geno))
#pairs(final[,7+1:10], col=as.numeric(final$sex)+1)
#pairs(final[,7+1:11], col=geno)
#pairs(final[,18+1:10], col=geno)

#library(ggplot2)
#library(reshape2)
#qplot(x=Var1, y=Var2, data=melt(cor(final[,8:29], use="p")), fill=value, geom="tile") +
#  scale_fill_gradient2(limits=c(-1, 1))

write.csv(final, file = "mouse.csv",row.names=FALSE)
