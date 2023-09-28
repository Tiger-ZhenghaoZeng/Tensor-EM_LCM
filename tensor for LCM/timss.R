library(CDM)
data("data.timss11.G4.AUT.part")
aut_part <- data.timss11.G4.AUT.part$data
aut_dat <- aut_part[,18:64] # The 47 items
mice_impu <- mice(aut_dat, m=5, seed=521)
imp_data1 <- complete(mice_impu,1)
imp_data2 <- complete(mice_impu,2)
imp_data3 <- complete(mice_impu,3)
imp_data4 <- complete(mice_impu,4)
imp_data5 <- complete(mice_impu,5)
write.table(imp_data1, 'timss1', col.names=F, row.names=F)
write.table(imp_data2, 'timss2', col.names=F, row.names=F)
write.table(imp_data3, 'timss3', col.names=F, row.names=F)
write.table(imp_data4, 'timss4', col.names=F, row.names=F)
write.table(imp_data5, 'timss5', col.names=F, row.names=F)


