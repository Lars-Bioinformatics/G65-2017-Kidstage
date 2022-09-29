args = commandArgs(trailingOnly=TRUE)
# Expecting 4 args: seqz_bin file, normal name, tumor name, output dir
# print(args)             

library(sequenza)
seqzdata = sequenza.extract(args[1], chromosome.list = paste0("chr", c(1:22,"X","Y")))
CP = sequenza.fit(seqzdata)

id = paste0("Sequenza_",args[3],"_vs_",args[2])
out.dir = args[4]
sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = id, out.dir=out.dir)