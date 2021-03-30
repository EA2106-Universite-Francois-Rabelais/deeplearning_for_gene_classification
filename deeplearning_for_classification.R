library(h2o)
library(data.table)

#load annotation
annot<-fread("cro_ap_asm_v2.gene_models.gff3", header = F)
annot.genes<-annot[V3=="gene"]
annot.genes[,GeneID:=sapply(V9, function(x)gsub("ID=CRO_", "CRO_T", strsplit(x, ";")[[1]][1]))]
annot.genes[,Functions:=sapply(V9, function(x)strsplit(x, "Note=")[[1]][2])]
annot.genes<-annot.genes[order(GeneID)]

#nonMIA are found among hypothetical proteins
hypothetical<-sort(table(annot.genes$Functions), decreasing=T)[1:5]

#MIA are identified by blast against referene datasets
besthits<-fread("MIA_besthists_cro_ap_asm_v2.transcripts.fasta", header=F)

#load expression matrix
exp.table<-fread("190812_expression.tpm.cro", header=F)
exp.table<-exp.table[order(V1)]

#label genes
#MIA genes
exp.table[, MIA:=ifelse(V1 %in% besthits$V2, "MIA", "others")]
exp.table[, MIAupstream:=ifelse(V1 %in% besthits[1:23]$V2, "MIA", "others")]
exp.table[, MIAupstreamvalid:=ifelse(V1 %in% besthits[24:nrow(besthits)]$V2, "MIA", "others")]
hypothetical.prots<-annot.genes[Functions=="hypothetical%20protein", GeneID]
#one nonMIA out of 2 will be kept
hypothetical.prots<-hypothetical.prots[seq(1, length(hypothetical.prots), 2)]
exp.table[V1 %in% hypothetical.prots, MIA:="nonMIA"]
exp.table[V1 %in% hypothetical.prots, MIAupstream:="nonMIA"]
exp.table[V1 %in% hypothetical.prots, MIAupstreamvalid:="nonMIA"]

#final counts
table(exp.table$MIA)
table(exp.table$MIAupstream)

#initialize cluster
h2o.init()

#export data to cluster
full<-as.h2o(exp.table[, c(2:86)])
full[,85]<-as.factor(full[, 85])

#prepare training dataset
training.with.upstream<-as.h2o(exp.table[MIAupstream %in% c("MIA", "nonMIA"), c(2:85, 87)])
training.with.upstream[,85]<-as.factor(training.with.upstream[,85])

#prepare validation dataset
validation.with.upstream<-as.h2o(exp.table[MIAupstreamvalid %in% c("MIA", "nonMIA"), c(2:85, 88)])
validation.with.upstream[,85]<-as.factor(validation.with.upstream[,85])
colnames(validation.with.upstream)[85]<-"MIAupstream"

#the training dataset have 23 known MIA
h2o.table(training.with.upstream[,85])

#prediction with upstream
mod1<-h2o.deeplearning(y=85, x=1:84, 
                 training_frame = training.with.upstream, 
                 validation_frame = validation.with.upstream, epochs=100, 
                 hidden=c(8), nfolds=10, seed=10, reproducible = T)
plot(mod1)
predictions<-h2o.predict(mod1, full[,-85])
h2o.table(predictions[,"predict"], full[,85])

##the best one: highest improvement on MIA
mod2<-h2o.deeplearning(y=85, x=1:84, 
                       training_frame = training.with.upstream, 
                       validation_frame = validation.with.upstream, epochs=100, 
                       hidden=c(16), nfolds=10, seed=10, reproducible = T)
plot(mod2)
predictions<-h2o.predict(mod2, full[,-85])
h2o.table(predictions[,"predict"], full[,85])
predicted.MIA<-unlist(exp.table[which(as.vector(as.data.frame(predictions)[,"predict"])=="MIA"),1])

#large increase in FP
mod3<-h2o.deeplearning(y=85, x=1:84, 
                       training_frame = training.with.upstream, 
                       validation_frame = validation.with.upstream, epochs=100, 
                       hidden=c(64), nfolds=10, seed=10, reproducible = T)
plot(mod3)
predictions<-h2o.predict(mod3, full[,-85])
h2o.table(predictions[,"predict"], full[,85])

#introducing l2
mod4<-h2o.deeplearning(y=85, x=1:84, 
                       training_frame = training.with.upstream, 
                       validation_frame = validation.with.upstream, epochs=100, 
                       hidden=c(200), nfolds=10, seed=10, reproducible = T, l2=1e-05)
plot(mod4)
predictions<-h2o.predict(mod4, full[,-85])
h2o.table(predictions[,"predict"], full[,85])

#introducing l2 46 TP, 865 FP
mod5<-h2o.deeplearning(y=85, x=1:84, 
                       training_frame = training.with.upstream, 
                       validation_frame = validation.with.upstream, epochs=100, 
                       hidden=c(200), nfolds=10, seed=10, reproducible = T, l2=1e-05, input_dropout_ratio = 0.2)
plot(mod5)
predictions<-h2o.predict(mod5, full[,-85])
h2o.table(predictions[,"predict"], full[,85])

#introducing l2 44 TP 483 FP
mod6<-h2o.deeplearning(y=85, x=1:84, 
                       training_frame = training.with.upstream, 
                       validation_frame = validation.with.upstream, epochs=100, 
                       hidden=c(400), nfolds=10, seed=10, reproducible = T, l2=1e-05, input_dropout_ratio = 0.2)
plot(mod6)
predictions<-h2o.predict(mod6, full[,-85])
h2o.table(predictions[,"predict"], full[,85])
