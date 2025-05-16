args<-as.numeric(commandArgs(TRUE))
block_index<-args[1]

library(Logica)
trait = "LDL"
ancestry_1 = "EUR"
ancestry_2 = "EAS"
library(data.table)
z1<-fread(file.path("/net/fantasia/home/borang/MALGC/pipeline_example/LDL_Block",paste0(trait , "_",ancestry_1,"_block_",block_index ,".txt")))
z2<-fread(file.path("/net/fantasia/home/borang/MALGC/pipeline_example/LDL_Block",paste0(trait , "_",ancestry_2,"_block_",block_index ,".txt")))
load(file.path("/net/fantasia/home/borang/MALGC/pipeline_example/LD_ref",paste0("LD_",block_index ,".RData")))
load(file.path("/net/fantasia/home/borang/MALGC/pipeline_example/LDL_Block",paste0(trait,"_intercept.RData")))

preprocess_data_file<-preprocess_data(z1,z2,ancestry_1_cov,ancestry_2_cov)
logica_res<-run_Logica(preprocess_data_file$pop1$z,preprocess_data_file$pop2$z, preprocess_data_file$pop1$R, preprocess_data_file$pop2$R, median(preprocess_data_file$pop1$z$N), 
                        median(preprocess_data_file$pop2$z$N), 
                        z1_intercept = ancestry_1_lder_intercept, z2_intercept = ancestry_2_lder_intercept, fix_intercept = TRUE)
write.table(logica_res,file.path("/net/fantasia/home/borang/MALGC/pipeline_example/LDL_Block",paste0("Logica_res_",trait , "_block_",block_index ,".txt")),col.names=T,row.names=F,quote=F)

