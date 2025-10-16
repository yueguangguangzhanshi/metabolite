args <- commandArgs(T)

data_path<-args[1];species<-args[2];kegg<-args[3];compare<-args[4];path_prefix<-args[5]

dir.create(paste0(data_path,"/KEGG_enrichment"),showWarnings = F)

library(reticulate)

kegg_script<-paste0('/home/saitoasuka/tools/kobas/scripts/run_kobas.py ', '-i ', data_path, '/', compare, '_diff.fasta', ' -t fasta:pro -s ', kegg, ' -N 12 -d K -o ', data_path, '/KEGG_enrichment/', compare, '_KEGG_enrichment.txt', ' -k ', path_prefix, '/db/kobas -b ', path_prefix, '/db/KEGG/', species, '_KEGG_annot.txt')

write_lines(kegg_script,paste0(data_path,"/kegg_script.sh"))

setwd(paste0(data_path,"/KEGG_enrichment"))

system(kegg_script,wait = T)
system(paste0("bash ",data_path,"/kegg_script.sh"),wait = T)

setwd("../")