library(cronR)

f="/home/saitoasuka/shiny-server/BioAnalysis/etc/dep/bg_script_metabolite.r"
cmd <- cron_rscript(f)
# cmd<-str_replace(cmd,"/usr/local/lib/R/bin/Rscript","/usr/lib/R/bin/Rscript")
# cmd <- cron_rscript(f,rscript_log =NULL)
cron_add(command = cmd, frequency = 'minutely', id = 'metabonomics', description = 'metabonomics', tags = c('lab', 'xyz'),ask=F)

cronR::cron_ls()

cron_rm(id="metabonomics",ask=FALSE)
cronR::cron_clear(ask=FALSE)

