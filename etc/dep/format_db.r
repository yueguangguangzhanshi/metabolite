#visirot db

visitor <- data.frame(
	"ip" = "0.0.0.0",
	"fingerprint" = "abcdefg",
	"time" = Sys.time()
)

visitordb <- dbConnect(RSQLite::SQLite(), "/home/saitoasuka/shiny-server/fluid/data/visitors.db")

dbWriteTable(visitordb, "Visitors", visitor, overwrite=TRUE)

dbDisconnect(visitordb)

#temp<- dbGetQuery(visitordb, paste("SELECT * FROM Visitors", sep=""))


#create user id database

usr <- data.frame(
	"uid" = 1,
	"username" = "admin",
	"password" = digest("12345", algo="sha1"),
	"email" = 'test@omicsolution.com',
	"createtime" = Sys.time(),
	"privilege" = 0,
	"sessionid" = "",
	"logintime"=as.character(Sys.time())
)

usrdb <- dbConnect(RSQLite::SQLite(), "/home/saitoasuka/shiny-server/usrdb.db", flags = SQLITE_RWC)

#usrdb_table <- dbGetQuery(usrdb, "SELECT * FROM User")

dbWriteTable(usrdb, "User", usr)

dbDisconnect(usrdb)

species <- data.frame(
	"name" = c("Homo_sapiens", "Mus_musculus", "Rattus_norvegicus","Arabidopsis_thaliana", paste0('dummy_species', 1:30)),
	"ch_names" = c("智人", "小鼠", "大鼠","拟南芥",paste0('测试物种',1:30)),
	"sci_names" = c("Homo sapiens", "Mus musculus", "Rattus norvegicus","Arabidopsis thaliana",paste0('dummy_species',1:30)),
	"othernames" = c("human", "mouse", "rat","",paste0('dummy_species',1:30)),
	"species_L0" = c("Eukaryotes", "Eukaryotes", "Eukaryotes","Eukaryotes", rep('dummy',30)),
	"species_L1" = c("Animals", "Animals", "Animals","Plants", rep('dummy',30)),
	"species_L2" = c("Vertebrates", "Vertebrates", "Vertebrates","Eudicots", rep('dummy',30)),
	"species_L3" = c("Mammals", "Mammals", "Mammals","Mustard family", rep('dummy',30)),
	"dbtype" = c("Uniprot", "Uniprot", "Uniprot", "Uniprot", rep('dummy',30)),
	"kegg" = c("hsa", "mmu", "rno", "ath", rep('dummy',30)),
	"taxid" = c("9606", "10090", "10116", "3702", rep('dummy',30)),
	"STRING.blast" = c("9606", "10090", "10116", "3702", rep('dummy',30)),
	"kegg.blast" = c("hsa", "mmu", "rno", "ath", rep('dummy',30)),
	"proteomeid" = c("UP000005640", "UP000000589", "UP000002494", "", rep('dummy',30)),
	"COG" = c("KOG", "KOG", "KOG", "KOG", rep('dummy',30)),
	"privilege" = c(99, 99, 99, 99, rep(99,30)),
	"updatetime" = c(rep(Sys.time(),4), rep(Sys.time(),30))
)

speciesdb <- dbConnect(RSQLite::SQLite(), "/home/saitoasuka/shiny-server/project/db/speciesdb.db")

#species <- dbGetQuery(speciesdb, "SELECT * FROM Species")

dbWriteTable(speciesdb, "Species", species, overwrite = T)

dbDisconnect(speciesdb)

#progress database
#0:排队 1:完成 2:中止 3:错误
projects_id <- list.files(paste("/home/saitoasuka/shiny-server/project/data/usrdata/",sep=""))

projects.list <- data.frame(
	'id' = project_id,
	'name' = rep('Project', length(projects_id)),
	'user' = rep('admin', length(projects_id)),
	'status' = rep(0, length(projects_id)),
	'time' = rep(as.character(Sys.time()), length(projects_id)),
	stringsAsFactors = FALSE
)

#projects.list <- data.frame()

projects_list <- dbConnect(RSQLite::SQLite(), "/home/saitoasuka/shiny-server/project/db/projects_list.db", flags = SQLITE_RW)

#project_status_table <- dbGetQuery(projects_list, "SELECT * FROM Project")

#project_status_table <- project_status_table[0, ]

#dbWriteTable(projects_list, "Project", project_status_table, overwrite = T)


dbWriteTable(projects_list, "Project", projects.list)

dbDisconnect(projects_list)

# msg-database subject, content, sender, typ, rcvr, rdtag, rcpt_group,  valid_from, valid_thru
#type: 1. 公告 2. 提醒 3.消息
msg.list <- data.frame(
	'msg_index' = 1:10, 
	'subject' = c('测试公告title1', '测试公告title2', '测试通知title3', '测试通知title4', '测试通知title4', '测试通知title4', '测试通知title4', '测试通知title4', '测试通知title4', '测试通知title4'),
	'content' = c('测试公告content1', '测试公告content2', '测试通知content3', '测试通知content4', '测试通知content4', '测试通知content4', '测试通知content4', '测试通知content4', '测试公告content4', '测试公告content4'),
	'sender' = c('system', 'system', 'system', 'system', 'system', 'system', 'system', 'system', 'system', 'system'),
	'typ' = c('1', '1', '2', '2', '2', '1', '1', '1', '3', '3'),
	'rcvr' = c('', '', 'admin', 'admin', 'admin', '', '', '', 'admin', 'admin'),
	'rdtag' = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
	'rcpt_group' = c('100', '100', '100', '100', '100', '100', '100', '100', '100', '100'),
	'valid_from' = c(1594905253, 1594905253, 1594905253, 1594905253, 1594905253, 1594905253, 1594905253, 1594905253, 1594905253, 1594905253),
	'valid_thru' = c(1752289509, 1752289509, 1752289509, 1752289509, 1752289509, 1752289509, 1752289509, 1752289509, 1752289509, 1752289509),
	stringsAsFactors = FALSE
)

msg_list <- dbConnect(RSQLite::SQLite(), "/home/saitoasuka/shiny-server/project/db/msg.db", flags = SQLITE_RW)

dbWriteTable(msg_list, "Message", msg.list, overwrite = T)

dbDisconnect(msg_list)


