# dataframe that holds usernames, passwords and other user data
user_base_raw <- tibble::tibble(
  user = c("liuhebin_guest","liuhebin","liuxiaohui","liuxiaohui_guest","lishunxiang","peter","yingweiwu","yating","huanglin","shenchengpin","qiankun","wanjingjing","yangzhuoni","yangzhuoni_guest","libin","qiziheng","Jane","tongbingrun","wangxiaoqing","wangxufei","yixiaotian","chenzhaozheng","bianxiaoyu","hedandan","lirongxin","mengshuhong","fuqin","chenxueying","liuconghui","xiangruiqi","zhouyutong","dukairui", "guozc2025"),
  password = c("liuhebin000","liuhebin000","liuxiaohui001","liuxiaohui001", "lishunxiang002","peter003","yingweiwu004","yating005","huanglin006","shenchengpin007","qiankun008","wanjingjing009","yangzhuoni010","yangzhuoni010","libin011","qiziheng012","Jane013","tongbingrun014","wangxiaoqing015","wangxufei016","yixiaotian017","chenzhaozheng018","bianxiaoyu019","hedandan020","lirongxin021","mengshuhong022","fuqin023","chenxueying024","liuconghui024","xiangruiqi025","zhouyutong026","dukairui027", "guozc202509"),
  permissions = c("standard","admin", "admin","standard", "standard","standard","standard","standard","standard","admin","admin","admin","admin","standard","standard","standard","admin","standard","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin"),
  name = c("guest","admin","刘晓慧","刘晓慧_Guest","李顺祥","Peter","应伟武","雅婷","黄琳","沈诚频","钱坤","万晶晶","杨卓妮","杨卓妮_Guest","李斌","祁子恒","Jane","仝冰润","王晓庆","王徐飞","易晓田","陈召政","边小禹","何丹丹","李荣鑫","蒙书红","付勤","陈雪莹","刘聪慧","向瑞奇","周昱彤","dukairui","郭振昌")
)

# write.csv(user_base_raw[-c(1:2),], paste0(Sys.getenv("HOME"), "/shiny-server/BioAnalysis/db/user_base.rds"))
write.csv(user_base_raw, paste0(Sys.getenv("HOME"), "/shiny-server/BioAnalysis/db/user_base.csv"))

user_base <- tibble::tibble(
  user = c("liuhebin_guest","liuhebin","liuxiaohui","liuxiaohui_guest","lishunxiang","peter","yingweiwu","yating","huanglin","shenchengpin","qiankun","wanjingjing","yangzhuoni","yangzhuoni_guest","libin","qiziheng","Jane","tongbingrun","wangxiaoqing","wangxufei","yixiaotian","chenzhaozheng","bianxiaoyu","hedandan","lirongxin","mengshuhong","fuqin","chenxueying","liuconghui","xiangruiqi","zhouyutong","dukairui", "guozc2025"),
  password = purrr::map_chr(c("liuhebin000","liuhebin000","liuxiaohui001","liuxiaohui001", "lishunxiang002","peter003","yingweiwu004","yating005","huanglin006","shenchengpin007","qiankun008","wanjingjing009","yangzhuoni010","yangzhuoni010","libin011","qiziheng012","Jane013","tongbingrun014","wangxiaoqing015","wangxufei016","yixiaotian017","chenzhaozheng018","bianxiaoyu019","hedandan020","lirongxin021","mengshuhong022","fuqin023","chenxueying024","liuconghui024","xiangruiqi025","zhouyutong026","dukairui027", "guozc202509"), sodium::password_store),
  permissions = c("standard","admin", "admin", "standard","standard","standard","standard","standard","standard","admin","admin","admin","admin","standard","standard","standard","admin","standard","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin","admin"),
  name = c("guest","admin","刘晓慧","刘晓慧_Guest","李顺祥","Peter","应伟武","雅婷","黄琳","沈诚频","钱坤","万晶晶","杨卓妮","杨卓妮_Guest","李斌","祁子恒","Jane","仝冰润","王晓庆","王徐飞","易晓田","陈召政","边小禹","何丹丹","李荣鑫","蒙书红","付勤","陈雪莹","刘聪慧","向瑞奇","周昱彤","dukairui","郭振昌")
)

saveRDS(user_base, paste0(Sys.getenv("HOME"), "/shiny-server/BioAnalysis/db/user_base.rds"))
# saveRDS(user_base, paste0(Sys.getenv("HOME"), "/mnt/os2/BioAnalysis/db/user_base.rds"))
