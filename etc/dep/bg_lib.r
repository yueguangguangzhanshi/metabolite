get.pca.plot <- function(quant.data,group,name){
  pc.cr<-prcomp(t(quant.data))
  p <- ggord::ggord(pc.cr, group, obslab=F, arrow=NULL, txt=NULL, ellipse =T, alpha=0.5, poly=F, coord_fix=F, facet = F, size=4, repel=T) +
    ggrepel::geom_text_repel(aes(label=lab),size=3,color="black",show.legend=F) +
    theme(plot.title = element_text(hjust = 0.5, size=12), text = element_text(size = 8)) +
    labs(title="PCA Plot",fill=NULL,colour=NULL,shape=NULL)
  ggsave(name,p)
}
