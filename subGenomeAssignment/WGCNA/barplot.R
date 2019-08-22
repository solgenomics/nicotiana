library(ggplot2)
df_root = data.frame(
  gene = rep(c('unbiased','S-biased','T-biased'),each=3),
  module = rep(c('same module','similar module','divergent module'),3),
  percentage = c(41.6, 27.9, 30.5, 
                 23.2, 27.4, 49.4, 
                 23.2, 30.9, 45.9)
)

df_leaf = data.frame(
  gene = rep(c('unbiased','S-biased','T-biased'),each=3),
  module = rep(c('same module','similar module','divergent module'),3),
  percentage = c(44.8, 29.1, 26.2, 
                 18.1, 31.4, 50.5, 
                 24.8, 25.2, 50.0)
)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

p_root <- ggplot(df_root, aes(x = gene, y = percentage))+
  geom_col(aes(fill = module), width = 0.7)+
  coord_flip()+
  scale_fill_grey()

p_leaf <- ggplot(df_leaf, aes(x = gene, y = percentage))+
  geom_col(aes(fill = module), width = 0.7)+
  coord_flip()+
  scale_fill_grey()

ggsave("homeologous gene in co-expression network.leaf.png",plot=p_leaf,
       path='C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA/round4')

ggsave("homeologous gene in co-expression network.root.png",plot=p_root,
       path='C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA/round3')

