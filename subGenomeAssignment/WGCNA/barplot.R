library(ggplot2)


df_root_08 = data.frame(
  gene = rep(c('unbiased(4367)','S-biased(967)','T-biased(1115)'),each=3),
  module = rep(c('same module','similar module','divergent module'),3),
  percentage = c(34.46, 37.07, 28.46, 
                 20.0, 33.4, 46.6, 
                 23.1, 33.5, 43.4)
)
df_root_08$module = factor(df_root_08$module, 
                           levels=c("divergent module","similar module","same module"))



df_leaf_08 = data.frame(
  gene = rep(c('unbiased(4225)','S-biased(866)','T-biased(987)'),each=3),
  module = rep(c('same module','similar module','divergent module'),3),
  percentage = c(37.4, 38.0, 24.7, 
                 18.7, 32.9, 48.4, 
                 18.7, 34.8, 46.5)
)

df_leaf_08$module = factor(df_leaf_08$module, 
                           levels=c("divergent module","similar module","same module"))


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

p_leaf_08 <- ggplot(df_leaf_08, aes(x = gene, y = percentage))+
  geom_col(aes(fill = module), width = 0.7)+
  coord_flip()+
  scale_fill_grey()

p_root_08 <- ggplot(df_root_08, aes(x = gene, y = percentage))+
  geom_col(aes(fill = module), width = 0.7)+
  coord_flip()+
  scale_fill_grey()


ggsave("homeologous gene in co-expression network.leaf.80%%.png",plot=p_leaf_08,
       path='C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA')


ggsave("homeologous gene in co-expression network.root.80%%.png",plot=p_root_08,
       path='C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA')
