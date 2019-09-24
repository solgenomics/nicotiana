library(ggplot2)
df_root = data.frame(
  gene = rep(c('unbiased','S-biased','T-biased'),each=3),
  module = rep(c('same module','similar module','divergent module'),3),
  percentage = c(41.6, 27.9, 30.5, 
                 23.2, 27.4, 49.4, 
                 23.2, 30.9, 45.9)
)

df_root$module = factor(df_root$module, 
                           levels=c("divergent module","similar module","same module"))


df_root_08 = data.frame(
  gene = rep(c('unbiased(4367)','S-biased(1115)','T-biased(967)'),each=3),
  module = rep(c('same module','similar module','divergent module'),3),
  percentage = c(41.6, 27.9, 30.5, 
                 26.2, 30, 43.8, 
                 29.5, 29.4, 41.1)
)
df_root_08$module = factor(df_root_08$module, 
                           levels=c("divergent module","similar module","same module"))

df_leaf = data.frame(
  gene = rep(c('unbiased','S-biased','T-biased'),each=3),
  module = rep(c('same module','similar module','divergent module'),3),
  percentage = c(44.8, 29.1, 26.2, 
                 18.1, 31.4, 50.5, 
                 24.8, 25.2, 50.0)
)
df_leaf$module = factor(df_leaf$module, 
                           levels=c("divergent module","similar module","same module"))



df_leaf_08 = data.frame(
  gene = rep(c('unbiased(4225)','S-biased(987)','T-biased(866)'),each=3),
  module = rep(c('same module','similar module','divergent module'),3),
  percentage = c(44.8, 29.1, 26.2, 
                 23.6, 28.6, 47.8, 
                 25.6, 29.3, 45.1)
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

ggsave("homeologous gene in co-expression network.leaf.png",plot=p_leaf,
       path='C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA/round4')

ggsave("homeologous gene in co-expression network.leaf.80%%.png",plot=p_leaf_08,
       path='C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA/round4')


ggsave("homeologous gene in co-expression network.root.png",plot=p_root,
       path='C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA/round3')

ggsave("homeologous gene in co-expression network.root.80%%.png",plot=p_root_08,
       path='C:/Users/10453/source/repos/SGN/nicotiana/subGenomeAssignment/WGCNA/round3')
