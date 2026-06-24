#Author: Juan Zegarra
#Date: 24/01/2026

library(dplyr)
library(ggplot2)

#kmer count file is a tsv file with 5 col
#sequence\tcoun\tsat\tlocation\ttotal
#normalization is made by dividing the count of each sat kmer by the total number of kmers
kmer_count_dmer$normalized <- kmer_count_dmer$count / kmer_count_dmer$total 
correct <- kmer_count_dmer[,c(3,4,6)]
correct <- correct %>% pivot_wider(values_from = normalized, names_from = location)
correct$log_ratio <- log2(correct$SER/correct$ANG)

#plot log-fold change between kmers
ggplot(total_log, aes(x = logfold, y = sat, fill = sat)) + 
  geom_bar(stat = "identity", width = 0.9)+
  theme(panel.background = element_rect(fill = "#F5F5F5", colour = "black", linewidth = 1),
        panel.spacing = unit(0.8, "cm"),
        panel.grid = element_blank(),
        axis.text.y =  element_blank(),
        axis.text.x = element_text(size = 13, face = "bold", colour = "black"), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 15, colour = "black", face = "bold"), 
        axis.title.x = element_text(size = 16, face = "bold", colour = "black")) + 
  scale_fill_viridis_d()+ 
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  labs(title = element_blank(),
       x = "Logfold change(Log2)", 
       y = element_blank()) +
  scale_x_continuous(breaks = seq(-10, 16, by = 2), limits = c(-10, 15))

#plot heatmap of kmer count 

kmer_count_dmer$normalized <- kmer_count_dmer$count / kmer_count_dmer$total
func_kmer <- kmer_count_dmer[,c(3,4,6)]
func_kmer$normalized_count <- log2(func_kmer$normalized_count)

ggplot(func_kmer, aes(x=sat, y = location, fill = normalized)) + geom_tile(colour = "black",
                                                                           linewidth = 0.6,
                                                                           height = 0.5) + scale_fill_viridis_c(option = "plasma") + 
  labs(fill = "Normalized kmer count"
  ) + 
  theme(
    axis.title = element_blank(), 
    axis.text = element_text(face = "bold", 
                             size = 12, color = "black"), 
    panel.background = element_blank(), 
    legend.title = element_text(face = "bold", size = 12), 
    legend.text = element_text(face = "bold", size = 11)) + scale_x_discrete(position = "top") + 
  coord_fixed()


#normalize mapping depth based on single copy genes
#utilizes both a tsv file of coverage of satellites and a tsv file of average single copy gene coverage in each genome
sat_cov<- coverage.stats[,c(1,2,8)]
normalized_d<- sat_cov %>% left_join(cov_single_copy, by = "location") %>% mutate(depth_norm=meandepth/mean_single_cov)
final_cov<- normalized_d[,c(1,2,5)]
final_cov <- final_cov %>% pivot_wider(names_from = location, values_from = depth_norm)
final_cov$logfold <- log2(final_cov$SER/final_cov$ANG)


#plot heatmap of normalized coverage
ggplot(normalized_d, aes(x=name, y = location, fill = depth_norm)) + 
  geom_tile(colour = "black", linewidth = 0.6, height = 0.5) + 
  scale_fill_viridis_c(option = "plasma") + 
  labs(fill = "Log2 of \nnormalized depth") +
  theme(axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 12, color = "black"),
        panel.background = element_blank(),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "bold", size = 11)) +
  scale_x_discrete(position = "top") + coord_fixed()


