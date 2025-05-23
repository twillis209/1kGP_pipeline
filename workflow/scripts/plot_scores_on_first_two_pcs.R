library(data.table)
library(ggplot2)

theme_set(theme_bw())

scores <- fread(snakemake@input$scores)
ped <- fread(snakemake@input$ped)

scores[ped, on = .(`#IID` = SampleID), `:=` (pop = i.Population, superpop = i.Superpopulation)]

ggplot(scores, aes(x = PC1_AVG, y = PC2_AVG)) +
  geom_point(aes(color = superpop), size = 0.5) +
  labs(x = "PC1", y = "PC2", color = "Superpopulation") +
  theme(legend.position = "right") +
  ggtitle("Scores on first two principal components")

ggsave(snakemake@output[[1]], width = 6, height = 4)
