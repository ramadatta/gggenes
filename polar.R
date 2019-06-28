library(tidyverse)
load_all()
example_genes <- example_genes %>%
  filter(molecule %in% c("Genome1", "Genome2"))

ggplot(example_genes, aes(
  xmin = start,
  xmax = end,
  y = molecule,
  fill = gene,
  forward = direction
)) +
  geom_gene_arrow() +
  coord_polar()
