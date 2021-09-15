
reaction <- c(rep('HMR_3215', 3), rep('HMR_3763',3),
              rep('HMR_3215', 3), rep('HMR_3763',3))
cell_type <-  c(rep('breast', 1), rep('brain',1), rep('lung', 1),
                rep('breast', 1), rep('brain', 1), rep('lung', 1),
                rep('breast', 1), rep('brain',1), rep('lung', 1),
                rep('breast', 1), rep('brain', 1), rep('lung', 1))
min_max <- c(rep('min_value', 6), rep('max_value', 6))
reaction_cell <- paste(reaction, cell_type)
value <- c(48.5491, 32.1392,40.5472,
           48.5491, 32.1392 ,40.5472,
           169.0969, 117.555,139.2393,
           169.0969, 117.555,139.2393)
valin_leucine <- data.frame(reaction, cell_type, min_max, reaction_cell, value)

valin_leucine %>% 
  ggplot(aes(x = value, y = reaction_cell)) +
  geom_line(aes(color = cell_type)) +
  geom_point(aes(color =cell_type)) 

fatty_oxidation <- read.csv("~/Documents/GitHub/metastasis/data/plotting/fatty_acid_oxidation.csv", 
                            row.names = NULL, header = TRUE, sep =";")
fatty_oxidation$reaction_cell <- paste(fatty_oxidation$reaction, fatty_oxidation$cell_type)
fatty_oxidation %>% 
  ggplot(aes(x = value, y = reaction_cell)) +
  geom_line(aes(color = cell_type)) +
  geom_point(aes(color =cell_type)) 

fatty_oxidation$value <- gsub(",", ".", fatty_oxidation$value)

fatty_oxidation_low <- fatty_oxidation %>% filter(value<=1)
fatty_oxidation_low %>% 
  ggplot(aes(x = value, y = reaction_cell)) +
  geom_line(aes(color = cell_type)) +
  geom_point(aes(color =cell_type)) + scale_x_log10()

folat_metabolism <- read.csv("~/Documents/GitHub/metastasis/data/plotting/folat_metabolism.csv", 
                            row.names = NULL, header = TRUE, sep =";")
folat_metabolism$value <- gsub(",", ".", folat_metabolism$value) %>% as.numeric()
folat_metabolism$reaction_cell <- paste(folat_metabolism$reaction, folat_metabolism$cell_type)

folat_metabolism %>%
  ggplot(aes(x = value, y = reaction_cell)) +
  geom_line(aes(color = cell_type)) +
  geom_point(aes(color =cell_type))

folat_metabolism <- read.csv("~/Documents/GitHub/metastasis/data/plotting/pyrimidne_ metabolism.csv", 
                             row.names = NULL, header = TRUE, sep =";")
folat_metabolism$value <- gsub(",", ".", folat_metabolism$value) %>% as.numeric()
folat_metabolism$reaction_cell <- paste(folat_metabolism$reaction, folat_metabolism$cell_type)

folat_metabolism %>%
  ggplot(aes(x = value, y = reaction_cell)) +
  geom_line(aes(color = cell_type)) +
  geom_point(aes(color =cell_type))
  

folat_metabolism <- read.csv("~/Documents/GitHub/metastasis/data/plotting/glutamine_metabolism.csv", 
                             row.names = NULL, header = TRUE, sep =";")
folat_metabolism$value <- gsub(",", ".", folat_metabolism$value) %>% as.numeric()
folat_metabolism$reaction_cell <- paste(folat_metabolism$reaction, folat_metabolism$cell_type)

folat_metabolism %>%
  ggplot(aes(x = value, y = reaction_cell)) +
  geom_line(aes(color = cell_type)) +
  geom_point(aes(color =cell_type))


dumbbell_plot <- function(path, titletext){
  df <- read.csv(path, row.names = NULL, header = TRUE, sep =";")
  df$value <- gsub(",", ".", df$value) %>% as.numeric()
  df$reaction_cell <- paste(df$reaction, df$cell_type)
  
  plot <- df %>% 
    ggplot(aes(x = value, y = reaction_cell)) +
    geom_line(aes(color = cell_type)) +
    geom_point(aes(color =cell_type)) +
    labs(y = "reaction", x = "flux", title = titletext)
  return(plot)
}

folat_metabolism <- read.csv("~/Documents/GitHub/metastasis/data/plotting/oxidative_phosphorylation.csv", 
                             row.names = NULL, header = TRUE, sep =";")
folat_metabolism$value <- gsub(",", ".", folat_metabolism$value) %>% as.numeric()
folat_metabolism$reaction_cell <- paste(folat_metabolism$reaction, folat_metabolism$cell_type)

folat_metabolism %>%
  ggplot(aes(x = value, y = reaction_cell)) +
  geom_line(aes(color = cell_type)) +
  geom_point(aes(color =cell_type))

dumbbell_plot("~/Documents/GitHub/metastasis/data/plotting/oxidative_phosphorylation.csv",
              "Oxidative Phosporylation")
dumbbell_plot("~/Documents/GitHub/metastasis/data/plotting/glutamine_metabolism.csv",
              "Glutamine Metabolism")


