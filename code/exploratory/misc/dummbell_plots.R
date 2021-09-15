library(dumbbell)
reaction <- c('HMR_0735', 'HMR_0753', 'HMR_0758', 'RE2675C', 'CRMte',
             'HMR_0735', 'HMR_0753', 'HMR_0758', 'RE2675C', 'CRMte',
             'HMR_0735', 'HMR_0753', 'HMR_0758', 'RE2675C', 'CRMte')
min_value <- c(0.44205,-383.0356, 181.1662, 176.3169, 4.8063,
              0.3959, -309.2845, 119.125, 114.7219, 4.3021,
              0.41919, -350.1657, 155.7064, 150.8956, 4.5552)
              
max_value <- c(0.6261, -176.3169, 388.7539, 383.0356, 6.8036,
              0.56027, -114.7219, 314.5242, 309.2845, 6.0883,
              0.60707, -150.8956, 355.3548, 350.1657, 6.5968)
min_max <- 


cell_type <- c('breast', 'breast', 'breast', 'breast', 'breast',
              'brain', 'brain', 'brain', 'brain', 'brain',
              'lung', 'lung', 'lung', 'lung', 'lung')
sphingolipid_df <- data.frame(reaction, min_value, max_value, cell_type)

dummbbell(xdf = sphingolipid_df, id="reaction", lab1 = "min_value", lab2 = "max_value")

value <- c(0.44205,-383.0356, 181.1662, 176.3169, 4.8063,
               0.3959, -309.2845, 119.125, 114.7219, 4.3021,
               0.41919, -350.1657, 155.7064, 150.8956, 4.5552,
               0.6261, -176.3169, 388.7539, 383.0356, 6.8036,
               0.56027, -114.7219, 314.5242, 309.2845, 6.0883,
               0.60707, -150.8956, 355.3548, 350.1657, 6.5968)
reaction <- c('HMR_0735', 'HMR_0753', 'HMR_0758', 'RE2675C', 'CRMte',
              'HMR_0735', 'HMR_0753', 'HMR_0758', 'RE2675C', 'CRMte',
              'HMR_0735', 'HMR_0753', 'HMR_0758', 'RE2675C', 'CRMte',
              'HMR_0735', 'HMR_0753', 'HMR_0758', 'RE2675C', 'CRMte',
              'HMR_0735', 'HMR_0753', 'HMR_0758', 'RE2675C', 'CRMte',
              'HMR_0735', 'HMR_0753', 'HMR_0758', 'RE2675C', 'CRMte')
min_max <-  c(rep('min_value', 15), rep('max_value', 15))
cell_type <- c(rep('breast', 5), rep('brain', 5), rep('lung', 5),
               rep('breast', 5), rep('brain', 5), rep('lung', 5))

df %>% 
  ggplot(aes(x = value, y = comdined)) +
  geom_line() +
  geom_point(aes(color =reaction)) 

dumbbell(v1 = sphingolipid_df$min_value, v2= sphingolipid_df$max_value, 
         group = sphingolipid_df$cell_type, text =TRUE, labels = sphingolipid_df$reaction,
         segments =TRUE)
df = data.frame(reaction, value, min_max, cell_type)
dumbbell(xdf = sphingolipid_df, id = "reaction", key = "cell_type", column1 = 'min_value', column2='max_value')






