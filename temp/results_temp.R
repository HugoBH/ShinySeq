 library(tidyverse)

 
df <- read.csv("data/Library_Genotypes_v3.forwardreadsbees_counts.csv", ) 
ef <- df %>% as_tibble
ef$Sample
loci <- ef %>% 
  filter(!str_detect(Sample, regex("negative", ignore_case = TRUE))) %>% 
  select(-c(Raw.Reads, On.Target.Reads, X.On.Target, X.GT, IFI)) 
  
  df %>% select(starts_with("X"))
  # Return locus columns (those starting with 'X')
  a <- grep("^X", names(df), value = TRUE)
  length(a)
  
  colnames(df)
  
  df <- loci %>% select(-Sample)  
  data.frame(
    Sample = colnames(df),
    Mean_ReadDepth = apply(df, 2, mean),
    Variance_ReadDepth = apply(df, 2, var),
    SD_ReadDepth = apply(df, 2, sd),
    SE_ReadDepth = apply(df, 2, sd) / sqrt(ncol(loci))
  )
  
  mean(unlist(df %>% select(starts_with("X"))))
  
#SampleGenotypes  
loci$SampleSuccess = rowSums(loci > 10) / ncol(loci)
loci <- loci %>% select(Sample, SampleSuccess)

sum(loci$SampleSuccess >= 0.8)

ggplot(loci, aes(x = SampleSuccess)) +
  geom_histogram()

#SAMPLES WITH AT LEAST 90% GENOTYPE

loci$SampleSuccess[> 0.8] 


#average locus success
mean(rowSums(loci > 10) / ncol(loci))
mean(colSums(loci > 10) / nrow(loci))


threshold <- input$threshold
success <- rowMeans(df[loci] > threshold)
df$SuccessRate <- rowMeans(df[loci] > threshold)
df[, c("Sample", "SuccessRate")]


89, 0, 84, 90, 89, 81