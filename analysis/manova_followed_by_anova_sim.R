library(tidyverse)
library(mvtnorm)
set.seed(1234)
manova_fit <- function(alpha=0.05) {
  n= 100
  sigma <- diag(x=1,nrow=3)
  mu <- c(0,1,-2)
  group1 <- data.frame(rmvnorm(n,mu,sigma))
  group1$group <-  rep("1",n)
  colnames(group1) <- c('x1','x2','x3','group')
  
  group2 <- data.frame(rmvnorm(n,mu,sigma))
  group2$group <-  rep("2",n)
  colnames(group2) <- c('x1','x2','x3','group')
  
  group3 <- data.frame(rmvnorm(n,mu,sigma))
  group3$group <-  rep("3",n)
  colnames(group3) <- c('x1','x2','x3','group')
  
  group_dat <- data.frame(rbind(group1,group2,group3))
  
  man_fit <- manova(cbind(x1,x2,x3) ~ group, data = group_dat)
  man_sum <- summary(man_fit)
  p_val_man <- man_sum$stats["group","Pr(>F)"]
  
  if(p_val_man < alpha){
    
    aov1 <- aov(x1 ~ group, data = group_dat)
    aov2 <- aov(x2 ~ group, data= group_dat)
    aov3 <- aov(x3 ~ group, data = group_dat)
    
    p_val1 <- summary(aov1)[[1]][["Pr(>F)"]][1]
    p_val2 <- summary(aov2)[[1]][["Pr(>F)"]][1]
    p_val3 <- summary(aov3)[[1]][["Pr(>F)"]][1]
  }
  else{
    p_val1 <- -9
    p_val2 <- -9
    p_val3 <- -9
  }
  
  return_dat <- list(
    p_val_man=p_val_man,
    p_val1=p_val1,
    p_val2=p_val2,
    p_val3=p_val3
  )
  
  return(return_dat)
}


k = 1000

temp <- replicate(k,manova_fit(),
                  simplify = FALSE)
results_df <- do.call(rbind, lapply(temp, as.data.frame))
saveRDS(results_df, "analysis/results_df.rds")

## Histogram of MANOVA p-vales
man_pval_dist <- ggplot(data = results_df, aes(x=p_val_man)) +
  geom_histogram(binwidth= 0.1, fill = "#c4bf33", color = "#c4bf45", alpha = 0.75) +
  geom_vline(xintercept = 0.05, color = "#ff3319", linetype = "dashed", linewidth = 1) +
  ggtitle("Distribution of MANOVA P-Values") +
  labs(x = "P-values") +
  theme_minimal()
saveRDS(man_pval_dist, "analysis/man_pval_dist.rds")

## Univarate ANOVA hists
uni_dat <- cbind(results_df$p_val1[results_df$p_val1>0], replicate(length(results_df$p_val1[results_df$p_val1>0]), 1))
temp1 <- cbind(results_df$p_val2[results_df$p_val2>0], replicate(length(results_df$p_val2[results_df$p_val2>0]), 2))
temp2 <- cbind(results_df$p_val3[results_df$p_val3>0], replicate(length(results_df$p_val3[results_df$p_val3>0]), 3))
uni_dat <- data.frame(rbind(uni_dat,temp1,temp2))
colnames(uni_dat) <- c("p_value","group")
uni_dat <- uni_dat |>
  mutate(group = factor(group,
                        levels = c(1,2,3),
                        labels = c("X1","X2", "X3")))

anova_pvals_dist <- ggplot(data = uni_dat, aes(x= p_value, fill = group)) +
  geom_histogram(binwidth= 0.1, fill = "#c4bf33", color = "#c4bf45", alpha = 0.75) +
  facet_wrap(~group, ) +
  geom_vline(xintercept = 0.05, color = "#ff3319", linetype = "dashed", linewidth = 1) +
  ggtitle("Histogram of P-Values from three Univarate ANOVAs") +
  labs(x = "P-values") +
  theme_minimal()

saveRDS(anova_pvals_dist, "analysis/anova_pvals_dist.rds")

sum(results_df$p_val_man[results_df$p_val_man>0] < 0.05) / length(results_df$p_val_man[results_df$p_val_man>0])
sum(results_df$p_val1[results_df$p_val1>0] < 0.05) / length(results_df$p_val1[results_df$p_val1>0])
sum(results_df$p_val2[results_df$p_val2>0] < 0.05) / length(results_df$p_val2[results_df$p_val1>0])
sum(results_df$p_val3[results_df$p_val3>0] < 0.05) / length(results_df$p_val3[results_df$p_val1>0])

sum(results_df$p_val_man[results_df$p_val_man>0] < 0.05) 
sum(results_df$p_val1[results_df$p_val1>0] < 0.05)
sum(results_df$p_val2[results_df$p_val2>0] < 0.05)
sum(results_df$p_val3[results_df$p_val3>0] < 0.05)
# heavily right skewed, clearly not normal, and very large type one errors

#### But what if bonferroni after manova #######
sum(results_df$p_val1[results_df$p_val1>0] < 0.0125) / length(results_df$p_val1[results_df$p_val1>0])
sum(results_df$p_val2[results_df$p_val2>0] < 0.0125) / length(results_df$p_val2[results_df$p_val1>0])
sum(results_df$p_val3[results_df$p_val3>0] < 0.0125) / length(results_df$p_val3[results_df$p_val1>0])


