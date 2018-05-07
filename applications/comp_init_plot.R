#!/usr/bin/Rscript
#  applications/comp_init_plot.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.07.2018

## Plot the results of comp_init_loop.R
load('data/comp_init.RData')
require(ggplot2)
require(gridExtra)

#Plot the two boxplots
ggdf <- as.data.frame(cbind(1, rand_times, rand_costs))
pt <- ggplot(ggdf, aes(x = V1, y= rand_times)) +
    geom_boxplot() + 
    geom_hline(yintercept = ldapca_time) + 
    theme_light()
pc <- ggplot(ggdf, aes(x = V1, y= rand_costs)) +
    geom_boxplot() + 
    geom_hline(yintercept = ldapca_cost) + 
    theme_light()

a <- arrangeGrob(pt, pc, nrow = 1)
ggsave('images/comp_init.pdf', a)
