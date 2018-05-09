#!/usr/bin/Rscript
#  applications/comp_init_plot.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.07.2018

## Plot the results of comp_init_loop.R
load('data/comp_init.RData')
require(ggplot2)
require(gridExtra)

ggdf$Cost <- as.numeric(as.character(ggdf$Cost))
ggdf$Time <- as.numeric(as.character(ggdf$Time))

#Plot the two boxplots
rownames(ggdf) <- NULL
pt <- ggplot(ggdf, aes(x = Type, y = Time)) +
    geom_boxplot() + 
    theme_light()
pc <- ggplot(ggdf, aes(x = Type, y = Cost)) +
    geom_boxplot() + 
    theme_light()

a <- arrangeGrob(pt, pc, nrow = 1)
ggsave('images/comp_init.pdf', a)
