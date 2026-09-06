################################################################################
### 08 DEPTH PLOT WITH 2 EQ IN ALL LAYERS #####
################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fatelibrary(tidyverse)

# A script to create a plots from results of 08 python script

# loading libraries
library(tidyverse)
library(this.path)
library(ggnewscale)
library(ggbeeswarm)
################################################################################
# function to load datatable, edit it and generate plot
load_and_plot <- function(depths.a,col,chain){
  results_pc = data.frame(depth = depths.a$depth, 
                        pc = depths.a$Pc_val, 
                        n4 = depths.a$n4,
                        eq = depths.a$equation,
                        group = rep(chain,nrow(depths.a)))

  # edit formating
  ggtheme <- function(){
    theme(
      text = element_text(size = 20),
      axis.text = element_text(size = 20,colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 20),
    
    )
  }
  results_pc1 <- results_pc %>%
    dplyr::filter(pc > 0.7) %>%
    mutate(pc_rounded = round(pc, 3)) %>%
    group_by(group, pc_rounded,depth, eq) %>%
    slice_max(pc, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-pc_rounded)

  ggplot(aes(x=group,y=pc,fill=group),data = results_pc1) +
    theme_light() +
    geom_violin(aes(fill = group),alpha = 0.3, 
              quantiles = c(0.5), quantile.linetype = "solid", quantile.linewidth = 2) +
    geom_beeswarm(data = results_pc1 %>% dplyr::filter(group == chain), 
                aes(color=depth), size=4, cex=2) +
    scale_color_gradient(low = col[2], high = col[3]) +
    stat_summary(fun = median, geom = "crossbar", 
               width = 0.5, color = "black", linewidth = 0.8, fatten = 0) +
    ggtheme() +
    stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 3)), 
               vjust = -0.5,hjust = -0.5, size = 8, fontface = "bold")+
    scale_fill_manual(values=col[1]) +
    theme(legend.title = element_blank()) +
    guides(fill = "none")
  ggsave(paste("08_",chain,"_layers_pc_plot.svg",sep=""),width = 16,height = 16, units = "cm")
  ggsave(paste("08_",chain,"_layers_pc_plot.png",sep=""),width = 16,height = 16, units = "cm")
}