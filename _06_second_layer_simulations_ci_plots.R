######################################################################################
### 06 CONFIDENCE INTERVALS IN SIMULATIONS PLOTS ###
######################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate

# This script works with the 2nd layer (sequences present in exactly two samples)

# It loads a table created by previous 06 python script and creates a plot

# load libraries
library(ggrepel)
library(tidyverse)
library(this.path)

# edit formating
ggtheme <- function(){
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 20,colour = "black"),
    axis.title = element_text(size = 20),
    legend.title = element_blank()
    
  )
  
}


# generate the plot
# input:
  # col1 = confidence interval lower value
  # col2 = confidence interval upper value
  # col3 = real mean in changed simulation
  # col4 = expected mean
  # val = true experimental value
  # color = color palette
  # data = full table
first_plot <- function(col1,col2,col3,col4,val,color,data){
  data2 <- data.frame(pc = data$true_pc,
                      value = data[[col1]],group = col1)
  data3 <- data.frame(pc = data$true_pc,value = data[[col2]],group = col2)
  data4 <- data.frame(pc = data$true_pc,value = data[[col3]],group = col3)
  data5 <- data.frame(pc = data$true_pc,value = data[[col4]],group = col4)
  data_plot = rbind(data2,data3,data4,data5)
  
  p <- ggplot(aes(x=value,y=pc, color = group, label = value),data = data_plot)+
    theme_light() +
    geom_point(stat = "identity",
               position = "jitter") +
    geom_line(aes(x = val)) +
    scale_color_manual(values = c("#535050","#535050",color,"red"))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    ggtheme()+
    theme(panel.grid.minor = element_line(color = "gray90", linewidth = 0.25))
  return(data_plot)
}
