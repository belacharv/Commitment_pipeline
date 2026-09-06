################################################################################
### 05 PARAMETER RECOVERY PLOTS IN 2ND LAYER ###
################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate

# This R script serves just to visualize the result from 05 step in python

# loading libraries
library("tidyverse")
library(ggpubr)
library(this.path)

# load data and edit the table so it's representable in ggplot
# input:
  # filename = filename with results from 05 python step in analysis
  # true_n8 = true number of CD8 biased sequences
  # cd4cd4 = number of sequences in group CD8_0 hence CD4+CD4
# output:
  # table suitable for ggplot visual interpretation
load_data <- function(filename, true_n8,cd4cd4){
  data_a <- read_csv(filename) %>%
    mutate(error_n8 = abs(cd8_cd8 - true_n8))
  # commitment Pc source table
  n_fit_data_a <- data.frame(id = data_a$sim_id, count = data_a$fitted_Pc) %>%
    mutate(bias = "CD4",type = "a_simulated")
  # add observed value
  n_fit_data_a <- rbind(n_fit_data_a, 
                      tibble(id = 100,count=data_a$true_Pc[1],
                             bias = "CD4",type="observed"))

  # get true_cd4 (number of sequences truly biased towards CD4 lineage)
  n4_fit_data_a <- data.frame(id = data_a$sim_id, count = data_a$cd4_cd4) %>%
  mutate(bias = "CD4",type = "a_simulated")

  # Add true cd4_cd4 as a separate row
  n4_fit_data_a <- rbind(n4_fit_data_a,
                       tibble(id = 100,count = cd4cd4,
                              bias = "CD4",type="observed"))
  return(list(n_fit_data_a,n4_fit_data_a))
}

# edit formating of the ggplot
ggtheme <- function(){
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = "none"
  )
}

plots_function <- function(datatables, color, obs_pc,obs_cd4cd4, chain){
# Pc plot
  p1 <- ggplot(aes(x = bias, y = count), data = datatables[[1]]) +
  #geom_boxplot(fill = "grey70", alpha = 0.4) +
  geom_dotplot(binaxis = 'y',
               dotsize = 1.3,
               stackdir = 'center',
               aes(fill = type)) +
  labs(title = paste('Pc param recov ',chain,sep=""), y = 'Pc') +
  theme_light() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_hline(color = color,yintercept = obs_pc,linewidth = 1.2) +
  ggtheme() +
  ggplot2::annotate("text", x = 1.45, y = obs_pc, label = obs_pc, 
                    color = color, vjust = -1, size = 7) +
  scale_fill_manual(values = c("a_simulated" = "grey80", "observed" = color))
  #theme(axis.x = element_blank())

  p1
  #setwd(mainDir)
  ggsave(filename = paste("plots/05_plot1_",chain,"_parameter_recovery.png",sep=""), width = 12, height = 12, units = "cm")
  ggsave(filename = paste("plots/05_plot1_",chain,"_parameter_recovery.svg",sep=""), width = 12, height = 12, units = "cm")

# true_cd4 plot
  p3 <- ggplot(aes(x = bias, y = count), data = datatables[[2]]) +
  #geom_boxplot(fill = "grey70", alpha = 0.4) +
  geom_dotplot(binaxis = 'y',
               dotsize = 1.3,
               stackdir = 'center',
               aes(fill = type)) +
  labs(title = paste('cd4cd4 param recov ',chain,sep=""), y = 'CD4CD4 count') +
  theme_light() +
  theme(legend.title = element_blank(),
        legend.text = element_blank()) +
  ggtheme() +
  geom_hline(yintercept = obs_cd4cd4, color = color,linewidth=1.2) +
  ggplot2::annotate("text", x = 1.45, y = obs_cd4cd4, label = obs_cd4cd4, 
                    color = color, vjust = -1, size = 7) +
  scale_fill_manual(values = c("a_simulated" = "grey80", "observed" = color))

  p3
  ggsave(filename = paste("plots/05_plot2_",chain,"_parameter_recovery.png",sep=""), width = 12, height = 12, units = "cm")
  ggsave(filename = paste("plots/05_plot2_",chain,"_parameter_recovery.svg",sep=""), width = 12, height = 12, units = "cm")
}
