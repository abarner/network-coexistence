
#### Source simulation functions ####

rm(list = ls())

source("foodwebs/intertidal_foodweb/intertidal_foodweb_model_modified.R")
source("foodwebs/intertidal_foodweb/do_intertidal_partition.R")


#### Load packages ####

library(tidyverse); library(grid)

#### load functions for plotting ####

# multiple plot function
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#### run overall partition, low, high supplies ####

# have 6 scenarios to run
# 100 x each

# want to run across full larval scenarios (redo this)
larval_scenarios_for_full_run <- tribble (
  ~scenario, ~species, ~supply_level,
  1, "B", "high",
  1, "C", "high",
  1, "L", "high",
  1, "P", "high",
  2, "B", "low",
  2, "C", "low",
  2, "L", "low",
  2, "P", "low",
  3, "B", "high",
  3, "C", "low",
  3, "L", "low",
  3, "P", "low",
  4, "B", "low",
  4, "C", "high",
  4, "L", "low",
  4, "P", "low",
  5, "B", "low",
  5, "C", "low",
  5, "L", "high",
  5, "P", "low",
  6, "B", "low",
  6, "C", "low",
  6, "L", "low",
  6, "P", "high"
)

larval_supply_full <- tribble (
  ~species, ~mean_recruit, ~variance_recruit, ~mean_supply_level, ~variance_supply_level,
  "B", 90000, sqrt(4.6*10^9), "high", "low",
  "B", 90000, sqrt(3.24*10^10), "high", "med",
  "B", 90000, sqrt(9.9*10^10), "high", "high",
  
  "B", 50000, sqrt(1.41*10^10), "med", "low", # not sure why variance higher than "med"
  "B", 50000, sqrt(1*10^10), "med", "med",
  "B", 50000, sqrt(3*10^10), "med", "high",
  
  "B", 6000, sqrt(2.025*10^7), "low", "low",
  "B", 6000, sqrt(1.44*10^8), "low", "med",
  "B", 6000, sqrt(4.41*10^8), "low", "high",
  
  "C", 70000, sqrt(2.75*10^9), "high", "low",
  "C", 70000, sqrt(1.96*10^10), "high", "med",
  "C", 70000, sqrt(6*10^10), "high", "high",
  
  "C", 30000, sqrt(5.1*10^8), "med", "low",
  "C", 30000, sqrt(3.6*10^9), "med", "med",
  "C", 30000, sqrt(1.1*10^10), "med", "high",
  
  "C", 6000, sqrt(2.025*10^7), "low", "low",
  "C", 6000, sqrt(1.44*10^8), "low", "med",
  "C", 6000, sqrt(4.41*10^8), "low", "high",
  
  "L", 3000, sqrt(3.8*10^6), "high", "low", 
  "L", 3000, sqrt(2.8*10^7), "high", "med", 
  "L", 3000, sqrt(8.1*10^7), "high", "high", 
  
  "L", 2400, sqrt(2.4*10^6), "med", "low",
  "L", 2400, sqrt(1.7*10^7), "med", "med",
  "L", 2400, sqrt(5.2*10^7), "med", "high",
  
  "L", 200, sqrt(16900), "low", "low",
  "L", 200, sqrt(1.2*10^5), "low", "med",
  "L", 200, sqrt(3.6*10^5), "low", "high",
  
  "P", 6873, sqrt(3.02*10^7), "high", "low",
  "P", 6873, sqrt(1.2*10^8), "high", "med",
  "P", 6873, sqrt(1.4*10^8), "high", "high",
  
  "P", 3800, sqrt(9.2*10^6), "med", "low",
  "P", 3800, sqrt(3.7*10^7), "med", "med",
  "P", 3800, sqrt(4.2*10^7), "med", "high",
  
  "P", 727, sqrt(3.4*10^5), "low", "low",
  "P", 727, sqrt(1.3*10^6), "low", "med",
  "P", 727, sqrt(1.5*10^6), "low", "high",
)

larval_scenarios_for_full_run %>%
  rename(mean_supply_level = "supply_level") %>%
  mutate(variance_supply_level = "low") %>%
  left_join(larval_supply_full) %>%
  select(-mean_supply_level, -variance_supply_level) %>%
  gather(key = variable, value = value, mean_recruit, variance_recruit) %>%
  unite(temp, species, variable) %>%
  spread(temp, value) -> larval_scenarios_input # need it to be named this




### main text run ####

# do in batches so nothing lost if there's a failure
csv_names <- c("001_100", "101_200", "201_300",
               "301_400", "401_500")

for (i in 1:5) {
  print(paste0("i = ", i))
  sim_output_list <- vector(mode = "list", length = 6)
  names(sim_output_list) <- c("high", "low", "balanus_high", "chthamalus_high",
                              "limpets_high", "pisaster_high")
  
  for (k in 1:length(sim_output_list)) {
    print(paste0("k = ", k))
    sim_output_list[[k]] <- do.larval.supply.simulation(k = k, n_sim = 100)
  }
  
  sim_output_list %>%
    bind_rows(.id = "larval_scenario") %>%
    write_csv(path = paste0("final_larval_maintext_results_", csv_names[i], ".csv"))
}

#### read in output ####

output_files <- list.files(path = "foodwebs/intertidal_foodweb/", 
                           pattern = "final_larval_maintext_results_", 
                           full.names = TRUE)
sim_output_list_all <- lapply(X = output_files, FUN = read_csv)
names(sim_output_list_all) <- paste0("rep_", csv_names)

# get data frame of all results
sim_output_list_all %>%
  bind_rows(.id = "rep") %>%
  # need to clear up unique ids
  mutate(rep_num = as.numeric(str_sub(rep, start = -7, end = -5)) - 1) %>%
  mutate(simulation_loop = simulation_loop + rep_num) %>%
  select(-rep_num, -rep) -> sim_output_df
  
sim_output_df %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs))


#### test plot for main text (no panels) ####

# get raw points
sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  mutate(mean_cs = coexistence_strength) %>%
  mutate(species = fct_recode(factor(species),
                              `Limpets` = "limpets",
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli")) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"))) -> sim_output_jitter


sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  mutate(species = fct_recode(factor(species),
                              `Limpets` = "limpets",
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli")) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = mean_cs)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition), size = .25, color = "black") + 
  geom_jitter(data = sim_output_jitter, 
              aes(x = coexistence_partition, y = mean_cs), 
              alpha = .25, size = .1) + 
  geom_errorbar(aes(ymin = mean_cs - sd_cs, ymax = mean_cs + sd_cs), color = "black", width = 0.2) +
  facet_grid(larval_scenario ~ species, scales = "free") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  scale_color_manual(values=fad) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(face = "italic")) + 
  scale_x_discrete(labels=xlab) + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare")

### maintext as boxplot instead #### 


## test
png(filename = "intertidal_boxplot.png", width = 7, height = 6, units = "in", res = 300)
sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  mutate(species = fct_recode(factor(species),
                              `Limpets` = "limpets",
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli")) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = coexistence_strength)) +
  geom_boxplot(aes(fill = coexistence_partition), outlier.shape = NA, 
               notch = TRUE, coef = 1) + 
  facet_grid(larval_scenario ~ species, scales = "free") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  scale_color_manual(values=fad) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(face = "italic")) + 
  scale_x_discrete(labels=xlab) + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare") #+
  ylim(-10, 5)
dev.off()


## 


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols<-c(1,6,3,8,2)
cbPalette <- cbbPalette[cols]

fad<-cbPalette
xlab=c(expression("r"[i]-"r"[r]) ,
       expression(Delta[i]^0),
       expression(Delta[i]^P),
       expression(Delta[i]^E),
       expression(Delta[i]^{E*P}))

tibble(larval_scenario = rep(c("low", "high"), times = 2)) %>%
  mutate(coexistence_partition = "r_bar") %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  # add a little space for labels
  mutate(coexistence_strength = rep(5)) -> range_supply
# panel labels
tibble(species = c("balanus_glandula", "balanus_glandula", 
                   "chthamalus_dalli", "chthamalus_dalli", 
                   "limpets", "limpets"),
       larval_scenario = rep(c("low", "high"), times = 3),
       panel_label = c("A)","D)","B)","E)","C)","F)")) %>%
  mutate(coexistence_partition = "r_bar") %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(mean_cs = max(range_supply$coexistence_strength)) %>%
  mutate(species = fct_recode(factor(species),
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli",
                              Limpets = "limpets")) -> label_supply

# left panel (b glandula)
sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  filter(species == "balanus_glandula") %>%
  mutate(species = fct_recode(factor(species),
                              `Limpets` = "limpets",
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli")) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = coexistence_strength)) +
  geom_boxplot(aes(fill = coexistence_partition), outlier.shape = NA, 
               notch = TRUE, coef = 1) + 
  facet_grid(larval_scenario ~ species, scales = "free") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  theme_bw() +
  scale_x_discrete(labels=xlab) + 
  ylim(-10, 5) +
  geom_text(data = filter(label_supply, species == "Balanus glandula"), 
            aes(x = coexistence_partition, y = 5, 
                label = panel_label),
            position = position_nudge(x = -0.25, y = 0.5))+
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(face = "italic", size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) + 
  xlab("") +
  ylab("Growth rate when rare") -> a

# middle panel (c dalli)
sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  filter(species == "chthamalus_dalli") %>%
  mutate(species = fct_recode(factor(species),
                              `Limpets` = "limpets",
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli")) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = coexistence_strength)) +
  geom_boxplot(aes(fill = coexistence_partition), outlier.shape = NA, 
               notch = TRUE, coef = 1) +                           
  facet_grid(larval_scenario ~ species) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  geom_text(data = filter(label_supply, species == "Chthamalus dalli"), 
            aes(x = coexistence_partition, y = 10, 
                label = panel_label),
            position = position_nudge(x = -0.25, y = 0.5))+
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(face = "italic", size = 14),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14)) + 
  scale_x_discrete(labels=xlab) + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare") +
  ylim(-10, 5) -> b

# right panel (limpets)
sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  filter(species == "limpets") %>%
  mutate(species = fct_recode(factor(species),
                              `Limpets` = "limpets",
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli")) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = coexistence_strength)) +
  geom_boxplot(aes(fill = coexistence_partition), outlier.shape = NA, 
               notch = TRUE, coef = 1) + 
  facet_grid(larval_scenario ~ species) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  geom_text(data = filter(label_supply, species == "Limpets"), 
            aes(x = coexistence_partition, y = 10, 
                label = panel_label),
            position = position_nudge(x = -0.25, y = 0.5))+
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.text = element_text( size = 14)) + 
  scale_x_discrete(labels=xlab) + 
  xlab("") +
  ylim(-10, 5) -> c

png(filename = "intertidal_maintext.png", width = 8, height = 6, units = "in", res = 300)
multiplot(a, b, c,
          layout = matrix(c(1,2,3,1,2,3), nrow=2, byrow=TRUE))
dev.off()








sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  mutate(species = fct_recode(factor(species),
                              `Limpets` = "limpets",
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli")) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = coexistence_strength)) +
  geom_boxplot(aes(fill = coexistence_partition), outlier.shape = NA, 
               notch = TRUE, coef = 1) + 
  facet_grid(larval_scenario ~ species, scales = "free") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  scale_color_manual(values=fad) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(face = "italic")) + 
  scale_x_discrete(labels=xlab) + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare") #+
ylim(-10, 5)




#### how many runs very positive or very negative ? ####

sim_output_df %>%
  filter(coexistence_partition == "r_bar") %>%
  mutate(cs_bin = ifelse(coexistence_strength > 10, "greater_10",
                         ifelse(coexistence_strength < -10, "negative_10", "fine"))) %>%
  group_by(cs_bin) %>%
  summarise(n_group = n() / 9000) # < 3%


#### test plot for main text (2 panels) ####

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols<-c(1,6,3,8,2)
cbPalette <- cbbPalette[cols]

fad<-cbPalette
xlab=c(expression("r"[i]-"r"[r]) ,
       expression(Delta[i]^0),
       expression(Delta[i]^P),
       expression(Delta[i]^E),
       expression(Delta[i]^{E*P}))

# left panel (barnacles)
sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  filter(species != "limpets") %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  mutate(species = fct_recode(factor(species),
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli")) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = mean_cs)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition), size = .25, color = "black") + 
  geom_errorbar(aes(ymin = mean_cs - se_cs, ymax = mean_cs + se_cs), color = "black", width = 0.2) +
  facet_grid(larval_scenario ~ species, scales = "free") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(face = "italic")) + 
  scale_x_discrete(labels=xlab) + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare") -> a

# right panel (limpets)
# want to add text over ones that we can't see
sim_limpets_label <- tibble(larval_scenario = c("low", "low", "high"),
                            coexistence_partition = c("delta_0", "delta_p", "delta_c"),
                            to_label = "yes",
                            label_hjust = c(0.75, 0.25, 0.5),
                            limpet_label = c("+", "+", "-"))

sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  filter(species == "limpets") %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  left_join(sim_limpets_label) %>%
  # mutate(limpet_label = ifelse(is.na(to_label), "",
  #                              paste0(round(mean_cs, 3), "±", 
  #                                     format(se_cs, scientific=TRUE, digits = 1)))) %>%
  mutate(species = "Limpets") %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = mean_cs)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition), size = .25, color = "black") + 
  geom_errorbar(aes(ymin = mean_cs - se_cs, ymax = mean_cs + se_cs), 
                color = "black", width = 0.2) +
  # geom_text(aes(label=limpet_label, hjust = label_hjust),
  #           position=position_dodge(width=0.9), vjust=-1.5, size = 3) +
  facet_grid(larval_scenario ~ species, scales = "free") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank()) + 
  scale_x_discrete(labels=xlab) + 
  xlab("") -> b

png(filename = "intertidal_maintext.png", width = 9, height = 6, units = "in", res = 300)
multiplot(a, b, 
          layout = matrix(c(1,1,2,1,1,2), nrow=2, byrow=TRUE))
dev.off()


#### main text ####


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols<-c(1,6,3,8,2)
cbPalette <- cbbPalette[cols]

fad<-cbPalette
xlab=c(expression("r"[i]-"r"[r]) ,
       expression(Delta[i]^0),
       expression(Delta[i]^P),
       expression(Delta[i]^E),
       expression(Delta[i]^{E*P}))

# get max by supply
sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(mean_max = mean_cs + se_cs,
         mean_min = mean_cs - se_cs) %>%
  ungroup() %>%
  group_by(larval_scenario) %>%
  summarise(max_by = max(mean_max),
            min_by = min(mean_min)) %>%
  gather(max_by, min_by, key = key, value = mean_cs) %>%
  select(-key) %>%
  mutate(coexistence_partition = "r_bar") %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  # add a little space for labels
  mutate(mean_cs = ifelse(mean_cs == max(mean_cs), 
                          mean_cs + 0.6, mean_cs)) -> range_supply
# panel labels
tibble(species = c("balanus_glandula", "balanus_glandula", 
                   "chthamalus_dalli", "chthamalus_dalli", 
                   "limpets", "limpets"),
       larval_scenario = rep(c("low", "high"), times = 3),
       panel_label = c("A)","D)","B)","E)","C)","F)")) %>%
  mutate(coexistence_partition = "r_bar") %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(mean_cs = max(range_supply$mean_cs)) %>%
  mutate(species = fct_recode(factor(species),
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli",
                              Limpets = "limpets")) -> label_supply

# left panel (b glandula)
sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  filter(species == "balanus_glandula") %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  mutate(species = fct_recode(factor(species),
                              `Balanus glandula` = "balanus_glandula")) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = mean_cs)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition), size = .25, color = "black") + 
  #geom_errorbar(aes(ymin = mean_cs - se_cs, ymax = mean_cs + se_cs), color = "black", width = 0.2) +
  geom_errorbar(aes(ymin = mean_cs - sd_cs, ymax = mean_cs + sd_cs), color = "black", width = 0.2) +
  # geom_point(data = range_supply, aes(x = coexistence_partition,
  #                                     y = mean_cs+0.1), color = "white", size = 0) +
  facet_grid(larval_scenario ~ species) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  geom_text(data = filter(label_supply, species == "Balanus glandula"), 
            aes(x = coexistence_partition, y = 10, 
                label = panel_label),
            position = position_nudge(x = -0.25, y = 0.5))+
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(face = "italic", size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) + 
  scale_x_discrete(labels=xlab) + 
  xlab("") +
  ylab("Growth rate when rare") -> a

# middle panel (c dalli)
sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  filter(species == "chthamalus_dalli") %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  mutate(species = fct_recode(factor(species),
                              `Chthamalus dalli` = "chthamalus_dalli")) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = mean_cs)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition), size = .25, color = "black") + 
  #geom_errorbar(aes(ymin = mean_cs - se_cs, ymax = mean_cs + se_cs), color = "black", width = 0.2) +
  geom_errorbar(aes(ymin = mean_cs - sd_cs, ymax = mean_cs + sd_cs), color = "black", width = 0.2) +
  # geom_point(data = range_supply, aes(x = coexistence_partition,
  #                                     y = mean_cs+0.1), color = "white", size = 0) +
  facet_grid(larval_scenario ~ species) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  geom_text(data = filter(label_supply, species == "Chthamalus dalli"), 
            aes(x = coexistence_partition, y = 10, 
                label = panel_label),
            position = position_nudge(x = -0.25, y = 0.5))+
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(face = "italic", size = 14),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14)) + 
  scale_x_discrete(labels=xlab) + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare") -> b

# right panel (limpets)
sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  filter(species == "limpets") %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  left_join(sim_limpets_label) %>%
  mutate(species = "Limpets") %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "high"
                                             )),
                                      `All low` = "low",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = mean_cs)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition), size = .25, color = "black") + 
  # geom_errorbar(aes(ymin = mean_cs - se_cs, ymax = mean_cs + se_cs), 
  #               color = "black", width = 0.2) +
  geom_errorbar(aes(ymin = mean_cs - sd_cs, ymax = mean_cs + sd_cs), 
                color = "black", width = 0.2) +
  # geom_point(data = range_supply, aes(x = coexistence_partition,
  #                                     y = mean_cs), color = "white") +
  facet_grid(larval_scenario ~ species) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  geom_text(data = filter(label_supply, species == "Limpets"), 
            aes(x = coexistence_partition, y = 10, 
                label = panel_label),
            position = position_nudge(x = -0.25, y = 0.5))+
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.text = element_text( size = 14)) + 
  scale_x_discrete(labels=xlab) + 
  xlab("") -> c

png(filename = "intertidal_maintext.png", width = 9, height = 6, units = "in", res = 300)
multiplot(a, b, c,
          layout = matrix(c(1,2,3,1,2,3), nrow=2, byrow=TRUE))
dev.off()


#### plot full paneled figure for supplement - all larval supply scenario ####

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols<-c(1,6,3,8,2)
cbPalette <- cbbPalette[cols]

fad<-cbPalette
xlab=c(expression("r"[i]-"r"[r]) ,
       expression(Delta[i]^0),
       expression(Delta[i]^P),
       expression(Delta[i]^E),
       expression(Delta[i]^{E*P}))

png(filename = "intertidal_supplement.png", width = 12, height = 6, units = "in", res = 300)
sim_output_df %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  mutate(species = fct_recode(factor(species),
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli", 
                              Limpets = "limpets")) %>%
  mutate(larval_scenario = fct_recode(factor(larval_scenario,
                                             levels = c(
                                               "low", "balanus_high",
                                               "chthamalus_high", "limpets_high",
                                               "pisaster_high", "high"
                                             )),
                                      `All low` = "low",
                                      `Balanus high` = "balanus_high",
                                      `Chthamalus high` = "chthamalus_high",
                                      `Limpets high` = "limpets_high", 
                                      `Pisaster high` = "pisaster_high",
                                      `All high` = "high")) %>%
  mutate(coexistence_partition = factor(coexistence_partition,
                                        levels = c(
                                          "r_bar",
                                          "delta_0", 
                                          "delta_p",
                                          "delta_c",
                                          "delta_cp"
                                        ))) %>% 
  ggplot(aes(x = coexistence_partition, y = mean_cs)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition)) + 
  #geom_errorbar(aes(ymin = mean_cs - se_cs, ymax = mean_cs + se_cs), color = "black", width = 0.2) +
  geom_errorbar(aes(ymin = mean_cs - sd_cs, ymax = mean_cs + sd_cs), color = "black", width = 0.2) +
  facet_grid(species ~ larval_scenario, scales = "free") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_discrete(labels=xlab) + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare")
dev.off()


#### run & plot time series for supplement ####

# use the high & low scenarios from above

n_sim <- 500

test_ab_list <- vector(mode = "list", length = 2)
names(test_ab_list) <- c("high", "low")

for (k in 1:2) {
  
  simulation_loop_output_tmp <- vector(mode = "list", length = n_sim)
  
  B.mean_input <- larval_scenarios_input$B_mean_recruit[k]
  B.stdev_input <- larval_scenarios_input$B_variance_recruit[k]
  C.mean_input <- larval_scenarios_input$C_mean_recruit[k]
  C.stdev_input <- larval_scenarios_input$C_variance_recruit[k]
  L.mean_input <- larval_scenarios_input$L_mean_recruit[k]
  L.stdev_input <- larval_scenarios_input$L_variance_recruit[k]
  P.mean_input <- larval_scenarios_input$P_mean_recruit[k]
  P.stdev_input <- larval_scenarios_input$P_variance_recruit[k]
  
  for (i in 1:n_sim) {
    print(i)
    simulation_loop_output_tmp[[i]] <- do.intertidal.simulation(years_set = 100,
                                                                B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                                C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                                L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                                P.mean = P.mean_input, P.stdev = P.stdev_input)
  }
  
  test_ab_list[[k]] <- simulation_loop_output_tmp
}

test_ab_list$high %>%
  map(top_n, n = 100, wt = time) %>%
  bind_rows(.id = "rep") %>%
  mutate(larval_supply = "High") -> test_ab_list_high
test_ab_list$low %>%
  map(top_n, n = 100, wt = time) %>%
  bind_rows(.id = "rep") %>%
  mutate(larval_supply = "Low") -> test_ab_list_low

png(filename = "intertidal_timeseries.png", width = 6, height = 8, units = "in", res = 300)
bind_rows(test_ab_list_high, test_ab_list_low) %>%
  gather(balanus_glandula:pisaster_ochraceus, key = "species", value = "abundance") %>%
  mutate(species = fct_recode(factor(species),
                              `Balanus glandula` = "balanus_glandula",
                              `Chthamalus dalli` = "chthamalus_dalli", 
                              Limpets = "limpets",
                              `Pisaster ochraceus` = "pisaster_ochraceus",
                              Whelks = "whelks")) %>%
  ggplot(aes(x = time, y = abundance)) +
  #geom_point(size = 2) +
  geom_line(aes(group = rep), alpha = 0.25) +
  facet_grid(species ~ larval_supply, scales = "free") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(x = "Time (months)", y = "Abundance")
dev.off()  

# store the data from the two runs
bind_rows(test_ab_list_high, test_ab_list_low) %>%
  write_csv("timeseries_final100_2larvalsupply.csv")








