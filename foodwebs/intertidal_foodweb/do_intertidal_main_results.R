
#### Source simulation functions ####

rm(list = ls())

source("foodwebs/intertidal_foodweb/intertidal_foodweb_model_modified.R")
source("foodwebs/intertidal_foodweb/do_intertidal_partition.R")


#### Load packages ####

library(tidyverse)


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

sim_output_list <- vector(mode = "list", length = 6)
names(sim_output_list) <- c("high", "low", "balanus_high", "chthamalus_high",
                            "limpets_high", "pisaster_high")

for (k in 1:length(sim_output_list)) {
  print(c("LOOP = ", k))
  sim_output_list[[k]] <- do.larval.supply.simulation(k = k, n_sim = 1)
}

sim_output_list %>%
  bind_rows(.id = "larval_scenario") %>%
  write_csv("final_larval_maintext_results.csv")

sim_output_list %>%
  bind_rows(.id = "larval_scenario") %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs))

# check to see if all the runs are within reasonable bounds
sim_output_list %>%
  bind_rows(.id = "larval_scenario") %>%
  filter(coexistence_strength > 10)

# want those specifically to run again one more time
sim_output_list %>%
  bind_rows(.id = "larval_scenario") %>%
  filter(coexistence_strength > 10) %>%
  distinct(larval_scenario, simulation_loop) %>%
  group_by(larval_scenario) %>%
  summarise(n_loops = n()) -> scenarios_to_run_again

# which scenario(s) exactly to drop?
sim_output_list %>%
  bind_rows(.id = "larval_scenario") %>%
  filter(coexistence_strength > 10) %>%
  distinct(larval_scenario, simulation_loop) %>%
  unite(larval_scenario, simulation_loop, col = "scenario_combn") %>%
  pull(scenario_combn) -> scenario_to_drop

scenario_num <- match(scenarios_to_run_again$larval_scenario, names(sim_output_list))
scenario_rep <- scenarios_to_run_again$n_loops

# take input from above and overwrite
larval_scenarios_input %>%
  filter(scenario %in% scenario_num) -> larval_scenarios_input

sim_output_list_redo <- vector(mode = "list", length = nrow(larval_scenarios_input))
names(sim_output_list_redo) <- scenarios_to_run_again$larval_scenario

for (l in 1:length(sim_output_list_redo)) {
  print(c("LOOP = ", l))
  sim_output_list_redo[[l]] <- do.larval.supply.simulation(k = l, n_sim = scenario_rep[l])
}

bind_rows(sim_output_list, .id = "larval_scenario") %>%
  unite(larval_scenario, simulation_loop, col = "scenario_combn", remove = FALSE) %>%
  filter(! scenario_combn %in% scenario_to_drop) %>%
  select(-scenario_combn) %>%
  bind_rows(bind_rows(sim_output_list_redo, .id = "larval_scenario")) %>%
  write_csv("final_larval_maintext_results_cleaned.csv")


#### make final df ####

# this is our final dataframe to work from:
bind_rows(sim_output_list, .id = "larval_scenario") %>%
  unite(larval_scenario, simulation_loop, col = "scenario_combn", remove = FALSE) %>%
  filter(! scenario_combn %in% scenario_to_drop) %>%
  select(-scenario_combn) %>%
  bind_rows(bind_rows(sim_output_list_redo, .id = "larval_scenario")) -> sim_output_df

  ## get summary stats of interest
    # sim_output_df %>%
    #   group_by(larval_scenario, coexistence_partition, species) %>%
    #   summarise(mean_cs = mean(coexistence_strength),
    #             sd_cs = sd(coexistence_strength),
    #             n_cs = n()) %>%
    #   mutate(se_cs = sd_cs/sqrt(n_cs))



#### plot for main text ####

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols<-c(1,6,3,8,2)
cbPalette <- cbbPalette[cols]

fad<-cbPalette
xlab=c(expression("r"[i]-"r"[r]) ,
       expression(Delta[i]^0),
       expression(Delta[i]^P),
       expression(Delta[i]^E),
       expression(Delta[i]^{E*P}))

sim_output_df %>%
  filter(larval_scenario %in% c("low", "high")) %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  mutate(species = fct_recode(factor(species),
                              `Balanus` = "balanus_glandula",
                              `Chthamalus` = "chthamalus_dalli", 
                              Limpets = "limpets")) %>%
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
  geom_bar(stat = "identity", aes(fill = coexistence_partition)) + 
  geom_errorbar(aes(ymin = mean_cs - se_cs, ymax = mean_cs + se_cs), color = "black", width = 0.2) +
  facet_grid(larval_scenario ~ species, scales = "free") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_discrete(labels=xlab) + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare")


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

sim_output_df %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  mutate(species = fct_recode(factor(species),
                              `Balanus` = "balanus_glandula",
                              `Chthamalus` = "chthamalus_dalli", 
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
  geom_errorbar(aes(ymin = mean_cs - se_cs, ymax = mean_cs + se_cs), color = "black", width = 0.2) +
  facet_grid(species ~ larval_scenario, scales = "free") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_discrete(labels=xlab) + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare")



#### results vs. larval supply ####

# results above are for covarying means & sd. in other words, low mean larval supply & low
# variation in larval supply, vs. high mean larval supply & high variation in larval supply

# run test with l-l, l-h, h-h, and h-l

larval_scenarios_for_supply <- tribble (
  ~scenario, ~species, ~mean_supply_level, ~variance_supply_level,
  1, "B", "low", "low",
  1, "C", "low", "low",
  1, "L", "low", "low",
  1, "P", "low", "low",
  2, "B", "low", "high",
  2, "C", "low", "high",
  2, "L", "low","high",
  2, "P", "low","high",
  3, "B", "high", "high",
  3, "C", "high", "high",
  3, "L", "high", "high",
  3, "P", "high", "high",
  4, "B", "high", "low",
  4, "C", "high", "low",
  4, "L", "high","low",
  4, "P", "high", "low"
)


larval_scenarios_for_supply %>%
  left_join(larval_supply_full) %>%
  select(-mean_supply_level, -variance_supply_level) %>%
  gather(key = variable, value = value, mean_recruit, variance_recruit) %>%
  unite(temp, species, variable) %>%
  spread(temp, value) -> larval_scenarios_input # need it to be named this

sim_output_list_var <- vector(mode = "list", length = 4)
names(sim_output_list_var) <- c("low-low", "low-high", "high-high", "high-low")

for (l in 1:length(sim_output_list_var)) {
  print(c("LOOP = ", l))
  sim_output_list_var[[l]] <- do.larval.supply.simulation(k = l, n_sim = 10)
}

# sim_output_list_var %>%
#   bind_rows(.id = "larval_scenario") %>%
#   write_csv("final_larval_maintext_results_varvar_n3.csv")

sim_output_list_var %>%
  bind_rows(.id = "larval_scenario") %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  mutate(species = fct_recode(factor(species),
                              `Balanus` = "balanus_glandula",
                              `Chthamalus` = "chthamalus_dalli", 
                              Limpets = "limpets")) %>%
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
  geom_errorbar(aes(ymin = mean_cs - se_cs, ymax = mean_cs + se_cs), color = "black", width = 0.2) +
  facet_grid(species ~ larval_scenario, scales = "free") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=fad) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_discrete(labels=xlab) + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare")

# larval supply doesn't change things


#### plot abundances/ free space ####

n_sim <- 3

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
    simulation_loop_output_tmp[[i]] <- do.intertidal.simulation(years_set = 100,
                                                                B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                                C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                                L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                                P.mean = P.mean_input, P.stdev = P.stdev_input)
  }
  
  test_ab_list[[k]] <- simulation_loop_output_tmp
}


lapply(X = test_ab_list$high, FUN = function(x) tail(x$free_space))
lapply(X = test_ab_list$low, FUN = function(x) tail(x$free_space))

lapply(X = test_ab_list$high, FUN = function(x) tail(x$chthamalus_dalli))
lapply(X = test_ab_list$high, FUN = function(x) tail(x$balanus_glandula))
lapply(X = test_ab_list$high, FUN = function(x) tail(x$whelks))
lapply(X = test_ab_list$high, FUN = function(x) tail(x$limpets))
lapply(X = test_ab_list$high, FUN = function(x) tail(x$pisaster_ochraceus))

lapply(X = test_ab_list$high, FUN = function(x) max(tail(x$whelks, n=200)))
lapply(X = test_ab_list$high, FUN = function(x) mean(tail(x$chthamalus_dalli, n = 500)))
lapply(X = test_ab_list$high, FUN = function(x) mean(tail(x$balanus_glandula, n = 500)))

lapply(X = test_ab_list$high, FUN = function(x) x$chthamalus_dalli)
lapply(X = test_ab_list$high, FUN = function(x) x$balanus_glandula)



#### why NA? ####


