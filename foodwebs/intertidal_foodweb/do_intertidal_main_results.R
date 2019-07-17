



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

larval_scenarios_for_full_run %>%
  left_join(larval_supply) %>%
  # select(-supply_level, -variance_recruit) %>%
  # mutate(supply_level = "med") %>%
  # left_join(select(larval_supply, -mean_recruit)) %>%
  select(-supply_level) %>%
  gather(key = variable, value = value, mean_recruit, variance_recruit) %>%
  unite(temp, species, variable) %>%
  spread(temp, value) -> larval_scenarios_input # need it to be named this

sim_output_list <- vector(mode = "list", length = 6)
names(sim_output_list) <- c("high", "low", "balanus_high", "chthamalus_high", 
                            "limpets_high", "pisaster_high")

for (l in 1:length(sim_output_list)) {
  print(c("LOOP = ", l))
  sim_output_list[[l]] <- do.larval.supply.simulation(k = l, n_sim = 100)
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

# this is written to work for only one scenario running again
# would need to rewrite slightly to work for multiple
scenario_num <- which(names(sim_output_list) %in% scenarios_to_run_again$larval_scenario)
scenario_rep <- scenarios_to_run_again$n_loops

# take input from above and overwrite
larval_scenarios_input %>%
  filter(scenario == scenario_num) -> larval_scenarios_input

sim_output_list_redo <- vector(mode = "list", length = nrow(larval_scenarios_input))
names(sim_output_list_redo) <- scenarios_to_run_again$larval_scenario

for (l in 1:length(sim_output_list_redo)) {
  print(c("LOOP = ", l))
  sim_output_list_redo[[l]] <- do.larval.supply.simulation(k = l, n_sim = scenario_rep)
}

bind_rows(bind_rows(sim_output_list_redo, .id = "larval_scenario"),
          bind_rows(sim_output_list, .id = "larval_scenario")) %>%
  unite(larval_scenario, simulation_loop, col = "scenario_combn", remove = FALSE) %>%
  filter(! scenario_combn %in% scenario_to_drop) %>%
  write_csv("final_larval_maintext_results.csv")

# this is our final dataframe to work from:
bind_rows(bind_rows(sim_output_list_redo, .id = "larval_scenario"),
          bind_rows(sim_output_list, .id = "larval_scenario")) %>%
  unite(larval_scenario, simulation_loop, col = "scenario_combn", remove = FALSE) %>%
  filter(! scenario_combn %in% scenario_to_drop) %>%
  select(-scenario_combn) -> sim_output_df

# get summary stats of interest
sim_output_df %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs))

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



