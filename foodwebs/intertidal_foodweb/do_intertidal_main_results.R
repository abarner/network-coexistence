



#### run overall partition, low, high supplies ####

# have 4 scenarios to run
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
  select(-supply_level, -variance_recruit) %>%
  mutate(supply_level = "med") %>%
  left_join(select(larval_supply, -mean_recruit)) %>%
  select(-supply_level) %>%
  gather(key = variable, value = value, mean_recruit, variance_recruit) %>%
  unite(temp, species, variable) %>%
  spread(temp, value) -> larval_scenarios_input # need it to be named this


sim_output_list <- vector(mode = "list", length = 4)
names(sim_output_list) <- c("high_mean_high_var",
                      "high_mean_low_var",
                      "low_mean_high_var",
                      "low_mean_low_var")

# note - loop may fail but seems to be due to a time out somehow
for (l in 1:4) {
  print(c("LOOP = ", l))
  sim_output_list[[l]] <- do.larval.supply.simulation(k = l, n_sim = 100)
}

sim_output_list %>%
  bind_rows(.id = "larval_scenario") %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  filter(species == "limpets" & coexistence_partition == "r_bar")



