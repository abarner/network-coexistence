
#### Source simulation functions ####

source("foodwebs/intertidal_foodweb/intertidal_foodweb_model_modified.R")
source("foodwebs/intertidal_foodweb/do_intertidal_partition.R")

#### Run main scenario ####

fd_results <- do.intertidal.simulation()

# plot.new()
png(filename = "forde_and_doak_new_dynamics.png", width = 10, height = 7, units = "in", res = 300)
fd_results %>%
  gather(balanus_glandula:free_space, key = "species",
         value = "abundance") %>%
  mutate(neg_value = ifelse(abundance < 0, "yes", "no")) %>%
  ggplot(aes(x = time, y = abundance)) +
  geom_point(size = 2, aes(col = neg_value)) +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  theme(legend.position = "none")
dev.off()

#### Run settlement variation scenarios ####

# because settlement rates were roughly estimated, want to run several scenarios under
# fixed larval supply. 
# which scenarios are of interest? 
# 1. current values
# 2. everyone goes up an order of magnitude
# 3. everyone goes down an order of magnitude
# 4. limpets > barnacles

settlement_scenarios <- tribble(
  ~scenario, ~settlement.B, ~settlement.C, ~settlement.L,
  1, .002 * 30 * 24, .002 * 30 * 24, .0002 * 30 * 24,
  2, .01 * 30 * 24, .01 * 30 * 24, .001 * 30 * 24,
  3, .0002 * 30 * 24, .0002 * 30 * 24, .00002 * 30 * 24,
  4, .0002 * 30 * 24, .0002 * 30 * 24, .002 * 30 * 24
)

fd_results_settle_list <- vector(mode = "list", length = 4)
names(fd_results_settle_list) <- c("main text, forde & doak values",
                                   "all increased one order of magnitude",
                                   "all decreased one order of magnitude",
                                   "forde & doak values, but limpets > barnacles")

for (i in 1:4) {
  fd_results_settle_list[[i]] <- do.intertidal.simulation(
    settlement.B = settlement_scenarios$settlement.B[i],
    settlement.C = settlement_scenarios$settlement.C[i],
    settlement.L = settlement_scenarios$settlement.L[i]
  )
}

png(filename = "forde_and_doak_larvalsupply_variation.png", width = 12, height = 7, units = "in", res = 300)
fd_results_settle_list %>%
  bind_rows(.id = "scenario") %>%
  # scenario 2 causes system to -> unrealistic free space and barnacle densitites
  filter(scenario != "settle_2") %>%
  gather(balanus_glandula:free_space, key = "species",
         value = "abundance") %>%
  ggplot(aes(x = time, y = abundance, color = scenario)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line() +
  facet_wrap(~species, scales = "free")
dev.off()

#### Run recruitment variation scenarios ####

# want to loop across multiple levels of recruitment variation
# for each species, from Table 1 in F & D, use "low" variance option each time

# see: https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
# for how rlnorm function is implemented

larval_supply <- tribble (
  ~species, ~mean_recruit, ~variance_recruit, ~supply_level,
  "B", 90000, sqrt(4.6*10^9), "high",
  "B", 50000, sqrt(1.41*10^10), "med",
  "B", 6000, sqrt(2.025*10^7), "low",
  "C", 70000, sqrt(2.75*10^9), "high",
  "C", 30000, sqrt(5.1*10^8), "med",
  "C", 6000, sqrt(2.025*10^7), "low",
  "L", 3000, sqrt(3.8*10^6), "high",
  "L", 2400, sqrt(2.4*10^6), "med",
  "L", 200, sqrt(169000), "low",
  "P", 6873, sqrt(3.02*10^7), "high",
  "P", 3800, sqrt(9.2*10^6), "med",
  "P", 727, sqrt(3.4*10^5), "low"
)

# which scenarios do we want to run? 
# B high, C high, L high, everyone high
larval_scenarios_table1 <- tribble (
  ~scenario, ~species, ~supply_level,
  1, "B", "high",
  1, "C", "low",
  1, "L", "low",
  1, "P", "low",
  2, "B", "low",
  2, "C", "high",
  2, "L", "low",
  2, "P", "low",
  3, "B", "low",
  3, "C", "low",
  3, "L", "high",
  3, "P", "low",
  4, "B", "high",
  4, "C", "high",
  4, "L", "high",
  4, "P", "high"
)

larval_scenarios_table1 %>%
  left_join(larval_supply) %>%
  select(-supply_level) %>%
  gather(key = variable, value = value, mean_recruit, variance_recruit) %>%
  unite(temp, species, variable) %>%
  spread(temp, value) -> larval_scenarios

fd_results_larval_list <- vector(mode = "list", length = 4)
names(fd_results_larval_list) <- c("B_high", "C_high", "L_high", "all_high")

for (i in 1:4) {
  fd_results_larval_list[[i]] <- do.intertidal.simulation(
    B.mean = larval_scenarios$B_mean_recruit[i],
    B.stdev = larval_scenarios$B_variance_recruit[i],
    C.mean = larval_scenarios$C_mean_recruit[i],
    C.stdev = larval_scenarios$C_variance_recruit[i],
    L.mean = larval_scenarios$L_mean_recruit[i],
    L.stdev = larval_scenarios$L_variance_recruit[i],
    P.mean = larval_scenarios$P_mean_recruit[i],
    P.stdev = larval_scenarios$P_variance_recruit[i]
  )
}

png(filename = "forde_and_doak_recruitment_variation.png", width = 12, height = 7, units = "in", res = 300)
fd_results_larval_list %>%
  bind_rows(.id = "scenario") %>%
  gather(balanus_glandula:free_space, key = "species",
         value = "abundance") %>%
  ggplot(aes(x = time, y = abundance, color = scenario)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line() +
  facet_wrap(~species, scales = "free")
dev.off()



#### Partition coexistence #####

# will ultimately want to run this n=100 times

# 1. run model to equilibrium, get low density growth rates for invader/resident combinations
fd_part_1 <- do.intertidal.rbar()

# 2. Partition coexistence

# Have four varying coexistence mechanisms:
# variation in: balanus, chthamalus, limpet, and pisaster recruitment

# Start just by asking, what happens when competitor recruitment variation removed 
# vs. predator variation removed?

# a. Set all varying parameters to long term average, run steps 1-2 again
# b. Set only "variation in competitor recruitment" = 0, run steps 1-2 again
# c. Set only "variation in predator population size" = 0, run steps 1-2 again

# a. Set all varying parameters to long term average, run steps 1-2 again

# get long term averages:
fd_part_1_var_C <- mean(fd_part_1$results_1$larvae.C)
fd_part_1_var_B <- mean(fd_part_1$results_1$larvae.B)
fd_part_1_var_P <- mean(fd_part_1$results_1$larvae.P)
fd_part_1_var_L <- mean(fd_part_1$results_1$larvae.L)

# run model to equilibrium and get long term and low density growth rates
fd_part_2a <- do.intertidal.all.average(var_P_input = fd_part_1_var_P, 
                                        var_B_input = fd_part_1_var_B, 
                                        var_C_input = fd_part_1_var_C, 
                                        var_L_input = fd_part_1_var_L)

# set variation in competitor recruitment to average (constant)
fd_part_2b <- do.intertidal.all.average(var_P_input = NULL, 
                                        var_B_input = fd_part_1_var_B, 
                                        var_C_input = fd_part_1_var_C, 
                                        var_L_input = fd_part_1_var_L)

# set variation in predator recruitment to average(constant)
fd_part_2c <- do.intertidal.all.average(var_P_input = fd_part_1_var_P, 
                                        var_B_input = NULL, 
                                        var_C_input = NULL, 
                                        var_L_input = NULL)

list(r_bar = fd_part_1$r_bar_result,
     delta_0 = fd_part_2a$r_bar_result,
     delta_p = fd_part_2b$r_bar_result,
     delta_c = fd_part_2c$r_bar_result) %>%
  bind_rows(.id = "id") -> fd_part_df

fd_part_df  %>%
  unite(col = "sp_pair", species_1, species_2, sep = "-", remove = FALSE) %>%
  select(-r_bar_species_2, -species_2) %>%
  rename(species = species_1) %>%
  spread(key = "id", value = r_bar_species_1) %>%
  mutate(delta_cp = r_bar - (delta_0 + delta_c + delta_p)) -> fd_part_df_1
fd_part_df  %>%
  unite(col = "sp_pair", species_1, species_2, sep = "-", remove = FALSE) %>%
  select(-r_bar_species_1, -species_1) %>%
  rename(species = species_2) %>%
  spread(key = "id", value = r_bar_species_2) %>%
  mutate(delta_cp = r_bar - (delta_0 + delta_c + delta_p)) %>%
  bind_rows(fd_part_df_1) -> fd_part_total

fd_part_total %>%
  gather(delta_0:delta_cp, key = "coexistence_partition", value = "coexistence_strength") %>%
  mutate(coexistence_partition = factor(coexistence_partition, levels = c("r_bar", 
                                                                           "delta_0", "delta_c", 
                                                                           "delta_p", "delta_cp"))) %>%
  ggplot(aes(x = coexistence_partition, y = coexistence_strength, color = coexistence_partition)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition)) + 
  facet_grid( ~ species) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(labels=xlab) +
  theme(legend.position = "none")






#### pick up here ####


# consider predator removal completely
# get long term averages:
fd_part_1_var_C <- mean(fd_part_1$results_1$larvae.C)
fd_part_1_var_B <- mean(fd_part_1$results_1$larvae.B)
fd_part_1_var_L <- mean(fd_part_1$results_1$larvae.L)
fd_part_1_P_avg <- mean(fd_part_1$results_1$pisaster_ochraceus)
fd_part_1_W_avg <- mean(fd_part_1$results_1$whelks)

# run model to equilibrium and get long term and low density growth rates
fd_part_3a <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                        var_C_input = fd_part_1_var_C, 
                                        var_L_input = fd_part_1_var_L)

# set variation in competitor recruitment to average (constant)
fd_part_3b <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                        var_C_input = fd_part_1_var_C, 
                                        var_L_input = fd_part_1_var_L)

# set variation in predator ABUNDANCE to average(constant)
fd_part_3c <- do.intertidal.predator.removal(var_B_input = NULL, 
                                        var_C_input = NULL, 
                                        var_L_input = NULL)








