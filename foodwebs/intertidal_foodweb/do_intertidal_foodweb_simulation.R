
#### Source simulation functions ####

source("foodwebs/intertidal_foodweb/intertidal_foodweb_model_modified.R")

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

# 1. run model to equilibrium

fd_results_1 <- do.intertidal.simulation(years_set = 500)

# 2. low density growth rate calculation

# balanus
fd_results_ldr_b_absent <- do.intertidal.simulation(years_set = 500, B_1 = 0)
fd_results_ldr_b_invade <- do.intertidal.simulation(years_set = 500, 
  B_1 = 1,
  C_1 = fd_results_ldr_b_absent$chthamalus_dalli[nrow(fd_results_ldr_b_absent)],
  L_1 = fd_results_ldr_b_absent$limpets[nrow(fd_results_ldr_b_absent)],
  W_1 = fd_results_ldr_b_absent$whelks[nrow(fd_results_ldr_b_absent)],
  P_1 = fd_results_ldr_b_absent$pisaster_ochraceus[nrow(fd_results_ldr_b_absent)],
  total_1 = fd_results_ldr_b_absent$free_space[nrow(fd_results_ldr_b_absent)]
)
# low density growthrates





