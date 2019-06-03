
## why aren't limpets coexisting in our simulations?


#### show the problem ####

fd_results <- do.intertidal.simulation()

# they look ok when you just run the simulation:
fd_results %>%
  gather(balanus_glandula:free_space, key = "species",
         value = "abundance") %>%
  mutate(neg_value = ifelse(abundance < 0, "yes", "no")) %>%
  ggplot(aes(x = time, y = abundance)) +
  geom_point(size = 2, aes(col = neg_value)) +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  theme(legend.position = "none")


# but when you run the full coexistence partition...
fd_part_1 <- do.intertidal.rbar()

# consider predator removal completely
fd_part_1_var_C <- mean(fd_part_1$results_1$larvae.C)
fd_part_1_var_B <- mean(fd_part_1$results_1$larvae.B)
fd_part_1_var_L <- mean(fd_part_1$results_1$larvae.L)
fd_part_1_var_P <- mean(fd_part_1$results_1$larvae.P)
fd_part_1_P_avg <- mean(fd_part_1$results_1$pisaster_ochraceus)
fd_part_1_W_avg <- mean(fd_part_1$results_1$whelks)

# run model to equilibrium and get long term and low density growth rates
fd_part_3a <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                             var_C_input = fd_part_1_var_C, 
                                             var_L_input = fd_part_1_var_L,
                                             var_P_input = fd_part_1_var_P,
                                             P_avg_input = fd_part_1_P_avg,
                                             W_avg_input = fd_part_1_W_avg)

# set variation in competitor recruitment to average (constant)
fd_part_3b <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                             var_C_input = fd_part_1_var_C, 
                                             var_L_input = fd_part_1_var_L,
                                             var_P_input = NULL,
                                             P_avg_input = NULL,
                                             W_avg_input = NULL)

# set variation in predator ABUNDANCE to average(constant)
fd_part_3c <- do.intertidal.predator.removal(var_B_input = NULL, 
                                             var_C_input = NULL, 
                                             var_L_input = NULL,
                                             var_P_input = fd_part_1_var_P,
                                             P_avg_input = fd_part_1_P_avg,
                                             W_avg_input = fd_part_1_W_avg)

list(r_bar = fd_part_1$r_bar_result,
     delta_0 = fd_part_3a$r_bar_result,
     delta_p = fd_part_3b$r_bar_result,
     delta_c = fd_part_3c$r_bar_result) %>%
  bind_rows(.id = "id") -> fd_part_3_df

fd_part_3_df  %>%
  spread(key = "id", value = r_bar) %>%
  mutate(delta_cp = r_bar - (delta_0 + delta_c + delta_p)) -> fd_part_3_total_unmodified

fd_part_3_total_unmodified %>%
  gather(delta_0:delta_cp, key = "coexistence_partition", value = "coexistence_strength") %>%
  mutate(coexistence_partition = factor(coexistence_partition, levels = c("r_bar", 
                                                                          "delta_0", "delta_c", 
                                                                          "delta_p", "delta_cp"))) %>%
  ggplot(aes(x = coexistence_partition, y = coexistence_strength, color = coexistence_partition)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition)) + 
  facet_grid( ~ species) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")


## limpets also do not coexist even if their larval supply rates are higher,
## though the coexistence strength weakens from negative to more neutral

#### barnacle survival rate and limpet coexistence ####

## some looking through old notes seems to indicate that increasing the survival
## of balanus glandula adults negatively impacted limpets

## so i decreased the survival rates of BG from 0.8 to 0.75 for both
## adults and recruits

## run again:

fd_part_1 <- do.intertidal.rbar()

# consider predator removal completely
fd_part_1_var_C <- mean(fd_part_1$results_1$larvae.C)
fd_part_1_var_B <- mean(fd_part_1$results_1$larvae.B)
fd_part_1_var_L <- mean(fd_part_1$results_1$larvae.L)
fd_part_1_var_P <- mean(fd_part_1$results_1$larvae.P)
fd_part_1_P_avg <- mean(fd_part_1$results_1$pisaster_ochraceus)
fd_part_1_W_avg <- mean(fd_part_1$results_1$whelks)

# run model to equilibrium and get long term and low density growth rates
fd_part_3a <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                             var_C_input = fd_part_1_var_C, 
                                             var_L_input = fd_part_1_var_L,
                                             var_P_input = fd_part_1_var_P,
                                             P_avg_input = fd_part_1_P_avg,
                                             W_avg_input = fd_part_1_W_avg)

# set variation in competitor recruitment to average (constant)
fd_part_3b <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                             var_C_input = fd_part_1_var_C, 
                                             var_L_input = fd_part_1_var_L,
                                             var_P_input = NULL,
                                             P_avg_input = NULL,
                                             W_avg_input = NULL)

# set variation in predator ABUNDANCE to average(constant)
fd_part_3c <- do.intertidal.predator.removal(var_B_input = NULL, 
                                             var_C_input = NULL, 
                                             var_L_input = NULL,
                                             var_P_input = fd_part_1_var_P,
                                             P_avg_input = fd_part_1_P_avg,
                                             W_avg_input = fd_part_1_W_avg)

list(r_bar = fd_part_1$r_bar_result,
     delta_0 = fd_part_3a$r_bar_result,
     delta_p = fd_part_3b$r_bar_result,
     delta_c = fd_part_3c$r_bar_result) %>%
  bind_rows(.id = "id") -> fd_part_3_df

fd_part_3_df  %>%
  spread(key = "id", value = r_bar) %>%
  mutate(delta_cp = r_bar - (delta_0 + delta_c + delta_p)) -> fd_part_3_total_lowerbarn

fd_part_3_total_lowerbarn %>%
  gather(delta_0:delta_cp, key = "coexistence_partition", value = "coexistence_strength") %>%
  mutate(coexistence_partition = factor(coexistence_partition, levels = c("r_bar", 
                                                                          "delta_0", "delta_c", 
                                                                          "delta_p", "delta_cp"))) %>%
  ggplot(aes(x = coexistence_partition, y = coexistence_strength, color = coexistence_partition)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition)) + 
  facet_grid( ~ species) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")

## still not quite right

## reduce further again
## this time make BG/CD +/- 0.1 of 0.7 (the original survival rates)
fd_part_1 <- do.intertidal.rbar()

# consider predator removal completely
fd_part_1_var_C <- mean(fd_part_1$results_1$larvae.C)
fd_part_1_var_B <- mean(fd_part_1$results_1$larvae.B)
fd_part_1_var_L <- mean(fd_part_1$results_1$larvae.L)
fd_part_1_var_P <- mean(fd_part_1$results_1$larvae.P)
fd_part_1_P_avg <- mean(fd_part_1$results_1$pisaster_ochraceus)
fd_part_1_W_avg <- mean(fd_part_1$results_1$whelks)

# run model to equilibrium and get long term and low density growth rates
fd_part_3a <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                             var_C_input = fd_part_1_var_C, 
                                             var_L_input = fd_part_1_var_L,
                                             var_P_input = fd_part_1_var_P,
                                             P_avg_input = fd_part_1_P_avg,
                                             W_avg_input = fd_part_1_W_avg)

# set variation in competitor recruitment to average (constant)
fd_part_3b <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                             var_C_input = fd_part_1_var_C, 
                                             var_L_input = fd_part_1_var_L,
                                             var_P_input = NULL,
                                             P_avg_input = NULL,
                                             W_avg_input = NULL)

# set variation in predator ABUNDANCE to average(constant)
fd_part_3c <- do.intertidal.predator.removal(var_B_input = NULL, 
                                             var_C_input = NULL, 
                                             var_L_input = NULL,
                                             var_P_input = fd_part_1_var_P,
                                             P_avg_input = fd_part_1_P_avg,
                                             W_avg_input = fd_part_1_W_avg)

list(r_bar = fd_part_1$r_bar_result,
     delta_0 = fd_part_3a$r_bar_result,
     delta_p = fd_part_3b$r_bar_result,
     delta_c = fd_part_3c$r_bar_result) %>%
  bind_rows(.id = "id") -> fd_part_3_df

fd_part_3_df  %>%
  spread(key = "id", value = r_bar) %>%
  mutate(delta_cp = r_bar - (delta_0 + delta_c + delta_p)) -> fd_part_3_total_lowerbothbarn

fd_part_3_total_lowerbothbarn %>%
  gather(delta_0:delta_cp, key = "coexistence_partition", value = "coexistence_strength") %>%
  mutate(coexistence_partition = factor(coexistence_partition, levels = c("r_bar", 
                                                                          "delta_0", "delta_c", 
                                                                          "delta_p", "delta_cp"))) %>%
  ggplot(aes(x = coexistence_partition, y = coexistence_strength, color = coexistence_partition)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition)) + 
  facet_grid( ~ species) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")

## still no limpet coexistence, but barnacle coexistence dynamics changed a lot


## just run BG slightly higher than CD
fd_part_1 <- do.intertidal.rbar()

# consider predator removal completely
fd_part_1_var_C <- mean(fd_part_1$results_1$larvae.C)
fd_part_1_var_B <- mean(fd_part_1$results_1$larvae.B)
fd_part_1_var_L <- mean(fd_part_1$results_1$larvae.L)
fd_part_1_var_P <- mean(fd_part_1$results_1$larvae.P)
fd_part_1_P_avg <- mean(fd_part_1$results_1$pisaster_ochraceus)
fd_part_1_W_avg <- mean(fd_part_1$results_1$whelks)

# run model to equilibrium and get long term and low density growth rates
fd_part_3a <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                             var_C_input = fd_part_1_var_C, 
                                             var_L_input = fd_part_1_var_L,
                                             var_P_input = fd_part_1_var_P,
                                             P_avg_input = fd_part_1_P_avg,
                                             W_avg_input = fd_part_1_W_avg)

# set variation in competitor recruitment to average (constant)
fd_part_3b <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                             var_C_input = fd_part_1_var_C, 
                                             var_L_input = fd_part_1_var_L,
                                             var_P_input = NULL,
                                             P_avg_input = NULL,
                                             W_avg_input = NULL)

# set variation in predator ABUNDANCE to average(constant)
fd_part_3c <- do.intertidal.predator.removal(var_B_input = NULL, 
                                             var_C_input = NULL, 
                                             var_L_input = NULL,
                                             var_P_input = fd_part_1_var_P,
                                             P_avg_input = fd_part_1_P_avg,
                                             W_avg_input = fd_part_1_W_avg)

list(r_bar = fd_part_1$r_bar_result,
     delta_0 = fd_part_3a$r_bar_result,
     delta_p = fd_part_3b$r_bar_result,
     delta_c = fd_part_3c$r_bar_result) %>%
  bind_rows(.id = "id") -> fd_part_3_df

fd_part_3_df  %>%
  spread(key = "id", value = r_bar) %>%
  mutate(delta_cp = r_bar - (delta_0 + delta_c + delta_p)) -> fd_part_3_total_bghigher

fd_part_3_total_bghigher %>%
  gather(delta_0:delta_cp, key = "coexistence_partition", value = "coexistence_strength") %>%
  mutate(coexistence_partition = factor(coexistence_partition, levels = c("r_bar", 
                                                                          "delta_0", "delta_c", 
                                                                          "delta_p", "delta_cp"))) %>%
  ggplot(aes(x = coexistence_partition, y = coexistence_strength, color = coexistence_partition)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition)) + 
  facet_grid( ~ species) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")

## still no fix


## double check that it worked at 0.7
fd_part_1 <- do.intertidal.rbar()

# consider predator removal completely
fd_part_1_var_C <- mean(fd_part_1$results_1$larvae.C)
fd_part_1_var_B <- mean(fd_part_1$results_1$larvae.B)
fd_part_1_var_L <- mean(fd_part_1$results_1$larvae.L)
fd_part_1_var_P <- mean(fd_part_1$results_1$larvae.P)
fd_part_1_P_avg <- mean(fd_part_1$results_1$pisaster_ochraceus)
fd_part_1_W_avg <- mean(fd_part_1$results_1$whelks)

# run model to equilibrium and get long term and low density growth rates
fd_part_3a <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                             var_C_input = fd_part_1_var_C, 
                                             var_L_input = fd_part_1_var_L,
                                             var_P_input = fd_part_1_var_P,
                                             P_avg_input = fd_part_1_P_avg,
                                             W_avg_input = fd_part_1_W_avg)

# set variation in competitor recruitment to average (constant)
fd_part_3b <- do.intertidal.predator.removal(var_B_input = fd_part_1_var_B, 
                                             var_C_input = fd_part_1_var_C, 
                                             var_L_input = fd_part_1_var_L,
                                             var_P_input = NULL,
                                             P_avg_input = NULL,
                                             W_avg_input = NULL)

# set variation in predator ABUNDANCE to average(constant)
fd_part_3c <- do.intertidal.predator.removal(var_B_input = NULL, 
                                             var_C_input = NULL, 
                                             var_L_input = NULL,
                                             var_P_input = fd_part_1_var_P,
                                             P_avg_input = fd_part_1_P_avg,
                                             W_avg_input = fd_part_1_W_avg)

list(r_bar = fd_part_1$r_bar_result,
     delta_0 = fd_part_3a$r_bar_result,
     delta_p = fd_part_3b$r_bar_result,
     delta_c = fd_part_3c$r_bar_result) %>%
  bind_rows(.id = "id") -> fd_part_3_df

fd_part_3_df  %>%
  spread(key = "id", value = r_bar) %>%
  mutate(delta_cp = r_bar - (delta_0 + delta_c + delta_p)) -> fd_part_3_total_barnequal

fd_part_3_total_barnequal %>%
  gather(delta_0:delta_cp, key = "coexistence_partition", value = "coexistence_strength") %>%
  mutate(coexistence_partition = factor(coexistence_partition, levels = c("r_bar", 
                                                                          "delta_0", "delta_c", 
                                                                          "delta_p", "delta_cp"))) %>%
  ggplot(aes(x = coexistence_partition, y = coexistence_strength, color = coexistence_partition)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition)) + 
  facet_grid( ~ species) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")

## ok no it's not barnacle survival rates that are causing the problem 


#### Test the effect of lower variance on limpet survival ####

# from FD paper: when variance in larval supply decreases, then limpet density increased

# try scenario: all high, but changed variance

larval_supply_full <- tribble (
  ~species, ~mean_recruit, ~variance_recruit, ~mean_supply_level, ~variance_supply_level,
  "B", 90000, sqrt(4.6*10^9), "high", "low",
  "B", 90000, sqrt(3.24*10^10), "high", "med",
  "B", 90000, sqrt(4.6*10^10), "high", "high",
  
  "B", 50000, sqrt(1.41*10^10), "med", "low",
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
  
  "L", 200, sqrt(169000), "low", "low",
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

tribble (
  ~scenario, ~species, ~variance_supply_level,
  1, "B", "high",
  1, "C", "high",
  1, "L", "high",
  1, "P", "high",
  2, "B", "low",
  2, "C", "low",
  2, "L", "low",
  2, "P", "low",
  3, "B", "med",
  3, "C", "med",
  3, "L", "med",
  3, "P", "med"
) -> larval_scenarios_variance_run

larval_scenarios_variance_run %>%
  mutate(mean_supply_level = "high") %>%
  left_join(larval_supply_full) %>%
  select(-mean_supply_level, -variance_supply_level) %>%
  gather(key = variable, value = value, mean_recruit, variance_recruit) %>%
  unite(temp, species, variable) %>%
  spread(temp, value) -> larval_scenarios_input

test_var_list <- vector(mode = "list", length = 3)
names(test_var_list) <- c("high", "low", "med")

# note - loop may fail but seems to be due to a time out somehow
for (l in 1:3) {
  print(c("LOOP = ", l))
  test_var_list[[l]] <- do.larval.supply.simulation(k = l, n_sim = 3)
}

png("high_mean_supply_variation_in_variance.png", width = 6, height = 6, units = "in", res = 300)
test_var_list %>%
  bind_rows(.id = "larval_scenario") %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  ggplot(aes(x = coexistence_partition, y = mean_cs)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition)) + 
  geom_errorbar(aes(ymin = mean_cs - se_cs, ymax = mean_cs + se_cs), color = "black", width = 0.2) +
  facet_grid(species ~ larval_scenario, scales = "free") +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare")
dev.off()

## does NOT improve coexistence for limpets
## higher variance seems to be more of a problem

## what about with low overall supply?
larval_scenarios_variance_run %>%
  mutate(mean_supply_level = "low") %>%
  left_join(larval_supply_full) %>%
  select(-mean_supply_level, -variance_supply_level) %>%
  gather(key = variable, value = value, mean_recruit, variance_recruit) %>%
  unite(temp, species, variable) %>%
  spread(temp, value) -> larval_scenarios_input

test_var_list_2 <- vector(mode = "list", length = 3)
names(test_var_list_2) <- c("high", "low", "med")

for (l in 1:3) {
  print(c("LOOP = ", l))
  test_var_list_2[[l]] <- do.larval.supply.simulation(k = l, n_sim = 3)
}

test_var_list_2 %>%
  bind_rows(.id = "larval_scenario") %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs)) %>%
  ungroup() %>%
  ggplot(aes(x = coexistence_partition, y = mean_cs)) +
  geom_bar(stat = "identity", aes(fill = coexistence_partition)) + 
  geom_errorbar(aes(ymin = mean_cs - se_cs, ymax = mean_cs + se_cs), color = "black", width = 0.2) +
  facet_grid(species ~ larval_scenario, scales = "free") +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none") + 
  xlab("Mechanistic partitioning") +
  ylab("Growth rate when rare")

### still no coexistence



#### set limpet supply rates to REALLY high? ####














