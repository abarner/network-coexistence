
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


####













