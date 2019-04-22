
### function to run the first step of coexistence partitioning

do.intertidal.rbar <- function() {
  
  # 1. run to eqilibrium
  
  fd_results_1 <- do.intertidal.simulation(years_set = 500)
  
  # 2. low density growth rate calculation
  
  # balanus
  # run without balanus adults or recruits -> equilibrium
  fd_results_ldr_b_absent <- do.intertidal.simulation(years_set = 500, B_1 = 0)
  # set adult balanus population size very small (~1/1000 of usual size) and run using 
  # densities from last time step
  fd_results_ldr_b_invade <- 
    do.intertidal.simulation(years_set = 500, 
                             B_1 = 1,
                             C_1 = fd_results_ldr_b_absent$chthamalus_dalli[nrow(fd_results_ldr_b_absent)],
                             L_1 = fd_results_ldr_b_absent$limpets[nrow(fd_results_ldr_b_absent)],
                             W_1 = fd_results_ldr_b_absent$whelks[nrow(fd_results_ldr_b_absent)],
                             P_1 = fd_results_ldr_b_absent$pisaster_ochraceus[nrow(fd_results_ldr_b_absent)],
                             total_1 = fd_results_ldr_b_absent$free_space[nrow(fd_results_ldr_b_absent)]
    )
  # calculate low density growth rates
  gr_b_invade <- do.growth.rates(results = fd_results_ldr_b_invade, col_nums = c(2:6))
  
  
  # chthamalus
  fd_results_ldr_c_absent <- do.intertidal.simulation(years_set = 500, C_1 = 0)
  fd_results_ldr_c_invade <- 
    do.intertidal.simulation(years_set = 500, 
                             B_1 = fd_results_ldr_c_absent$balanus_glandula[nrow(fd_results_ldr_c_absent)],
                             C_1 = 1,
                             L_1 = fd_results_ldr_c_absent$limpets[nrow(fd_results_ldr_c_absent)],
                             W_1 = fd_results_ldr_c_absent$whelks[nrow(fd_results_ldr_c_absent)],
                             P_1 = fd_results_ldr_c_absent$pisaster_ochraceus[nrow(fd_results_ldr_c_absent)],
                             total_1 = fd_results_ldr_c_absent$free_space[nrow(fd_results_ldr_c_absent)]
    )
  # low density growthrates
  gr_c_invade <- do.growth.rates(results = fd_results_ldr_c_invade, col_nums = c(2:6))
  
  
  # limpets
  fd_results_ldr_l_absent <- do.intertidal.simulation(years_set = 500, L_1 = 0)
  fd_results_ldr_l_invade <- 
    do.intertidal.simulation(years_set = 500, 
                             B_1 = fd_results_ldr_l_absent$balanus_glandula[nrow(fd_results_ldr_l_absent)],
                             C_1 = fd_results_ldr_l_absent$chthamalus_dalli[nrow(fd_results_ldr_l_absent)],
                             L_1 = 1,
                             W_1 = fd_results_ldr_l_absent$whelks[nrow(fd_results_ldr_l_absent)],
                             P_1 = fd_results_ldr_l_absent$pisaster_ochraceus[nrow(fd_results_ldr_l_absent)],
                             total_1 = fd_results_ldr_l_absent$free_space[nrow(fd_results_ldr_l_absent)]
    )
  # low density growthrates
  gr_l_invade <- do.growth.rates(results = fd_results_ldr_l_invade, col_nums = c(2:6))
  
  
  # whelks
  fd_results_ldr_w_absent <- do.intertidal.simulation(years_set = 500, W_1 = 0)
  fd_results_ldr_w_invade <- 
    do.intertidal.simulation(years_set = 500, 
                             B_1 = fd_results_ldr_w_absent$balanus_glandula[nrow(fd_results_ldr_w_absent)],
                             C_1 = fd_results_ldr_w_absent$chthamalus_dalli[nrow(fd_results_ldr_w_absent)],
                             L_1 = fd_results_ldr_w_absent$limpets[nrow(fd_results_ldr_w_absent)],
                             W_1 = 1,
                             P_1 = fd_results_ldr_w_absent$pisaster_ochraceus[nrow(fd_results_ldr_w_absent)],
                             total_1 = fd_results_ldr_w_absent$free_space[nrow(fd_results_ldr_w_absent)]
    )
  # low density growthrates
  gr_w_invade <- do.growth.rates(results = fd_results_ldr_w_invade, col_nums = c(2:6))
  
  
  # sea stars (pisaster)
  fd_results_ldr_p_absent <- do.intertidal.simulation(years_set = 500, P_1 = 0)
  fd_results_ldr_p_invade <- 
    do.intertidal.simulation(years_set = 500, 
                             B_1 = fd_results_ldr_p_absent$balanus_glandula[nrow(fd_results_ldr_p_absent)],
                             C_1 = fd_results_ldr_p_absent$chthamalus_dalli[nrow(fd_results_ldr_p_absent)],
                             L_1 = fd_results_ldr_p_absent$limpets[nrow(fd_results_ldr_p_absent)],
                             W_1 = fd_results_ldr_p_absent$whelks[nrow(fd_results_ldr_p_absent)],
                             P_1 = 1,
                             total_1 = fd_results_ldr_p_absent$free_space[nrow(fd_results_ldr_p_absent)]
    )
  # low density growthrates
  gr_p_invade <- do.growth.rates(results = fd_results_ldr_p_invade, col_nums = c(2:6))
  
  ## pairwise coexistence strengths
  
  # balanus & chthamalus
  r_bar_bc_b <- mean(gr_b_invade[, "balanus_glandula"], na.rm=TRUE) - mean(gr_b_invade[, "chthamalus_dalli"], na.rm=TRUE)
  r_bar_bc_c <- mean(gr_c_invade[, "chthamalus_dalli"], na.rm=TRUE) - mean(gr_c_invade[, "balanus_glandula"], na.rm=TRUE)
  #r_bar_bc_b; r_bar_bc_c
  
  # balanus & limpets
  r_bar_bl_b <- mean(gr_b_invade[, "balanus_glandula"], na.rm=TRUE) - mean(gr_b_invade[, "limpets"], na.rm=TRUE)
  r_bar_bl_l <- mean(gr_l_invade[, "limpets"], na.rm=TRUE) - mean(gr_l_invade[, "balanus_glandula"], na.rm=TRUE)
  #r_bar_bl_b; r_bar_bl_l
  
  # chthamalus & limpets
  r_bar_cl_c <- mean(gr_c_invade[, "chthamalus_dalli"], na.rm=TRUE) - mean(gr_c_invade[, "limpets"], na.rm=TRUE)
  r_bar_cl_l <- mean(gr_l_invade[, "limpets"], na.rm=TRUE) - mean(gr_l_invade[, "chthamalus_dalli"], na.rm=TRUE)
  #r_bar_cl_c; r_bar_cl_l
  return(list(r_bar_result = tibble(r_bar_species_1 = c(r_bar_bc_b, r_bar_bl_b, r_bar_cl_c),
                                    r_bar_species_2 = c(r_bar_bc_c, r_bar_bl_l, r_bar_cl_l),
                                    species_1 = c("balanus_glandula", "balanus_glandula", "chthamalus_dalli"),
                                    species_2 = c("chthamalus_dalli", "limpets", "limpets")),
              var_average = tibble(var_b_average = mean(fd_results_1$larvae.B),
                                   var_c_average = mean(fd_results_1$larvae.C),
                                   var_l_average = mean(fd_results_1$larvae.L),
                                   var_p_average = mean(fd_results_1$larvae.P)),
              # also return results in case want to look at each
              # results from overall run
              results_1 = fd_results_1,
              # results from absence/invasion of each species
              results_2_b_invade = bind_rows(fd_results_ldr_b_absent, fd_results_ldr_b_invade),
              results_2_c_invade = bind_rows(fd_results_ldr_c_absent, fd_results_ldr_c_invade),
              results_2_l_invade = bind_rows(fd_results_ldr_l_absent, fd_results_ldr_l_invade),
              results_2_p_invade = bind_rows(fd_results_ldr_p_absent, fd_results_ldr_p_invade)
              
              )
         )
}



### function to run the first step of coexistence partitioning

do.intertidal.all.average <- function(var_P, var_B, var_C, var_L) {
  
  # low density growth rate calculation
  
  # balanus
  # run without balanus adults or recruits -> equilibrium
  fd_results_ldr_b_absent <- do.intertidal.simulation(years_set = 500, B_1 = 0,
                                                      var_P = var_P, var_B = var_B, var_C = var_C, var_L = var_L)
  # set adult balanus population size very small (~1/1000 of usual size) and run using 
  # densities from last time step
  fd_results_ldr_b_invade <- 
    do.intertidal.simulation(years_set = 500, 
                             B_1 = 1,
                             C_1 = fd_results_ldr_b_absent$chthamalus_dalli[nrow(fd_results_ldr_b_absent)],
                             L_1 = fd_results_ldr_b_absent$limpets[nrow(fd_results_ldr_b_absent)],
                             W_1 = fd_results_ldr_b_absent$whelks[nrow(fd_results_ldr_b_absent)],
                             P_1 = fd_results_ldr_b_absent$pisaster_ochraceus[nrow(fd_results_ldr_b_absent)],
                             total_1 = fd_results_ldr_b_absent$free_space[nrow(fd_results_ldr_b_absent)],
                             var_P = var_P, var_B = var_B, var_C = var_C, var_L = var_L
    )
  # calculate low density growth rates
  gr_b_invade <- do.growth.rates(results = fd_results_ldr_b_invade, col_nums = c(2:6))
  
  
  # chthamalus
  fd_results_ldr_c_absent <- do.intertidal.simulation(years_set = 500, C_1 = 0,
                                                      var_P = var_P, var_B = var_B, var_C = var_C, var_L = var_L)
  fd_results_ldr_c_invade <- 
    do.intertidal.simulation(years_set = 500, 
                             B_1 = fd_results_ldr_c_absent$balanus_glandula[nrow(fd_results_ldr_c_absent)],
                             C_1 = 1,
                             L_1 = fd_results_ldr_c_absent$limpets[nrow(fd_results_ldr_c_absent)],
                             W_1 = fd_results_ldr_c_absent$whelks[nrow(fd_results_ldr_c_absent)],
                             P_1 = fd_results_ldr_c_absent$pisaster_ochraceus[nrow(fd_results_ldr_c_absent)],
                             total_1 = fd_results_ldr_c_absent$free_space[nrow(fd_results_ldr_c_absent)],
                             var_P = var_P, var_B = var_B, var_C = var_C, var_L = var_L
    )
  # low density growthrates
  gr_c_invade <- do.growth.rates(results = fd_results_ldr_c_invade, col_nums = c(2:6))
  
  
  # limpets
  fd_results_ldr_l_absent <- do.intertidal.simulation(years_set = 500, L_1 = 0,
                                                      var_P = var_P, var_B = var_B, var_C = var_C, var_L = var_L)
  fd_results_ldr_l_invade <- 
    do.intertidal.simulation(years_set = 500, 
                             B_1 = fd_results_ldr_l_absent$balanus_glandula[nrow(fd_results_ldr_l_absent)],
                             C_1 = fd_results_ldr_l_absent$chthamalus_dalli[nrow(fd_results_ldr_l_absent)],
                             L_1 = 1,
                             W_1 = fd_results_ldr_l_absent$whelks[nrow(fd_results_ldr_l_absent)],
                             P_1 = fd_results_ldr_l_absent$pisaster_ochraceus[nrow(fd_results_ldr_l_absent)],
                             total_1 = fd_results_ldr_l_absent$free_space[nrow(fd_results_ldr_l_absent)],
                             var_P = var_P, var_B = var_B, var_C = var_C, var_L = var_L
    )
  # low density growthrates
  gr_l_invade <- do.growth.rates(results = fd_results_ldr_l_invade, col_nums = c(2:6))
  
  
  # whelks
  fd_results_ldr_w_absent <- do.intertidal.simulation(years_set = 500, W_1 = 0,
                                                      var_P = var_P, var_B = var_B, var_C = var_C, var_L = var_L)
  fd_results_ldr_w_invade <- 
    do.intertidal.simulation(years_set = 500, 
                             B_1 = fd_results_ldr_w_absent$balanus_glandula[nrow(fd_results_ldr_w_absent)],
                             C_1 = fd_results_ldr_w_absent$chthamalus_dalli[nrow(fd_results_ldr_w_absent)],
                             L_1 = fd_results_ldr_w_absent$limpets[nrow(fd_results_ldr_w_absent)],
                             W_1 = 1,
                             P_1 = fd_results_ldr_w_absent$pisaster_ochraceus[nrow(fd_results_ldr_w_absent)],
                             total_1 = fd_results_ldr_w_absent$free_space[nrow(fd_results_ldr_w_absent)],
                             var_P = var_P, var_B = var_B, var_C = var_C, var_L = var_L
    )
  # low density growthrates
  gr_w_invade <- do.growth.rates(results = fd_results_ldr_w_invade, col_nums = c(2:6))
  
  
  # sea stars (pisaster)
  fd_results_ldr_p_absent <- do.intertidal.simulation(years_set = 500, P_1 = 0,
                                                      var_P = var_P, var_B = var_B, var_C = var_C, var_L = var_L)
  fd_results_ldr_p_invade <- 
    do.intertidal.simulation(years_set = 500, 
                             B_1 = fd_results_ldr_p_absent$balanus_glandula[nrow(fd_results_ldr_p_absent)],
                             C_1 = fd_results_ldr_p_absent$chthamalus_dalli[nrow(fd_results_ldr_p_absent)],
                             L_1 = fd_results_ldr_p_absent$limpets[nrow(fd_results_ldr_p_absent)],
                             W_1 = fd_results_ldr_p_absent$whelks[nrow(fd_results_ldr_p_absent)],
                             P_1 = 1,
                             total_1 = fd_results_ldr_p_absent$free_space[nrow(fd_results_ldr_p_absent)],
                             var_P = var_P, var_B = var_B, var_C = var_C, var_L = var_L
    )
  # low density growthrates
  gr_p_invade <- do.growth.rates(results = fd_results_ldr_p_invade, col_nums = c(2:6))
  
  ## pairwise coexistence strengths
  
  # balanus & chthamalus
  r_bar_bc_b <- mean(gr_b_invade[, "balanus_glandula"], na.rm=TRUE) - mean(gr_b_invade[, "chthamalus_dalli"], na.rm=TRUE)
  r_bar_bc_c <- mean(gr_c_invade[, "chthamalus_dalli"], na.rm=TRUE) - mean(gr_c_invade[, "balanus_glandula"], na.rm=TRUE)
  #r_bar_bc_b; r_bar_bc_c
  
  # balanus & limpets
  r_bar_bl_b <- mean(gr_b_invade[, "balanus_glandula"], na.rm=TRUE) - mean(gr_b_invade[, "limpets"], na.rm=TRUE)
  r_bar_bl_l <- mean(gr_l_invade[, "limpets"], na.rm=TRUE) - mean(gr_l_invade[, "balanus_glandula"], na.rm=TRUE)
  #r_bar_bl_b; r_bar_bl_l
  
  # chthamalus & limpets
  r_bar_cl_c <- mean(gr_c_invade[, "chthamalus_dalli"], na.rm=TRUE) - mean(gr_c_invade[, "limpets"], na.rm=TRUE)
  r_bar_cl_l <- mean(gr_l_invade[, "limpets"], na.rm=TRUE) - mean(gr_l_invade[, "chthamalus_dalli"], na.rm=TRUE)
  #r_bar_cl_c; r_bar_cl_l
  return(list(r_bar_result = tibble(r_bar_species_1 = c(r_bar_bc_b, r_bar_bl_b, r_bar_cl_c),
                                    r_bar_species_2 = c(r_bar_bc_c, r_bar_bl_l, r_bar_cl_l),
                                    species_1 = c("balanus_glandula", "balanus_glandula", "chthamalus_dalli"),
                                    species_2 = c("chthamalus_dalli", "limpets", "limpets")),
              var_average = tibble(var_b_average = mean(fd_results_1$larvae.B),
                                   var_c_average = mean(fd_results_1$larvae.C),
                                   var_l_average = mean(fd_results_1$larvae.L),
                                   var_p_average = mean(fd_results_1$larvae.P)),
              # also return results in case want to look at each
              # results from overall run
              results_1 = fd_results_1,
              # results from absence/invasion of each species
              results_2_b_invade = bind_rows(fd_results_ldr_b_absent, fd_results_ldr_b_invade),
              results_2_c_invade = bind_rows(fd_results_ldr_c_absent, fd_results_ldr_c_invade),
              results_2_l_invade = bind_rows(fd_results_ldr_l_absent, fd_results_ldr_l_invade),
              results_2_p_invade = bind_rows(fd_results_ldr_p_absent, fd_results_ldr_p_invade)
              
  )
  )
}


### 









