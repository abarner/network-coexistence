
#### Source simulation functions ####

source("foodwebs/intertidal_foodweb/intertidal_foodweb_model_modified.R")

#### Run main scenario ####

fd_results <- do.intertidal.simulation()

# plot.new()
# png(filename = "forde_and_doak_new_dynamics.png", width = 10, height = 7, units = "in", res = 300)
fd_results %>%
  gather(balanus_glandula:free_space, key = "species",
         value = "abundance") %>%
  mutate(neg_value = ifelse(abundance < 0, "yes", "no")) %>%
  ggplot(aes(x = time, y = abundance)) +
  geom_point(size = 2, aes(col = neg_value)) +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  theme(legend.position = "none")
#dev.off()



#### Run multiple scenarios ####

## Settlement variation

# because settlement rates were roughly estimated, want to run several scenarios under
# fixed larval supply. 
# which scenarios are of interest? 
# 1. current values
# 2. everyone goes up an order of magnitude
# 3. everyone goes down an order of magnitude
# 4. limpets > barnacles

settlement.B <- .002 * 30 * 24 # to convert into per month
settlement.C <- .002 * 30 * 24
settlement.L <- .0002 * 30 * 24

settlement_scenarios <- tribble(
  ~scenario, ~settlement.B, ~settlement.C, ~settlement.L,
  1, .002 * 30 * 24, .002 * 30 * 24, .0002 * 30 * 24,
  2, .01 * 30 * 24, .01 * 30 * 24, .001 * 30 * 24,
  3, .0002 * 30 * 24, .0002 * 30 * 24, .00002 * 30 * 24,
  4, .0002 * 30 * 24, .0002 * 30 * 24, .002 * 30 * 24
)


## Recruitment variation

# want to loop across multiple levels of recruitment variation
# for each species, from Table 1 in F & D, use "low" variance option each time
location_function <- function(x_mean, x_stdev) log(x_mean^2 / sqrt(x_stdev^2 + x_mean^2))
shape_function <- function(x_mean, x_stdev) sqrt(log(1 + (x_stdev ^2 / x_mean^2)))

larval_supply <- tribble (
  ~species, ~mean_recruit, ~variance_recruit, ~supply_level,
  "B", 90000, 4.6*10^9, "high",
  "B", 50000, 1.41*10^10, "med",
  "B", 6000, 2.025*10^7, "low",
  "C", 70000, 2.75*10^9, "high",
  "C", 30000, 5.1*10^8, "med",
  "C", 6000, 2.025*10^7, "low",
  "L", 3000, 3.8*10^6, "high",
  "L", 2400, 2.4*10^6, "med",
  "L", 200, 169000, "low",
  "P", 6873, 3.02*10^7, "high",
  "P", 3800, 9.2*10^6, "med",
  "P", 727, 3.4*10^5, "low"
)

# which scenarios do we want to run? 
# B high, C high
larval_scenarios_table1 <- tribble (
  ~scenario, ~species, ~supply_level,
  1, "B", "high",
  1, "C", "low",
  1, "L", "low",
  1, "P", "low",
  2, "B", "low",
  2, "C", "high",
  2, "L", "low",
  2, "P", "low"
)

larval_scenarios_table1 %>%
  left_join(larval_supply) %>%
  mutate(location = location_function(mean_recruit, variance_recruit),
         shape = shape_function(mean_recruit, variance_recruit)) %>%
  select(-mean_recruit, -variance_recruit, - supply_level) %>%
  gather(key = variable, value = value, location, shape) %>%
  unite(temp, species, variable) %>%
  spread(temp, value) -> larval_scenarios













