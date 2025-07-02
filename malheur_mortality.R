##################################################
######## R Code for Johnston et al. 2025, ########
###### Significant mortality of old trees ########
### across a dry forest landscape, Oregon, USA ###
##################################################

#Contact jamesjoh@uoregon.edu with comments or suggestions for improvement.
#Note that some of this code may take several minutes or more to run.  
#The graphical outputs are optimized for my MacBook Pro and some of the graphical workflow is a little idiosyncratic.  You will probably want to customize the graphical output for a good fit with your workflow.  

##################################################
## Load data, create some columns and summaries ##
##################################################

#Clear the workspace
rm(list=ls())

#Load libraries needed for this section of code
library(tidyverse)
library(broom)

#Load focal tree data
focal_trees <- read.csv("focal_trees.csv")
head(focal_trees)

#Create age class column
focal_trees$age_class <- ifelse(focal_trees$age <150, "<150", ifelse(focal_trees$age >=150 & focal_trees$age <300, "150-300", "300+"))

#Add column indicating whether the site experienced fire during the observation period
fire_site_list <- c("cny", "dca", "nma", "rey")
focal_trees$fire <- ifelse(focal_trees$site1 %in% fire_site_list, "Y", "N")
table(focal_trees$fire, focal_trees$site1)

#Summary of trees by status:  alive (A), dead (D) or not relocated(U)
status_count <- focal_trees %>%
  group_by(status) %>%
  tally() %>%
  arrange(-n)
status_count
sum(status_count$n)

#Reduce data to relocated trees
focal_trees <- focal_trees %>% filter(status=="A" | status=="D")
nrow(focal_trees)

#Count of trees by age class and sample strategy—whether the tree is in the objective sample of species and age or not
sample_count <- focal_trees %>%
  group_by(objective, age_class) %>%
  tally() %>%
  arrange(objective, -n)
sample_count
sum(sample_count$n)

#Counts of species
table(focal_trees$spp)

#Summary of age 
ages <- focal_trees %>%
  #filter(objective=="Y")  %>%
  summarise(min=min(age, na.rm = TRUE),
            max=max(age, na.rm = TRUE),
            mean=mean(age, na.rm = TRUE)) %>%
  arrange(mean)
ages

#Summary of size
sizes <- focal_trees %>%
  summarise(min=min(dbh_cm, na.rm = TRUE),
            max=max(dbh_cm, na.rm = TRUE),
            mean=mean(dbh_cm, na.rm = TRUE)) %>%
  arrange(mean)
sizes

#Graph age~species
#Define species names and colors
species_names <- c(
  "PICO" = "Lodgepole",
  "JUOC" = "Juniper",
  "LAOC" = "Larch",
  "PSME" = "Douglas-fir",
  "ABGR" = "Grand fir",
  "PIPO" = "Ponderosa"
)
species_colors <- c(
  "Lodgepole" = "gray44",
  "Juniper" = "cadetblue",
  "Larch" = "yellow3",
  "Douglas-fir" = "darkolivegreen3",
  "Grand fir" = "darkgreen",
  "Ponderosa" = "chocolate"
)

#Filter and relabel species
plot_data <- focal_trees %>%
  filter(!is.na(age), !is.na(dbh_cm), spp %in% names(species_names)) %>%
  mutate(species = factor(species_names[spp], levels = species_names))

#Fit regression models by species
reg_results <- plot_data %>%
  group_by(species) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(dbh_cm ~ age, data = .x)),
    coefs = map(model, coef),
    stats = map(model, glance),
    label = map2_chr(coefs, stats, ~ {
      intercept <- round(.x[1], 1)
      slope <- round(.x[2], 3)
      r2 <- round(.y$r.squared, 2)
      glue::glue("y = {intercept} + {slope}x\nR² = {r2}")
    }),
    label_x = 25,
    label_y = 165
  ) %>%
  select(species, label, label_x, label_y)

#Merge label data back into main plot data
plot_data_labeled <- left_join(plot_data, reg_results, by = "species")

#Plot age - species relationship
age_dbh_plot <- ggplot(plot_data_labeled, aes(x = age, y = dbh_cm, color = species)) +
  geom_point(alpha = 0.6, size = 1.25) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
  geom_text(aes(x = label_x, y = label_y, label = label), 
            hjust = 0, vjust = 1, size = 3.5, family = "Helvetica", color = "black") +
  facet_wrap(~ species, scales = "fixed") +
  scale_color_manual(values = species_colors, name = "") +
  scale_x_continuous("Age (years)", expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous("DBH (cm)", limits = c(0, 175), breaks = c(50, 100, 150)) +
  theme_bw(base_size = 14, base_family = "Helvetica") +
  theme(
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "plain"),
    legend.position = "bottom",
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )
ggsave(
  filename = "age_dbh_plot.png",
  plot = age_dbh_plot,
  path = "/Users/james/Desktop",
  width = 9,
  height = 6.5,
  units = "in",
  dpi = 600
)

#Calculate mortality by site
mortality_summary1 <- focal_trees %>%
  group_by(site1) %>%
  summarise(
    total_trees = n(),
    dead_trees = sum(status == "D"),
    mortality_prop = dead_trees / total_trees,
    .groups = "drop"
  ) %>%
  arrange(mortality_prop)
data.frame(mortality_summary1)

#Calculate mortality site and age class
mortality_summary2 <- focal_trees %>%
  group_by(site1, age_class) %>%
  summarise(
    total_trees = n(),
    dead_trees = sum(status == "D"),
    mortality_prop = dead_trees / total_trees,
    .groups = "drop"
  ) %>%
  arrange(age_class, mortality_prop)
data.frame(mortality_summary2)

#Recode status="A" to "Alive and status="D" to "Dead"... this is mostly for the benefit of graphics to come
focal_trees <- focal_trees %>%
  filter(status %in% c("A", "D")) %>%
  mutate(status = recode(status, "A" = "Alive", "D" = "Dead"))


##############################################
## Summarize mortality by site for Figure 1 ##
##############################################

#Define age class mapping
age_map <- c("<150" = "Young", "150-300" = "Old", "300+" = "Very Old")

#Site-level and all-sites summaries
status_site <- bind_rows(
#By site
focal_trees %>%
  filter(age_class %in% names(age_map)) %>%
  mutate(age_cat = age_map[age_class]) %>%
  group_by(site1, age_cat, status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(site1, age_cat) %>%
   mutate(prop = n / sum(n)) %>%
  ungroup(),
#All sites combined
focal_trees %>%
  filter(age_class %in% names(age_map)) %>%
  mutate(age_cat = age_map[age_class]) %>%
  group_by(age_cat, status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(age_cat) %>%
  mutate(prop = n / sum(n), site1 = "All sites") %>%
  ungroup()
) %>%
mutate(
  site1 = if_else(tolower(site1) == "all sites", "All sites", toupper(site1)),
  age_class = factor(recode(age_cat, "Young" = "<150", "Old" = "150-300", "Very Old" = "300+"), levels = c("<150", "150-300", "300+")),age_cat = factor(age_cat, levels = c("Very Old", "Old", "Young")),
    site1 = factor(site1, levels = c('DEE', 'JUG', 'THO', 'REY', 'DCA', 'CNY', 'STI', 'NMA', 'CRA', 'MAL', 'MYR', 'All sites'))
  ) 
data.frame(status_site)

#Make side graphic portion of Figure 1
site_mortality <- ggplot(status_site, aes(x = age_class, y = prop, fill = status)) +
  geom_col(position = "stack", width = 0.95) +
  scale_fill_manual(name = "", values = c('Dead' = 'red3', 'Alive' = 'darkolivegreen4')) +
  scale_y_continuous("", position = "left", expand = c(0, 0), breaks = c(0, 0.5), limits = c(0, 1)) +
  scale_x_discrete("", expand = c(0, 0)) +
  facet_grid(site1 ~ .) +
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text.y.left = element_text(angle = 90),
    panel.spacing.y = unit(0.25, "lines"),
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
#site_mortality
ggsave(path = "/Users/james/Desktop", filename = "site_mortality.png", width = 2, height = 6.85, units = "in", dpi = 500)


##########################################
## Make Figure 2 age structure graphics ##
##########################################

#Another library for this section
library(ggpubr)

#Filter to the objective sample of age and species
mal_subplots <- focal_trees %>% filter(objective=="Y")

#Reorder species factors
focal_trees$spp <- factor(focal_trees$spp, levels = c('PIMO', 'PIEN', 'CADE', 'TSME', 'PICO', 'JUOC', 'LAOC', 'PSME', 'ABGC', 'ABGR', 'PIPO'))
mal_subplots$spp <- factor(mal_subplots$spp, levels = c('PIMO', 'PIEN', 'CADE', 'TSME', 'PICO', 'JUOC', 'LAOC', 'PSME', 'ABGC', 'ABGR', 'PIPO'))

#Age histogram 
all_age_structure <- ggplot() +
  geom_vline(aes(xintercept = 150), linewidth=.25, linetype="dashed") +
  geom_vline(aes(xintercept = 300), linewidth=.25, linetype="dashed") +
  geom_histogram(data=focal_trees, aes(x=age, fill=spp), binwidth=3) +
  scale_fill_manual(name = "Species", values = c('PICO' = 'gray44', 'JUOC' = 'cadetblue', 'LAOC' = 'yellow3', 'PSME' = 'darkolivegreen3', 'ABGR' = 'darkgreen', 'PIPO' = 'chocolate'), labels = c("Lodgepole (n=22)", "Juniper (n=31)", "Larch (n=37)", "Douglas-fir (n=215)", "Grand fir (n=296)", "Ponderosa (n=1,016)")) +
  annotate("text", x=75, y=46.7, label= "n=668", size=3.75) + 
  annotate("text", x=225, y=46.7, label= "n=527", size=3.75) + 
  annotate("text", x=400, y=46.7, label= "n=422", size=3.75) + 
  scale_x_continuous("Age (years)", breaks = seq(0, 700, 100), limits=c(0, 700)) + 
  scale_y_continuous("Tree count", breaks = c(0, 20, 40), limits=c(0, 47)) + 
  ggtitle("A. All sampled trees (for all statistical models)") +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(legend.title=element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=90,hjust=1,vjust=0.25)) +
  theme(panel.grid.minor = element_line(linewidth = 0.1), panel.grid.major = element_line(linewidth = .3)) +
  theme(strip.text.x = element_blank(), strip.background = element_rect(colour="white", fill="white"), 
        legend.position = "inside",
        legend.position.inside = c(.94,.825),
        legend.justification = c(.87,.725),
        legend.box.margin = margin(5, l = 5, unit = "mm"),
        legend.background = element_rect(
          fill = alpha("grey95", .2), linewidth =  0.1,
          linetype = "solid", color = "grey40"
        )) +
  theme(plot.title = element_text(size = 14, hjust = .5, vjust = -7.3))

#Objective sample age structure
sub_age_structure <- ggplot() +
  geom_vline(aes(xintercept = 150), linewidth=.25, linetype="dashed") +
  geom_vline(aes(xintercept = 300), linewidth=.25, linetype="dashed") +
  geom_histogram(data=mal_subplots, aes(x=age, fill=spp), binwidth=3) +
  scale_fill_manual(name = "Species", values = c('PICO' = 'gray44', 'JUOC' = 'cadetblue', 'LAOC' = 'yellow3', 'PSME' = 'darkolivegreen3', 'ABGR' = 'darkgreen', 'PIPO' = 'chocolate'), labels = c("Lodgepole (n=15)", "Juniper (n=21)", "Larch (n=4)", "Douglas-fir (n=89)", "Grand fir (n=160)", "Ponderosa (n=194) ")) +
  annotate("text", x=75, y=46.7, label= "n=369", size=3.75) + 
  annotate("text", x=225, y=46.7, label= "n=75", size=3.75) + 
  annotate("text", x=400, y=46.7, label= "n=39", size=3.75) + 
  scale_x_continuous("Age (years)", breaks = seq(0, 700, 100), limits=c(0, 700)) + 
  scale_y_continuous("Tree count", breaks = c(0, 20, 40), limits=c(0, 47)) + 
  ggtitle("B. Unbiased sample of species and age class (for simulations)") +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=90,hjust=1,vjust=0.25)) +
  theme(panel.grid.minor = element_line(linewidth = 0.1), panel.grid.major = element_line(linewidth = .3)) +
  theme(strip.text.x = element_blank(), strip.background = element_rect(colour="white", fill="white"), 
        legend.position = "inside",
        legend.position.inside = c(.94,.825),
        legend.justification = c(.87,.725),
        legend.box.margin = margin(5, l = 5, unit = "mm"),
        legend.background = element_rect(
          fill = alpha("grey95", .2), linewidth =  0.1,
          linetype = "solid", color = "grey40"
        )) +
  theme(plot.title = element_text(size = 14, hjust = .5, vjust = -7.3))

require(grid)   # for the textGrob() function
quartz(width=7.5, height=7)
age_structure_plot <- ggarrange(all_age_structure + rremove("ylab") + rremove("xlab"), sub_age_structure + rremove("ylab") + rremove("xlab"), heights = c(1, 1), ncol = 1, nrow = 2, align = "v") 
annotate_figure(age_structure_plot, left = textGrob("Tree counts", rot = 90, vjust = 1, gp = gpar(cex = 1.3)), bottom = textGrob("Age (years)", gp = gpar(cex = 1.3)))
ggsave(path="/Users/james/Desktop", "age_structure_plot.png", width = 7.5, height = 7, units = "in", dpi=600)


#########################################
#### Calculate annualized mortality #####
#########################################

#Note:  Running this section may take up to 30 seconds (calculation of CIs is computationally intensive)

#Another library for this section
library(scales)

#Function for bootstrapping 95% CIs for annualized mortality
bootstrap_mortality_ci <- function(data, n_boot = 1000) {
  reps <- map(1:n_boot, ~ {
    resampled <- sample_frac(data, replace = TRUE)
    alive <- sum(resampled$status == "Alive")
    total <- nrow(resampled)
    if (total == 0) return(NA_real_)
    1 - (alive / total)^(1 / 10)
  })
  q <- quantile(unlist(reps), probs = c(0.025, 0.975), na.rm = TRUE)
  tibble(mort_rate_lower = q[[1]], mort_rate_upper = q[[2]])
}

#Function for calculating annualized mortality and CIs by specified groupings
summarize_mortality <- function(data, group_vars, rename_map = NULL, recodes = list()) {
  # Check for required columns
  missing_cols <- setdiff(group_vars, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing column(s) in data: ", paste(missing_cols, collapse = ", "))
  }
  alive_summary <- data %>%
    filter(status == "Alive") %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(final_census = n_distinct(sample_id), .groups = "drop")
  all_summary <- data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(initial_census = n_distinct(sample_id), .groups = "drop") %>%
    left_join(alive_summary, by = group_vars) %>%
    mutate(
      final_census = replace_na(final_census, 0),
      annualized_mortality = 1 - (final_census / initial_census)^(1 / 10)
    )
  ci_summary <- data %>%
    group_by(across(all_of(group_vars))) %>%
    group_split() %>%
    map_dfr(~ {
      meta <- .x %>% select(all_of(group_vars)) %>% distinct()
      ci <- bootstrap_mortality_ci(.x)
      bind_cols(meta, ci)
    })
  joined <- all_summary %>%
    left_join(ci_summary, by = group_vars) %>%
    mutate(
      annualized_mortality = percent(annualized_mortality),
      mort_rate_lower = percent(mort_rate_lower),
      mort_rate_upper = percent(mort_rate_upper)
    )
  for (v in names(recodes)) {
    joined[[v]] <- recode_factor(joined[[v]], !!!recodes[[v]])
  }
  if (!is.null(rename_map)) {
    joined <- rename(joined, !!!rename_map)
  }
  joined
}

#Levels for ordering
fire_levels <- c("Unburned", "Burned")
age_levels <- c("<150", "150-300", "300+")

#Annualized mortality by species
spp_all_ci <- summarize_mortality(
  focal_trees, "spp",
  rename_map = c(
    Species = "spp",
    `Initial census` = "initial_census",
    `Final census` = "final_census",
    `Mortality yr` = "annualized_mortality",
    `Lower CI` = "mort_rate_lower",
    `Upper CI` = "mort_rate_upper"
  ),
  recodes = list(
    spp = c(PIPO = "P. ponderosa", ABGR = "A grandis", PSME = "P. menziesii", LAOC = "L. occidentalis", PICO = "P. contorta", JUOC = "J. occidentalis")
  )
) %>%
  arrange(Species)
spp_all_ci

#Annualized mortality by age class
stopifnot("age_class" %in% names(focal_trees))
ageclass_all_ci <- summarize_mortality(
  focal_trees, "age_class",
  rename_map = c(
    `Age class (years old)` = "age_class",
    `Initial census` = "initial_census",
    `Final census` = "final_census",
    `Mortality yr` = "annualized_mortality",
    `Lower CI` = "mort_rate_lower",
    `Upper CI` = "mort_rate_upper"
  )
) %>%
  mutate(`Age class (years old)` = factor(`Age class (years old)`, levels = age_levels)) %>%
  arrange(`Age class (years old)`)
ageclass_all_ci

#Annualized mortality by age class and fire status
ageclass_spp_all_ci <- summarize_mortality(
  focal_trees, c("fire", "age_class"),
  rename_map = c(
    `Fire status` = "fire",
    `Age class (years old)` = "age_class",
    `Initial census` = "initial_census",
    `Final census` = "final_census",
    `Mortality yr` = "annualized_mortality",
    `Lower CI` = "mort_rate_lower",
    `Upper CI` = "mort_rate_upper"
  ),
  recodes = list(fire = c(N = "Unburned", Y = "Burned"))
) %>%
  mutate(
    `Fire status` = factor(`Fire status`, levels = fire_levels),
    `Age class (years old)` = factor(`Age class (years old)`, levels = age_levels)
  ) %>%
  arrange(`Fire status`, `Age class (years old)`)
ageclass_spp_all_ci

#Annualized mortality by fire status, age class and species
fireetc_all_ci <- summarize_mortality(
  focal_trees, c("fire", "spp", "age_class"),
  rename_map = c(
    `Fire status` = "fire",
    Species = "spp",
    `Age class (years old)` = "age_class",
    `Initial census` = "initial_census",
    `Final census` = "final_census",
    `Mortality yr` = "annualized_mortality",
    `Lower CI` = "mort_rate_lower",
    `Upper CI` = "mort_rate_upper"
  ),
  recodes = list(
    spp = c(PIPO = "P. ponderosa", ABGR = "A grandis", PSME = "P. menziesii", LAOC = "L. occidentalis", PICO = "P. contorta", JUOC = "J. occidentalis"),
    fire = c(N = "Unburned", Y = "Burned")
  )
) %>%
  mutate(
    `Fire status` = factor(`Fire status`, levels = fire_levels),
    `Age class (years old)` = factor(`Age class (years old)`, levels = age_levels)
  ) %>%
  arrange(`Fire status`, `Age class (years old)`, Species)
fireetc_all_ci

#Annualized mortality by fire status
fire_all_ci <- summarize_mortality(
  focal_trees, "fire",
  rename_map = c(
    `Fire status` = "fire",
    `Initial census` = "initial_census",
    `Final census` = "final_census",
    `Mortality yr` = "annualized_mortality",
    `Lower CI` = "mort_rate_lower",
    `Upper CI` = "mort_rate_upper"
  ),
  recodes = list(fire = c(N = "Unburned", Y = "Burned"))
) %>%
  mutate(`Fire status` = factor(`Fire status`, levels = fire_levels)) %>%
  arrange(`Fire status`)
fire_all_ci

#Annualized mortality all trees
initial_census <- nrow(focal_trees)
final_census <- sum(focal_trees$status == "Alive")
overall_rate <- 1 - (final_census / initial_census)^(1 / 10)
overall_ci <- bootstrap_mortality_ci(focal_trees)
overall_mortality_final <- tibble(
  `Initial census` = initial_census,
  `Final census` = final_census,
  `Mortality yr` = percent(overall_rate, accuracy = 0.01),
  `Lower CI` = percent(overall_ci$mort_rate_lower, accuracy = 0.01),
  `Upper CI` = percent(overall_ci$mort_rate_upper, accuracy = 0.01)
)
overall_mortality_final

#Summaries
print(spp_all_ci, n = 100, width = Inf)
print(ageclass_all_ci, n = 100, width = Inf)
print(ageclass_spp_all_ci, n = 100, width = Inf)
print(fireetc_all_ci, n = 100, width = Inf)
print(overall_mortality_final, n = 100, width = Inf)
print(fire_all_ci, n = 100, width = Inf)

## Combine and export summaries as one .csv ##

#Library for this portion
library(readr)

#Define the target column structure
standard_cols <- c("Fire status", "Age class (years old)", "Species", "Initial census", "Final census", "Mortality yr", "Lower CI", "Upper CI")

#Column standardization function
standardize <- function(df) {
  df <- df %>% mutate(across(everything(), as.character))
  for (col in standard_cols) {
    if (!col %in% names(df)) {
      df[[col]] <- NA_character_
    }
  }
  df %>% select(all_of(standard_cols))
}

#Apply to all data frames
combined_tbl <- bind_rows(
  standardize(spp_all_ci),
  standardize(ageclass_all_ci),
  standardize(ageclass_spp_all_ci),
  standardize(fireetc_all_ci),
  standardize(overall_mortality_final),
  standardize(fire_all_ci)
)

#Write to .csv (you'll have to do just a little manual clean up of this document to produce Table 2)
write_csv(combined_tbl, "~/Desktop/mortality_summary_combined.csv")

#We could test these summaries manually to make sure the code above is working as expected.  Manual calculations of subsets of data is also useful for some specific calculations, for instance mortality of ponderosa pine of different age classes and other statistics reported in the text.  
#Some test data
test_data <- focal_trees %>% filter(spp=="PIPO" & age>100  & fire=="N")
#test_data <- focal_trees %>% filter(spp=="PIPO" & age_class=="300+")
#test_data <- focal_trees %>% filter(spp=="PIPO" & age_class=="<150" & fire=="N")
#test_data <- focal_trees %>% filter(spp=="ABGR")
#test_data <- focal_trees %>% filter(age_class=="300+")
nrow(test_data)
table(test_data$status)
#Run the manual version of Eq. 6 from Sheil et al. on the nrow call (initial census) and the final census from the dead column of the table
manual_annualized_mortality=1-(564/741)^(1/10)
percent(manual_annualized_mortality, accuracy = 0.01)


########################################
###### Prop tests to evaluate ##########
## mortality by species and age class ##
########################################

#Another library for this section
library(freqtables)

#Filter to sites that experienced fire and didn't experience fire
fire_site_list <- c("cny", "dca", "nma", "rey")
focal_trees_fire <- focal_trees %>% filter(site1 %in% fire_site_list)
focal_trees_nofire <- focal_trees %>% filter(!site1 %in% fire_site_list)

#Min, max, and mean mortality/survivorship by sites that didn't experience fire
#Alive vs. dead by site
status_site_nofire <- focal_trees_nofire %>%
  group_by(site1, status) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))
data.frame(status_site_nofire)
#Combined
nofire_site_summary <- status_site_nofire  %>%
  group_by(status) %>%
  summarize(min_prop=min(prop),
            max_prop=max(prop),
            mean_prop=mean(prop))
nofire_site_summary

#Min, max, and mean mortality/survivorship by sites that experienced fire
#Alive vs. dead by site
status_site_fire <- focal_trees_fire %>%
  group_by(site1, status) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))
data.frame(status_site_fire)
#Combined
fire_site_summary <- status_site_fire  %>%
  group_by(status) %>%
  summarize(min_prop=min(prop),
            max_prop=max(prop),
            mean_prop=mean(prop))
fire_site_summary

#Function to compute summary stats (min, max, mean proportion dead/alive by group)
site_status_summary <- function(data, group_var) {
  data %>%
    group_by({{group_var}}, status) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(prop = n / sum(n)) %>%
    group_by(status) %>%
    summarise(min_prop = min(prop), max_prop = max(prop), mean_prop = mean(prop), .groups = "drop")
}

#Function to create mortality table, prop test, and pairwise Fisher's tests
mortality_tests <- function(data, group_var, label = "Group") {
  tab <- data %>%
    group_by({{group_var}}, status) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
    mutate(Total = Alive + Dead)
  #Prop test
  death_vector <- tab$Dead
  total_vector <- tab$Total
  names(death_vector) <- tab[[as_label(enquo(group_var))]]
  names(total_vector) <- tab[[as_label(enquo(group_var))]]
  prop_test <- prop.test(death_vector, total_vector)
  #Pairwise Fisher tests
  pairs <- combn(tab[[as_label(enquo(group_var))]], 2, simplify = FALSE)
  fisher_results <- map_dfr(pairs, function(pair) {
    g1 <- tab %>% filter(!!sym(as_label(enquo(group_var))) == pair[1])
    g2 <- tab %>% filter(!!sym(as_label(enquo(group_var))) == pair[2])
    m <- matrix(c(g1$Dead, g1$Alive, g2$Dead, g2$Alive), nrow = 2, byrow = TRUE)
    test <- fisher.test(m)
    tibble(
      !!paste0(label, "1") := pair[1],
      !!paste0(label, "2") := pair[2],
      p_value = test$p.value
    )
  }) %>%
    mutate(p_adj_holm = p.adjust(p_value, method = "holm"),
           p_adj_bonf = p.adjust(p_value, method = "bonferroni"))
  
  list(summary_table = tab, prop_test = prop_test, fisher = fisher_results)
}

#Function to fit quadratic GLM for mortality by age... this isn't reported in the paper, see below.
fit_quadratic_glm <- function(data) {
  data %>%
    mutate(mortality = if_else(status == "Dead", 1, 0)) %>%
    glm(mortality ~ age + I(age^2), family = binomial, data = .)
}

#Run summaries and tests
nofire_summary <- site_status_summary(focal_trees_nofire, site1)
fire_summary   <- site_status_summary(focal_trees_fire, site1)

all_spp_tests       <- mortality_tests(focal_trees, spp, label = "Species")
unburned_spp_tests  <- mortality_tests(focal_trees_nofire, spp, label = "Species")
burned_spp_tests    <- mortality_tests(focal_trees_fire, spp, label = "Species")

all_age_tests       <- mortality_tests(focal_trees, age_class, label = "group")
unburned_age_tests  <- mortality_tests(focal_trees_nofire, age_class, label = "group")
burned_age_tests    <- mortality_tests(focal_trees_fire, age_class, label = "group")

#Run GLMs... this is bonus material.  We didn't report on the GLM results in the paper, but they are consistent with findings of prop tests and with the Random Forest models
glm_all       <- fit_quadratic_glm(focal_trees)
glm_unburned  <- fit_quadratic_glm(focal_trees_nofire)
glm_burned    <- fit_quadratic_glm(focal_trees_fire)

#View results
nofire_summary 
fire_summary  

all_spp_tests       
unburned_spp_tests  
burned_spp_tests    

all_age_tests       
unburned_age_tests  
burned_age_tests    

glm_all       
glm_unburned  
glm_burned    


##################################################
## Proportion mortality CIs and Figures 5 and 6 ##
##################################################

#Function to generate summary tables using freq_table function
generate_summary <- function(data, row_var, col_var, fire_label) {
  data %>%
    freq_table({{row_var}}, {{col_var}}) %>%
    select(row_cat, col_cat, n, percent_row, lcl_row, ucl_row) %>%
    drop_na() %>%
    mutate(fire = fire_label)
}

#Manually adjust the JUOC row
add_juoc_row <- function(df) {
  juoc_rows <- df %>% filter(row_cat == "JUOC")
  #If "Dead" is missing, add it with n = 0 but NA CIs so no bar/error shown
  if (!"Dead" %in% juoc_rows$col_cat) {
    juoc_rows <- bind_rows(
      juoc_rows,
      tibble(
        row_cat = "JUOC", col_cat = "Dead", n = 0,
        percent_row = 0, lcl_row = NA_real_, ucl_row = NA_real_,
        fire = "B. Burned sites"
      )
    )
  }
  #If "Alive" is missing (defensive)
  if (!"Alive" %in% juoc_rows$col_cat) {
    juoc_rows <- bind_rows(
      juoc_rows,
      tibble(
        row_cat = "JUOC", col_cat = "Alive", n = 1,
        percent_row = 100, lcl_row = 100, ucl_row = 100,
        fire = "B. Burned sites"
      )
    )
  }
  #Combine with rest of dataset
  df %>% filter(row_cat != "JUOC") %>% bind_rows(juoc_rows)
}

#Create species proportion mortality summary tables
status_species_fire <- generate_summary(focal_trees_fire, spp, status, "B. Burned sites") %>%
  add_juoc_row()
status_species_nofire <- generate_summary(focal_trees_nofire, spp, status, "A. Unburned sites")
status_species <- bind_rows(status_species_nofire, status_species_fire) %>%
  mutate(
    row_cat = recode_factor(
      row_cat,
      PIPO = "Ponderosa", ABGR = "Grand fir",
      PSME = "Douglas-fir", LAOC = "Larch",
      PICO = "Lodgepole", JUOC = "Juniper"
    ),
    fire = factor(fire, levels = c("A. Unburned sites", "B. Burned sites")),
    col_cat = factor(col_cat, levels = c("Alive", "Dead"))  # << this line is new
  )

#Create age class proportion mortality summary tables
status_age_nofire <- generate_summary(focal_trees_nofire, age_class, status, "A. Unburned sites")
status_age_fire <- generate_summary(focal_trees_fire, age_class, status, "B. Burned sites")
status_age <- bind_rows(status_age_nofire, status_age_fire) %>%
  mutate(
    row_cat = factor(row_cat, levels = c("<150", "150-300", "300+")),
    fire = factor(fire, levels = c("A. Unburned sites", "B. Burned sites"))
  )

#Species barplot
status_species_plot <- ggplot(status_species, aes(
  x = row_cat, y = percent_row / 100,
  ymin = lcl_row / 100, ymax = ucl_row / 100,
  fill = col_cat
)) +
  geom_col(position = "dodge") +
  geom_errorbar(position = position_dodge(0.9), width = 0.15, linewidth = 0.3) +
  geom_text(
    aes(y = ucl_row / 100 + 0.025, label = n),
    position = position_dodge(0.85),
    size = 2.1, color = "black"
  ) +
  scale_fill_manual(
    name = "", values = c('Dead' = 'red3', 'Alive' = 'darkolivegreen4'),
    labels = c("Lived", "Died")
  ) +
  facet_wrap(~fire) +
  scale_y_continuous("Proportion") +
  scale_x_discrete("") +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(
    strip.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major = element_line(colour = "grey50", linewidth = 0.05),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "right"
  )
ggsave("/Users/james/Desktop/status_species_plot.png", status_species_plot, width = 7, height = 4, units = "in", dpi = 500)

#Age class barplot
status_age_plot <- ggplot(status_age, aes(
  x = row_cat, y = percent_row / 100,
  ymin = lcl_row / 100, ymax = ucl_row / 100,
  fill = col_cat
)) +
  geom_col(position = "dodge") +
  geom_errorbar(position = position_dodge(0.9), width = 0.15, linewidth = 0.3) +
  geom_text(
    aes(y = ucl_row / 100 + 0.025, label = n),
    position = position_dodge(0.85),
    size = 2.1, color = "black"
  ) +
  scale_fill_manual(
    name = "", values = c('Dead' = 'red3', 'Alive' = 'darkolivegreen4'),
    labels = c("Lived", "Died")
  ) +
  facet_wrap(~fire) +
  scale_y_continuous("Proportion") +
  scale_x_discrete("Age class (years)") +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(
    strip.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major = element_line(colour = "grey50", linewidth = 0.05),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    legend.position = "right"
  )
ggsave("/Users/james/Desktop/status_age_plot.png", status_age_plot, width = 7, height = 4, units = "in", dpi = 500)


##################################################
## Variable selection for  Random Forest models ##
##################################################

#Libraries for this section
library(Boruta)

##Load data, create some columns, and clean data

#Let's clear the workspace and start from scratch to keep the environment from becoming to cluttered.
rm(list=ls())

#Load libraries needed for this section of code
library(tidyverse)

#Load focal tree data
focal_trees <- read.csv("focal_trees.csv")
head(focal_trees)

#Reduce to relocated trees for all future analysis
focal_trees <- focal_trees %>% filter(status=="A" | status=="D")
nrow(focal_trees)

#Add column indicating whether the site experienced fire during the observation period
fire_site_list <- c("cny", "dca", "nma", "rey")
focal_trees$fire <- ifelse(focal_trees$site1 %in% fire_site_list, "Y", "N")

#Add a binary column for survival
focal_trees$survived <- ifelse(focal_trees$status=="A", 1, 0)
focal_trees$status <- NULL

#Any columns still have NAs?  Boruta won't run on NAs
names(which(colSums(is.na(focal_trees)) > 0))
#Filter to NA rows in a BAI column
has_nas <- focal_trees[is.na(focal_trees$bai50_total),]
has_nas$sample_id
#Drop the NA rows and columns we know we won't use and rename that file "all_trees".
all_trees <- focal_trees %>% tidyr::drop_na()
all_trees <- all_trees %>% select(-latitude, -longitude, -sample_id, -objective)

#Create files for trees in unburned and burned sites
#Filter to unburned sites
burned_sites <- c("cny", "dca", "nma", "rey")
unburned_trees <- all_trees %>% filter(!site1 %in% burned_sites)
burned_trees <- all_trees %>% filter(site1 %in% burned_sites)
table(unburned_trees$site1)
table(burned_trees$site1)

##Unburned sites

#Variable selection with Boruta for trees in unburned sites
set.seed(234)
unburned_bor <- Boruta(survived~., data = unburned_trees, doTrace = 2)

#Plot unburned sites Boruta results 
par(mar = c(9, 3, 3, 3))
plot(unburned_bor, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(unburned_bor$ImpHistory),function(i)
  unburned_bor$ImpHistory[is.finite(unburned_bor$ImpHistory[,i]),i])
names(lz) <- colnames(unburned_bor$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(unburned_bor$ImpHistory), cex.axis = 0.9)

#Get importance values for unburned site trees
getSelectedAttributes(unburned_bor, withTentative = F)
unburned_imp_df <- attStats(unburned_bor)
unburned_imp_df <- unburned_imp_df[order(-unburned_imp_df$meanImp),]
unburned_imp_df

#The next step is to examine correlation charts and stepwise eliminate variables that are correlated with more explanatory variables and rerun Boruta until you're left with high performing non-correlated variables.  
corrs_test <- unburned_trees %>% select(-site1, -site2, -spp, -pvt, -forest_group, -survived, -fire)
PerformanceAnalytics::chart.Correlation(corrs_test, histogram=TRUE, pch=19)

##Burned sites

#Variable selection with Boruta for trees in burned sites
set.seed(234)
burned_bor <- Boruta(survived~., data = burned_trees, doTrace = 2)

#Plot burned sites Boruta results 
par(mar = c(9, 3, 3, 3))
plot(burned_bor, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(burned_bor$ImpHistory),function(i)
  burned_bor$ImpHistory[is.finite(burned_bor$ImpHistory[,i]),i])
names(lz) <- colnames(burned_bor$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(burned_bor$ImpHistory), cex.axis = 0.9)

#Get importance values for burned site trees
getSelectedAttributes(burned_bor, withTentative = F)
burned_imp_df <- attStats(burned_bor)
burned_imp_df <- burned_imp_df[order(-burned_imp_df$meanImp),]
burned_imp_df

#The next step is to examine correlation charts and stepwise eliminate variables that are correlated with more explanatory variables and rerun Boruta until you're left with high performing non-correlated variables.  
corrs_test <- burned_trees %>% select(-site1, -site2, -spp, -pvt, -forest_group, -survived, -fire)
PerformanceAnalytics::chart.Correlation(corrs_test, histogram=TRUE, pch=19)


##########################################
### Model mortality with random forest ###
##########################################

#Let's clear the workspace and start from scratch to keep the environment from becoming to cluttered.
rm(list=ls())

#Reload libraries from above and load some new libraries for this section
library(tidyverse)
library(tidymodels)
library(ranger)
library(cowplot)
library(pdp)
library(patchwork)
library(grid)

##Load data, create some columns, and clean data

#Load libraries needed for this section of code
library(tidyverse)

#Load focal tree data
focal_trees <- read.csv("focal_trees.csv")
head(focal_trees)

#Reduce to relocated trees
focal_trees <- focal_trees %>% filter(status=="A" | status=="D")
nrow(focal_trees)

#Create binary survival response variable
focal_trees$survived <- ifelse(focal_trees$status=="A", 1, 0)
focal_trees$status <- NULL

#Add column indicating whether the site experienced fire during the observation period
fire_site_list <- c("cny", "dca", "nma", "rey")
focal_trees$fire <- ifelse(focal_trees$site1 %in% fire_site_list, "Y", "N")

#Select best variables, i.e., the most predictive non-correlated (>0.65 or <-0.65) variables from Boruta procedures above
all_trees <- focal_trees %>% select(tree_ht_m, site2, bai5_scaled_mean, temp_mean_c, precip_mm, modis_npp, cwd_1981_2010, slope, basal15_m2ha, trees15_ha, forest_group, pvt, age, spp, asw_mm, survived)

#Any columns still have NAs?  
names(which(colSums(is.na(all_trees)) > 0))
#Filter to NA rows in a BAI column 
has_nas <- all_trees[is.na(all_trees$bai50_total),]
has_nas$sample_id
#Drop the NA rows
all_trees <- all_trees %>% tidyr::drop_na()

#Remove sample ID, which we don't need anymore, convert character strings to factors, and view file
all_trees$sample_id <- NULL
all_trees <- all_trees %>%
  dplyr::mutate_if(is.character, factor)
all_trees$survived <- as.factor(all_trees$survived)
str(all_trees)
head(all_trees)
nrow(all_trees)

#Filter to burned sites and unburned sites... with respect to the site variable, note that we are going to model plot or transect within sites from this point forward.
burned_sites <- c("cny1", "cny2", "cny3", "cny5", "D3", "G3", "G4", "dca1", "dca2", "dca3", "rey1", "rey2", "rey3", "nma1", "nma2", "nma3")
unburned <- all_trees %>% filter(!site2 %in% burned_sites)
burned <- all_trees %>% filter(site2 %in% burned_sites)
table(unburned$site2)
table(burned$site2)

## Random Forest model with tuned hyper-parameters— unburned sites... this may take several minutes to run

#Split into training and testing sets
set.seed(123)
unburned_split  <- initial_split(unburned, prop = 0.7, strata = "survived")
unburned_train  <- training(unburned_split)
unburned_test   <- testing(unburned_split)
nrow(unburned_train)
nrow(unburned_test)

#Get number of features (the number of predictors)
unburned_n_features <- length(setdiff(names(unburned_train), "survived"))

#Train a random forest model with default settings
unburned_rf_default <- ranger(
  survived ~ ., 
  data = unburned_train,
  mtry = floor(unburned_n_features / 3),
  respect.unordered.factors = "order",
  seed = 123
)

#Get OOB RMSE
unburned_default_rmse <- sqrt(unburned_rf_default$prediction.error)
unburned_default_rmse

#Create hyperparameter grid
unburned_hyper_grid <- expand.grid(
  mtry = floor(unburned_n_features * c(.05, .1, .15, .2, .25, .3, .333, .35, .4, .45)),
  min.node.size = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.1, .2, .3, .4, .5, .6, .7, .8, .9, 1 ),                       
  rmse = NA                                               
)

#Full cartesian grid search (last line exports OOB error)
for(i in seq_len(nrow(unburned_hyper_grid))) {
  unburned_fit <- ranger(
    formula         = survived ~ ., 
    data            = unburned_train, 
    num.trees       = unburned_n_features * 10,
    mtry            = unburned_hyper_grid$mtry[i],
    min.node.size   = unburned_hyper_grid$min.node.size[i],
    replace         = unburned_hyper_grid$replace[i],
    sample.fraction = unburned_hyper_grid$sample.fraction[i],
    verbose         = FALSE,
    seed            = 234,
    respect.unordered.factors = 'order',
  )
  unburned_hyper_grid$rmse[i] <- sqrt(unburned_fit$prediction.error)
}

#Have a look at the top 10 models
unburned_hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (unburned_default_rmse - rmse) / unburned_default_rmse * 100) %>%
  head(10)

#Re-run model with impurity-based variable importance and hyperparameters from the best model... The hyperparameters mtry, min.node.size and sample.fraction determine the degree of randomness, and should be tuned (Probst, Philipp, Marvin N Wright, and Anne-Laure Boulesteix. 2019. “Hyperparameters and Tuning Strategies for Random Forest.” Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery. Wiley Online Library, e1301.). The sample.fraction parameter specifies the fraction of observations to be used in each tree. Smaller fractions lead to greater diversity, and thus less correlated trees which often is desirable.  Note the probability = TRUE syntax for probability of outcome.
unburned_rf <- ranger(
  formula = survived ~ ., 
  data = unburned_train, 
  probability = TRUE,
  num.trees = 2000,
  mtry =  6,
  min.node.size = 2,
  sample.fraction = .4,
  replace = FALSE,
  importance = "impurity",
  #importance = "impurity_corrected",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed = 234
)

#Examine variable importance
unburned_vi <- data.frame(unburned_rf$variable.importance)
unburned_vi$variable <- rownames(unburned_vi)
colnames(unburned_vi) <- c("importance","variable")
rownames(unburned_vi) <- NULL
unburned_vi <- unburned_vi[order(-unburned_vi$importance),]
unburned_vi

#Predict for test set (note the funny syntax:all_test)$predictions for probabilities to put this result in terms of probabilities) 
predict_unburned <- cbind(unburned_test, predict(unburned_rf, data = unburned_test)$predictions)
colnames(predict_unburned)[which(names(predict_unburned) == "0")] <- "predict_0"
colnames(predict_unburned)[which(names(predict_unburned) == "1")] <- "predict_1"
head(predict_unburned)

#Evaluate accuracy using test data
#Kappa statistic:  < than 0.2 = poor agreement; 0.2–0.4 = fair agreement; 0.41–0.6 = moderate agreement; 0.61–0.8 = substantial agreement; > than 0.8 = great agreement (Landis, J. R., & Koch, G. G. (1977). An application of hierarchical kappa-type statistics in the assessment of majority agreement among multiple observers. Biometrics, 363-374)
predict_unburned$predicted_outcome <- ifelse(predict_unburned$predict_1 >= 0.50, "A", "D")
predict_unburned$actual_outcome <- ifelse(predict_unburned$survived=="1", "A", "D")
predict_unburned$predicted_outcome <- as.factor(predict_unburned$predicted_outcome)
predict_unburned$actual_outcome <- as.factor(predict_unburned$actual_outcome)
unburned_confusion <- caret::confusionMatrix(data=predict_unburned$predicted_outcome, reference = predict_unburned$actual_outcome)
unburned_confusion

## Random Forest model with tuned hyper-parameters for burned site trees 

#Split into training and testing sets
set.seed(234)
burned_split  <- initial_split(burned, prop = 0.7, strata = "survived")
burned_train  <- training(burned_split)
burned_test   <- testing(burned_split)
nrow(burned_train)
nrow(burned_test)

#Get number of features (the number of predictors)
burned_n_features <- length(setdiff(names(burned_train), "survived"))

#Train a random forest model with default settings
burned_rf_default <- ranger(
  survived ~ ., 
  data = burned_train,
  mtry = floor(burned_n_features / 3),
  respect.unordered.factors = "order",
  seed = 234
)

#Get OOB RMSE
burned_default_rmse <- sqrt(burned_rf_default$prediction.error)
burned_default_rmse

#Create hyperparameter grid
burned_hyper_grid <- expand.grid(
  mtry = floor(burned_n_features * c(.05, .1, .15, .2, .25, .3, .333, .35, .4, .45)),
  min.node.size = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.1, .2, .3, .4, .5, .6, .7, .8, .9, 1 ),                       
  rmse = NA                                               
)

#Full cartesian grid search (last line exports OOB error)
for(i in seq_len(nrow(burned_hyper_grid))) {
  burned_fit <- ranger(
    formula         = survived ~ ., 
    data            = burned_train, 
    num.trees       = burned_n_features * 10,
    mtry            = burned_hyper_grid$mtry[i],
    min.node.size   = burned_hyper_grid$min.node.size[i],
    replace         = burned_hyper_grid$replace[i],
    sample.fraction = burned_hyper_grid$sample.fraction[i],
    verbose         = FALSE,
    seed            = 123,
    respect.unordered.factors = 'order',
  )
  burned_hyper_grid$rmse[i] <- sqrt(burned_fit$prediction.error)
}

#Have a look at the top 10 models
burned_hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (burned_default_rmse - rmse) / burned_default_rmse * 100) %>%
  head(10)

#Re-run model with impurity-based variable importance and hyperparameters from the best model... The hyperparameters mtry, min.node.size and sample.fraction determine the degree of randomness, and should be tuned (Probst, Philipp, Marvin N Wright, and Anne-Laure Boulesteix. 2019. “Hyperparameters and Tuning Strategies for Random Forest.” Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery. Wiley Online Library, e1301.). The sample.fraction parameter specifies the fraction of observations to be used in each tree. Smaller fractions lead to greater diversity, and thus less correlated trees which often is desirable.  Note the probability = TRUE syntax for probability of outcome.
burned_rf <- ranger(
  formula = survived ~ ., 
  data = burned_train, 
  probability = TRUE,
  num.trees = 2000,
  mtry =  5,
  min.node.size = 4,
  sample.fraction = 1,
  replace = TRUE,
  importance = "impurity",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 567
)

#Variable importance
burned_vi <- data.frame(burned_rf$variable.importance)
burned_vi$variable <- rownames(burned_vi)
colnames(burned_vi) <- c("importance","variable")
rownames(burned_vi) <- NULL
burned_vi <- burned_vi[order(-burned_vi$importance),]
burned_vi

#Predict for test set (note the funny syntax:all_test)$predictions for probabilities to put this in terms of probabilities) 
predict_burned <- cbind(burned_test, predict(burned_rf, data = burned_test)$predictions)
colnames(predict_burned)[which(names(predict_burned) == "0")] <- "predict_0"
colnames(predict_burned)[which(names(predict_burned) == "1")] <- "predict_1"
head(predict_burned)

#Evaluate accuracy using test data
#Kappa statistic:  < than 0.2 = poor agreement; 0.2–0.4 = fair agreement; 0.41–0.6 = moderate agreement; 0.61–0.8 = substantial agreement; > than 0.8 = great agreement (Landis, J. R., & Koch, G. G. (1977). An application of hierarchical kappa-type statistics in the assessment of majority agreement among multiple observers. Biometrics, 363-374)
predict_burned$predicted_outcome <- ifelse(predict_burned$predict_1 >= 0.50, "A", "D")
predict_burned$actual_outcome <- ifelse(predict_burned$survived=="1", "A", "D")
predict_burned$predicted_outcome <- as.factor(predict_burned$predicted_outcome)
predict_burned$actual_outcome <- as.factor(predict_burned$actual_outcome)
burned_confusion <- caret::confusionMatrix(data=predict_burned$predicted_outcome, reference = predict_burned$actual_outcome)
burned_confusion

#Revisit correlations if desired
#corrs_test <- mal_focal %>% select(age, tree_ht)
#quartz(h=7, w=10)
#PerformanceAnalytics::chart.Correlation(corrs_test, histogram=TRUE, pch=19)

## Variable importance plot... we're going to build a bunch of plots and combine them as a last step

#Unburned variable importance plot
unburned_vi
#Rename the variables we'll plot... we're doing this manually for the purposes of a graphic, so take care the order of these variables matches the order of the unburned_vi "variable" column
unburned_vi$new_variable_name <- c("BAI", "Tree height", "Age", "Slope", "Site strata", "CWD", "Species", "Precip", "NPP", "Mean temp", "Basal area", "Density",  "PVT", "ASW", "Forest type")
unburned_vi_plot <- unburned_vi %>%
  arrange(importance) %>%
  mutate(new_variable_name=factor(new_variable_name, levels=new_variable_name)) %>%
  ggplot() +
  geom_col(aes(x=importance, y=new_variable_name, alpha=importance), fill = "midnightblue", width = .75) +
  scale_x_continuous(
    limits = c(0, 30),
    breaks = seq(0, 30, by = 10), 
    expand = c(0, 1), 
    position = "top"
  ) +
  scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
  scale_alpha(range = c(.35, 1)) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
    axis.ticks.length = unit(0, "mm"),
    axis.title = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(family = "Arial", size = 14)
  ) +
  geom_text(
    data = subset(unburned_vi, importance > 40),
    aes(0, y = new_variable_name, label = new_variable_name
    ),
    hjust = 0,
    nudge_x = 0.3,
    colour = "white",
    family = "Arial",
    size = 5.5
  ) + 
  geom_text(
    data = subset(unburned_vi, importance > 20 & importance < 40),
    aes(0, y = new_variable_name, label = new_variable_name
    ),
    hjust = 0,
    nudge_x = 0.3,
    colour = "white",
    family = "Arial",
    size = 4.75
  ) +
  geom_text(
    data = subset(unburned_vi, importance < 20 & importance >8),
    aes(0, y = new_variable_name, label = new_variable_name),
    hjust = 0,
    nudge_x = 0.3,
    colour = "white",
    family = "Arial",
    size = 3.75
  ) +
  geom_text(
    data = subset(unburned_vi, importance <8 & importance >2.6),
    aes(0, y = new_variable_name, label = new_variable_name),
    hjust = 0,
    nudge_x = 0.3,
    colour = "white",
    family = "Arial",
    size = 2.25
  ) +
  shadowtext::geom_shadowtext(
    data = subset(unburned_vi, importance <=2.6),
    aes(importance, y = new_variable_name, label = new_variable_name),
    hjust = 0,
    nudge_x = 0.3,
    colour = "midnightblue",
    bg.colour = "white",
    bg.r = 0.2,
    family = "Arial",
    size = 2.25) + 
  annotate(geom="text", x=17.5, y=6.85, label="Kappa=0.26", color="black", size=3.5) +
  ggtitle("A. Unburned stands variable importance") +
  theme(plot.title = element_text(hjust = 0.5, size=14)) +
  theme(legend.position="none", plot.margin=unit(c(.2,1,.2,1),"cm")) 
#unburned_vi_plot
#ggsave(path="/Users/james/Desktop", "vi_plot.png", width = 6, height = 8, units = "in", dpi=500)

#Extract the confusion matrix values as data.frame
unburned_cm_d <- as.data.frame(unburned_confusion$table)
#Confusion matrix statistics as data.frame
unburned_cm_st <-data.frame(unburned_confusion$overall)
#Round the values
unburned_cm_st$cm.overall <- round(unburned_cm_st$unburned_confusion.overall,3)
#Percentage
unburned_cm_p <- as.data.frame(prop.table(unburned_confusion$table))
unburned_cm_d$Perc <- round(unburned_cm_p$Freq*100,1)
unburned_cm_d$Prediction <- ifelse(unburned_cm_d$Prediction=="A", "Lived", "Died")
unburned_cm_d$Reference <- ifelse(unburned_cm_d$Reference=="A", "Lived", "Died")
unburned_confusion
unburned_cm_d

#Plot the matrix
unburned_cm_plot <-  ggplot(data = unburned_cm_d, aes(x = Prediction, y =  Reference, fill = Reference, alpha=Perc)) +
  geom_tile(color="black") +
  scale_fill_manual(name = "", values = c('Died' = 'red3', 'Lived' = 'darkolivegreen4')) +
  scale_alpha(range = c(.75, 1)) +
  geom_text(aes(label=paste0(Freq," ","(",Perc,"%",")")), color = 'black', size = 3.5) +
  theme_bw(base_size=11, base_family = "Arial") +
  theme(legend.position="none") 
#unburned_cm_plot

#Combine variable importance and confusion matrix plots
#quartz(width=5, height=7.75)
unburned_vi_confusion <- ggdraw() +
  draw_plot(unburned_vi_plot) +
  draw_plot(unburned_cm_plot,
            height = 0.38,
            width = .65,
            x = 0.35,
            y = .0005
  )
#unburned_vi_confusion

#Burned variable importance plot
burned_vi
#Rename the variables we'll plot... we're doing this manually for the purposes of a graphic, so take care the order of these variables matches the order of the unburned_vi "variable" column
burned_vi$new_variable_name <- c("Site strata", "Tree height", "BAI", "Slope", "Age", "Mean temp", "Precip", "PVT", "CWD", "Density", "Basal area", "ASW", "NPP", "Species", "Forest type")
burned_vi
burned_vi_plot <- burned_vi %>%
  arrange(importance) %>%
  mutate(new_variable_name=factor(new_variable_name, levels=new_variable_name)) %>%
  ggplot() +
  geom_col(aes(x=importance, y=new_variable_name, alpha=importance), fill = "midnightblue", width = .75) +
  scale_x_continuous(
    limits = c(0, 30),
    breaks = seq(0, 30, by = 10), 
    expand = c(0, 1), 
    position = "top"
  ) +
  scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
  scale_alpha(range = c(.35, 1)) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
    axis.ticks.length = unit(0, "mm"),
    axis.title = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(family = "Arial", size = 14)
  ) +
  geom_text(
    data = subset(burned_vi, importance > 20),
    aes(0, y = new_variable_name, label = new_variable_name
    ),
    hjust = 0,
    nudge_x = 0.3,
    colour = "white",
    family = "Arial",
    size = 4.75
  ) +
  geom_text(
    data = subset(burned_vi, importance < 20 & importance >7),
    aes(0, y = new_variable_name, label = new_variable_name),
    hjust = 0,
    nudge_x = 0.3,
    colour = "white",
    family = "Arial",
    size = 3.75
  ) +
  geom_text(
    data = subset(burned_vi, importance <7 & importance >2.75),
    aes(0, y = new_variable_name, label = new_variable_name),
    hjust = 0,
    nudge_x = 0.3,
    colour = "white",
    family = "Arial",
    size = 2.25
  ) +
  shadowtext::geom_shadowtext(
    data = subset(burned_vi, importance <=2.75),
    aes(importance, y = new_variable_name, label = new_variable_name),
    hjust = 0,
    nudge_x = 0.3,
    colour = "midnightblue",
    bg.colour = "white",
    bg.r = 0.2,
    family = "Arial",
    size = 2.25) + 
  annotate(geom="text", x=17.5, y=6.85, label="Kappa=0.61", color="black", size=3.5) +
  ggtitle("B. Burned stands variable importance") +
  theme(plot.title = element_text(hjust = 0.5, size=14)) +
  theme(legend.position="none", plot.margin=unit(c(.2,1,.2,1),"cm")) 
#burned_vi_plot

#Extract the confusion matrix values as data.frame
burned_cm_d <- as.data.frame(burned_confusion$table)
#Confusion matrix statistics as data.frame
burned_cm_st <-data.frame(burned_confusion$overall)
#Round the values
burned_cm_st$cm.overall <- round(burned_cm_st$burned_confusion.overall,3)
#Percentage
burned_cm_p <- as.data.frame(prop.table(burned_confusion$table))
burned_cm_d$Perc <- round(burned_cm_p$Freq*100,1)
burned_cm_d$Prediction <- ifelse(burned_cm_d$Prediction=="A", "Lived", "Died")
burned_cm_d$Reference <- ifelse(burned_cm_d$Reference=="A", "Lived", "Died")
burned_confusion
burned_cm_d

#Plot the matrix
burned_cm_plot <-  ggplot(data = burned_cm_d, aes(x = Prediction , y =  Reference, fill = Reference, alpha=Perc)) +
  geom_tile(color="black") +
  scale_fill_manual(name = "", values = c('Died' = 'red3', 'Lived' = 'darkolivegreen4')) +
  scale_alpha(range = c(.75, 1)) +
  geom_text(aes(label=paste0(Freq," ","(",Perc,"%",")")), color = 'black', size = 3.5) +
  theme_bw(base_size=11, base_family = "Arial") +
  theme(legend.position="none")
#burned_cm_plot

#Combine variable importance and confusion matrix plots
burned_vi_confusion <- ggdraw() +
  draw_plot(burned_vi_plot) +
  draw_plot(burned_cm_plot,
            height = 0.38,
            width = .65,
            x = 0.35,
            y = .0005
  )

#Combine the two final graphics
final_rf_model_graphic <- ggpubr::ggarrange(unburned_vi_confusion, burned_vi_confusion, ncol = 2, nrow = 1)
#final_rf_model_graphic
ggsave(path="/Users/james/Desktop", "final_rf_model_graphic.png", width = 9.15, height = 7, units = "in", dpi=500)


## Partial dependence plots

#Function to cretae partial dependence plots
build_pd_plot <- function(data, rf_model, var, xlab, xbreaks, xlim) {
  pd <- partial(rf_model, pred.var = var, prob = TRUE, rug = TRUE)
  ggplot() +
    geom_line(data = pd, aes(x = !!sym(var), y = yhat), color = "darkred", linewidth = 0.5) +
    geom_smooth(data = pd, aes(x = !!sym(var), y = yhat), color = "midnightblue", fill = "midnightblue", linewidth = 1) +
    geom_rug(data = data, aes(x = !!sym(var)), color = "grey10", alpha = 0.75, linewidth = 0.1) +
    scale_y_continuous("", breaks = c(0.2, 0.4, 0.6, 0.8), limits = c(0.1, 0.81)) +
    scale_x_continuous(xlab, breaks = xbreaks, limits = xlim) +
    theme_bw(base_size = 10, base_family = "Arial") +
    theme(plot.margin = margin(2, 2, 2, 2))
}

#Create partial dependence plots with specified variables
unburned_bai_pd_plot    <- build_pd_plot(unburned, unburned_rf, "bai5_scaled_mean",    "BAI (scaled)",     c(-1.5, 0, 1.5, 3),    c(-1.8, 3.5))
burned_bai_pd_plot      <- build_pd_plot(burned,   burned_rf,   "bai5_scaled_mean",    "BAI (scaled)",     c(-1.5, 0, 1.5, 3),    c(-1.8, 3.5))
unburned_height_pd_plot <- build_pd_plot(unburned, unburned_rf, "tree_ht_m",          "Tree height (m)",  c(20, 40),            c(0, 55))
burned_height_pd_plot   <- build_pd_plot(burned,   burned_rf,   "tree_ht_m",          "Tree height (m)",  c(20, 40),            c(0, 55))
unburned_age_pd_plot    <- build_pd_plot(unburned, unburned_rf, "age",                "Age (years)",      c(100, 300, 500),     c(16, 662))
burned_age_pd_plot      <- build_pd_plot(burned,   burned_rf,   "age",                "Age (years)",      c(100, 300, 500),     c(16, 662))
unburned_slope_pd_plot  <- build_pd_plot(unburned, unburned_rf, "slope",              "Slope (%)",        c(10, 20, 30),        c(1, 36))
burned_slope_pd_plot    <- build_pd_plot(burned,   burned_rf,   "slope",              "Slope (%)",        c(10, 20, 30),        c(1, 36))
unburned_cwd_pd_plot    <- build_pd_plot(unburned, unburned_rf, "cwd_1981_2010",      "CWD",              c(100, 150, 200),     c(93, 233))
burned_cwd_pd_plot      <- build_pd_plot(burned,   burned_rf,   "cwd_1981_2010",      "CWD",              c(100, 150, 200),     c(93, 233))

#Add column headers
unburned_bai_pd_plot <- unburned_bai_pd_plot + ggtitle("Unburned
stands") + theme(plot.title = element_text(hjust = 0.5, size = 11))
burned_bai_pd_plot   <- burned_bai_pd_plot   + ggtitle("Burned
stands")   + theme(plot.title = element_text(hjust = 0.5, size = 11))

#Remove redundant y-axis labels from right column
strip_y <- function(p) p + theme(axis.title.y = element_blank())
burned_bai_pd_plot      <- strip_y(burned_bai_pd_plot)
burned_height_pd_plot   <- strip_y(burned_height_pd_plot)
burned_age_pd_plot      <- strip_y(burned_age_pd_plot)
burned_slope_pd_plot    <- strip_y(burned_slope_pd_plot)
burned_cwd_pd_plot      <- strip_y(burned_cwd_pd_plot)

#Combine panels
row1 <- unburned_bai_pd_plot    + burned_bai_pd_plot
row2 <- unburned_height_pd_plot + burned_height_pd_plot
row3 <- unburned_age_pd_plot    + burned_age_pd_plot
row4 <- unburned_slope_pd_plot  + burned_slope_pd_plot
row5 <- unburned_cwd_pd_plot    + burned_cwd_pd_plot
main_plot <- row1 / row2 / row3 / row4 / row5 +
  plot_layout(ncol = 1)

#Save
png("/Users/james/Desktop/partial_dependence_final_ultratight.png", width = 3.72, height = 6.75, units = "in", res = 600)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(0.035, 0.965), "npc"))))
grid.text("Probability of tree mortality", x = unit(0.021, "npc"), rot = 90, gp = gpar(fontsize = 12, fontfamily = "Arial"))
print(main_plot, vp = viewport(layout.pos.col = 2))
dev.off()

## The interaction graphic (Figure 9)... this may take several minutes to run

#Compute partial dependence for unburned BAI and age
unburned_pd_bai_age <- partial(unburned_rf, pred.var = c("bai5_scaled_mean", "age"), prob = TRUE)
head(unburned_pd_bai_age)

#Interpolate the partial dependence values
unburned_bai_age_dens <- akima::interp(x = unburned_pd_bai_age$bai5_scaled_mean, y = unburned_pd_bai_age$age, z = unburned_pd_bai_age$yhat)
#Color palette (100 colors)
col.pal<-colorRampPalette(c("darkolivegreen", "red3"))
colors<-col.pal(100)
#Height of facets
z.facet.center <- (unburned_bai_age_dens$z[-1, -1] + unburned_bai_age_dens$z[-1, -ncol(unburned_bai_age_dens$z)] + unburned_bai_age_dens$z[-nrow(unburned_bai_age_dens$z), -1] + unburned_bai_age_dens$z[-nrow(unburned_bai_age_dens$z), -ncol(unburned_bai_age_dens$z)])/4
#Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)
#Plot BAI and age
quartz(h=7, w=12)
unburned_bai_age_persp <- persp(x = unburned_bai_age_dens$x, y = unburned_bai_age_dens$y, z = unburned_bai_age_dens$z, theta=45, phi=25, ticktype="detailed", xlab="\nBAI (scaled)", ylab="\nTree age (years)", zlab="\nProb. Mortality", expand=2/3, shade=0.25, cex.lab=1.35, cex.axis=1.25, col=colors[z.facet.range])
unburned_bai_age_persp


########################################
##### Simulate mortality over time #####
########################################

#Let's clear the workspace so we don't get conflicts with previous code
rm(list=ls())

#Reload the libraries we'll need
library(tidyverse)
library(ranger)

##Prepare data 

#Load focal tree data
focal_trees <- read.csv("focal_trees.csv")
head(focal_trees)

#Reduce to relocated trees
focal_trees <- focal_trees %>% filter(status=="A" | status=="D")
nrow(focal_trees)

#Create survived categorical column
focal_trees$survived <- ifelse(focal_trees$status=="A", 1, 0)

#Change column names and make site a factor
focal_trees <- focal_trees %>%
  rename("bai" = bai5_scaled_mean,
         "tree_ht" = tree_ht_m)
focal_trees$site2 <- as.factor(focal_trees$site2)

#Predict heights with GAM
library(mgcv)
gam_ht <- gam(tree_ht ~ 
                s(age) + 
                spp +
                pvt +
                s(site2, bs = 're'),                
              data=focal_trees)

#Reduce to relocated trees
mal_relocated <- focal_trees %>% filter(status=="A" | status=="D")

#Select best variables for random forest, i.e., the most predictive non-correlated (>0.65 or <-0.65) variables from Boruta.
all_trees <- mal_relocated %>% select(site2, tree_ht, bai, temp_mean_c, precip_mm, modis_npp, cwd_1981_2010, slope, basal15_m2ha, trees15_ha, forest_group, pvt, age, spp, asw_mm, survived)

#Any columns still have NAs?  
names(which(colSums(is.na(all_trees)) > 0))
#Filter to NA rows in a BAI column
has_nas <- all_trees[is.na(all_trees$bai50_total),]
has_nas$sample_id
#Drop the NA rows
all_trees <- all_trees %>% tidyr::drop_na()

#Remove sample ID, which we don't need anymore, and convert character strings to factors, and view file
all_trees$sample_id <- NULL
all_trees <- all_trees %>%
  dplyr::mutate_if(is.character, factor)
all_trees$survived <- as.factor(all_trees$survived)
str(all_trees)
head(all_trees)
nrow(all_trees)

#Change '81-'10 CWD to cwd (that's the one we'll use for the models)
colnames(all_trees)[colnames(all_trees) == 'cwd_1981_2010'] <- 'cwd'

#Filter to burned sites and unburned sites
burned_sites <- c("cny1", "cny2", "cny3", "cny5", "D3", "G3", "G4", "dca1", "dca2", "dca3", "rey1", "rey2", "rey3", "nma1", "nma2", "nma3")
unburned <- all_trees %>% filter(!site2 %in% burned_sites)
burned <- all_trees %>% filter(site2 %in% burned_sites)
table(unburned$site2)
table(burned$site2)

## Rerun random forest models... this time we'll use all the data instead of with-holding some data for testing (which we already did)

## First, the unburned model

#Get number of features (the number of predictors)
unburned_n_features <- length(setdiff(names(unburned), "survived"))

#Train a random forest model with default settings
unburned_rf_default <- ranger(
  survived ~ ., 
  data = unburned,
  mtry = floor(unburned_n_features / 3),
  respect.unordered.factors = "order",
  seed = 123
)

#Get OOB RMSE
unburned_default_rmse <- sqrt(unburned_rf_default$prediction.error)
unburned_default_rmse

#Create hyperparameter grid
unburned_hyper_grid <- expand.grid(
  mtry = floor(unburned_n_features * c(.05, .1, .15, .2, .25, .3, .333, .35, .4, .45)),
  min.node.size = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
  replace = c(TRUE, FALSE),                             
  sample.fraction = c(.1, .2, .3, .4, .5, .6, .7, .8, .9, 1 ),                       
  rmse = NA                                               
)

#Full cartesian grid search (last line exports OOB error)
for(i in seq_len(nrow(unburned_hyper_grid))) {
  unburned_fit <- ranger(
    formula         = survived ~ ., 
    data            = unburned, 
    num.trees       = unburned_n_features * 10,
    mtry            = unburned_hyper_grid$mtry[i],
    min.node.size   = unburned_hyper_grid$min.node.size[i],
    replace         = unburned_hyper_grid$replace[i],
    sample.fraction = unburned_hyper_grid$sample.fraction[i],
    verbose         = FALSE,
    seed            = 234,
    respect.unordered.factors = 'order',
  )
  unburned_hyper_grid$rmse[i] <- sqrt(unburned_fit$prediction.error)
}

#Have a look at the top 10 models
unburned_hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (unburned_default_rmse - rmse) / unburned_default_rmse * 100) %>%
  head(10)

#Re-run model with selected hyper-parameters
unburned_rf <- ranger(
  formula = survived ~ ., 
  data = unburned, 
  probability = TRUE,
  num.trees = 2000,
  mtry =  6,
  min.node.size = 3,
  sample.fraction = .6,
  replace = FALSE,
  importance = "impurity",
  #importance = "impurity_corrected",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed = 567
)

#Variable importance
unburned_vi <- data.frame(unburned_rf$variable.importance)
unburned_vi$variable <- rownames(unburned_vi)
colnames(unburned_vi) <- c("importance","variable")
rownames(unburned_vi) <- NULL
unburned_vi <- unburned_vi[order(-unburned_vi$importance),]

## Next, the burned model

#Get number of features (the number of predictors)
burned_n_features <- length(setdiff(names(burned), "survived"))

#Train a random forest model with default settings
burned_rf_default <- ranger(
  survived ~ ., 
  data = burned,
  mtry = floor(burned_n_features / 3),
  respect.unordered.factors = "order",
  seed = 234
)

#Get OOB RMSE
burned_default_rmse <- sqrt(burned_rf_default$prediction.error)
burned_default_rmse

#Create hyperparameter grid
burned_hyper_grid <- expand.grid(
  mtry = floor(burned_n_features * c(.05, .1, .15, .2, .25, .3, .333, .35, .4, .45)),
  min.node.size = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.1, .2, .3, .4, .5, .6, .7, .8, .9, 1 ),                       
  rmse = NA                                               
)

#Full cartesian grid search (last line exports OOB error)
for(i in seq_len(nrow(burned_hyper_grid))) {
  burned_fit <- ranger(
    formula         = survived ~ ., 
    data            = burned, 
    num.trees       = burned_n_features * 10,
    mtry            = burned_hyper_grid$mtry[i],
    min.node.size   = burned_hyper_grid$min.node.size[i],
    replace         = burned_hyper_grid$replace[i],
    sample.fraction = burned_hyper_grid$sample.fraction[i],
    verbose         = FALSE,
    seed            = 123,
    respect.unordered.factors = 'order',
  )
  burned_hyper_grid$rmse[i] <- sqrt(burned_fit$prediction.error)
}

#Have a look at the top 10 models
burned_hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (burned_default_rmse - rmse) / burned_default_rmse * 100) %>%
  head(10)

#Re-run model with selected hyper-parameters
burned_rf <- ranger(
  formula = survived ~ ., 
  data = burned, 
  probability = TRUE,
  num.trees = 2000,
  mtry =  3,
  min.node.size = 2,
  sample.fraction = .7,
  replace = FALSE,
  importance = "impurity",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 123
)

#Variable importance
burned_vi <- data.frame(burned_rf$variable.importance)
burned_vi$variable <- rownames(burned_vi)
colnames(burned_vi) <- c("importance","variable")
rownames(burned_vi) <- NULL
burned_vi <- burned_vi[order(-burned_vi$importance),]

#View variable importance
unburned_vi
burned_vi

## Run simulations

#Reduce to the objective sample of age structure and species composition, then filter to the columns we need for prediction (from the RF model)
mal_subplots <- mal_relocated %>% filter(objective=="Y")
mal_subplots <- mal_subplots %>% select(site2, sample_id, bai, tree_ht, age, slope, cwd_1981_2010, cwd_2041_2070, spp, temp_mean_c, precip_mm, trees15_ha, basal15_m2ha, modis_npp, pvt, asw_mm, forest_group)

#The primary/secondary simulation function (for applying mortality predictions across multiple decades)
primary_simulation_function <- function(data_for_sim, unburned_rf, burned_rf, prob_unburned, primary_cwd, secondary_cwd, bai_adj, start_decade, end_decade)
{
  initial_sim_function <- function(data_for_sim, unburned_rf, burned_rf) {
    primary_df <- data_for_sim 
    primary_df$cwd <-primary_df[[primary_cwd]]
    primary_df$site1 <- substr(primary_df$site2, 1, 3)
    site1 <- unique(primary_df$site1)
    if_burn <- sample(0:1, 10, replace=T, 
                      prob=c(prob_unburned, 1-prob_unburned))
    site1_names <- data.frame(cbind(site1, if_burn))
    primary_df <- primary_df %>% left_join(site1_names, 
                                           by="site1")
    for_unburned <- primary_df %>% filter(if_burn=="0")
    for_burned <- primary_df %>% filter(if_burn=="1")
    if (nrow(for_unburned) > 1) {
      unburned_predicts <- 
        cbind(for_unburned, 
              predict(unburned_rf, for_unburned))
    } else {
      unburned_predicts <- for_unburned
    }
    if (nrow(for_burned) > 1) {
      burned_predicts <- cbind(for_burned, 
                               predict(burned_rf, for_burned))
    } else {
      burned_predicts <- for_burned
    }
    initial_predicts <- 
      dplyr::bind_rows(unburned_predicts, burned_predicts)
    initial_predicts$rand_num <- 
      runif(nrow(initial_predicts), min = 0, max = 1)
    initial_predicts$mortality <- 
      ifelse(initial_predicts$rand_num > 
               initial_predicts$X1, "Died", "Survived") 
    initial_predicts <- initial_predicts %>% 
      filter(mortality=="Survived")
    initial_predicts$decade <- rep(start_decade, 
                                   nrow(initial_predicts))
    initial_predicts[ ,c('if_burn', 'X0', 'X1', 
                         'rand_num', 'mortality')] <- list(NULL)
    initial_predicts
  }
  secondary_sim_function <- function(data_for_sim, unburned_rf, burned_rf){
    max_decade <- max(data_for_sim$decade)
    this_sim <- data_for_sim %>% filter(decade==max_decade)
    this_sim$cwd <-this_sim[[secondary_cwd]]
    this_sim$est_ht1 <- 
      predict.gam(gam_ht, newdata = this_sim)
    this_sim$age <- this_sim$age+10
    this_sim$decade <- this_sim$decade+10
    this_sim$est_ht2 <- 
      predict.gam(gam_ht, newdata = this_sim)
    this_sim$ht_diff <- this_sim$est_ht2 - this_sim$est_ht1
    this_sim$ht_diff <- this_sim$ht_diff+.5
    this_sim$tree_ht <- this_sim$tree_ht + this_sim$ht_diff
    this_sim$bai <- ifelse(
      this_sim$bai >=0,
      this_sim$bai +
        (this_sim$bai*bai_adj), this_sim$bai + 
        (this_sim$bai*-bai_adj))
    site1 <- unique(this_sim$site1)
    if_burn <- sample(0:1, 10, replace=T, 
                      prob=c(prob_unburned, 1-prob_unburned))
    site1_names <- data.frame(cbind(site1, if_burn))
    for_rf <- this_sim %>% left_join(site1_names, by="site1")
    for_unburned <- for_rf %>% 
      filter(if_burn=="0")
    for_burned <- for_rf %>% 
      filter(if_burn=="1")
    if (nrow(for_unburned) > 1) {
      unburned_predicts <- cbind(for_unburned, 
                                 predict(unburned_rf, for_unburned))
    } else {
      unburned_predicts <- for_unburned
    }
    if (nrow(for_burned) > 1) {
      burned_predicts <- cbind(for_burned, 
                               predict(burned_rf, for_burned))
    } else {
      burned_predicts <- for_burned
    }
    secondary_predicts <- dplyr::bind_rows(
      unburned_predicts, burned_predicts)
    secondary_predicts$rand_num <- runif(nrow
                                         (secondary_predicts), min = 0, max = 1)
    secondary_predicts$mortality <- ifelse(secondary_predicts$rand_num > secondary_predicts$X1, "Died", "Survived") 
    secondary_predicts <- secondary_predicts %>% 
      filter(mortality=="Survived")
    secondary_predicts[ ,c('if_burn', 'X0', 'X1',
                           'rand_num', 'mortality', 'est_ht1', 'est_ht2', 
                           'ht_diff')] <- list(NULL)
    final_results <- rbind(data_for_sim, secondary_predicts)
    final_results
  }
  if(!"decade" %in% colnames(data_for_sim)) {
    result <- initial_sim_function(data_for_sim, unburned_rf, burned_rf)
  } 
  if(!end_decade %in% result$decade) {
    result <- secondary_sim_function(result, unburned_rf, burned_rf)
  } 
  repeat {
    result <- secondary_sim_function(result, unburned_rf, burned_rf)
    if(end_decade %in% result$decade) break
  }
  return(result)
} 

#The tertiary simulation function (for applying manually specified probabilities of mortality following initial runs with the primary/secondary simulation function)
tertiary_simulation_function <- function(data_for_sim, prob_young, prob_old, secondary_cwd, bai_adj, start_decade, end_decade)
{
  tertiary_sim_function <- function(data_for_sim, prob_young, prob_old){
    max_decade <- max(data_for_sim$decade)
    this_sim <- data_for_sim %>% filter(decade==max_decade)
    this_sim$cwd <-this_sim[[secondary_cwd]]
    this_sim$est_ht1 <- 
      predict.gam(gam_ht, newdata = this_sim)
    this_sim$age <- this_sim$age+10
    this_sim$decade <- this_sim$decade+10
    this_sim$est_ht2 <- 
      predict.gam(gam_ht, newdata = this_sim)
    this_sim$ht_diff <- this_sim$est_ht2 - this_sim$est_ht1
    this_sim$ht_diff <- this_sim$ht_diff+.5
    this_sim$tree_ht <- this_sim$tree_ht + this_sim$ht_diff
    this_sim$bai <- ifelse(
      this_sim$bai >=0,
      this_sim$bai +
        (this_sim$bai*bai_adj), this_sim$bai + 
        (this_sim$bai*-bai_adj))
    this_sim$X1 <- ifelse(this_sim$age < 150, .97, .99)
    this_sim$rand_num <- runif(nrow
                               (this_sim), min = 0, max = 1)
    this_sim$mortality <- ifelse(this_sim$rand_num > 
                                   this_sim$X1, "Died", "Survived") 
    this_sim <- this_sim %>% 
      filter(mortality=="Survived")
    this_sim[ ,c('X1', 'rand_num', 'mortality', 'est_ht1',
                 'est_ht2', 'ht_diff')] <- list(NULL)
    final_results <- rbind(data_for_sim, this_sim)
  }
  if(!end_decade %in% data_for_sim$decade) {
    result <- tertiary_sim_function(data_for_sim, prob_young, prob_old)
  } 
  repeat {
    result <- tertiary_sim_function(result, prob_young, prob_old)
    if(end_decade %in% result$decade) break
  }
  return(result)
} 

#Run simulations.  Warning:  Ten simulations or less (see the "reps" line) will take a minute or so  A thousand simulations could take an hour or so  Ten thousand simulations for multiple scenarios might take the better part of a day.

#Run simulation scenario #1:  Persistent trend:  Realistic CWD, same fire, same BAI
reps <- 1000
scenario1 <- map_df(1:reps, ~primary_simulation_function(data_for_sim=mal_subplots, unburned_rf=unburned_rf, burned_rf=burned_rf, prob_unburned=.63, primary_cwd="cwd_1981_2010", secondary_cwd="cwd_2041_2070", bai_adj=0, start_decade=2020, end_decade=2080), .id="sim_no") 
scenario1$scenario <- "Scenario 1"
head(scenario1)
nrow(scenario1)
table(scenario1$decade)

#Run simulation scenario #2:  Worsening trend:  Realistic CWD, more fire, same BAI
reps <- 1000
scenario2 <- map_df(1:reps, ~primary_simulation_function(data_for_sim=mal_subplots, unburned_rf=unburned_rf, burned_rf=burned_rf, prob_unburned=.33, primary_cwd="cwd_1981_2010", secondary_cwd="cwd_2041_2070", bai_adj=0, start_decade=2020, end_decade=2080), .id="sim_no") 
scenario2$scenario <- "Scenario 2"
head(scenario2)
nrow(scenario2)
table(scenario2$decade)

#Run simulation scenario #3:  Improving trend:  Realistic CWD, less fire, augmented BAI
reps <- 1000
scenario3 <- map_df(1:reps, ~primary_simulation_function(data_for_sim=mal_subplots, unburned_rf=unburned_rf, burned_rf=burned_rf, prob_unburned=.75, primary_cwd="cwd_1981_2010", secondary_cwd="cwd_2041_2070", bai_adj=.1, start_decade=2020, end_decade=2080), .id="sim_no") 
scenario3$scenario <- "Scenario 3"
head(scenario3)
nrow(scenario3)
table(scenario3$decade)

#Bind the first three scenarios
first_batch_scenarios <- do.call("rbind", list(scenario1, scenario2, scenario3))

#Run simulation scenario #4:  Disturbance restoration: Realistic CWD, more fire, same BAI for three decades, and then vastly improved mortality probabilities (run primary simulation function for three decades and then switch to tertiary simulation with mortality probabilities manually specified for four decades):  
reps <- 1000
scenario4 <- map_df(1:reps, ~primary_simulation_function(data_for_sim=mal_subplots, unburned_rf=unburned_rf, burned_rf=burned_rf, prob_unburned=.5, primary_cwd="cwd_1981_2010", secondary_cwd="cwd_2041_2070", bai_adj=.0, start_decade=2020, end_decade=2040), .id="sim_no") 
scenario4$scenario <- "Scenario 4"
head(scenario4)
nrow(scenario4)
table(scenario4$decade)

#Run simulation scenario #5:  Recovery from anomaly: Realistic CWD, same fire, same BAI for just one decade, and then vastly improved mortality probabilities (run primary simulation function for one decade and then switch to tertiary simulation with mortality probabilities manually specified for seven decades):  
reps <- 1000
scenario5 <- map_df(1:reps, ~primary_simulation_function(data_for_sim=mal_subplots, unburned_rf=unburned_rf, burned_rf=burned_rf, prob_unburned=.63, primary_cwd="cwd_1981_2010", secondary_cwd="cwd_2041_2070", bai_adj=0, start_decade=2020, end_decade=2020), .id="sim_no") 
scenario5$scenario <- "Scenario 5"
scenario5 <- scenario5  %>% filter(decade=="2020")
head(scenario5)
nrow(scenario5)
table(scenario5$decade)

#Feed the results of scenario 4 tertiary simulation function (manually specified probabilities)
reps <- 1
scenario4$scenario <- NULL
scenario4.1 <- map_df(1:reps, ~tertiary_simulation_function(data_for_sim=scenario4, prob_young=.85, prob_old=.90, secondary_cwd="cwd_2041_2070", bai_adj=.1, start_decade=2050, end_decade=2080)) 
scenario4.1$scenario <- "Scenario 4"
head(scenario4.1)
nrow(scenario4.1)
table(scenario4.1$decade)

#Feed the results of the truncated simulation to the tertiary simulation function (manually specified probabilities)
reps <- 1
scenario5$scenario <- NULL
scenario5.1 <- map_df(1:reps, ~tertiary_simulation_function(data_for_sim=scenario5, prob_young=.85, prob_old=.90, secondary_cwd="cwd_2041_2070", bai_adj=0, start_decade=2030, end_decade=2080)) 
scenario5.1$scenario <- "Scenario 5"
head(scenario5.1)
nrow(scenario5.1)
table(scenario5.1$decade)

#Bind the first batch of scenarios to the second batch of scenarios
all_scenarios <- do.call("rbind", list(first_batch_scenarios, scenario4.1, scenario5.1))

#Add age class columns to simulated data
all_scenarios$age_class1 <- ifelse(all_scenarios$age <150, "<150", ifelse(all_scenarios$age >=150 & all_scenarios$age <250, "150-250", ifelse(all_scenarios$age >=250 & all_scenarios$age <350, "250-350", ifelse(all_scenarios$age >=350 & all_scenarios$age <450, "350-450", "450+"))))
all_scenarios$age_class2 <- ifelse(all_scenarios$age <150, "<150", ifelse(all_scenarios$age >=150 & all_scenarios$age <300, "150-300", "300+"))
all_scenarios$age_class3 <- ifelse(all_scenarios$age <150, "young", "old")

#Add age class columns to real data and decade column
mal_subplots$age_class1 <- ifelse(mal_subplots$age <150, "<150", ifelse(mal_subplots$age >=150 & mal_subplots$age <250, "150-250", ifelse(mal_subplots$age >=250 & mal_subplots$age <350, "250-350", ifelse(mal_subplots$age >=350 & mal_subplots$age <450, "350-450", "450+"))))
mal_subplots$age_class2 <- ifelse(mal_subplots$age <150, "<150", ifelse(mal_subplots$age >=150 & mal_subplots$age <300, "150-300", "300+"))
mal_subplots$age_class3 <- ifelse(mal_subplots$age <150, "young", "old")
mal_subplots$decade <- 2010

#Counts of trees by simulation run, age class, and decade, then filter out young trees and add type ("simulations") column
age_counts_sim <- all_scenarios %>%
  group_by(sim_no, decade, age_class2, scenario) %>%
  tally() %>%
  arrange(sim_no)
age_counts_sim <- age_counts_sim %>% filter(!age_class2=="<150")
age_counts_sim$type <- "simulated"
age_counts_sim

#Counts of trees for real data by age class, and decade, then filter out young trees and add type ("real") and ("dummy") sim_no column, then add scenarios (just a duplicate)
age_counts_real <- mal_subplots %>%
  group_by(decade, age_class2) %>%
  tally() 
age_counts_real <- age_counts_real %>% filter(!age_class2=="<150")
age_counts_real$type <- "real"
age_counts_real$sim_no <- "0"
scenario <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5")
age_counts_real <- crossing(age_counts_real, scenario)
age_counts_real

#Bind real and simulated
age_counts_all <- rbind(age_counts_real, age_counts_sim)
nrow(age_counts_all)
head(age_counts_all,10)

#Recode the decade
age_counts_all <- age_counts_all %>% mutate(decade=recode(decade, '2010'='Decade 0', '2020'='Decade 1', '2030'='Decade 2', '2040'='Decade 3', '2050'='Decade 4', '2060'='Decade 5', '2070'='Decade 6', '2080'='Decade 7'))

#Create columns that put tree counts in terms of proportion of the baseline
old_baseline <- age_counts_all %>% 
  filter(age_class2=="150-300" & decade=="Decade 0") %>% 
  select(n) 
very_old_baseline <- age_counts_all %>% 
  filter(age_class2=="300+" & decade=="Decade 0") %>% 
  select(n) 
age_counts_all$prop <- ifelse(age_counts_all$age_class2=="150-300", age_counts_all$n/mean(old_baseline$n), age_counts_all$n/mean(very_old_baseline$n))
print(age_counts_all, n=20)

#Split into different age class dataframes 
age_counts_old <- age_counts_all %>% filter(age_class2=="150-300")
age_counts_very_old <- age_counts_all %>% filter(age_class2=="300+")

#Mean number of trees per decade and scenario
library(rcompanion)
old_means <- groupwiseMean(n ~ type + decade + scenario,
                           data   = age_counts_old,
                           conf   = 0.95,
                           digits = 3)
very_old_means <- groupwiseMean(n ~ type + decade + scenario,
                                data   = age_counts_very_old,
                                conf   = 0.95,
                                digits = 3)

#Summary of simulation tree counts
decade0_old_means <- old_means %>% filter(decade=="Decade 0") %>% select(type, decade, scenario, Mean, Trad.lower, Trad.upper)
decade0_old_means
decade1_old_means <- old_means %>% filter(decade=="Decade 1") %>% select(type, decade, scenario, Mean, Trad.lower, Trad.upper)
decade1_old_means
decade7_old_means <- old_means %>% filter(decade=="Decade 7") %>% select(type, decade, scenario, Mean, Trad.lower, Trad.upper)
decade7_old_means
decade0_very_old_means <- very_old_means %>% filter(decade=="Decade 0") %>% select(type, decade, scenario, Mean, Trad.lower, Trad.upper)
decade0_very_old_means
decade1_very_old_means <- very_old_means %>% filter(decade=="Decade 1") %>% select(type, decade, scenario, Mean, Trad.lower, Trad.upper)
decade1_very_old_means
decade7_very_old_means <- very_old_means %>% filter(decade=="Decade 7") %>% select(type, decade, scenario, Mean, Trad.lower, Trad.upper)
decade1_very_old_means

#Graph old trees scenarios (proportions)
old_graph <- ggplot(age_counts_old, aes(x = decade, y = prop, color = scenario)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, jitter.height=.1, dodge.width = 0.8), shape=21, alpha=.075, size=.1) +
  geom_boxplot(fill=NA, position = position_dodge(width = 0.8), outlier.shape = NA, outlier.size = 0, width=.7) +
  geom_hline(yintercept = 0.676, linetype="dashed", linewidth=.25, alpha=.5) +
  scale_x_discrete("") +
  scale_y_continuous("Proportion reference (Decade 0) tree counts", breaks = c(0, 1, 2, 3), limits=c(0, 3.45)) +
  scale_color_manual("", values = c('Scenario 1' = 'gold4', 'Scenario 2' = 'orangered4', 'Scenario 3' = 'darkmagenta', 'Scenario 4' = 'royalblue4', 'Scenario 5' = 'darkgreen'), labels = c("Peristent Trend", "Worsening Trend", "Improving Trend", "Disturbance Restoration", "Recovery from Anomaly")) +
  ggtitle("A. 150-300 year old trees") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw(base_size=15, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=45, hjust=.9, vjust=.85)) +
  theme(axis.title.x=element_blank()) +
  theme(strip.text.y.left = element_text(angle = 90), panel.spacing.y = unit(0.25, "lines"), strip.background = element_blank(), legend.background = element_rect(fill = "white", colour = "white")) +
  theme(plot.title = element_text(size = 15, hjust = 0.02, vjust=10))
#old_graph 

#Graph very old trees scenarios
very_old_graph <- ggplot(age_counts_very_old, aes(x = decade, y = prop, color = scenario)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, jitter.height=.1, dodge.width = 0.8), shape=21, alpha=.075, size=.1) +
  geom_boxplot(fill=NA, position = position_dodge(width = 0.8), outlier.shape = NA, outlier.size = 0, width=.7) +
  geom_hline(yintercept = 0.746, linetype="dashed", linewidth=.25, alpha=.5) +
  scale_x_discrete("") +
  scale_y_continuous("Proportion reference (Decade 0) tree counts", breaks = c(0, 1), limits=c(0, 1.3)) +
  scale_color_manual("", values = c('Scenario 1' = 'gold4', 'Scenario 2' = 'orangered4', 'Scenario 3' = 'darkmagenta', 'Scenario 4' = 'royalblue4', 'Scenario 5' = 'darkgreen'), labels = c("Peristent Trend", "Worsening Trend", "Improving Trend", "Disturbance Restoration", "Recovery from Anomaly")) +
  ggtitle("B. 300+ year old trees") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw(base_size=15, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=45, hjust=.9, vjust=.85)) +
  theme(strip.text.y.left = element_text(angle = 90), panel.spacing.y = unit(0.25, "lines"), strip.background = element_blank(), legend.background = element_rect(fill = "white", colour = "white")) +
  theme(plot.title = element_text(size = 15, hjust = 0.02, vjust=10))
#very_old_graph

#Combine
library(ggpubr)
quartz(width=7, height=8)
old_counts_plot <- ggarrange(old_graph + rremove("ylab"), very_old_graph+ rremove("ylab"), heights = c(1, 1), ncol = 1, nrow = 2, align = "v", common.legend = TRUE, legend="bottom")
#Annotate the figure by adding a y-axis common label
annotate_figure(old_counts_plot, left = text_grob("Proportion reference (Decade 0) tree counts", size = 13, family="Arial", rot = 90, vjust = .25, hjust = .3))

## Evaluate compositional shifts across simulations 

#Species composition of real data in Decade 0
real_comp <- focal_trees %>% filter(objective=="Y")
real_comp$age_class2 <- ifelse(real_comp$age <150, "<150", ifelse(real_comp$age >=150 & real_comp$age <300, "150-300", "300+"))
real_comp$decade <- 2010
real_comp <- real_comp %>% select(spp, age_class2, decade, status)
unique(real_comp$decade)
real_comp_decade0 <- real_comp %>% 
  group_by(age_class2, decade) %>% 
  count(spp) %>% 
  mutate(prop = n / sum(n))
real_comp_decade0

#Species composition of real data in Decade 1
real_comp_alive <- real_comp %>% filter(status=="A")
real_comp_alive$decade <- 2010
real_comp_decade1 <- real_comp_alive %>% 
  group_by(age_class2, decade) %>% 
  count(spp) %>% 
  mutate(prop = n / sum(n))
real_comp_decade1

#Counts of species by simulation number, scenario, decade, and age class
sims_comp <- all_scenarios %>% select(spp, age_class2, decade, sim_no, scenario)
unique(sims_comp$decade)
unique(sims_comp$scenario)
#unique(sims_comp$sim_no)
sims_comp <- sims_comp %>% 
  group_by(sim_no, decade, scenario, age_class2) %>% 
  count(spp) %>% 
  mutate(prop = n / sum(n)) %>% 
  arrange(decade, age_class2, spp)
print(sims_comp, n=14)

#Filter to a scenario and 150-300 and >300 trees
#Scenario 1 = Peristent Trend, Scenario 2 = Worsening Trend, Scenario 3 = Improving Trend, Scenario 4 = Disturbance Restoration, Scenario 5 = Recovery from Anomaly) +
sims_comp_scen <- sims_comp %>% filter(scenario=="Scenario 1")
sims_comp_scen_old <- sims_comp_scen %>% filter(!age_class2=="<150")

#Bind real and simulated data
real_comp_decade0$sim_no <- 0
real_comp_decade0$scenario <- "Real"
real_comp_decade0 <- real_comp_decade0 %>% select(scenario, decade, sim_no, age_class2, spp, prop)
sims_comp_scen_old <- sims_comp_scen_old %>% select(scenario, decade, sim_no, age_class2, spp, prop)
sims_comp_scen_old$sim_no <- as.numeric(sims_comp_scen_old$sim_no)
sims_comp_scen_old$sim_no <- as.numeric(sims_comp_scen_old$sim_no)
all_comp <- rbind(real_comp_decade0, sims_comp_scen_old)
all_comp

#Recode the decade
all_comp <- all_comp %>% mutate(decade=recode(decade, '2010'='Decade 0', '2020'='Decade 1', '2030'='Decade 2', '2040'='Decade 3', '2050'='Decade 4', '2060'='Decade 5', '2070'='Decade 6', '2080'='Decade 7'))

#Split into different age class dataframes 
all_comp_old <- all_comp %>% filter(age_class2=="150-300")
all_comp_very_old <- all_comp %>% filter(age_class2=="300+")

#Actual composition in Decade 1
real_comp_decade1

#Average simulated composition by decade (old)
mean_all_comp_old <- all_comp_old %>% 
  group_by(decade, spp) %>% 
  summarize(mean_prop=mean(prop))
data.frame(mean_all_comp_old)

#Average simulated composition by decade (very old)
mean_all_comp_very_old <- all_comp_very_old %>% 
  group_by(decade, spp) %>% 
  summarize(mean_prop=mean(prop))
data.frame(mean_all_comp_very_old)

#Graph composition (proportions) of 150-300
#Reorder species factor
all_comp_old$spp <- factor(all_comp_old$spp, levels = c('PIPO', 'ABGR', 'PSME', 'LAOC', 'JUOC', 'PICO'))
old_comp_graph <- ggplot(all_comp_old, aes(x = decade, y = prop, color = spp)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.4, jitter.height=.1, dodge.width = 0.8), shape=21, alpha=.075, size=.1) +
  geom_boxplot(fill=NA, position = position_dodge(width = 0.8), outlier.shape = NA, outlier.size = 0, width=.7) +
  geom_hline(yintercept = 0.558, linetype="dashed", size=.25, alpha=.5, color="darkgreen") +
  geom_hline(yintercept = 0.0233, linetype="dashed", size=.25, alpha=.5, color="cadetblue") +
  geom_hline(yintercept = 0.0229, linetype="dashed", size=.25, alpha=.5, color="gray44") +
  geom_hline(yintercept = 0.349, linetype="dashed", size=.25, alpha=.5, color="chocolate") +
  geom_hline(yintercept = 0.0465, linetype="dashed", size=.25, alpha=.5, color="darkolivegreen3") +
  scale_x_discrete("") +
  scale_y_continuous("Proportion reference (Decade 0) tree composition", breaks = c(0, .5, 1), limits=c(0, 1)) +
  scale_color_manual(name = "", values = c('PIPO' = 'chocolate', 'ABGR' = 'darkgreen', 'PSME' = 'darkolivegreen3', 'LAOC' = 'yellow3', 'JUOC' = 'cadetblue', 'PICO' = 'gray44'), labels = c("Ponderosa", "Grand fir", "Douglas-fir", "Larch", "Juniper", "Lodgepole")) +
  ggtitle("A. 150-300 year old trees") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw(base_size=15, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=45, hjust=.9, vjust=.85)) +
  theme(axis.title.x=element_blank()) +
  theme(strip.text.y.left = element_text(angle = 90), panel.spacing.y = unit(0.25, "lines"), strip.background = element_blank(), legend.background = element_rect(fill = "white", colour = "white")) +
  theme(plot.title = element_text(size = 15, hjust = 0.02, vjust=10))
#old_graph 

#Graph composition (proportions) of 300+
all_comp_very_old$spp <- factor(all_comp_very_old$spp, levels = c('PIPO', 'ABGR', 'PSME', 'LAOC', 'JUOC', 'PICO'))
very_old_comp_graph <- ggplot(all_comp_very_old, aes(x = decade, y = prop, color = spp)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height=.1, dodge.width = 0.8), shape=21, alpha=.075, size=.1) +
  geom_boxplot(fill=NA, position = position_dodge(width = 0.8), outlier.shape = NA, outlier.size = 0, width=.7) +
  geom_hline(yintercept = 0.0357, linetype="dashed", size=.25, alpha=.5, color="yellow3") +
  geom_hline(yintercept = 0.893, linetype="dashed", size=.25, alpha=.5, color="chocolate") +
  geom_hline(yintercept = 0.0714, linetype="dashed", size=.25, alpha=.5, color="darkolivegreen3") +
  scale_x_discrete("") +
  scale_y_continuous("Proportion reference (Decade 0) tree composition", breaks = c(0, .5, 1), limits=c(0, 1)) +
  scale_color_manual(name = "", values = c('PIPO' = 'chocolate', 'ABGR' = 'darkgreen', 'PSME' = 'darkolivegreen3', 'LAOC' = 'yellow3', 'JUOC' = 'cadetblue', 'PICO' = 'gray44'), labels = c("Ponderosa", "Grand fir", "Douglas-fir", "Larch", "Juniper", "Lodgepole")) +
  ggtitle("B. 300+ year old trees") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw(base_size=15, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=45, hjust=.9, vjust=.85)) +
  theme(strip.text.y.left = element_text(angle = 90), panel.spacing.y = unit(0.25, "lines"), strip.background = element_blank(), legend.background = element_rect(fill = "white", colour = "white")) +
  theme(plot.title = element_text(size = 15, hjust = 0.02, vjust=10))
#very_old_graph

#Combine
library(ggpubr)
quartz(width=7, height=8)
figure <- ggarrange(old_comp_graph + rremove("ylab"), very_old_comp_graph + rremove("ylab"), heights = c(1, 1), ncol = 1, nrow = 2, align = "v", common.legend = TRUE, legend="bottom")
#Annotate the figure by adding a y-axis common label
annotate_figure(figure, left = text_grob("Proportion reference (Decade 0) tree counts", size = 13, family="Arial", rot = 90, vjust = .25, hjust = .315))
#ggsave(path="/Users/james/Desktop", "combine_graphs.png", width = 7, height = 8, units = "in", dpi=500)


###########################################
#### Make Figure 12 to evalute climate ####
########## and pine butterfly #############
###########################################

#Clear the workspace
rm(list=ls())

#Reload library
library(tidyverse)

#Read PRISM climate data
prism_climate <- read.csv("prism_climate.csv")

#Total area defoliated by pine butterfly from Forest Service ADS data (available from:  https://www.fs.usda.gov/science-technology/data-tools-products/fhp-mapping-reporting/detection-surveys)
pine_butterfly_area_by_year <- structure(list(
  year = c(2009L, 2010L, 2011L, 2012L),
  total_area_ha = c(1618, 9717, 101022, 36852)
), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -4L))

#Reduce to March through August
climate_growing <- prism_climate %>%
  filter(month %in% 3:9)

#Means across site by year and month
climate_summary <- climate_growing %>%
  group_by(year, month) %>%
  summarise(
    mean_ppt = mean(ppt, na.rm = TRUE),
    mean_tmax = mean(tmax, na.rm = TRUE),
    mean_maxvpd = mean(vpdmax, na.rm = TRUE),
    .groups = "drop"
  )
climate_summary
nrow(climate_summary)

##Bar chart of VPD and precipitation anomalies

#Calculate long-term monthly means
monthly_means <- climate_summary %>%
  group_by(month) %>%
  summarise(
    mean_tmax_month = mean(mean_tmax, na.rm = TRUE),
    mean_ppt_month = mean(mean_ppt, na.rm = TRUE),
    mean_maxvpd_month = mean(mean_maxvpd, na.rm = TRUE),
    .groups = "drop"
  )

#Join and compute anomalies
climate_anomaly <- climate_summary %>%
  left_join(monthly_means, by = "month") %>%
  mutate(
    deviation_tmax = mean_tmax - mean_tmax_month,
    deviation_ppt = mean_ppt - mean_ppt_month,
    deviation_maxvpd = mean_maxvpd - mean_maxvpd_month
  )

#Filter months for plotting
plot_data <- climate_anomaly %>%
  filter(month %in% 6:8, year >= 1973, year <= 2023) %>%
  mutate(
    month_label = factor(month, levels = 6:8, labels = month.name[6:8])
  )

#Expand across months in plot_data
months_used <- unique(plot_data$month_label)
severity_overlay <- expand.grid(
  year = pine_butterfly_area_by_year$year,
  month_label = months_used
) %>%
  left_join(pine_butterfly_area_by_year, by = "year") %>%
  mutate(
    xmin = year - 0.5,
    xmax = year + 0.5,
    ymin = -Inf,
    ymax = Inf
  )

#Plot
vpd_anomaly_graph <- ggplot(plot_data, aes(x = year)) +
  geom_rect(
    data = severity_overlay,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = total_area_ha),
    fill = "brown4",
    inherit.aes = FALSE
  ) +
  scale_alpha(
    range = c(0.05, 0.3),
    name = "Total Defoliation (ha)"
  ) +
  geom_col(aes(y = deviation_maxvpd, fill = deviation_maxvpd > 0), width = 0.95, show.legend = FALSE) +
  geom_line(aes(y = deviation_ppt / 5), color = "blue4", linewidth = 0.45) +
  #geom_smooth(aes(y = deviation_ppt / 5), color = "blue4", linewidth = 0.25, se = FALSE, span = 0.2) +
  geom_point(aes(y = deviation_ppt / 5), color = "blue4", size = 0.25) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
  geom_vline(xintercept = c(2012, 2023), linetype = "dashed", color = "black", linewidth = 0.6) +
  facet_wrap(~ month_label, ncol = 1) +
  scale_fill_manual(
    values = c("TRUE" = "darkgoldenrod3", "FALSE" = "cornflowerblue")
  ) +
  scale_x_continuous(
    breaks = seq(1980, 2020, by = 5),
    limits = c(1980, 2024),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(-15.2, 9),
    name = "Max. VPD Anomaly (hPA)",
    sec.axis = sec_axis(~ . * 5, name = "PPT Anomaly (mm)\n")
  ) +
  labs(x = NULL) +
  theme_minimal(base_family = "Arial", base_size = 14) +
  theme(
    strip.text = element_text(face = "plain"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3),
    panel.grid.minor = element_line(linewidth = 0.15),
    legend.position = "bottom"  # <--- Move alpha legend to bottom
  )
ggsave(
  filename = "~/Desktop/vpd_anomaly_graph.png",
  width = 8,
  height = 7.6,
  dpi = 500,
  units = "in",
  bg = "white"
)

#Filter baseline period (1924–2011) for June–August and compute monthly climatology
monthly_means_jja <- prism_climate %>%
  filter(year >= 1924, year <= 2011, month %in% 6:8) %>%
  group_by(month) %>%
  summarise(
    mean_vpd = mean(vpdmax, na.rm = TRUE),
    sd_vpd   = sd(vpdmax, na.rm = TRUE),
    mean_ppt = mean(ppt, na.rm = TRUE),
    sd_ppt   = sd(ppt, na.rm = TRUE),
    .groups = "drop"
  )

#Join monthly means to full dataset and compute anomalies
climate_anomalies_jja <- prism_climate %>%
  filter(month %in% 6:8) %>%  # Only JJA months for all years
  left_join(monthly_means_jja, by = "month") %>%
  mutate(
    anomaly_vpd = vpdmax - mean_vpd,
    anomaly_ppt = ppt - mean_ppt,
    z_vpd = anomaly_vpd / sd_vpd,
    z_ppt = anomaly_ppt / sd_ppt
  )

#Summarize 2012–2023 JJA anomalies
anomaly_summary_jja <- climate_anomalies_jja %>%
  filter(year >= 2012, year <= 2023) %>%
  summarise(
    mean_vpd_anomaly     = mean(anomaly_vpd, na.rm = TRUE),
    mean_ppt_anomaly     = mean(anomaly_ppt, na.rm = TRUE),
    mean_vpd_zscore      = mean(z_vpd, na.rm = TRUE),
    mean_ppt_zscore      = mean(z_ppt, na.rm = TRUE),
    max_vpd_zscore       = max(z_vpd, na.rm = TRUE),
    min_ppt_zscore       = min(z_ppt, na.rm = TRUE),
    n_records            = n()
  )
data.frame(anomaly_summary_jja)

