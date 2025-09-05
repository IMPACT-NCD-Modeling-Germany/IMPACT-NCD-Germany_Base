
##################################################################################################################
#----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------    Semaglutide Modelling    ------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#-----------------------------------    Scenario 1: Individual-level CEA     ------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
##################################################################################################################
scenario_0_fn <- function(sp) {

  sp$pop[, c("bmi_delta", "sbp_delta", "tchol_delta") := 0]

}

#scenario_0_fn(sp)

### Scenario 1 - Cost-effectiveness analysis of semaglutide 2.4 mg for weight management in German adult population

scenario_1_fn <- function(sp) {

  # Set scenario variables #
  ###########################################################################################################
  # Semag scenario variables:
  # 1. Baseline year: 2024
  # 2. Intervention year: 2025 (years of drug roll-out)
  # 3. simulation horizon: 20 years (the end of simulation could be 2043)
  # 4. Eligibility criteria: BMI >= 35
  # 5. Treatment uptake_rates:
  #    a) For CEA: 100 % uptake for the eligible patients at the year of drug roll-out (2025)
  #    b) For Budget impact: time-varying
  # 6. Efficacy input trajectory (BMI, SBP, total cholesterol) (*uncertainty)
  ###########################################################################################################
  # Efficacy for bmi at 1st year after treatment
  tbl <- read_fst("./inputs/other_parameters/bmi_1y_samples.fst", as.data.table = TRUE)
  bmi_1y <- as.numeric(tbl[mc == sp$mc_aggr, "bmi_1y"])

  # Efficacy for bmi at 2nd year after treatment
  tbl <- read_fst("./inputs/other_parameters/bmi_2y_samples.fst", as.data.table = TRUE)
  bmi_2y <- as.numeric(tbl[mc == sp$mc_aggr, "bmi_2y"])

  # Efficacy for bmi at 3rd year after treatment
  tbl <- read_fst("./inputs/other_parameters/bmi_3y_samples.fst", as.data.table = TRUE)
  bmi_3y <- as.numeric(tbl[mc == sp$mc_aggr, "bmi_3y"])

  # Efficacy for bmi at 4th year after treatment
  tbl <- read_fst("./inputs/other_parameters/bmi_4y_samples.fst", as.data.table = TRUE)
  bmi_4y <- as.numeric(tbl[mc == sp$mc_aggr, "bmi_4y"])

  # Efficacy for SBP at 1st year after treatment
  tbl <- read_fst("./inputs/other_parameters/sbp_1y_samples.fst", as.data.table = TRUE)
  sbp_1y <- as.numeric(tbl[mc == sp$mc_aggr, "sbp_1y"])

  # Efficacy for SBP at 2nd year of treatment
  tbl <- read_fst("./inputs/other_parameters/sbp_2y_samples.fst", as.data.table = TRUE)
  sbp_2y <- as.numeric(tbl[mc == sp$mc_aggr, "sbp_2y"])

  # Efficacy for total cholesterol at 1st year of treatment
  tbl <- read_fst("./inputs/other_parameters/tchol_1y_samples.fst", as.data.table = TRUE)
  tchol_1y <- as.numeric(tbl[mc == sp$mc_aggr, "tchol_1y"])

  # Efficacy for total cholesterol at 2nd year of treatment
  tbl <- read_fst("./inputs/other_parameters/tchol_2y_samples.fst", as.data.table = TRUE)
  tchol_2y <- as.numeric(tbl[mc == sp$mc_aggr, "tchol_2y"])
  ###########################################################################################################

  # a variable for how long it has been since the baseline (2024)
  baseline_yr <- 24 ## it's 2024, but all years should be relative to year 2000
  sp$pop[, ys_rollout := year - baseline_yr] # baseline_yr == 2024 / year == 24 --> baseline year

  # a variable for uptake rates based on years after drug roll-out
  # * For CEA: 100 % uptake for the eligible patients at the first year of drug roll-out (year == 25 --> rollout year)
  # * This uptake rate is independent of individual simulation history, it only depends on years after roll-out
  sp$pop[, uptake_rate := ifelse(ys_rollout == 1, 1, 0)] # Set the uptake_rate to be 100% and only possible in year 2025 in CEA

  # a variable for defining individual's eligibility for accessing Semaglutide by BMI
  sp$pop[, eligible_bi := ifelse(bmi_curr_xps >= 35 & age <= 80, 1, 0)]

  # a variable for the uptake of the drug conditioning on meeting eligibility (person_year)
  sp$pop[, uptake_psyr := ifelse(eligible_bi ==1, rbinom(.N, size = 1, prob = uptake_rate), 0)]

  # a variable for the once-in-a-lifetime uptake during an individual's simulation years
  sp$pop[, uptake_one := 0]
  sp$pop[uptake_psyr == 1, uptake_one := as.integer(.I == .I[which.min(year)]), by = pid]
  ### In this testing, only two ids got the treatment: pid 156, 175

  # a variable for the year someone uptaking the drug
  sp$pop[, anchor_year:=year[uptake_one == 1], by = pid]

  # a variable for the year someone enters the synthetic population
  sp$pop[, entry_year:= min(year), by = pid]

  # a variable for theorectical treatment trajectory
  ### Step 0: Initialize the column
  sp$pop[, trtm_theo := NA_integer_]
  ### Step 1: Identify first year of treatment uptake
  sp$pop[uptake_one == 1, trtm_theo := fifelse(anchor_year > entry_year, 1,0)]
  ### Step 2: Fill the sequence in subsequent years
  sp$pop[, trtm_theo := fifelse(cumsum(!is.na(trtm_theo)) > 0,   # cumsum(): computes a cumulative count of non-NA entries:
                                                                 # fill in a sequence starting from the first non-NA value
                               seq_len(.N) - min(which(!is.na(trtm_theo))) + fifelse(anchor_year > entry_year, 1L, 0L),
                               NA_integer_),
        by = pid]
  ### Step 3: Replace NAs with 0
  sp$pop[, trtm_theo := fifelse(is.na(trtm_theo), 0, trtm_theo)]

  #scenario_1_fn(sp) # Works!
  #test_uptake <- sp$pop[pid %in% sp$pop[entry_year >= 25 & eligible_bi ==1, unique(pid)]]
  #test_uptake <- sp$pop[pid %in% sp$pop[uptake_one == 1, unique(pid)]]

  ########################################################################################
  ################################          BMI         ##################################
  ########################################################################################
  # a variable for the change in BMI after the uptake of semaglutide
  sp$pop[, `:=`(
    bmi_delta = fcase(
      ### All values until year 5, are % weight loss relative to baseline
      trtm_theo == 0, 0,       ### Baseline year, the year patients become eligible
      trtm_theo == 1, bmi_1y,  ### replace '-0.149' with 'bmi_1y', and link to the uncertainty file: mc + bmi_1y, (Jane today)
                               ### 1st year of treatment: ITT effects from STEP 1 (1 year f/u): 14.9% weight loss
      trtm_theo == 2, bmi_2y,  ### replace '-0.152' with 'bmi_2y', and link to the uncertainty file: mc + bmi_2y, (Jane today)
                               ### 2nd year of treatment: ITT effects from STEP 5 (2 year f/u): 15.2% weight loss
      trtm_theo == 3, bmi_3y,  ### STEP 1 extension: regain 2/3 of weight loss in 1st year treatment cessation
                               ###
      trtm_theo == 4, bmi_4y,  ### weight regain at 2nd year treatment cessation (uncertainty?)
      trtm_theo >= 5, 0        ### weight return to baseline weight by 3rd year treatment cessation (uncertainty?)
      ### The efficacy inputs are assigned until here, the rest we will revert to the original exposure column (bmi_curr_xps)
    ))
    ]

  # a variable for each individual's baseline weight (baseline: year == 24 / ys_rollout == 0)
  # ------------------------------------------ Aug 22, 2025, Jane ------------------------------------------------- #
  # *but, some individuals in the dataset does not have entry at baseline_yr, they came in later than year 24
  # *and, some individuals came in later than year 24, and because they have a chance to uptake drug begining year 25
  # *they might already uptake the drug at the first year, which means they won't have a year (anchor-1),
  # *so what shall we use as a baseline BMI value here...
  sp$pop[, baseline_bmi := {              # sp$pop[, ... , by = pid]: grouping the operation by pid
    anchor <- unique(anchor_year)         # So we can distinguish people who have entry at year 24 and who do not
    if (is.na(anchor)) {                  # If no anchor is found,
      bmi_curr_xps                        # the function simply returns the original bmi_curr_xps as-is
    } else {                              # If yes anchor is found
      bmi_curr_xps[year == (year[trtm_theo == 1][1] - 1)] # baseline_bmi will take the BMI at the year before trtm_theo==1
    }
  }, by = pid] #Works!

  ###############################################################################################################
  #----------------------------  Tutorial: if_else function inside data.table  ---------------------------------#
  ###############################################################################################################
  sp$pop[, bmi_shift := {             # sp$pop[, ... , by = pid]: grouping the operation by pid
                                      # { ... } operetion per group: inside {}, you can run any custom R code, it will
                                      # return a vector that gets assigned to 'bmi_shift' for the rows in that group
    anchor <- unique(anchor_year)     # This line finds a single anchor year (e.g., year 30) for each person.
    if (is.na(anchor)) {              # If no anchor is found (to avoid breaking your logic when a person doesn't uptake the drug)
      bmi_curr_xps                    # the function simply returns the original bmi_curr_xps as-is
    } else {                          # If yes anchor is found, You shift the bmi_curr_xps backward in time (type='lag')
      shift(bmi_curr_xps, type = "lag", n = 5) # by 5 years: 5 years is the duration of efficacy inputs (weight loss and reboud)
    }
  }, by = pid] #Works!
  ###############################################################################################################
  #--------------------------------------------  Tutorial: End  ------------------------------------------------#
  ###############################################################################################################

  # a variable for new weight after applying efficacy input during treatment influenced period
  sp$pop[, new_bmi := {
    anchor <- unique(anchor_year)                             # This line finds a single anchor year for each person.
    if (is.na(anchor)) {                                      # If no anchor is found,
      bmi_curr_xps                                            # the function simply returns the original bmi_curr_xps as-is
    } else {                                                  # If yes anchor is found,
      fcase(                                                  # We need to assign value to the' new_bmi' based on the time periods
        year < anchor,                       bmi_curr_xps,    # Before anchor year/uptake: original values (including baseline year)
        year >= anchor & year <= anchor + 4, baseline_bmi*(1 + bmi_delta),  # During the 5 years where treatment effects manifest
        year > anchor + 4,                   bmi_shift        # After rebound to baseline: revert to the original bmi_curr_xps
        )
    }
  }, by = pid]

  #scenario_1_fn(sp) # Works!
  #test_uptake <- sp$pop[pid %in% sp$pop[entry_year >= 25 & eligible_bi ==1, unique(pid)]]
  #test_uptake <- sp$pop[pid %in% sp$pop[uptake_one == 1, unique(pid)]]

  # Overwriting the original bmi exposure column with the newly created bmi column
  sp$pop[, bmi_curr_xps := new_bmi]

  ########################################################################################
  ################################          SBP         ##################################
  ########################################################################################
  # a variable for the change in SBP after the uptake of semaglutide
  sp$pop[, `:=`(
    sbp_delta = fcase(
      ### All values until year 5, are % weight loss relative to baseline
      trtm_theo == 0, 0,       ### Baseline year, the year patients become eligible
      trtm_theo == 1, sbp_1y,    ### 1st year of treatment: ITT effects from STEP 1 (1 year f/u): mean decrease of 6.2 mmHg from baseline
      trtm_theo == 2, sbp_2y,    ### 2nd year of treatment: ITT effects from STEP 5 (2 year f/u): mean decrease of 5.7 mmHg from baseline
      trtm_theo >= 3, 0        ### STEP 1 extension: assume mean SBP will revert to baseline level after one year without treatment
      ### The efficacy inputs are assigned until here, the rest we will revert to the original exposure column (sbp_curr_xps)
    ))
  ]

  # a variable for each individual's baseline SBP (baseline: year == 24 / ys_rollout == 0)
  # *but, some individuals in the dataset does not have entry at baseline_yr, they came in later than year 24
  sp$pop[, baseline_sbp := {                # sp$pop[, ... , by = pid]: grouping the operation by pid
    anchor <- unique(anchor_year)           # So we can distinguish people who have entry at year 24 and who do not
    if (is.na(anchor)) {                    # If no anchor is found,
      sbp_curr_xps                          # the function simply returns the original sbp_curr_xps as-is
    } else {                                # If yes anchor is found, baseline_sbp will take the sbp at year 24
      sbp_curr_xps[year == (year[trtm_theo == 1][1] - 1)]  # The baseline year will be (anchor - 1): one year before uptake
    }
  }, by = pid] #Works!

  ###############################################################################################################
  #----------------------------  Tutorial: if_else function inside data.table  ---------------------------------#
  ###############################################################################################################
  sp$pop[, sbp_shift := {                 # sp$pop[, ... , by = pid]: grouping the operation by pid
    # { ... } operetion per group: inside {}, you can run any custom R code, it will
    # return a vector that gets assigned to 'sbp_shift' for the rows in that group

    anchor <- unique(anchor_year)        # This line finds a single anchor year (e.g., year 30) for each person.
    if (is.na(anchor)) {                 # If no anchor is found, (This is useful to avoid breaking your logic when a person doesn't meet your anchor condition)
      sbp_curr_xps                       # the function simply returns the original sbp_curr_xps as-is
    } else {                             # If yes anchor is found, You shift the sbp_curr_xps backward in time (type='lag')
      shift(sbp_curr_xps, type = "lag", n = 3) # by 3 years: 3 years is the duration of efficacy inputs (decrease in SBP and rebound)
    }
  }, by = pid] #Works!

  #  A if_else function inside of :={} for creating a new variable
  #  if(condition) {
  #   what should happen if condition is TRUE
  #  } else {
  #   what should happen if condition is FALSE
  #  }
  ###############################################################################################################
  #--------------------------------------------  Tutorial: End  ------------------------------------------------#
  ###############################################################################################################

  # a variable for new weight after applying efficacy input during treatment influenced period
  sp$pop[, new_sbp := {
    anchor <- unique(anchor_year)                             # This line finds a single anchor year for each person.
    if (is.na(anchor)) {                                      # If no anchor is found,
      sbp_curr_xps                                            # the function simply returns the original sbp_curr_xps as-is
    } else {                                                  # If yes anchor is found,
      fcase(                                                  # We need to assign value to the 'new_sbp' based on the time periods
        year < anchor,                       sbp_curr_xps,    # Before anchor year/uptake: original values (including baseline year)
        year >= anchor & year <= anchor + 2, baseline_sbp + sbp_delta,  # During the 2 years where treatment effect manifest
        year >= anchor + 3,                  sbp_shift        # After sbp rebound to baseline: revert to the original sbp_curr_xps
      )
    }
  }, by = pid]

  # Overwriting the original sbp exposure column with the newly created sbp column
    sp$pop[, sbp_curr_xps := new_sbp]

  ########################################################################################
  ###############################          Tchol         #################################
  ########################################################################################
  # a variable for the change in SBP after the uptake of semaglutide
  sp$pop[, `:=`(
    tchol_delta = fcase(
      trtm_theo == 0, 0,        ### Baseline year, the year patients become eligible
      trtm_theo == 1, tchol_1y,    ### 1st year of treatment: ITT effects from STEP 1: mean decrease of 3% from baseline
      trtm_theo == 2, tchol_2y,   ### 2nd year of treatment: ITT effects from STEP 5: mean decrease of 3.3% from baseline
      trtm_theo >= 3, 0         ### STEP 1 extension: assume mean tchol will revert to baseline level after one year without treatment
      ### The efficacy inputs are assigned until here, the rest we will revert to the original exposure column (tchol_curr_xps)
    ))
  ]

  # a variable for each individual's baseline Tchol (baseline: year == 24 / ys_rollout == 0)
  # *but, some individuals in the dataset does not have entry at baseline_yr, they came in later than year 24
  sp$pop[, baseline_tchol := {                 # sp$pop[, ... , by = pid]: grouping the operation by pid
    anchor <- unique(anchor_year)           # So we can distinguish people who have entry at year 24 and who do not
    if (is.na(anchor)) {                    # If no anchor is found,
      tchol_curr_xps                        # the function simply returns the original tchol_curr_xps as-is
    } else {                                # If yes anchor is found, baseline_tchol will take the tchol at year 24
      tchol_curr_xps[year == (year[trtm_theo == 1][1] - 1)]  # The baseline year will be (anchor - 1): one year before uptake
    }
  }, by = pid] #Works!

  ###############################################################################################################
  #----------------------------  Tutorial: if_else function inside data.table  ---------------------------------#
  ###############################################################################################################
  sp$pop[, tchol_shift := {                 # sp$pop[, ... , by = pid]: grouping the operation by pid
    # { ... } operetion per group: inside {}, you can run any custom R code, it will
    # return a vector that gets assigned to 'tchol_shift' for the rows in that group

    anchor <- unique(anchor_year)        # This line finds a single anchor year (e.g., year 30) for each person.
    if (is.na(anchor)) {                 # If no anchor is found, (This is useful to avoid breaking your logic when a person doesn't meet your anchor condition)
      tchol_curr_xps                       # the function simply returns the original tchol_curr_xps as-is
    } else {                             # If yes anchor is found, You shift the tchol_curr_xps backward in time (type='lag')
      shift(tchol_curr_xps, type = "lag", n = 3) # by 3 years: 3 years is the duration of efficacy inputs (decrease in tchol and reboud)
    }
  }, by = pid] #Works!

  # test <- sp$pop[pid %in% c(156, 175)]
  ###############################################################################################################
  #--------------------------------------------  Tutorial: End  ------------------------------------------------#
  ###############################################################################################################

  # a variable for new weight after applying efficacy input during treatment influenced period
  sp$pop[, new_tchol := {
    anchor <- unique(anchor_year)                             # This line finds a single anchor year for each person.
    if (is.na(anchor)) {                                      # If no anchor is found,
      tchol_curr_xps                                          # the function simply returns the original tchol_curr_xps as-is
    } else {                                                  # If yes anchor is found,
      fcase(                                                  # We need to assign value to the 'new_tchol' based on the time periods
        year < anchor,                       tchol_curr_xps,  # Before anchor year/uptake: original values (including baseline year)
        year >= anchor & year <= anchor + 2, baseline_tchol*(1 + tchol_delta),  # During the 2 years where treatment effect manifest
        year >= anchor + 3,                  tchol_shift      # After tchol rebound to baseline: revert to the original tchol_curr_xps
      )
    }
  }, by = pid]

  # Overwriting the original tchol exposure column with the newly created tchol column
    sp$pop[, tchol_curr_xps := new_tchol]
  ### Jane Aug 2025, test scenario on sp$pop --> Looks like until here everything worked!!!

  ##################### Get rid of unnecessary variables ##########################
  # Delete unnecessary variables from synthpop #
  sp$pop[, c("rankstat_sbp", "rankstat_tchol", "ys_rollout", "uptake_rate",
             "eligible_bi", "uptake_psyr", "uptake_one", "anchor_year",
             "entry_year", "trtm_theo", "baseline_bmi", "bmi_shift", "new_bmi",
             "baseline_sbp", "sbp_shift", "new_sbp", "baseline_tchol",
             "tchol_shift", "new_tchol") := NULL]

}

  # Good to get rid of all the intermediate variables,which were created to get to the new exposure column (tchol_curr_xps)
  # The model runs faster without them
  # Keep the variables needed for further analysis

  #scenario_1_fn(sp) # Works!
  #test_uptake <- sp$pop[pid %in% sp$pop[entry_year >= 25 & eligible_bi ==1, unique(pid)]]
  #test_uptake <- sp$pop[pid %in% sp$pop[uptake_one == 1, unique(pid)]]

##################################################################################################################
#----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------    Semaglutide Modelling    ------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------   Scenario 2: Budget impact     --------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
##################################################################################################################
##################################################################################################################
####                                                                                                          ####
####           2. Budget impact scenario:                                                                     ####
####              - Open cohort. Baseline year: 2024, Drug roll-out year: 2025                                ####
####              - Uptake possible at any time point after roll-out                                          ####
####              - Eligibility criteria: BMI >= 35 (eligible_bi), assessed at each year                      ####
####              ############################################################################                ####
####              ** Time-varying uptake rates (by Year after roll-out):                                      ####
####                   Year 1-2 (2025-2026): 20%  - Moderate uptake for early adopters                        ####
####                   Year 3-5 (2027-2029): 40%  - Higher uptake as awareness increases                      ####
####                   Year >5  (>2030)    : 60%  - Potential plateau                                         ####
####              ** ??? This could be problematic, because:                                                  ####
####              ** !!! The overall prevalence in the population over a time period is different from        ####
####              ** !!! the probability of individual uptake each year (incidence)                           ####
####              ############################################################################                ####
####              Assumptions on the time-varying uptake rate:                                                ####
####              1) After 10 year, ~overall prevalence of uptake: 50-60%                                     ####
####              2) We will assume a 10% uptake rate every year                                              ####
####                 * This gives us a 58% overall uptake prevalence over 10 years, 75% over 20 ys            ####
####                 * This gives us a curve of declining numbers of people uptaking the drug over years      ####
####              ############################################################################                ####
####              - People can only uptake the treatment *once* during lifetime                               ####
####              - Treatment last for 2 years (ITT effect, no discontinuation applied to individuals)        ####
####              - Efficacy trajectory (weight loss--rebound--natural progression) applied identical for     ####
####                all eligible patients                                                                     ####
####              - Follow-up period: 5, 10, 20 year                                                          ####
####                                                                                                          ####
##################################################################################################################
##################################################################################################################
#### What do we need as input/function arguments in this scenario: ### Jane Aug 11
#### 1. lifecourse dataset, with original BMI, SBP, Cholesterol
####    Semag:    sppop, bmi_curr_xps, sbp_curr_xps, tchol_curr_xps
####
#### 2. Time horizon: initial year, intervention year, simulation horizon
####    Semag:    2013-2043, with 2024 as baseline year
####
#### 3. Policy input: Eligibility and uptake of treatment
####    Semag:    ys_rollout, eligible_bi, uptake_rate
####
#### 4. Efficacy input: change in BMI, SBP, Cholesterol from RCT
####    Semag:    efficacy input trajectory
####                  Treatment_y1,  xyz,     ### 1st year of treatment: ITT effects from RCT
####                  Treatment_y2,  xyz,     ### 2nd year of treatment: ITT effects from RCT
####                  Treatment_y3,  xyz,     ### Results from STEP 1 extension: waning of 2-year treatment effect
####                  ### The efficacy should be assigned until here, the rest we will revert to the original exposure column
####                  ******
####                  Uncertainty: consistent among estimates (with CI) from multiple RCTs,
####                               In each iteration, sample values from uniform distribution, and use these values
####                               as quantile [0. 1] to sample values from distributions of efficacy inputs
####                               ---> individual fixed effect (?)
####                  Draw random values from uniform distribution for as many MCs
####                  There will be one dataset / input file (?) for efficacy inputs,
####                  ---> In this dataset, each columns is for one efficacy inputs (for all exposures)
####                  One column for mc (1-1000): Col_MC  Col_Quant  Col_BMI_Y1 Col_BMI_Y2 Col_SBP_Y2 (?)
####                  ******
####
#### 5. lifecourse dataset, with modified exposure variables (replacing the original ones)
####    Semag:    bmi_curr_xps := new_bmi

##################################################################################################################

### Scenario 2 - Budget impact of semaglutide 2.4 mg for weight management in German adult population

scenario_2_fn <- function(sp) {

  # Efficacy for bmi at 1st year after treatment
  tbl <- read_fst("./inputs/other_parameters/bmi_1y_samples.fst", as.data.table = TRUE)
  bmi_1y <- as.numeric(tbl[mc == sp$mc_aggr, "bmi_1y"])

  # Efficacy for bmi at 2nd year after treatment
  tbl <- read_fst("./inputs/other_parameters/bmi_2y_samples.fst", as.data.table = TRUE)
  bmi_2y <- as.numeric(tbl[mc == sp$mc_aggr, "bmi_2y"])

  # Efficacy for bmi at 3rd year after treatment
  tbl <- read_fst("./inputs/other_parameters/bmi_3y_samples.fst", as.data.table = TRUE)
  bmi_3y <- as.numeric(tbl[mc == sp$mc_aggr, "bmi_3y"])

  # Efficacy for bmi at 4th year after treatment
  tbl <- read_fst("./inputs/other_parameters/bmi_4y_samples.fst", as.data.table = TRUE)
  bmi_4y <- as.numeric(tbl[mc == sp$mc_aggr, "bmi_4y"])

  # Efficacy for SBP at 1st year after treatment
  tbl <- read_fst("./inputs/other_parameters/sbp_1y_samples.fst", as.data.table = TRUE)
  sbp_1y <- as.numeric(tbl[mc == sp$mc_aggr, "sbp_1y"])

  # Efficacy for SBP at 2nd year of treatment
  tbl <- read_fst("./inputs/other_parameters/sbp_2y_samples.fst", as.data.table = TRUE)
  sbp_2y <- as.numeric(tbl[mc == sp$mc_aggr, "sbp_2y"])

  # Efficacy for total cholesterol at 1st year of treatment
  tbl <- read_fst("./inputs/other_parameters/tchol_1y_samples.fst", as.data.table = TRUE)
  tchol_1y <- as.numeric(tbl[mc == sp$mc_aggr, "tchol_1y"])

  # Efficacy for total cholesterol at 2nd year of treatment
  tbl <- read_fst("./inputs/other_parameters/tchol_2y_samples.fst", as.data.table = TRUE)
  tchol_2y <- as.numeric(tbl[mc == sp$mc_aggr, "tchol_2y"])
  ###########################################################################################################

  # a variable for how long it has been since the baseline (2024)
  baseline_yr <- 24 ## it's 2024, but all years should be relative to year 2000
  sp$pop[, ys_rollout := year - baseline_yr] # baseline_yr == 2024 / year == 24 --> baseline year

  # a variable for uptake rates based on years after drug roll-out
  # * For Budget impact, We assume a 10% uptake rate every year
  # * This uptake rate is independent of individual simulation history, it only depends on years after roll-out
  sp$pop[, uptake_rate := ifelse(ys_rollout >= 1, 0.1, 0)]

  # a variable for defining individual's eligibility for accessing Semaglutide by BMI
  sp$pop[, eligible_bi := ifelse(bmi_curr_xps >= 35, 1, 0)]

  # a variable for the uptake of the drug conditioning on meeting eligibility (person_year)
  sp$pop[, uptake_psyr := ifelse(eligible_bi ==1, rbinom(.N, size = 1, prob = uptake_rate), 0)]

  # a variable for the once-in-a-lifetime uptake during an individual's simulation years
  sp$pop[, uptake_one := 0]
  sp$pop[uptake_psyr == 1, uptake_one := as.integer(.I == .I[which.min(year)]), by = pid]
  ### In this testing, only two ids got the treatment: pid 156, 175

  # a variable for the year someone uptaking the drug
  sp$pop[, anchor_year:=year[uptake_one == 1], by = pid]

  # a variable for the year someone enters the synthetic population
  sp$pop[, entry_year:= min(year), by = pid]

  # a variable for theorectical treatment trajectory
  ### Step 0: Initialize the column
  sp$pop[, trtm_theo := NA_integer_]
  ### Step 1: Identify first year of treatment uptake
  sp$pop[uptake_one == 1, trtm_theo := fifelse(anchor_year > entry_year, 1,0)]
  ### Step 2: Fill the sequence in subsequent years
  sp$pop[, trtm_theo := fifelse(cumsum(!is.na(trtm_theo)) > 0,   # cumsum(): computes a cumulative count of non-NA entries:
                                # fill in a sequence starting from the first non-NA value
                                seq_len(.N) - min(which(!is.na(trtm_theo))) + fifelse(anchor_year > entry_year, 1L, 0L),
                                NA_integer_),
         by = pid]
  ### Step 3: Replace NAs with 0
  sp$pop[, trtm_theo := fifelse(is.na(trtm_theo), 0, trtm_theo)]

  ########################################################################################
  ################################          BMI         ##################################
  ########################################################################################
  # a variable for the change in BMI after the uptake of semaglutide
  sp$pop[, `:=`(
    bmi_delta = fcase(
      ### All values until year 5, are % weight loss relative to baseline
      trtm_theo == 0, 0,       ### Baseline year, the year patients become eligible
      trtm_theo == 1, bmi_1y,  ### 1st year of treatment: ITT effects from STEP 1 (1 year f/u): 14.9% weight loss (uncertainty?)
      trtm_theo == 2, bmi_2y,  ### 2nd year of treatment: ITT effects from STEP 5 (2 year f/u): 15.2% weight loss (uncertainty?)
      trtm_theo == 3, bmi_3y,  ### STEP 1 extension: regain 2/3 of weight loss in 1st year treatment cessation (uncertainty?)
      trtm_theo == 4, bmi_4y,  ### weight regain at 2nd year treatment cessation (uncertainty?)
      trtm_theo >= 5, 0        ### weight return to baseline weight by 3rd year treatment cessation (uncertainty?)
      ### The efficacy inputs are assigned until here, the rest we will revert to the original exposure column (bmi_curr_xps)
    ))
  ]

  # a variable for each individual's baseline weight (baseline: year == 24 / ys_rollout == 0)
  # ------------------------------------------ Aug 22, 2025, Jane ------------------------------------------------- #
  # *but, some individuals in the dataset does not have entry at baseline_yr, they came in later than year 24
  # *and, some individuals came in later than year 24, and because they have a chance to uptake drug begining year 25
  # *they might already uptake the drug at the first year, which means they won't have a year (anchor-1),
  # *so what shall we use as a baseline BMI value here...
  sp$pop[, baseline_bmi := {              # sp$pop[, ... , by = pid]: grouping the operation by pid
    anchor <- unique(anchor_year)         # So we can distinguish people who have entry at year 24 and who do not
    if (is.na(anchor)) {                  # If no anchor is found,
      bmi_curr_xps                        # the function simply returns the original bmi_curr_xps as-is
    } else {                              # If yes anchor is found
      bmi_curr_xps[year == (year[trtm_theo == 1][1] - 1)] # baseline_bmi will take the BMI at the year before trtm_theo==1
    }
  }, by = pid] #Works!

  ###############################################################################################################
  #----------------------------  Tutorial: if_else function inside data.table  ---------------------------------#
  ###############################################################################################################
  sp$pop[, bmi_shift := {             # sp$pop[, ... , by = pid]: grouping the operation by pid
    # { ... } operetion per group: inside {}, you can run any custom R code, it will
    # return a vector that gets assigned to 'bmi_shift' for the rows in that group
    anchor <- unique(anchor_year)     # This line finds a single anchor year (e.g., year 30) for each person.
    if (is.na(anchor)) {              # If no anchor is found (to avoid breaking your logic when a person doesn't uptake the drug)
      bmi_curr_xps                    # the function simply returns the original bmi_curr_xps as-is
    } else {                          # If yes anchor is found, You shift the bmi_curr_xps backward in time (type='lag')
      shift(bmi_curr_xps, type = "lag", n = 5) # by 5 years: 5 years is the duration of efficacy inputs (weight loss and reboud)
    }
  }, by = pid] #Works!
  ###############################################################################################################
  #--------------------------------------------  Tutorial: End  ------------------------------------------------#
  ###############################################################################################################

  # a variable for new weight after applying efficacy input during treatment influenced period
  sp$pop[, new_bmi := {
    anchor <- unique(anchor_year)                             # This line finds a single anchor year for each person.
    if (is.na(anchor)) {                                      # If no anchor is found,
      bmi_curr_xps                                            # the function simply returns the original bmi_curr_xps as-is
    } else {                                                  # If yes anchor is found,
      fcase(                                                  # We need to assign value to the' new_bmi' based on the time periods
        year < anchor,                       bmi_curr_xps,    # Before anchor year/uptake: original values (including baseline year)
        year >= anchor & year <= anchor + 4, baseline_bmi*(1 + bmi_delta),  # During the 5 years where treatment effects manifest
        year > anchor + 4,                   bmi_shift        # After rebound to baseline: revert to the original bmi_curr_xps
      )
    }
  }, by = pid]

  # Overwriting the original bmi exposure column with the newly created bmi column
    sp$pop[, bmi_curr_xps := new_bmi]

  ########################################################################################
  ################################          SBP         ##################################
  ########################################################################################
  # a variable for the change in SBP after the uptake of semaglutide
  sp$pop[, `:=`(
    sbp_delta = fcase(
      ### All values until year 5, are % weight loss relative to baseline
      trtm_theo == 0, 0,       ### Baseline year, the year patients become eligible
      trtm_theo == 1, sbp_1y,    ### 1st year of treatment: ITT effects from STEP 1 (1 year f/u): mean decrease of 6.2 mmHg from baseline
      trtm_theo == 2, sbp_2y,    ### 2nd year of treatment: ITT effects from STEP 5 (2 year f/u): mean decrease of 5.7 mmHg from baseline
      trtm_theo >= 3, 0        ### STEP 1 extension: assume mean SBP will revert to baseline level after one year without treatment
      ### The efficacy inputs are assigned until here, the rest we will revert to the original exposure column (sbp_curr_xps)
    ))
  ]

  # a variable for each individual's baseline SBP (baseline: year == 24 / ys_rollout == 0)
  # *but, some individuals in the dataset does not have entry at baseline_yr, they came in later than year 24
  sp$pop[, baseline_sbp := {                # sp$pop[, ... , by = pid]: grouping the operation by pid
    anchor <- unique(anchor_year)           # So we can distinguish people who have entry at year 24 and who do not
    if (is.na(anchor)) {                    # If no anchor is found,
      sbp_curr_xps                          # the function simply returns the original sbp_curr_xps as-is
    } else {                                # If yes anchor is found, baseline_sbp will take the sbp at year 24
      sbp_curr_xps[year == (year[trtm_theo == 1][1] - 1)]  # The baseline year will be (anchor - 1): one year before uptake
    }
  }, by = pid] #Works!

  ###############################################################################################################
  #----------------------------  Tutorial: if_else function inside data.table  ---------------------------------#
  ###############################################################################################################
  sp$pop[, sbp_shift := {                 # sp$pop[, ... , by = pid]: grouping the operation by pid
    # { ... } operetion per group: inside {}, you can run any custom R code, it will
    # return a vector that gets assigned to 'sbp_shift' for the rows in that group

    anchor <- unique(anchor_year)        # This line finds a single anchor year (e.g., year 30) for each person.
    if (is.na(anchor)) {                 # If no anchor is found, (This is useful to avoid breaking your logic when a person doesn't meet your anchor condition)
      sbp_curr_xps                       # the function simply returns the original sbp_curr_xps as-is
    } else {                             # If yes anchor is found, You shift the sbp_curr_xps backward in time (type='lag')
      shift(sbp_curr_xps, type = "lag", n = 3) # by 3 years: 3 years is the duration of efficacy inputs (decrease in SBP and rebound)
    }
  }, by = pid] #Works!

  #  A if_else function inside of :={} for creating a new variable
  #  if(condition) {
  #   what should happen if condition is TRUE
  #  } else {
  #   what should happen if condition is FALSE
  #  }
  ###############################################################################################################
  #--------------------------------------------  Tutorial: End  ------------------------------------------------#
  ###############################################################################################################

  # a variable for new weight after applying efficacy input during treatment influenced period
  sp$pop[, new_sbp := {
    anchor <- unique(anchor_year)                             # This line finds a single anchor year for each person.
    if (is.na(anchor)) {                                      # If no anchor is found,
      sbp_curr_xps                                            # the function simply returns the original sbp_curr_xps as-is
    } else {                                                  # If yes anchor is found,
      fcase(                                                  # We need to assign value to the 'new_sbp' based on the time periods
        year < anchor,                       sbp_curr_xps,    # Before anchor year/uptake: original values (including baseline year)
        year >= anchor & year <= anchor + 2, baseline_sbp + sbp_delta,  # During the 2 years where treatment effect manifest
        year >= anchor + 3,                  sbp_shift        # After sbp rebound to baseline: revert to the original sbp_curr_xps
      )
    }
  }, by = pid]

  # Overwriting the original sbp exposure column with the newly created sbp column
    sp$pop[, sbp_curr_xps := new_sbp]

  ########################################################################################
  ###############################          Tchol         #################################
  ########################################################################################
  # a variable for the change in SBP after the uptake of semaglutide
  sp$pop[, `:=`(
    tchol_delta = fcase(
      trtm_theo == 0, 0,        ### Baseline year, the year patients become eligible
      trtm_theo == 1, tchol_1y,    ### 1st year of treatment: ITT effects from STEP 1: mean decrease of 3% from baseline
      trtm_theo == 2, tchol_2y,   ### 2nd year of treatment: ITT effects from STEP 5: mean decrease of 3.3% from baseline
      trtm_theo >= 3, 0         ### STEP 1 extension: assume mean tchol will revert to baseline level after one year without treatment
      ### The efficacy inputs are assigned until here, the rest we will revert to the original exposure column (tchol_curr_xps)
    ))
  ]

  # a variable for each individual's baseline Tchol (baseline: year == 24 / ys_rollout == 0)
  # *but, some individuals in the dataset does not have entry at baseline_yr, they came in later than year 24
  sp$pop[, baseline_tchol := {                 # sp$pop[, ... , by = pid]: grouping the operation by pid
    anchor <- unique(anchor_year)           # So we can distinguish people who have entry at year 24 and who do not
    if (is.na(anchor)) {                    # If no anchor is found,
      tchol_curr_xps                        # the function simply returns the original tchol_curr_xps as-is
    } else {                                # If yes anchor is found, baseline_tchol will take the tchol at year 24
      tchol_curr_xps[year == (year[trtm_theo == 1][1] - 1)]  # The baseline year will be (anchor - 1): one year before uptake
    }
  }, by = pid] #Works!

  ###############################################################################################################
  #----------------------------  Tutorial: if_else function inside data.table  ---------------------------------#
  ###############################################################################################################
  sp$pop[, tchol_shift := {                 # sp$pop[, ... , by = pid]: grouping the operation by pid
    # { ... } operetion per group: inside {}, you can run any custom R code, it will
    # return a vector that gets assigned to 'tchol_shift' for the rows in that group

    anchor <- unique(anchor_year)        # This line finds a single anchor year (e.g., year 30) for each person.
    if (is.na(anchor)) {                 # If no anchor is found, (This is useful to avoid breaking your logic when a person doesn't meet your anchor condition)
      tchol_curr_xps                       # the function simply returns the original tchol_curr_xps as-is
    } else {                             # If yes anchor is found, You shift the tchol_curr_xps backward in time (type='lag')
      shift(tchol_curr_xps, type = "lag", n = 3) # by 3 years: 3 years is the duration of efficacy inputs (decrease in tchol and reboud)
    }
  }, by = pid] #Works!

  # test <- sp$pop[pid %in% c(156, 175)]
  ###############################################################################################################
  #--------------------------------------------  Tutorial: End  ------------------------------------------------#
  ###############################################################################################################

  # a variable for new weight after applying efficacy input during treatment influenced period
  sp$pop[, new_tchol := {
    anchor <- unique(anchor_year)                             # This line finds a single anchor year for each person.
    if (is.na(anchor)) {                                      # If no anchor is found,
      tchol_curr_xps                                          # the function simply returns the original tchol_curr_xps as-is
    } else {                                                  # If yes anchor is found,
      fcase(                                                  # We need to assign value to the 'new_tchol' based on the time periods
        year < anchor,                       tchol_curr_xps,  # Before anchor year/uptake: original values (including baseline year)
        year >= anchor & year <= anchor + 2, baseline_tchol*(1 + tchol_delta),  # During the 2 years where treatment effect manifest
        year >= anchor + 3,                  tchol_shift      # After tchol rebound to baseline: revert to the original tchol_curr_xps
      )
    }
  }, by = pid]

  # Overwriting the original tchol exposure column with the newly created tchol column
    sp$pop[, tchol_curr_xps := new_tchol]
  ### Jane Aug 2025, test scenario on sp$pop --> Looks like until here everything worked!!!

  ##################### Get rid of unnecessary variables ##########################
  # Delete unnecessary variables from synthpop #
  sp$pop[, c("rankstat_sbp", "rankstat_tchol", "ys_rollout", "uptake_rate",
             "eligible_bi", "uptake_psyr", "uptake_one", "anchor_year",
             "entry_year", "trtm_theo", "baseline_bmi", "bmi_shift", "new_bmi",
             "baseline_sbp", "sbp_shift", "new_sbp", "baseline_tchol",
             "tchol_shift", "new_tchol") := NULL]

}

# scenario_2_fn(sp) # Works! Aug 25, 2025, Jane
#test_uptake <- sp$pop[pid %in% sp$pop[entry_year >= 25 & eligible_bi ==1, unique(pid)]]
#test_uptake <- sp$pop[pid %in% sp$pop[uptake_one == 1, unique(pid)]]
#test_uptake <- sp$pop[pid %in% sp$pop[entry_year == anchor_year, unique(pid)]]

