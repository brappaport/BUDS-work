

## Load packages
library(dplyr) # for data manipulation
library(psych) # for calculating descriptives, autocorrelation, RMSSD, etc.
library(lme4) # for multilevel models
library(lmerTest) # for getting p-values from lmer models
library(sjPlot) # for quick and dirty plots
library(here)

select <- dplyr::select
dat_all <- read.csv(file=here("./data/FRT_ERP_EMA_NF_Prepped_wcov_revLPP_30_60_90_Acc_timestamps.csv"))

#### -------------------- General data cleaning -------------------- 

dat_all <- dat_all %>% 
  ## convert to factors
  mutate(sub = factor(sub),
         time_of_dayF = factor(time_of_day, levels = 1:4, 
                               labels = c('morning', 'afternoon', 'evening', 'night')),
         gender = factor(gender, levels = c(1,2), labels = c('male', 'female')),
         site = factor(site, levels = c(1, 2), labels = c('CU', 'NU')),
         group = factor(group, levels = c(1, 3, 2))) %>% 
  # convert to date/time format
  mutate(dt_start = as.POSIXct(strptime(paste(date_start, time_start), format = "%m/%d/%Y %I:%M:%S %p")),
         dt_end = as.POSIXct(strptime(paste(date_end, time_end), format = "%m/%d/%Y %I:%M:%S %p")),
         secsToComplete = as.numeric(difftime(dt_end, dt_start, units = "secs")),
         minsToComplete = as.numeric(difftime(dt_end, dt_start, units = "mins")))

## make time-varying (daily or within-day) variables
dat_all <- dat_all %>% 
  group_by(sub, day) %>% 
  mutate(nSurveys_daily = sum(!is.na(date_end)),
         nConsecPairs_daily = sum(!is.na(date_end) & !is.na(lag(date_end))),
         ## center within day (within people; .cwd = centered within day)
         time_of_day.cwd = time_of_day - mean(time_of_day, na.rm = TRUE),
         ## calculate amount of time between surveys (excluding overnight lags and missing surveys)
         SecsSinceLastCompleted = as.numeric(difftime(dt_end, lag(zoo::na.locf(dt_end, na.rm = FALSE)), units = "secs")),
         MinsSinceLastCompleted = SecsSinceLastCompleted / 60) %>% 
  ungroup()

## Create participant-level variables for EMA completion rates
dat_all <- dat_all %>% 
  group_by(sub) %>% 
  mutate(nSurveysTotal = sum(!is.na(date_end)),
         ComplianceRate = nSurveysTotal / 28) %>% 
  ungroup()

## Center EMA variables
# Person-mean centering (more generally, "centering within cluster")
# .pm = person mean
# .cwp = centered within person
dat_all <- dat_all %>% 
  group_by(sub) %>%
  mutate(
    across(.cols = c('happy_EMA', 'angry_EMA', 'sad_EMA'),
           .fns = list('pm' = ~ mean(., na.rm = TRUE), # create person-level means
                       'psd' = ~ sd(., na.rm = TRUE), # create person-level SDs (variability)
                       'cwp' = ~ . - mean(., na.rm = TRUE)), # create person-centered variables
           .names = '{.col}.{.fn}')
  ) %>% 
  ungroup()

## Make lagged and leaded variables for instability and inertia analyses
dat_all <- dat_all %>%
  group_by(sub, day) %>% # note: grouping by participant *and day* excludes overnight lags
  mutate(
    across(.cols = c('happy_EMA', 'angry_EMA', 'sad_EMA', # uncentered EMA variables
                     paste0(c('happy_EMA', 'angry_EMA', 'sad_EMA'), '.cwp')), # person-centered EMA variables
           .fns = list('lead' = ~ lead(., n = 1L), # create leaded affect (at time t+1)
                       'lag' = ~ lag(., n = 1L)), # create lagged affect (at time t-1)
           .names = '{.fn}.{.col}')
  ) %>% 
  ungroup()

#### -------------------- Data preparation for instability analyses (adjusting for variability in time between surveys) -------------------- 
  
## note: code adapted from Sarah Sperry
## note: this implements the method from Jahng et al. (2008) [https://doi.org/10.1037/a0014173]
  
## Create df wih no rows of NAs
dat_noNA <- subset(dat_all, !is.na(dt_end))

## Create lagged or adjusted time variables
dat_noNA <- dat_noNA %>% 
  # group by subject and day so it excludes (a) lags spanning 2 subjects, and (b) overnight lags by setting the lagged variable to NA
  group_by(sub, day) %>%
  mutate(
    ## Create lagged time variables
    lag.dt_start = lag(dt_start),
    lag.dt_end = lag(dt_end)) %>%
  mutate(
    ## Time difference between surveys
    Timedif_start = as.numeric(difftime(dt_start, lag.dt_start, units = "mins")),
    Timedif_end = as.numeric(difftime(dt_end, lag.dt_end, units = "mins")),
    ## Create variable representing time lag between each row (time t) and the next observation (time t+1)
    lead.Timedif_start = lead(Timedif_start),
    lead.Timedif_end = lead(Timedif_end)
  ) %>% 
  ungroup()

## Exclude time lags > 1.5 SD above the sample mean
dat_noNA_1.5SD <- dat_noNA %>% 
  filter(Timedif_start <= mean(Timedif_start, na.rm = TRUE) + 1.5*sd(Timedif_start, na.rm = TRUE) | is.na(Timedif_start))

## How many pairs of surveys are there within the same day?
dat_noNA_1.5SD %>% with(sum(!is.na(dt_start) & !is.na(lag.dt_start)))

## Find the median time lag between surveys
median.value_start <- median(dat_noNA_1.5SD$Timedif_start, na.rm = TRUE) 
## Divide each time lag by the median
dat_noNA_1.5SD$Timedif_start.adj <- dat_noNA_1.5SD$Timedif_start / median.value_start

## Calculate the WASD (which is needed to calculate lambda) by creating temporary ('temp_') variables that will get deleted later
dat_noNA_1.5SD <- dat_noNA_1.5SD %>% 
  group_by(sub, day) %>% 
  mutate(
    # Calculate the successive difference
    across(.cols = c('happy_EMA', 'angry_EMA', 'sad_EMA'),
           .fns = ~ lead(.) - ., 
           .names = 'temp_{.col}.SD'),
    # Divide successive difference by (time diff/median)
    across(.cols = paste0('temp_', c('happy_EMA', 'angry_EMA', 'sad_EMA'), '.SD'),
           .fns = ~ . / Timedif_start.adj, 
           .names = '{.col}W')
    ) %>% 
  ## change the ".SDW" suffix to ".WSD"
  rename_with(.cols = ends_with('.SDW'),
              .fn = ~ gsub('.SDW', '.WSD', .x)) %>% 
  mutate(
    # Take the absolute value
    across(.cols = paste0('temp_', c('happy_EMA', 'angry_EMA', 'sad_EMA'), '.WSD'),
           .fns = ~ abs(.), 
           .names = '{.col}A')
    ) %>% 
  ## change the ".WSDA" suffix to ".WASD"
  rename_with(.cols = ends_with('.WSDA'),
              .fn = ~ gsub('.WSDA', '.WASD', .x)) %>% 
  ungroup()

## Listwise deletion b/c some surveys are only missing some of the EMA items and smooth.spline doesn't handle NAs
newdata_happy <- na.omit(dat_noNA_1.5SD[,c('temp_happy_EMA.WASD', 'Timedif_start.adj')])
newdata_angry <- na.omit(dat_noNA_1.5SD[,c('temp_angry_EMA.WASD', 'Timedif_start.adj')])
newdata_sad <- na.omit(dat_noNA_1.5SD[,c('temp_sad_EMA.WASD', 'Timedif_start.adj')])

## Calculate lambda for each EMA item separately
spline.modelHappy <- smooth.spline(x = newdata_happy$Timedif_start.adj, y = newdata_happy$temp_happy_EMA.WASD)
lambdaHappy <- spline.modelHappy$lambda

spline.modelAngry <- smooth.spline(x = newdata_angry$Timedif_start.adj, y = newdata_angry$temp_angry_EMA.WASD)
lambdaAngry <- spline.modelAngry$lambda

spline.modelSad <- smooth.spline(x = newdata_sad$Timedif_start.adj, y = newdata_sad$temp_sad_EMA.WASD)
lambdaSad <- spline.modelSad$lambda

## Create ADJUSTED successive difference variables incorporating the lambda values
dat_noNA_1.5SD <- dat_noNA_1.5SD %>% 
  group_by(sub, day) %>% 
  mutate(
    # Take the successive difference
    happy_EMA.SD = lead.happy_EMA - happy_EMA,
    angry_EMA.SD = lead.angry_EMA - angry_EMA,
    sad_EMA.SD = lead.sad_EMA - sad_EMA,
    # divide successive difference by (time diff/median)*lamda
    happy_EMA.WSD = happy_EMA.SD / (Timedif_start.adj^lambdaHappy),
    angry_EMA.WSD = angry_EMA.SD / (Timedif_start.adj^lambdaAngry),
    sad_EMA.WSD = sad_EMA.SD / (Timedif_start.adj^lambdaSad),
    # Take the absolute value
    happy_EMA.WASD = abs(happy_EMA.WSD),
    angry_EMA.WASD = abs(angry_EMA.WSD),
    sad_EMA.WASD = abs(sad_EMA.WSD),
    # Square the weighted absolute successive difference
    happy_EMA.WSSD = happy_EMA.WASD^2,
    angry_EMA.WSSD = angry_EMA.WASD^2,
    sad_EMA.WSSD = sad_EMA.WASD^2
  ) %>% 
  ungroup()

## remove the 'temporary' variables that were used to calculate lambda
dat_noNA_1.5SD <- dat_noNA_1.5SD %>% 
  select(-starts_with('temp_'))



#### -------------------- Descriptives --------------------

## Create a dataframe with 1 row per person
dat_all_1rowPerPerson <- dat_all %>% 
  filter(!duplicated(sub)) %>% # filter so there is 1 row per person
  select(-c(happy_EMA, angry_EMA, sad_EMA, # remove some of the variables that vary within-person (should remove them all to avoid confusion)
            starts_with('lead.'), starts_with('lag.'), ends_with('.cwp'))) 


## Analyzable N
dat_all_1rowPerPerson %>% with(n_distinct(sub))

## calculate the average number of EMA surveys per day for each person
dat_all_1rowPerPerson <- dat_all_1rowPerPerson %>% 
  full_join(dat_all %>% 
              distinct(sub, day, .keep_all = TRUE) %>% # filter so there's one row per person-day
              group_by(sub) %>% 
              summarize(nSurveys_daily.pm = mean(nSurveys_daily, na.rm = TRUE),
                        nConsecPairs.psum = sum(nConsecPairs_daily, na.rm = TRUE),
                        nConsecPairs_daily.pm = mean(nConsecPairs_daily, na.rm = TRUE)))

## table of participant characteristics (demographics, etc.)
# tableone::CreateTableOne(vars = c('gender', 'age', 'nSurveysTotal', 'ComplianceRate', 'nSurveys_daily.pm', 
#                                   'nConsecPairs.psum', 'nConsecPairs_daily.pm'),
#                          data = dat_1rowPerPerson,
#                          strata = 'site',
#                          addOverall = TRUE)


#### -------------------- Inertia analysis example --------------------

## note: beep = survey number (in this case, 1-28 because there are 4 EMA survey prompts/day for 7 days)
## note: this example includes the LPP to happy faces (LPPem_Happy_CT) as a level 2 moderator of the random autoregressive slope of happiness
## note: it is common to control for the linear effect of time (e.g., beep) because some EMA studies have reported linear time trends (e.g., https://pubmed.ncbi.nlm.nih.gov/25844974)

inertia_LPP_happy <- lmer(lead.happy_EMA ~ day + happy_EMA.cwp*LPPem_Happy_CT + (1+happy_EMA.cwp|sub), 
                          data = dat_all, REML = FALSE)
# plot_model(inertia_LPP_happy, type = 'pred', terms = c('happy_EMA.cwp', 'LPPem_Happy_CT')) # quick and dirty plot of the cross-level moderation


#### -------------------- Instability analysis example -------------------- 

instability_LPP_happy <- lmer(happy_EMA.WSSD ~ LPPem_Happy_CT + (1|sub), 
                              data = dat_noNA_1.5SD, REML = TRUE)


#### -------------------- Calculate inertia estimates for each person using BLUPs -------------------- 

## Note: this generates values called "best linear unbiased predictions" (BLUPs)
## BLUPs can be useful for plotting (see Figure 2 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7914176/ for an example), but you may need to partial out covariates first

blup_inertia_happy <- lmer(lead.happy_EMA ~ trial + happy_EMA.cwp + (1+happy_EMA.cwp|sub), 
                           data = dat_all, REML = TRUE)
blup_inertia_angry <- lmer(lead.angry_EMA ~ trial + angry_EMA.cwp + (1+angry_EMA.cwp|sub), 
                           data = dat_all, REML = TRUE)
blup_inertia_sad <- lmer(lead.sad_EMA.cwp ~ trial + sad_EMA.cwp + (1+sad_EMA.cwp|sub), 
                         data = dat_all, REML = TRUE)

## Merge these BLUPs into the participant-level data frame
dat_all_1rowPerPerson <- dat_all_1rowPerPerson %>% 
  full_join(data.frame(sub = row.names(coef(modAutocorBlup_happy)$sub),
                       inertiaBLUP_happy = coef(modAutocorBlup_happy)$sub$happy_EMA.cwp),
            by = "sub") %>% 
  full_join(data.frame(sub = row.names(coef(modAutocorBlup_angry)$sub),
                       inertiaBLUP_angry = coef(modAutocorBlup_angry)$sub$angry_EMA.cwp),
            by = "sub") %>% 
  full_join(data.frame(sub = row.names(coef(modAutocorBlup_sad)$sub),
                       inertiaBLUP_sad = coef(modAutocorBlup_sad)$sub$sad_EMA.cwp),
            by = "sub")


#### -------------------- Use the two-step approach to calculate autocorrelation and MSSD --------------------

# autoR_output <- autoR(dat_all[,c("happy_EMA","angry_EMA","sad_EMA")],
#                       group = dat_all$sub,
#                       lag = 1,
#                       use = "pairwise")
# 
# TwoStepDF <- data.frame(autoR = autoR_output$autoR,
#                         rmssd = autoR_output$rmssd,
#                         mssd = mssd(dat[,c("happy_EMA","angry_EMA","sad_EMA")],
#                                     group = dat_all$sub, lag = 1, na.rm = TRUE))
# TwoStepDF$sub <- row.names(TwoStepDF)
# 
# ## merge autocorrelations, RMSSDs, and MSSDs into the participant-level data frame
# dat_1rowPerPerson <- left_join(dat_1rowPerPerson,
#                                TwoStepDF)