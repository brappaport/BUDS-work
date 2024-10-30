# Run scripts for processing of all data for BUDS project
library(rmarkdown)
library(knitr)
library(here)

here <- here::here
here::i_am("BUDS_run_all_processing.R")

# Extract the R code from the RMarkdown file
purl(here("management/data01_BUDS.Rmd"), output = here("management/data01_BUDS.R"))
purl(here("management/data02_BUDS.Rmd"), output = here("management/data02_BUDS.R"))
purl(here("management/data03_BUDS.Rmd"), output = here("management/data03_BUDS.R"))
purl(here("management/data04_BUDS.Rmd"), output = here("management/data04_BUDS.R"))
purl(here("management/data05_BUDS_EMA.Rmd"), output = here("management/data05_BUDS_EMA.R"))
purl(here("management/data06_BUDS_ratings.Rmd"), output = here("management/data06_BUDS_ratings.R"))

# Source the extracted R scripts
rm(list=ls()); source("management/data01_BUDS.R"); rm(list=ls())
rm(list=ls()); source("management/data02_BUDS.R"); rm(list=ls())
rm(list=ls()); source("management/data03_BUDS.R"); rm(list=ls())
rm(list=ls()); source("management/data04_BUDS.R"); rm(list=ls())
rm(list=ls()); source("management/data05_BUDS_EMA.R"); rm(list=ls())
rm(list=ls()); source("management/data06_BUDS_ratings.R")

print("Done running all BUDS processing scripts!")
