#merge information relevant to creating line graphs that show trends of AMR
#prevelance per serovars. 
library(tidyverse)
library(lubridate)
library(stringr)
library(readr)

amr <- read_tsv(
  "combined_amr_with_sample.tsv",
  col_types = cols(
    Sample = col_character(),
    `Element symbol` = col_character(),
    Type = col_character(),
    Class = col_character()
  ),
  col_select = c("Sample", "Element symbol", "Type", "Class")
)

meta <- read_csv(
  "beef_isolates.csv",
  col_types = cols(
    Run            = col_character(),
    Serovar        = col_character(),
    `Create date`  = col_datetime(),  # readr handles ISO-like dttm; we’ll standardize TZ next
    `Collection date` = col_guess()   # can be numeric or date-like; we’ll normalize below
  ),
  col_select = c("Run", "Serovar", "Create date", "Collection date")
)
problems(amr)
problems(meta)