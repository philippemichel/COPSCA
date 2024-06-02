

#  ------------------------------------------------------------------------
#
# Title : Import COPSA
#    By : PhM
#  Date : 2024-05-31
#
#  ------------------------------------------------------------------------

library(tidyverse)
library(readODS)
library(labelled)
library(janitor)
library(lubridate)
library(baseph)

  #
importr <- function(nn){
  nna <- c("NA",""," ", "C", "D", "K","A")
  nnd <- paste0("datas/",nn, ".ods")
  nn <- read_ods(nnd, na = nna) |>
    clean_names() |>
    mutate_if(is.character, as.factor)
 bn <- read_ods(nnd, sheet=2)
 var_label(nn) <- bn$nom
  return(nn)
}

fich <- c("atcd","demog","finet","patho","resultats")
for (f in fich){
  assign(f, importr(f))
}

# Demog

demog <- demog |>
  mutate(naisdte = as.character(naisdte)) |>
  mutate(naisdte = paste0(str_sub(naisdte,1,6),"19",str_sub(naisdte,7,8))) |>
  mutate(incldte = dmy(as.character(incldte))) |>
  mutate(naisdte = dmy(naisdte)) |>
  mutate(age = as.numeric(incldte - naisdte)/365.25) |>
  relocate(age, .after = "naisdte") |>
  mutate(imc = bmiph(imc))

var_label(demog$age) <- "Âge"
var_label(demog$imc) <- "IMC"

# Patho

difdate <- function(dd1,hh1,dd2,hh2){
  zz1 <- as.numeric(dmy_hms(paste0(dd1," ",hh1)))
  zz2 <- as.numeric(dmy_hms(paste0(dd2," ",hh2)))
  dif <- (zz2 - zz1)/60
  return(dif)
}

patho <- patho |>
  mutate(evah3dte = as.character(evah3dte)) |>
  mutate(sympt_urg = difdate(symptodte, symptohr, urgencdte, urgenhr)) |>
  mutate(urg_h3 =difdate(urgencdte, urgenhr,evah3dte, evah3hr))

var_label(patho$sympt_urg) <- "Délai symptômes-urgences (min)"
var_label(patho$urg_h3) <- "Délai urgences-évaluation h3 (min)"

# Résultats

tt <- left_join(resultats, demog, by="subjid") |>
  mutate(sca3  = ifelse((tropoh3> 34.2 & sex == "Masculin") |(tropoh3> 15.6 & sex == "Féminin") ,"yes","no")) |>
  mutate(sca3 = as.factor(sca3)) |>
  mutate(sca0  = ifelse((tropoh0> 34.2 & sex == "Masculin") |(tropoh0> 15.6 & sex == "Féminin") ,"yes","no")) |>
  mutate(sca0 = as.factor(sca0))


#
save(atcd, demog, finet, patho, tt, file="datas/copsca.RData")
