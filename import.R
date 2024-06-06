

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
  nna <- c("NA",""," ", "C", "D", "K","A", "Non disponible", "no disponible")
  nnd <- paste0("data/",nn, ".ods")
  nn <- read_ods(nnd, na = nna) |>
    clean_names() |>
    mutate_if(is.character, as.factor)
 bn <- read_ods(nnd, sheet=2)
 var_label(nn) <- bn$nom
  return(nn)
}

fich <- c("atcd","demog","finet","patho","result")
for (f in fich){
  assign(f, importr(f))
}

# Demog

demog <- demog |>
  mutate(naisdte = as.character(naisdte)) |>
  mutate(incldte = dmy(as.character(incldte))) |>
  mutate(naisdte = dmy(naisdte)) |>
  mutate(age = as.numeric(incldte - naisdte)/365.25) |>
  relocate(age, .after = "naisdte") |>
  mutate(imc = bmiph(imc))

var_label(demog$age) <- "Âge"
var_label(demog$imc) <- "IMC"

# ATCD

atcd <- atcd |> 
  mutate(tabacon = fct_relevel(tabacon,
    "Aucun", "Actif", "Sevré"))

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

tt <- left_join(result, demog, by="subjid") |>
  mutate(sca3  = ifelse((tropoh3> 34.2 & sex == "Masculin") |(tropoh3> 15.6 & sex == "Féminin") ,"yes","no")) |>
  mutate(sca3 = as.factor(sca3)) |>
  mutate(sca3 = fct_relevel(sca3,"yes", "no")) |>
  drop_na(sca3) |>
  mutate(tp0 = ifelse((tropoh0 > 15.6 &
                         sex == "Féminin") |
                        (tropoh0 > 34.2 & sex == "Masculin"),
                      "yes",
                      "no"
  )) |>
  mutate(tp0 = as.factor(tp0)) |>
  mutate(tp0 = fct_relevel(tp0, "yes", "no")) |>
  mutate(cp0 = ifelse(copepth0 > 10, "yes", "no")) |>
  mutate(cp0 = as.factor(cp0)) |> 
  mutate(cp0 = fct_relevel(cp0, "yes", "no")) |> 
  mutate(tpcp0 = ifelse(tp0 == "yes" | cp0 == "yes", "yes", "no")) |> 
  mutate(tpcp0 = as.factor(tpcp0)) |> 
  mutate(tpcp0 = fct_relevel(tpcp0, "yes", "no"))
#
var_label(tt$tp0) <- "Troponine h0 anormale"
var_label(tt$cp0) <- "Copeptine h0 anormale"
var_label(tt$tpcp0) <- "Troponine ou copeptine h0 anormale"
var_label(tt$sca3) <- "SCA ST-"
#
# finet
zz1 <- as.numeric(dmy_hms(paste0(finet$sortiurgdte," ",finet$sortiurghr)))
zz2 <- as.numeric(dmy_hms(paste0(patho$urgencdte," ", patho$urgenhr)))
zz <- (zz1 - zz2)/3600
finet <- finet |> 
  mutate(duree_urg = zz)

#
save(atcd, demog, finet, patho, tt, file="data/copsca.RData")