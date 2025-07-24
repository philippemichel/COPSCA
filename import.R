#  ------------------------------------------------------------------------
#
# Title : Import COPSA
#    By : PhM
#  Date : 2024-05-31
#
#  ------------------------------------------------------------------------


importph <- function() {
  library(tidyverse)
  library(readODS)
  library(labelled)
  library(janitor)
  library(stringr)
  library(lubridate)
  library(baseph)

  #
  importr <- function(nn) {
    nna <- c("NA", "", " ", "C", "D", "K", "A", "Non disponible", "no disponible")
    nnd <- paste0("data/", nn, ".ods")
    nn <- read_ods(nnd, na = nna) |>
      clean_names() |>
      mutate_if(is.character, as.factor)
    bn <- read_ods(nnd, sheet = 2)
    var_label(nn) <- bn$bnom
    return(nn)
  }

  fich <- c("atcd", "demog", "finet", "patho", "result")
  for (f in fich) {
    assign(f, importr(f))
  }

  # Demog

  demog <- demog |>
    mutate(naisdte = dmy(str_replace(naisdte, "ND", "15"))) |>
    mutate(incldte = dmy(as.character(incldte))) |>
    mutate(age = as.numeric(incldte - naisdte) / 365.25) |>
    relocate(age, .after = "naisdte") |>
    mutate(imc = bmiph(imc)) |>
    mutate(imc = fct_recode(imc,
      "dénutrition/maigreur" = "dénutrition",
      "dénutrition/maigreur" = "maigreur"
    )) |>
    dplyr::select(-naisdte, -incldte)
  bn <- read_ods("data/demog.ods", sheet = 2)
  var_label(demog) <- bn$nom

  # ATCD

  atcd <- atcd |>
    mutate(tabacon = fct_relevel(
      tabacon,
      "Aucun", "Actif", "Sevré"
    )) |>
    mutate(tabacpa = cut(tabacpa,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(0, 1, 20, 40, 100),
      labels = c("0", "1-19", "20-39", "> 39")
    ))
  #
  bn <- read_ods("data/atcd.ods", sheet = 2)
  var_label(atcd) <- bn$nom
  tt1 <- left_join(demog, atcd, by = "subjid")


  # Patho

  difdate <- function(dd1, hh1, dd2, hh2) {
    zz1 <- as.numeric(dmy_hms(paste0(dd1, " ", hh1)))
    zz2 <- as.numeric(dmy_hms(paste0(dd2, " ", hh2)))
    dif <- (zz2 - zz1) / 60
    return(dif)
  }

  patho <- patho |>
    mutate(evah3dte = as.character(evah3dte)) |>
    mutate(sympt_urg = difdate(symptodte, symptohr, urgencdte, urgenhr)) |>
    mutate(urg_h3 = difdate(urgencdte, urgenhr, evah3dte, evah3hr)) |>
    dplyr::select(-c(symptodte, symptohr, urgencdte, urgenhr, urgencdte, urgenhr, evah3dte, evah3hr))
  bn <- read_ods("data/patho.ods", sheet = 2)
  var_label(patho) <- bn$nom
  tt2 <- left_join(tt1, patho, by = "subjid")
  # Résultats

  bn <- read_ods("data/result.ods", sheet = 2)
  var_label(result) <- bn$nom
  tt3 <- left_join(tt2, result, by = "subjid") |>
    mutate(sca3 = ifelse((tropoh3 > 34.2 & sex == "Masculin") | (tropoh3 > 15.6 & sex == "Féminin"), "Oui", "Non")) |>
    mutate(sca3 = as.factor(sca3)) |>
    mutate(sca3 = fct_relevel(sca3, "Oui", "Non")) |>
    drop_na(sca3) |>
    mutate(tp0 = ifelse((tropoh0 > 15.6 & sex == "Féminin") |
      (tropoh0 > 34.2 & sex == "Masculin"),
    "Oui", "Non"
    )) |>
    mutate(tp0 = as.factor(tp0)) |>
    mutate(tp0 = fct_relevel(tp0, "Oui", "Non")) |>
    mutate(cp0 = ifelse(copepth0 > 10, "Oui", "Non")) |>
    mutate(cp0 = as.factor(cp0)) |>
    mutate(cp0 = fct_relevel(cp0, "Oui", "Non")) |>
    mutate(tpcp0 = ifelse(tp0 == "Oui" | cp0 == "Oui", "Oui", "Non")) |>
    mutate(tpcp0 = as.factor(tpcp0)) |>
    mutate(tpcp0 = fct_relevel(tpcp0, "Oui", "Non"))
  #
  var_label(tt3$tp0) <- "Troponine h0 anormale"
  var_label(tt3$cp0) <- "Copeptine h0 anormale"
  var_label(tt3$tpcp0) <- "Troponine ou copeptine h0 anormale"
  var_label(tt3$sca3) <- "SCA ST-"
  #
  # finet
  zz1 <- as.numeric(dmy_hms(paste0(finet$sortiurgdte, " ", finet$sortiurghr)))
  zz2 <- as.numeric(dmy_hms(paste0(patho$urgencdte, " ", patho$urgenhr)))
  zz <- (zz1 - zz2) / 3600
  finet <- finet |>
    mutate(duree_urg = zz)

  #
  bn <- read_ods("data/finet.ods", sheet = 2)
  var_label(finet) <- bn$nom

  tt <- left_join(tt3, finet, by = "subjid") |>
    mutate(
      scanonston =
        fct_relevel(
          scanonston,
          "Oui", "Non"
        )
    )
  #
  save(atcd, demog, finet, patho, tt, file = "data/copsca.RData")
}


importph()
load("data/copsca.RData")
