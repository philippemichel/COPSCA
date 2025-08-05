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
    nna <- c("NA", "", " ", "C", "D", "K", "A", "Non disponible", "no disponible", "NK")
    nnd <- paste0("data/", nn, ".ods")
    nn <- read_ods(nnd, na = nna) |>
      clean_names() |>
      mutate_if(is.character, as.factor)
    bn <- read_ods(nnd, sheet = 2)
    var_label(nn) <- bn$nom
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
  var_label(demog$age) <- "Âge"
  var_label(demog$imc) <- "IMC"


  # ATCD

  atcd <- atcd |>
    mutate(tabacon = fct_relevel(
      tabacon,
      "Non", "Actif", "Sevré"
    )) |>
    mutate(tabacpa = cut(tabacpa,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(0, 1, 20, 40, 100),
      labels = c("0", "1-19", "20-39", "> 39")
    ))
  var_label(atcd$tabacpa) <- "Paquets-années"
  #
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
    mutate(urg_h3 = difdate(urgencdte, urgenhr, evah3dte, evah3hr))
  var_label(patho$sympt_urg) <- "Délai premiers symptômes/urgences"
  var_label(patho$urg_h3) <- "Délai arrivée urgence/évaluation à h3"

  tt2 <- left_join(tt1, patho, by = "subjid")
  # Résultats


  tt3 <- left_join(tt2, result, by = "subjid") |>
    mutate(sca3 = ifelse((tropoh3 > 34.2 & sex == "Masculin") | (tropoh3 > 15.6 & sex == "Féminin"), "Élevé", "Normal")) |>
    mutate(sca3 = as.factor(sca3)) |>
    mutate(tp0 = ifelse((tropoh0 > 15.6 & sex == "Féminin") |
      (tropoh0 > 34.2 & sex == "Masculin"),
    "Élevée", "Normale"
    )) |>
    mutate(cp0 = ifelse(copepth0 > 10, "Élevée", "Normale")) |>
    mutate(cp0 = as.factor(cp0)) |>
    mutate(tpcp0 = ifelse(tp0 == "Élevée" | cp0 == "Élevée", "Positif", "Négatif")) |>
    mutate(tpcp0 = as.factor(tpcp0)) |>
    mutate(tpcp0 = fct_relevel(tpcp0, "Positif", "Négatif"))
  #
  var_label(tt3$tp0) <- "Troponine h0"
  var_label(tt3$cp0) <- "Copeptine h0"
  var_label(tt3$tpcp0) <- "combinaison Troponine/Copeptine"
  var_label(tt3$sca3) <- "Troponine H3"
  #
  # finet

  #


  tt <- left_join(tt3, finet, by = "subjid") |>
    mutate(
      scanonston =
        fct_relevel(
          scanonston,
          "Oui", "Non"
        )
    ) |>
    mutate(duree_urg = difdate(urgencdte, urgenhr, sortiurgdte, sortiurghr)) |>
    mutate(doulprev = difdate(tt$symptodte, tt$symptohr, tt$prelh0dte, tt$prelh0hr)) |>
    mutate(ddoul = doulh3 - doulh0)
  dplyr::select(-ends_with(c("dtex", "hrx")))

  var_label(tt$duree_urg) <- "Temps passé aux urgences"
  var_label(tt$doulprev) <- "Délai premiers smptômes/prélèvement"
  var_label(tt$ddoul) <- "Évolution douleur H0-H3"

  #
  save(atcd, demog, finet, patho, tt, file = "data/copsca.RData")
}


importph()
load("data/copsca.RData")
