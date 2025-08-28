tbssb <- function(ll, sca = tt$scanonston) {
  taa <- table(ll, sca)
  naa <- sum(taa)
  aa <- epi.tests(taa)[1]$detail |>
    as_tibble() |>
    dplyr::slice(3, 4, 9, 10) |>
    mutate(res = paste0(npc(est), " % [", npc(lower), ",", npc(upper), "]")) |>
    dplyr::select(statistic, res) |>
    pivot_wider(names_from = statistic, values_from = res) |>
    mutate(n = naa) |>
    relocate(n, .before = se)
  tabx <- rbind(tabx, aa)
  return(tabx)
}
