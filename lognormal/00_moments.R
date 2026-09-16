# log-normal moment conditions

g_lognormal <- function(tet, obs){
  y <- as.matrix(obs)[ ,1]
  g_1 <- y - exp(tet[1] + tet[2]/2)
  g_2 <- y^2 - exp(2 * tet[1] + 2 * tet[2])
  g_3 <- y^3 - exp(3 * tet[1] + 9 * tet[2] / 2)
  cbind(g_1, g_2, g_3)
}



