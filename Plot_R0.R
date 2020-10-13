#' Plot transmission events by category.
#'
#'
#' @param dat tracking.inf as provided as output from the CWD model functions
#'
#' @return a plot of infection events by sex-based transmission category
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom forcats fct_recode fct_reorder
#' @examples
#' params <- list(fawn.an.sur = 0.6, juv.an.sur = 0.8, ad.an.f.sur = 0.95, 
#' ad.an.m.sur = 0.9, fawn.repro = 0, juv.repro = 0.6, ad.repro = 1, 
#' hunt.mort.fawn = 0.01, hunt.mort.juv.f = 0.1, hunt.mort.juv.m = 0.1,
#' hunt.mort.ad.f = 0.2, hunt.mort.ad.m = 0.2, ini.fawn.prev = 0.02,
#' ini.juv.prev = 0.03, ini.ad.f.prev = 0.04,  ini.ad.m.prev = 0.04,
#' n.age.cats = 12,  p = 0.43, env.foi = 0,  beta.ff = 0.08, 
#' gamma.mm = 1, gamma.mf = 1, gamma.fm = 1.5,
#' theta = 1, n0 = 2000, n.years = 15, rel.risk = 1.0)
#' 
#' out <- cwd_wiw_track_R0(params)
#' plot_R0(out$R0)
#' 
#' @export

plot_R0 <- function (dat) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(dat[[1]],
                        mapping = aes(x = 0, 
                                      y = value, 
                                      shape = name),
                        size = 2) +
    ggplot2::geom_linerange(data = dat[[2]],
                            mapping = aes(x = 0, 
                                          ymin = min, 
                                          ymax = max,
                                          linetype = name),
                            stat = "identity",
                            position = "identity",
                            inherit.aes = F) +
    ggplot2::scale_linetype_discrete(name  ="Approximation",
                                     labels="Growth Rate") +
    ggplot2::scale_shape_discrete(name  ="Approximation",
                                  labels=c("Next Generation","Final Prevalence"))+
    ggplot2::labs(y = "R0") +
    ggplot2::theme_classic()+
    ggplot2::theme(axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank())
  p
}
