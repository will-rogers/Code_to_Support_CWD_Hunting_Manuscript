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
# params <- list(fawn.an.sur = 0.6, juv.an.sur = 0.8, ad.an.f.sur = 0.95,
# ad.an.m.sur = 0.9, fawn.repro = 0, juv.repro = 0.6, ad.repro = 1,
# hunt.mort.fawn = 0.01, hunt.mort.juv.f = 0.1, hunt.mort.juv.m = 0.1,
# hunt.mort.ad.f = 0.2, hunt.mort.ad.m = 0.2, ini.fawn.prev = 0.02,
# ini.juv.prev = 0.03, ini.ad.f.prev = 0.04,  ini.ad.m.prev = 0.04,
# n.age.cats = 12,  p = 0.43, env.foi = 0,  beta.ff = mat$B.ff[3],
# gamma.mm = mat$G.mm[3], gamma.mf = mat$G.mf[3], gamma.fm = mat$G.fm[3],
# theta = 1, n0 = 2000, n.years = 15, rel.risk = 1.0)
# 
# outa <- cwd_wiw_track_R0(params)
# plot_track_sex_trans(outa$tracking.inf)
#' 
#' @export
plot_track_sex_trans <- function (dat, totals) {
  if(missing(dat)==TRUE) warning("missing data to plot")
  
  if(missing(totals) == TRUE) totals = FALSE
  
  if(totals == FALSE) {
  plot <- dat %>% 
    filter(category == "mm.cases"|
             category == "mf.cases"|
             category == "fm.cases"|
             category == "ff.cases") %>% 
    mutate(year = floor(year)) %>% 
    group_by(category, year) %>% 
    summarise(population = sum(population)) %>% 
    ggplot(aes(year, population, color = category)) +
    geom_line() +
    scale_color_discrete(name = "Infection Type",
                         labels = c("F-F", "F-M", "M-F", "M-M")) +
    labs(y = "Transmission Events",
         x = "Year") +
    theme_classic()
  return(plot)
  }
  if (totals == TRUE) {
    df <- dat %>% 
      mutate(group = ifelse(category == "fm.cases"|
                              category == "ff.cases", "Female", NA)) %>% 
      mutate(group = ifelse(category == "mm.cases"|
                              category == "mf.cases", "Male", group)) %>% 
      mutate(group = ifelse(category == "em.cases"|
                              category == "ef.cases", "Environment", group))
    
    plot <- df %>% 
      group_by(group, age, month) %>% 
      mutate(n = sum(population)) %>% 
      filter(group == "Female"|
               group == "Male"|
               group == "Environment") %>% 
      mutate(year = floor(year)) %>% 
      group_by(group, year) %>% 
      summarise(n = sum(n)) %>% 
      ggplot(aes(year, n, color = group)) +
      geom_line() +
      scale_color_discrete(name = "Infection Type",
                           labels = c("E-", "F-", "M-")) +
      labs(y = "Transmission Events",
           x = "Year") +
      theme_classic()
    return(plot)
  }
}

