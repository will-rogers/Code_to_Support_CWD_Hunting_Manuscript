R0_from_prev <- function(dat){
  N <- dat %>% 
    group_by(year) %>% 
    summarize(n = sum(population)) # total pop by time step
  I <- dat %>% 
    filter(str_detect(category, "I")) %>%
    group_by(year) %>% 
    summarize(n = sum(population)) # total infected by time step
  
  
  if(I$n[nrow(I)]/I$n[nrow(I)-1] >= 1.00001){
    message("WARNING: Prevalence-based R0 likely inaccurate as epidemic is still growing. Consider extending projection in time.")
  }
  
  # prevalence approximation of R0 based on max prevalence
  s.inf <- min((N$n-I$n)/N$n)
  R0.prev.approx <- log(s.inf)/(s.inf-1)
  return(R0.prev.approx)
}
R0_from_prev(out$counts)
