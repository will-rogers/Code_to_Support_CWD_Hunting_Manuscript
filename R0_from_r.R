dat<-counts.long
R0_from_r <- function(dat, p){
  I <- dat %>% 
    filter(str_detect(category, "I")) %>%
    group_by(year) %>% 
    summarize(n = sum(population))
  
  r <- max(log(na.omit(I[2:nrow(I),2]/I[1:(nrow(I)-1),2]))) #this is the max monthly exponential rate of growth
  
  ser <- (10/p) + c(-1,1) * sqrt(10/(p^2)) # these are a series of generation times that are reasonable given p
  R0.r.approx <- 1+(r*ser) #This comes from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1578275/
  return(R0.r.approx)
}
