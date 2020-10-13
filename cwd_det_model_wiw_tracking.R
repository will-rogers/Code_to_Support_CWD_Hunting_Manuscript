#' CWD who-infects-who deterministic model 
#' 
#' Deterministic monthly age and sex structured model with constant
#'  environmental transmission and dynamic direct transmission. 2x2 matrix of 
#'  transmission rates between males and females. 
#'  
#' @param params A list with the following parameters included: 
#' 
#' fawn.an.sur = annual fawn survival (scaler value between 0 and 1),  
#' 
#' juv.an.sur = annual juvenile survival (scaler value between 0 and 1),  
#' 
#' ad.an.f.sur = annual adult female survival (scaler value between 0 and 1),  
#' 
#' ad.an.m.sur = annual adult male survival (scaler value between 0 and 1),   
#' 
#' fawn.repro = fawn reproduction (scaler value >= 0),  
#' 
#' juv.repro = juvenile reproduction (scaler value >= 0),  
#' 
#' ad.repro = adult reproduction (scaler value >= 0),  
#' 
#' hunt.mort.fawn = percentage of fawns hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.juv.f = percentage of juvenile females hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.juv.m = percentage of juvenile males hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.ad.f = percentage of adult females hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.ad.m = percentage of adult males hunted (scaler value between 0 and 1),  
#' 
#' ini.fawn.prev = percentage of fawns infected at the start (scaler value between 0 and 1),  
#' 
#' ini.juv.prev = percentage of juveniles infected at the start (scaler value between 0 and 1),  
#' 
#' ini.ad.f.prev = percentage of adult females infected at the start (scaler value between 0 and 1),  
#' 
#' ini.ad.m.prev = percentage of adult males infected at the start (scaler value between 0 and 1),   
#' 
#' n.age.cats = number of age categories to monitor. (scaler value greater than 3).
#' The final age category includes all those of that age or greater. 
#' 
#' p = rate of movement in the infectious categories (scaler values between 0 and 1). 
#' See model documentation vignette for how this relates to disease induced mortality.   
#' 
#' env.foi = % of the population that is infected by the environment per month (scaler value between 0 and 1) #' 
#' beta.ff = female to female transmission coefficient (scaler value greater than 0),  
#' 
#' gamma.mm = relative increase in the male-male transmission coefficient compared 
#' to female-female transmission (scaler value greater than 0. A value of 1 indicates equal transmission)  
#' 
#' gamma.mf = relative increase in the male-female transmission coefficient compared 
#' to female-female transmission (scaler value greater than 0. A value of 1 indicates equal transmission)  
#' 
#' gamma.fm = relative increase in the female-male transmission coefficient compared 
#' to female-female transmission (scaler value greater than 0. A value of 1 indicates equal transmission) 
#' 
#' theta = effect of population size on transmission (1 = frequency dependence, 0 = density dependent).  
#' 
#' n0 = initial population size (scaler value greater than 0)
#' 
#' n.years = number of years to run the model (scaler value greater than 2),  
#' 
#' rel.risk = relative risk of infected individuals being hunted. A value of 1 
#' indicates no hunter preference for infected individuals  
#' 
#' 
#' @return A list with 3 dataframes is returned as output: 
#' 
#' 1. counts of the # of individuals in the susceptible and infectious 
#' categories by over time. 
#' 
#'  Columns include: 
#' 
#'  age (in years)
#' 
#'  month of simulation,
#' 
#'  population = number of individuals
#' 
#'  category: St.f = susceptible females, St.m = susceptible males, Ixt.f = 
#'  infectious females in the x category (1-10), Ixt.m = infectious males in the 
#'  x infectious category (1-10) 
#' 
#'  sex = female or male
#' 
#' disease = yes or no for susceptible or infectious
#' 
#'   
#' 2. deaths--how individuals died over time (hunting, natural or disease).
#' 
#'  Columns include: 
#'  
#'  age in years, 
#'  
#'  month of the simulation, 
#'  
#'  population = # of individuals, 
#'  
#'  category: Ht.f = hunted females, Ht.m = hunted males, Dt.f = natural 
#'  mortality females, Dt.m = natural mortality males, CWDt.f = disease mortality 
#'  females, CWDt.m = disease mortality males.  
#'    
#'  year = year of the simulation
#'  
#'  sex
#'  
#' 3. tracking.inf--infection tracting over time (sex-based contraction, sex-based transmission).
#' 
#'  Columns include: 
#'  
#'  age in years, 
#'  
#'  month of the simulation, 
#'  
#'  population = # of events, 
#'  
#'  category: new.f.inf = new female infections, 
#'  new.m.inf new male infections, mm.cases = male to male transmission events,
#'  mf.cases = male to female transmission events, 
#'  ff.cases = female to female transmission events,
#'  fm.cases = female to male transmission events,
#'  em.cases = environment to male transmission events, 
#'  ef.cases = environment to female transmission events.  
#'    
#'  year = year of the simulation
#'
#' @importFrom popbio stable.stage
#' @importFrom dplyr rename mutate
#' @importFrom reshape2 melt
#' @examples 
#' params <- list(fawn.an.sur = 0.6, juv.an.sur = 0.8, ad.an.f.sur = 0.95, 
#' ad.an.m.sur = 0.9, fawn.repro = 0, juv.repro = 0.6, ad.repro = 1, 
#' hunt.mort.fawn = 0.01, hunt.mort.juv.f = 0.1, hunt.mort.juv.m = 0.1,
#' hunt.mort.ad.f = 0.1, hunt.mort.ad.m = 0.2, ini.fawn.prev = 0.02,
#' ini.juv.prev = 0.03, ini.ad.f.prev = 0.04,  ini.ad.m.prev = 0.04,
#' n.age.cats = 12,  p = 0.43, env.foi = 0,  beta.ff = 0.06, 
#' gamma.mm = 2, gamma.mf = 2, gamma.fm = 1,
#' theta = 1, n0 = 2000, n.years = 10, rel.risk = 1.0)
#' 
#' out <- cwd_det_model_wiw_tracking(params)
#' 
#' plot_track_sex_trans(out$tracking.inf)
#' 
#' @export
cwd_det_model_wiw_tracking <- function (params) {
  for (v in 1:length(params)) assign(names(params)[v], params[[v]])
  if (exists("fawn.an.sur") == FALSE) {
    message("fawn survival is missing, using default value")
    fawn.an.sur <- 0.6
  }
  if (exists("juv.an.sur") == FALSE) {
    message("juvenile survival is missing, using default value")
    juv.an.sur <- 0.8
  }
  if (exists("ad.an.f.sur") == FALSE) {
    message("adult female survival is missing, using default value")
    ad.an.f.sur <- 0.95
  }
  if (exists("ad.an.m.sur") == FALSE) {
    message("adult male survival is missing, using default value")
    ad.an.m.sur <- 0.9
  }
  if (exists("fawn.repro") == FALSE) {
    message("fawn repro is missing, using default value")
    fawn.repro <- 0
  }
  if (exists("juv.repro") == FALSE) {
    message("juvenile repro is missing, using default value")
    juv.repro <- 0.6
  }
  if (exists("ad.repro") == FALSE) {
    message("adult repro is missing, using default value")
    ad.repro <- 1
  }
  if (exists("hunt.mort.fawn") == FALSE) {
    message("fawn hunting mortality is missing, using default value")
    hunt.mort.fawn <- 0.01
  }
  if (exists("hunt.mort.juv.f") == FALSE) {
    message("juv. female  hunting mortality is missing, using default value")
    hunt.mort.juv.f <- 0.1
  }
  if (exists("hunt.mort.juv.m") == FALSE) {
    message("juv. male hunting mortality is missing, using default value")
    hunt.mort.juv.m <- 0.1
  }
  if (exists("hunt.mort.ad.f") == FALSE) {
    message("adult female hunting mortality is missing, using default value")
    hunt.mort.ad.f <- 0.2
  }
  if (exists("hunt.mort.ad.m") == FALSE) {
    message("adult male hunting mortality is missing, using default value")
    hunt.mort.ad.m <- 0.2
  }
  if (exists("ini.fawn.prev") == FALSE) {
    message("initial fawn prevalence is missing, using default value")
    ini.fawn.prev <- 0.01
  }
  if (exists("ini.juv.prev") == FALSE) {
    message("initial juvenile prevalence is missing, using default value")
    ini.juv.prev <- 0.03
  }
  if (exists("ini.ad.f.prev") == FALSE) {
    message("initial adult female prevalence is missing, using default value")
    ini.ad.f.prev <- 0.04
  }
  if (exists("ini.ad.m.prev") == FALSE) {
    message("initial adult male prevalence is missing, using default value")
    ini.ad.m.prev <- 0.04
  }
  if (exists("n.age.cats") == FALSE) {
    message("# of age categories is missing, using default value")
    n.age.cats <- 12
  }
  if (exists("p") == FALSE) {
    message("disease mortality index p is missing, using default value")
    p <- 0.43
  }
  if (exists("beta.ff") == FALSE) {
    message("female transmission beta.ff is missing, using default value")
    beta.f <- 0.05
  }
  if (exists("theta") == FALSE) {
    message("theta is missing, using default value")
    theta <- 1
  }
  if (exists("n0") == FALSE) {
    message("initial population size n0 is missing, using default value")
    n0 <- 1000
  }
  if (exists("n.years") == FALSE) {
    message("n.years is missing, using default value")
    n.years <- 10
  }
  if (exists("rel.risk") == FALSE) {
    message("rel.risk is missing, using default value")
    rel.risk <- 1
  }
  if (exists("gamma.mm") == FALSE) {
    message("gamma.mm is missing, using default value")
    gamma.mm <- 2
  }
  if (exists("gamma.fm") == FALSE) {
    message("gamma.fm is missing, using default value")
    gamma.fm <- 1
  }
  if (exists("gamma.mf") == FALSE) {
    message("gamma.mf is missing, using default value")
    gamma.mf <- 1
  }
  if (fawn.an.sur <= 0) 
    warning("fawn survival must be positive")
  if (fawn.an.sur > 1) 
    warning("fawn survival must be <= 1")
  if (juv.an.sur <= 0) 
    warning("juvenile survival must be positive")
  if (juv.an.sur > 1) 
    warning("juvenile survival must be <= 1")
  if (ad.an.f.sur <= 0) 
    warning("adult female survival must be positive")
  if (ad.an.f.sur > 1) 
    warning("adult female survival must be <= 1")
  if (fawn.repro < 0) 
    warning("fawn.repro must be positive")
  if (juv.repro <= 0) 
    warning("juv.repro must be >= 0 ")
  if (ad.repro <= 0) 
    warning("ad.repro must be >= 0 ")
  if (hunt.mort.fawn < 0) 
    warning("hunt.mort.fawn must be >=0")
  if (hunt.mort.fawn > 1) 
    warning("hunt.mort.fawn must be < 1")
  if (hunt.mort.juv.f < 0) 
    warning("hunt.mort.juv.f must be >=0")
  if (hunt.mort.juv.f > 1) 
    warning("hunt.mort.juv.f must be < 1")
  if (hunt.mort.juv.m < 0) 
    warning("hunt.mort.juv.m must be >=0")
  if (hunt.mort.juv.m > 1) 
    warning("hunt.mort.juv.m must be < 1")
  if (hunt.mort.ad.f < 0) 
    warning("hunt.mort.ad.f must be >=0")
  if (hunt.mort.ad.f > 1) 
    warning("hunt.mort.ad.f must be < 1")
  if (hunt.mort.ad.m < 0) 
    warning("hunt.mort.ad.m must be >=0")
  if (hunt.mort.ad.m > 1) 
    warning("hunt.mort.ad.m must be < 1")
  if (ini.fawn.prev < 0) 
    warning("ini.fawn.prev must >=0")
  if (ini.fawn.prev > 1) 
    warning("ini.fawn.prev must be <= 1")
  if (ini.juv.prev < 0) 
    warning("ini.juv.prev must >=0")
  if (ini.juv.prev > 1) 
    warning("ini.juv.prev must be <= 1")
  if (ini.ad.f.prev < 0) 
    warning("ini.ad.f.prev must >=0")
  if (ini.ad.f.prev > 1) 
    warning("ini.ad.f.prev must be <= 1")
  if (ini.ad.m.prev < 0) 
    warning("ini.ad.m.prev must >=0")
  if (ini.ad.m.prev > 1) 
    warning("ini.ad.m.prev must be <= 1")
  if (n.age.cats < 3) 
    warning("n.age.cats must be 3 or more")
  if (p < 0) 
    warning("p must be between 0 and 1")
  if (p > 1) 
    warning("p must be between 0 and 1")
  if (env.foi < 0) 
    warning("env.foi must be between 0 and 1")
  if (env.foi > 1) 
    warning("env.foi must be between 0 and 1")
  if (beta.ff < 0) 
    warning("beta.ff cannot be negative")
  if (n0 <= 0) 
    warning("n0 must be positive")
  if (n.years <= 0) 
    warning("n.years must be positive")
  if (rel.risk <= 0) 
    warning("n.years must be positive")
  if (gamma.mm <= 0) 
    warning("repro.var must be positive")
  if (gamma.fm <= 0) 
    warning("fawn.sur.var must be positive")
  if (gamma.mf <= 0) 
    warning("sur.var must be positive")
  #Disease
  beta.mm <- beta.ff * gamma.mm #male-male transmission
  beta.mf <- beta.ff * gamma.mf #male-female transmission
  beta.fm <- beta.ff * gamma.fm #female-male transmission
  months <- seq(1, n.years * 12) #total months of sim
  hunt.mo <- rep(0, n.years * 12) #establishing when hunting occurs
  hunt.mo[months%%12 == 7] <- 1 #only making the 7th time period when hunting occurs
  fawn.sur <- fawn.an.sur^(1/12) #monthly fawn survival
  juv.sur <- juv.an.sur^(1/12) #"" juvenile
  ad.f.sur <- ad.an.f.sur^(1/12) #"" female
  ad.m.sur <- ad.an.m.sur^(1/12) #"" male
  ini.f.prev <- c(ini.fawn.prev, ini.juv.prev, rep(ini.ad.f.prev, 
                                                   (n.age.cats - 2))) # initial vector of female prevalence
  ini.m.prev <- c(ini.fawn.prev, ini.juv.prev, rep(ini.ad.m.prev, 
                                                   (n.age.cats - 2))) # initial vector of male prevalence
  Sur.f <- c(fawn.sur, juv.sur, rep(ad.f.sur, n.age.cats - 
                                      2)) # initial vector of female survival
  Sur.m <- c(fawn.sur, juv.sur, rep(ad.m.sur, n.age.cats - 
                                      2)) # initial vector of male survival
  M <- matrix(rep(0, n.age.cats * 2 * n.age.cats * 2), nrow = n.age.cats * 
                2) #making blank matrix
  M[row(M) == (col(M) + 1)] <- c(juv.an.sur * (1 - hunt.mort.juv.f), 
                                 rep(ad.an.f.sur * (1 - hunt.mort.ad.f), n.age.cats - 
                                       2), 0, c(juv.an.sur * (1 - hunt.mort.juv.m), rep(ad.an.m.sur * 
                                                                                          (1 - hunt.mort.ad.m), n.age.cats - 2))) # making off-diagonal
  M[n.age.cats, n.age.cats] <- ad.an.f.sur * (1 - hunt.mort.ad.f) #making survival in last age class non-zero
  M[n.age.cats * 2, n.age.cats * 2] <- ad.an.m.sur * (1 - hunt.mort.ad.m)  #same for males
  M[1, 1:n.age.cats] <- c(0, juv.repro, rep(ad.repro, n.age.cats - 
                                              2)) * 0.5 * fawn.an.sur * (1 - hunt.mort.fawn) #post-breeding fecundity for females only
  M[n.age.cats + 1, 1:n.age.cats] <- M[1, 1:n.age.cats] # copying and pasting, though no one inhabits applicable classes here
  tmp <- matrix(0, nrow = n.age.cats, ncol = n.years * 12) #making a temporary matrix 
  St.f <- tmp #for susceptible females
  St.m <- tmp #for susceptible males
  It.m <- array(rep(tmp), dim = c(n.age.cats, n.years * 12, 
                                  10)) # making this single matrix applicable for 10 disease categories for females
  It.f <- array(rep(tmp), dim = c(n.age.cats, n.years * 12, 
                                  10)) # making this single matrix applicable for 10 disease categories for males
  Ht.f <- tmp # for hunted females
  Ht.m <- tmp # for hunted males
  Dt.f <- tmp # for naturally killed females
  Dt.m <- tmp # for naturally killed males
  CWDt.f <- tmp # for cwd killed females
  CWDt.m <- tmp # for cwd killed females 
  
  ###################
  #### Amendment ####
  f.inf <- tmp
  m.inf <- tmp
  mm.cases <- tmp
  mf.cases <- tmp
  ff.cases <- tmp
  fm.cases <- tmp
  em.cases <- tmp
  ef.cases <- tmp
  ###################
  
  St.f[, 1] <- popbio::stable.stage(M)[1:n.age.cats] * n0 * 
    (1 - ini.f.prev) #filling the first matrix of females for susceptible and stable age distribution from M
  St.m[, 1] <- popbio::stable.stage(M)[(n.age.cats + 1):(n.age.cats * 
                                                           2)] * n0 * (1 - ini.m.prev) #filling the first matrix of males for susceptible and stable age distribution from M
  It.m[, 1, 1:10] <- popbio::stable.stage(M)[1:n.age.cats] * 
    n0/10 * ini.m.prev #filling the first matrix of females for infectious and stable age distribution from M
  It.f[, 1, 1:10] <- popbio::stable.stage(M)[(n.age.cats + 
                                                1):(n.age.cats * 2)] * n0/10 * ini.f.prev #filling the first matrix of males for infectious and stable age distribution from M
  for (t in 2:(n.years * 12)) {
    if (t%%12 == 2) { # when the month is 2
      St.f[2:(n.age.cats - 1), t] <- St.f[1:(n.age.cats - 
                                               2), t - 1] # individuals move to the next age class in the second time period
      St.f[n.age.cats, t] <- St.f[n.age.cats, t - 1] + 
        St.f[(n.age.cats - 1), t - 1] # individuals in the final age class remain plus those from 11
      St.m[2:(n.age.cats - 1), t] <- St.m[1:(n.age.cats - 
                                               2), t - 1] # ""
      St.m[n.age.cats, t] <- St.m[n.age.cats, t - 1] + 
        St.m[(n.age.cats - 1), t - 1] # ""
      It.f[2:(n.age.cats - 1), t, ] <- It.f[1:(n.age.cats - 
                                                 2), t - 1, ] # for infected, respective to disease state
      It.f[n.age.cats, t, ] <- It.f[n.age.cats, t - 1, 
                                    ] + It.f[(n.age.cats - 1), t - 1, ] #""
      It.m[2:(n.age.cats - 1), t, ] <- It.m[1:(n.age.cats - 
                                                 2), t - 1, ] #""
      It.m[n.age.cats, t, ] <- It.m[n.age.cats, t - 1, 
                                    ] + It.m[(n.age.cats - 1), t - 1, ] #""
      I_juv <- sum(It.f[2, t - 1, ]) #finding a total for infectious juveniles
      I_adults <- sum(It.f[3:n.age.cats, t - 1, ]) #finding a total for infectious adults
      St.f[1, t] <- ((St.f[2, t - 1] + I_juv) * juv.repro + 
                       (sum(St.f[3:n.age.cats, t - 1]) + I_adults) * 
                       ad.repro) * 0.5 #infected and susceptible adults and juveniles reproduce to make suscpeptible fawns, 50% are female
      St.m[1, t] <- ((St.f[2, t - 1] + I_juv) * juv.repro + 
                       (sum(St.f[3:n.age.cats, t - 1]) + I_adults) * 
                       ad.repro) * 0.5 #"", 50% are male
    }
    if (t%%12 != 2) { # if the month is not 2, just carry over classes to the next month
      St.f[, t] <- St.f[, t - 1] 
      St.m[, t] <- St.m[, t - 1]
      It.f[, t, ] <- It.f[, t - 1, ]
      It.m[, t, ] <- It.m[, t - 1, ]
    }
    St.f[, t] <- St.f[, t] * Sur.f # susceptible females survive at a monthly, age specific rate
    St.m[, t] <- St.m[, t] * Sur.m # susceptible males""
    It.f[, t, ] <- It.f[, t, ] * Sur.f #infected females ""
    It.m[, t, ] <- It.m[, t, ] * Sur.m #infected males ""
    Dt.f[, t] <- (St.f[, t] + rowSums(It.f[, t, ])) * (1 - Sur.f) # the monthly deaths of susceptible and infections females
    Dt.m[, t] <- (St.m[, t] + rowSums(It.m[, t, ])) * (1 - Sur.m) # "" males
    if (hunt.mo[t] == 1) { # in the 7th month, remember the setting for hunting month above^^^
      Iall.f <- rowSums(It.f[, t, ]) # all infectious females
      Iall.m <- rowSums(It.m[, t, ]) # all infectious males
      Nt.f <- St.f[, t] + Iall.f # total pop females
      Nt.m <- St.m[, t] + Iall.m # total pop males
      hunted.f <- Nt.f * c(hunt.mort.fawn, hunt.mort.juv.f, 
                           rep(hunt.mort.ad.f, n.age.cats - 2)) # hunted females by age class
      hunted.m <- Nt.m * c(hunt.mort.fawn, hunt.mort.juv.m, 
                           rep(hunt.mort.ad.m, n.age.cats - 2)) # hunted males by age class
      Ht.f[, t] <- hunted.f # filling in the matrix for females
      Ht.m[, t] <- hunted.m # for males
      hunted.i.f <- (rel.risk * Iall.f * hunted.f)/(St.f[, t] + rel.risk * Iall.f) # finding how many hunted cases are infected females
      hunted.i.m <- (rel.risk * Iall.m * hunted.m)/(St.m[, t] + rel.risk * Iall.m) # finding how many hunted cases are infected males
      hunted.i.f[which(is.na(hunted.i.f))] <- 0 # if a given scenario fails
      hunted.i.m[which(is.na(hunted.i.m))] <- 0 # ""
      St.f[, t] <- St.f[, t] - (hunted.f - hunted.i.f) # hunted population of susceptible females to carry forward
      St.m[, t] <- St.m[, t] - (hunted.m - hunted.i.m) # hunted population of susceptible males to carry forward
      It.f[, t, ] <- It.f[, t, ] * (1 - hunted.i.f/Iall.f) # hunted population of infected females to carry forward
      It.m[, t, ] <- It.m[, t, ] * (1 - hunted.i.m/Iall.m) # hunted population of infected males to carry forward
    }
    f.move <- It.f[, t, ] * p # establishing the monthly movement between disease classes by monthly rate p
    m.move <- It.m[, t, ] * p # establishing the monthly movement between disease classes by monthly rate p
    It.f[, t, 1] <- It.f[, t, 1] - f.move[, 1] # removing from the first disease class
    It.f[, t, 2:10] <- It.f[, t, 2:10] - f.move[, 2:10] + 
      f.move[, 1:9] # removing for 2-10 and adding from 1-9, for 2-10
    It.m[, t, 1] <- It.m[, t, 1] - m.move[, 1] # for males
    It.m[, t, 2:10] <- It.m[, t, 2:10] - m.move[, 2:10] + 
      m.move[, 1:9] # for females
    CWDt.f[, t] <- f.move[, 10] # cwd related mortality for progressing past last disease state
    CWDt.m[, t] <- m.move[, 10] # for males
    Iall <- sum(It.f[, t, ] + It.m[, t, ]) # all infectious males and females in all disease categories for a given month
    Nall <- sum(St.f[, t] + St.m[, t]) + Iall # total population
    cases.f <- St.f[, t] * (1 - exp(-(beta.ff * sum(It.f[, t, ])/Nall^theta + beta.mf * sum(It.m[, t, ])/Nall^theta))) # new cases of female infection
    cases.m <- St.m[, t] * (1 - exp(-(beta.mm * sum(It.m[, t, ])/Nall^theta + beta.fm * sum(It.f[, t, ])/Nall^theta))) # of males
    
    ###################
    #### Amendment ####
    mm.cases[, t] <- St.m[, t] * (1 - exp(-(beta.mm * sum(It.m[, t, ])/Nall^theta)))
    mf.cases[, t] <- St.f[, t] * (1 - exp(-(beta.mf * sum(It.m[, t, ])/Nall^theta)))
    ff.cases[, t] <- St.f[, t] * (1 - exp(-(beta.ff * sum(It.f[, t, ])/Nall^theta)))
    fm.cases[, t] <- St.m[, t] * (1 - exp(-(beta.fm * sum(It.f[, t, ])/Nall^theta)))
    em.cases[, t] <- envcases.m
    ef.cases[, t] <- envcases.f
    ###################
    
    St.f[, t] <- St.f[, t] - cases.f # removing newly infected females
    St.m[, t] <- St.m[, t] - cases.m # for males 
    It.f[, t, 1] <- It.f[, t, 1] + cases.f # adding newly infected females
    It.m[, t, 1] <- It.m[, t, 1] + cases.m # for males
    envcases.f <- St.f[, t] * env.foi # environmental cases of females
    envcases.m <- St.m[, t] * env.foi # of males
    
    ###################
    #### Amendment ####
    f.inf[, t] <- cases.f + envcases.f
    m.inf[, t] <- cases.m + envcases.m
    em.cases[, t] <- envcases.m
    ef.cases[, t] <- envcases.f
    ###################
    
    St.f[, t] <- St.f[, t] - envcases.f # removing environmentally infected females
    St.m[, t] <- St.m[, t] - envcases.m # males
    It.f[, t, 1] <- It.f[, t, 1] + envcases.f # adding environmentally infected females
    It.m[, t, 1] <- It.m[, t, 1] + envcases.m # males
  }
  counts <- list(St.f = St.f, St.m = St.m, I1t.f = It.f[, , 
                                                        1], I1t.m = It.m[, , 1], I2t.f = It.f[, , 2], I2t.m = It.m[, 
                                                                                                                   , 2], I3t.f = It.f[, , 3], I3t.m = It.m[, , 3], I4t.f = It.f[,, 4], 
                 I4t.m = It.m[, , 4], I5t.f = It.f[, , 5], I5t.m = It.m[, 
                                                                        , 5], I6t.f = It.f[, , 6], I6t.m = It.m[, , 6], I7t.f = It.f[, 
                                                                                                                                     , 7], I7t.m = It.m[, , 7], I8t.f = It.f[, , 8], I8t.m = It.m[, 
                                                                                                                                                                                                  , 8], I9t.f = It.f[, , 9], I9t.m = It.m[, , 9], I10t.f = It.f[, 
                                                                                                                                                                                                                                                                , 10], I10t.m = It.m[, , 10]) # tracking how many individuals are in susceptible and infectious classes
  deaths <- list(Ht.f = Ht.f, Ht.m = Ht.m, Dt.f = Dt.f, Dt.m = Dt.m, 
                 CWDt.f = CWDt.f, CWDt.m = CWDt.m) # tracking cause of death
  ###################
  #### Amendment ####
  tracking.inf <- list(new.f.inf = f.inf, 
                       new.m.inf = m.inf, 
                       mm.cases = mm.cases,
                       mf.cases = mf.cases,
                       ff.cases = ff.cases,
                       fm.cases = fm.cases,
                       em.cases = em.cases,
                       ef.cases = ef.cases)
  ###################
  counts.long <- reshape::melt(counts) %>% 
    dplyr::rename(age = X1, month = X2, population = value, category = L1) %>% 
    mutate(year = (month - 1)/12, sex = as.factor(str_sub(category, -1)), disease = "no") #compiling based on year, sex, and susceptible
  counts.long$disease[str_sub(counts.long$category, 1, 1) == "I"] <- "yes" #labeling infected population
  counts.long$disease <- as.factor(counts.long$disease) 
  deaths.long <- reshape::melt(deaths) %>% 
    dplyr::rename(age = X1, month = X2, population = value, category = L1) %>% 
    mutate(year = (month - 1)/12, sex = as.factor(str_sub(category, -1))) # again for death
  
  ###################
  #### Amendment ####
  tracking.inf.long <- reshape::melt(tracking.inf) %>% 
    dplyr::rename(age = X1, month = X2, population = value, category = L1) %>%
    mutate(year = (month - 1)/12)
  ###################
  
  output <- list(counts = counts.long, deaths = deaths.long, tracking.inf = tracking.inf.long)
}

