require(popbio)
require(dplyr)
require(reshape2)
require(stringr)
require(CWDsims)

###### Example ######

#realistic population
params <- list(fawn.an.sur = 0.354, juv.an.sur = 0.883, ad.an.f.sur = 0.972, 
               ad.an.m.sur = 0.968, fawn.repro = 0,  juv.repro = 0.198, 
               ad.repro = .928, hunt.mort.fawn = 0.01,  hunt.mort.juv.f = 0.1, 
               hunt.mort.juv.m = 0.1, hunt.mort.ad.f = 0.1,  hunt.mort.ad.m = 0.1, 
               ini.fawn.prev = 0.02, ini.juv.prev = 0.03,  ini.ad.f.prev = 0.04,  
               ini.ad.m.prev = 0.04, n.age.cats = 12,   p = 0.43,   env.foi.f = 0.001,env.foi.m = 0.01,  
               beta.ff = 0.045/2000,  gamma.mm = 1,  gamma.mf = 1,  gamma.fm = 2,
               theta = 0,  n0 = 2000,  n.years = 25, rel.risk = 1.0)

#high survival population
params <- list(fawn.an.sur = 0.6, juv.an.sur = 0.8, ad.an.f.sur = 0.95, 
               ad.an.m.sur = 0.9, fawn.repro = 0, juv.repro = 0.6, ad.repro = 1, 
               hunt.mort.fawn = 0.01, hunt.mort.juv.f = 0.1, hunt.mort.juv.m = 0.1,
               hunt.mort.ad.f = .1, hunt.mort.ad.m = .1, ini.fawn.prev = 0.02,
               ini.juv.prev = 0.03, ini.ad.f.prev = 0.04,  ini.ad.m.prev = 0.04,
               n.age.cats = 12,  p = 0.43, env.foi = 0.00,  beta.ff = 0.06/2000,
               gamma.mm = 2, gamma.mf = 1.5, gamma.fm = 1,
               theta = 0, n0 = 2000, n.years = 20, rel.risk = 1.0)
outc <- cwd_wiw_track_R0(params = params)
outc$R0
# plot_prev_time(outb$counts)
a <- outc$counts %>% 
  filter(month%%12 == 10) %>% 
  group_by(year, disease, sex) %>% 
  summarize(n = sum(population)) %>% 
  spread(key = disease, value = n) %>% 
  mutate(prev = yes/(no + yes))

ggplot(a, aes(x = year, y = prev, color = sex)) + geom_line(size = 1.5) + 
  ylim(0, 1) + ylab("Prevalence") + xlab("Year") + theme_light() + 
  theme(text = element_text(size = 18), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())

##### Function #####
cwd_wiw_track_R0 <- function (params) {
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
  if (exists("env.foi") == FALSE) {
    message("indirect transmission env.foi is missing, using default value")
    env.foi <- 0
  }
  if (exists("beta.ff") == FALSE) {
    message("female transmission beta.f is missing, using default value")
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
                                  10)) # making this single matrix applicable for 10 disease categories for males
  It.f <- array(rep(tmp), dim = c(n.age.cats, n.years * 12, 
                                  10)) # making this single matrix applicable for 10 disease categories for females
  Ht.f <- tmp # for hunted females
  Ht.m <- tmp # for hunted males
  Dt.f <- tmp # for naturally killed females
  Dt.m <- tmp # for naturally killed males
  CWDt.f <- tmp # for cwd killed females
  CWDt.m <- tmp # for cwd killed females 

  f.inf <- tmp
  m.inf <- tmp
  mm.cases <- tmp
  mf.cases <- tmp
  ff.cases <- tmp
  fm.cases <- tmp
  em.cases <- tmp
  ef.cases <- tmp

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
    St.f[, t] <- St.f[, t] - cases.f # removing newly infected females
    St.m[, t] <- St.m[, t] - cases.m # for males 
    It.f[, t, 1] <- It.f[, t, 1] + cases.f # adding newly infected females
    It.m[, t, 1] <- It.m[, t, 1] + cases.m # for males
    envcases.f <- St.f[, t] * env.foi.f # environmental cases of females
    envcases.m <- St.m[, t] * env.foi.m # of males
    St.f[, t] <- St.f[, t] - envcases.f # removing environmentally infected females
    St.m[, t] <- St.m[, t] - envcases.m # males
    It.f[, t, 1] <- It.f[, t, 1] + envcases.f # adding environmentally infected females
    It.m[, t, 1] <- It.m[, t, 1] + envcases.m # males
    
    f.inf[, t] <- cases.f + envcases.f # new cases for females
    m.inf[, t] <- cases.m + envcases.m # and males
    #which of these new infections, on average, came from which transmision type
    mm.cases[, t] <- St.m[, t] * (1 - exp(-(beta.mm * sum(It.m[, t, ])/Nall^theta))) 
    mf.cases[, t] <- St.f[, t] * (1 - exp(-(beta.mf * sum(It.m[, t, ])/Nall^theta)))
    ff.cases[, t] <- St.f[, t] * (1 - exp(-(beta.ff * sum(It.f[, t, ])/Nall^theta)))
    fm.cases[, t] <- St.m[, t] * (1 - exp(-(beta.fm * sum(It.f[, t, ])/Nall^theta)))
    em.cases[, t] <- envcases.m
    ef.cases[, t] <- envcases.f
  }
  counts <- list(St.f = St.f, St.m = St.m, I1t.f = It.f[, , 1], 
                 I1t.m = It.m[, , 1], I2t.f = It.f[, , 2], I2t.m = It.m[,, 2], 
                 I3t.f = It.f[, , 3], I3t.m = It.m[, , 3], I4t.f = It.f[,, 4], 
                 I4t.m = It.m[, , 4], I5t.f = It.f[, , 5], I5t.m = It.m[,, 5], 
                 I6t.f = It.f[, , 6], I6t.m = It.m[, , 6], I7t.f = It.f[,, 7], 
                 I7t.m = It.m[, , 7], I8t.f = It.f[, , 8], I8t.m = It.m[,, 8], 
                 I9t.f = It.f[, , 9], I9t.m = It.m[, , 9], I10t.f = It.f[,, 10], 
                 I10t.m = It.m[, , 10]) # tracking how many individuals are in susceptible and infectious classes
  
  deaths <- list(Ht.f = Ht.f, Ht.m = Ht.m, Dt.f = Dt.f, Dt.m = Dt.m, 
                 CWDt.f = CWDt.f, CWDt.m = CWDt.m) # tracking cause of death

  tracking.inf <- list(new.f.inf = f.inf, 
                       new.m.inf = m.inf, 
                       mm.cases = mm.cases,
                       mf.cases = mf.cases,
                       ff.cases = ff.cases,
                       fm.cases = fm.cases,
                       em.cases = em.cases,
                       ef.cases = ef.cases) # compiling

  counts.long <- reshape2::melt(counts) %>% 
    dplyr::rename(age = Var1, month = Var2, population = value, category = L1) %>% 
    dplyr::mutate(year = (month - 1)/12, sex = as.factor(str_sub(category, -1)), disease = "no") #compiling based on year, sex, and susceptible
  counts.long$disease[str_sub(counts.long$category, 1, 1) == "I"] <- "yes" #labeling infected population
  counts.long$disease <- as.factor(counts.long$disease) 
  deaths.long <- reshape2::melt(deaths) %>% 
    dplyr::rename(age = Var1, month = Var2, population = value, category = L1) %>% 
    dplyr::mutate(year = (month - 1)/12, sex = as.factor(str_sub(category, -1))) # again for death

  tracking.inf.long <- reshape2::melt(tracking.inf) %>% 
    dplyr::rename(age = Var1, month = Var2, population = value, category = L1) %>%
    dplyr::mutate(year = (month - 1)/12) #making into usable df
  
  N <- n0 # population
  
  ####################################
  ####################################
  ####################################
  # S <- rep(N/(n.age.cats*2), n.age.cats*2) #equal age distribution of disease free equilibrium
  S <- popbio::stable.stage(M)*N #stable age distribution of disease free equilibrium
  ####################################
  ####################################
  ####################################
  
  age <- rep(c(1:n.age.cats), 2) #age naming vector
  sex <- rep(c("f","m"), each = n.age.cats) # sex naming vector
  S. <- paste0("S", age, sex) # making a vector of variable names
  for (i in 1:length(age)){
    val <- S[i] #calling the value from stable DFE
    var <- S.[i] #calling the name from above
    assign(var, val) #creating a variable in the local environment
  }
  
  age <- rep(c(rep(c(1:n.age.cats), each = 10)),2) # for infectious classes
  sex <- rep(c("f","m"), each = n.age.cats*10) # for infectious classes
  cat <- rep(c(rep(c(1:10), n.age.cats)),2) # for infectious classes
  I <- paste0("I", age, sex, cat) # for infectious classes
  for (i in 1:length(age)){
    var <- I[i] # for infectious classes
    assign(var, 0) # zero because no infectious animals at dfe
  }
  
  A <- str2lang(paste(I[which(str_detect(I, "f") == T)], collapse = " + ")) # this creates a sum of all infectious female classes
  B <- str2lang(paste(I[which(str_detect(I, "m") == T)], collapse = " + ")) # this creates a sum of all infectious male classes
  
  F.mat <- rep(list(NA), length(age)) #creating a place to store a ton of expressions
  names(F.mat) <- paste0("F", age, sex, cat) #naming each list for my own ability to check
  
  for (i in 1:(2*10*n.age.cats)){ #for each infectious age class
    if (i == 1){
      F.mat[[i]] <- substitute(C * (1 - exp(-(((beta.ff*(A))+
                                                 (beta.mf*(B)))
                                              /(N^theta)))), # this equation comes from the model vignette
                                list(A = A, B = B, C = str2lang(S.[1]))) #if fawn
    }
    if (i == 11){
      F.mat[[i]] <- substitute(C * (1 - exp(-(((beta.ff*(A))+
                                                 (beta.mf*(B)))
                                              /(N^theta)))),
                               list(A = A, B = B, C = str2lang(S.[2]))) #if juvenile f
    }
    if (i %in% seq(21, 10*n.age.cats, by = 10)){
      for (j in 3:n.age.cats){
        F.mat[[i]] <- substitute(C * (1 - exp(-(((beta.ff*(A))+
                                                   (beta.mf*(B)))
                                                /(N^theta)))),
                                 list(A = A, B = B, C = str2lang(S.[j]))) # if adult female in 3:n.age.cat
      }
    }
    if (i == 10*n.age.cats+1){
      F.mat[[i]] <- substitute(C * (1 - exp(-(((beta.mm*(B))+
                                                 (beta.fm*(A)))
                                              /(N^theta)))),
                               list(A = A, B = B, C = str2lang(S.[1+n.age.cats]))) # for males
    }
    if (i == 10*n.age.cats+11){
      F.mat[[i]] <- substitute(C * (1 - exp(-(((beta.mm*(B))+
                                                 (beta.fm*(A)))
                                              /(N^theta)))),
                               list(A = A, B = B, C = str2lang(S.[2+n.age.cats]))) # for males
    }
    if (i %in% seq(10*n.age.cats + 21, 2*10*n.age.cats, by = 10)){
      for (j in (n.age.cats+3):(2*n.age.cats)) {
        F.mat[[i]] <- substitute(C * (1 - exp(-(((beta.mm*(B))+
                                                   (beta.fm*(A)))
                                                /(N^theta)))),
                                 list(A = A, B = B, C = str2lang(S.[j]))) # for males
      }
    }
    if (!i %in% seq(1, 2*10*n.age.cats, by = 10)) {F.mat[i] <- 0} # if it is not the first infectious class, there are no gains of infection, only transfers
  }
  
  Vm.mat <- rep(list(NA), length(age)) #another expression holder
  names(Vm.mat) <- paste0("Vm", age, sex, cat) #named
  
  hunt <- rep(0, length(I))
  hunt[1:10] <- "hunt.mort.fawn"
  hunt[11:20] <- "hunt.mort.juv.f"
  hunt[21:(n.age.cats*10)] <- "hunt.mort.ad.f"
  hunt[(n.age.cats*10+1):(n.age.cats*10+10)] <- "hunt.mort.fawn"
  hunt[(n.age.cats*10+11):(n.age.cats*10+20)] <- "hunt.mort.juv.m"
  hunt[(n.age.cats*10+21):(n.age.cats*10*2)] <- "hunt.mort.ad.m" #creating vector of appropriate hunting variables
  
  mort <- rep(0, length(I))
  mort[1:10] <- "fawn.sur"
  mort[11:20] <- "juv.sur"
  mort[21:(n.age.cats*10)] <- "ad.f.sur"
  mort[(n.age.cats*10+1):(n.age.cats*10+10)] <- "fawn.sur"
  mort[(n.age.cats*10+11):(n.age.cats*10+20)] <- "juv.sur"
  mort[(n.age.cats*10+21):(n.age.cats*10*2)] <- "ad.m.sur" #creating vector of appropriate natural mortality variables
  
  for (i in 1:(2*10*n.age.cats)){
    # this says that infectious age classes lose infection because of hunting, mortality, and transfers
    Vm.mat[[i]] <- substitute(A*((B/12)+(1-C) + p),
                              list(A = as.symbol(I[i]), 
                                   B = as.symbol(hunt[i]), 
                                   C = as.symbol(mort[i]))) 
  }
  
  Vp.mat <- rep(list(NA), length(age)) # another expression holder
  names(Vp.mat) <- paste0("Vp", age, sex, cat)  # names
  
  for (i in 1:(2*10*n.age.cats)){
    #there are no transfers into first disease class
    if (i %in% seq(1, (2*10*n.age.cats), by = 10)){
      Vp.mat[[i]] <- 0
    }
    #there are for disease classes 2:10
    if (!i %in% seq(1, (2*10*n.age.cats), by = 10)){
      Vp.mat[[i]] <- substitute(A*p, list(A = as.symbol(I[i-1]))) 
    }
  }
  
  V.mat <- rep(list(NA), length(age)) # another expression holder
  names(V.mat) <- paste0("V", age, sex, cat) # named
  
  for (i in 1:(2*10*n.age.cats)){
    V.mat[[i]] <- substitute(a - b, list(a = Vm.mat[[i]], b = Vp.mat[[i]])) # subtracting appropriate expression to create loss of infection
  }
  
  f.matrix <- matrix(0, nrow = (2*10*n.age.cats), ncol = (2*10*n.age.cats)) # a blank matrix to store evaluations of gains
  v.matrix <- matrix(0, nrow = (2*10*n.age.cats), ncol = (2*10*n.age.cats)) # a blank matrix to store evaluations of losses
  
  for (i in 1:(2*10*n.age.cats)){
    for (j in 1:(2*10*n.age.cats)){
      eq <- F.mat[[i]] #calling gain equations
      disease.class <- I[j] #calling infectious classes
      partial.d <- D(eq, disease.class) # creating a partial derivivative of the "q-th" eq with respect to the "r-th" infectious class
      f.matrix[i,j] <- eval(partial.d) # Evaluating the partial derivative at disease free equilibrium, specified above
    }
  }
  
  for (i in 1:(2*10*n.age.cats)){
    for (j in 1:(2*10*n.age.cats)){
      eq <- V.mat[[i]] #ditto but for losses
      disease.class <- I[j]
      partial.d <- D(eq, disease.class)
      v.matrix[i,j] <- eval(partial.d)
    }
  }
  R0 <- max(Re(eigen(f.matrix %*% solve(v.matrix))$values)) # this finds the largest eigen value of the matrix product of gains times inverse losses
  #NOTE: sometimes, this evaluates with complex numbers (e.g. +0i). Re() finds the real component of this complex number.
  
  N <- counts.long %>% 
    group_by(year) %>% 
    summarize(n = sum(population)) # total pop by time step
  I <- counts.long %>% 
    filter(str_detect(category, "I")) %>%
    group_by(year) %>% 
    summarize(n = sum(population)) # total infected by time step
  
  # prevalence approximation of R0 based on max prevalence
  s.inf <- min((N$n-I$n)/N$n)
  R0.prev.approx <- log(s.inf)/(s.inf-1)
  
  I1 <- rep(NA, (n.years*12)-1)
  I2 <- rep(NA, (n.years*12)-1)
  for (i in 2:(n.years*12)){
    I1[i] <- as.numeric(counts.long %>% 
      filter(year == unique(year)[i-1] &
               str_detect(category, "I")) %>% 
      summarize(n = sum(population))) # number of infected at kth time step
    I2[i] <- as.numeric(counts.long %>% 
      filter(year == unique(year)[i] &
               str_detect(category, "I")) %>% 
      summarize(n = sum(population))) # number of infected at lth time step
  }
  r <- max(log(na.omit(I2/I1))) #this is the max monthly exponential rate of growth
  ser <- seq(10, 23, length.out = 5) # these are a series of generation times that are reasonable given p
  R0.r.approx <- 1+(r*ser) #This comes from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1578275/
  
  R0.point <- data.frame(value = c(R0, R0.prev.approx),
                      name = c("Next Generation", "Prev. Approx."))
  
  R0.range <- data.frame(min = min(R0.r.approx), 
                         max = max(R0.r.approx),
                         name = c("r Approx."))
  
  output <- list(counts = counts.long, 
                 deaths = deaths.long, 
                 tracking.inf = tracking.inf.long, 
                 R0 = list(R0.point, R0.range)) #prints it out
}
