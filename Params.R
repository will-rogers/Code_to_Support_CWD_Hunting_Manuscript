# Default Parameters for Elk
############################



# from Reithel et al. 2007 and Lubow et al. 2010
params <- list(fawn.an.sur = 0.354, 
               juv.an.sur = 0.883, 
               ad.an.f.sur = 0.972, 
               ad.an.m.sur = 0.968, 
               fawn.repro = 0, 
               juv.repro = 0.198, 
               ad.repro = .928, 
               hunt.mort.fawn = 0, 
               hunt.mort.juv.f = 0, 
               hunt.mort.juv.m = 0,
               hunt.mort.ad.f = 0, 
               hunt.mort.ad.m = 0, 
               ini.fawn.prev = 0.02,
               ini.juv.prev = 0.03, 
               ini.ad.f.prev = 0.04,  
               ini.ad.m.prev = 0.04,
               n.age.cats = 12,  
               p = 0.43, 
               env.foi = 0,  
               beta.ff = 0.07, 
               gamma.mm = 1, 
               gamma.mf = 1, 
               gamma.fm = 1,
               theta = 1, 
               n0 = 2000, 
               n.years = 25, 
               rel.risk = 1.0)

cwd_det_model
for (v in 1:length(params)) assign(names(params)[v], params[[v]])
months <- seq(1, n.years * 12)
hunt.mo <- rep(0, n.years * 12)
hunt.mo[months%%12 == 7] <- 1
fawn.sur <- fawn.an.sur^(1/12)
juv.sur <- juv.an.sur^(1/12)
ad.f.sur <- ad.an.f.sur^(1/12)
ad.m.sur <- ad.an.m.sur^(1/12)
ini.f.prev <- c(ini.fawn.prev, ini.juv.prev, rep(ini.ad.f.prev, 
                                                 (n.age.cats - 2)))
ini.m.prev <- c(ini.fawn.prev, ini.juv.prev, rep(ini.ad.m.prev, 
                                                 (n.age.cats - 2)))
Sur.f <- c(fawn.sur, juv.sur, rep(ad.f.sur, n.age.cats - 
                                    2))
Sur.m <- c(fawn.sur, juv.sur, rep(ad.m.sur, n.age.cats - 
                                    2))
M <- matrix(rep(0, n.age.cats * 2 * n.age.cats * 2), nrow = n.age.cats * 
              2)
M[row(M) == (col(M) + 1)] <- c(juv.an.sur * (1 - hunt.mort.juv.f), 
                               rep(ad.an.f.sur * (1 - hunt.mort.ad.f), n.age.cats - 
                                     2), 0, c(juv.an.sur * (1 - hunt.mort.juv.m), rep(ad.an.m.sur * 
                                                                                        (1 - hunt.mort.ad.m), n.age.cats - 2)))
M[n.age.cats, n.age.cats] <- ad.an.f.sur * (1 - hunt.mort.ad.f)
M[n.age.cats * 2, n.age.cats * 2] <- ad.an.m.sur * (1 - hunt.mort.ad.m)
M[1, 1:n.age.cats] <- c(0, juv.repro, rep(ad.repro, n.age.cats - 
                                            2)) * 0.5 * fawn.an.sur * (1 - hunt.mort.fawn)
M[n.age.cats + 1, 1:n.age.cats] <- M[1, 1:n.age.cats]
tmp <- matrix(0, nrow = n.age.cats, ncol = n.years * 12)
St.f <- tmp
St.m <- tmp
It.m <- array(rep(tmp), dim = c(n.age.cats, n.years * 12, 
                                10))
It.f <- array(rep(tmp), dim = c(n.age.cats, n.years * 12, 
                                10))
Ht.f <- tmp
Ht.m <- tmp
Dt.f <- tmp
Dt.m <- tmp
CWDt.f <- tmp
CWDt.m <- tmp
St.f[, 1] <- popbio::stable.stage(M)[1:n.age.cats] * n0 * 
  (1 - ini.f.prev)
St.m[, 1] <- popbio::stable.stage(M)[(n.age.cats + 1):(n.age.cats * 
                                                         2)] * n0 * (1 - ini.m.prev)
It.m[, 1, 1:10] <- popbio::stable.stage(M)[1:n.age.cats] * 
  n0/10 * ini.m.prev
It.f[, 1, 1:10] <- popbio::stable.stage(M)[(n.age.cats + 
                                              1):(n.age.cats * 2)] * n0/10 * ini.f.prev
for (t in 2:(n.years * 12)) {
  if (t%%12 == 2) {
    St.f[2:(n.age.cats - 1), t] <- St.f[1:(n.age.cats - 
                                             2), t - 1]
    St.f[n.age.cats, t] <- St.f[n.age.cats, t - 1] + 
      St.f[(n.age.cats - 1), t - 1]
    St.m[2:(n.age.cats - 1), t] <- St.m[1:(n.age.cats - 
                                             2), t - 1]
    St.m[n.age.cats, t] <- St.m[n.age.cats, t - 1] + 
      St.m[(n.age.cats - 1), t - 1]
    It.f[2:(n.age.cats - 1), t, ] <- It.f[1:(n.age.cats - 
                                               2), t - 1, ]
    It.f[n.age.cats, t, ] <- It.f[n.age.cats, t - 1, 
                                  ] + It.f[(n.age.cats - 1), t - 1, ]
    It.m[2:(n.age.cats - 1), t, ] <- It.m[1:(n.age.cats - 
                                               2), t - 1, ]
    It.m[n.age.cats, t, ] <- It.m[n.age.cats, t - 1, 
                                  ] + It.m[(n.age.cats - 1), t - 1, ]
    I_juv <- sum(It.f[2, t - 1, ])
    I_adults <- sum(It.f[3:n.age.cats, t - 1, ])
    St.f[1, t] <- ((St.f[2, t - 1] + I_juv) * juv.repro + 
                     (sum(St.f[3:n.age.cats, t - 1]) + I_adults) * 
                     ad.repro) * 0.5
    St.m[1, t] <- ((St.f[2, t - 1] + I_juv) * juv.repro + 
                     (sum(St.f[3:n.age.cats, t - 1]) + I_adults) * 
                     ad.repro) * 0.5
  }
  if (t%%12 != 2) {
    St.f[, t] <- St.f[, t - 1]
    St.m[, t] <- St.m[, t - 1]
    It.f[, t, ] <- It.f[, t - 1, ]
    It.m[, t, ] <- It.m[, t - 1, ]
  }
  St.f[, t] <- St.f[, t] * Sur.f
  St.m[, t] <- St.m[, t] * Sur.m
  It.f[, t, ] <- It.f[, t, ] * Sur.f
  It.m[, t, ] <- It.m[, t, ] * Sur.m
  Dt.f[, t] <- (St.f[, t] + rowSums(It.f[, t, ])) * (1 - 
                                                       Sur.f)
  Dt.m[, t] <- (St.m[, t] + rowSums(It.m[, t, ])) * (1 - 
                                                       Sur.m)
  if (hunt.mo[t] == 1) {
    Iall.f <- rowSums(It.f[, t, ])
    Iall.m <- rowSums(It.m[, t, ])
    Nt.f <- St.f[, t] + Iall.f
    Nt.m <- St.m[, t] + Iall.m
    hunted.f <- Nt.f * c(hunt.mort.fawn, hunt.mort.juv.f, 
                         rep(hunt.mort.ad.f, n.age.cats - 2))
    hunted.m <- Nt.m * c(hunt.mort.fawn, hunt.mort.juv.m, 
                         rep(hunt.mort.ad.m, n.age.cats - 2))
    Ht.f[, t] <- hunted.f
    Ht.m[, t] <- hunted.m
    hunted.i.f <- (rel.risk * Iall.f * hunted.f)/(St.f[, 
                                                       t] + rel.risk * Iall.f)
    hunted.i.m <- (rel.risk * Iall.m * hunted.m)/(St.m[, 
                                                       t] + rel.risk * Iall.m)
    hunted.i.f[which(is.na(hunted.i.f))] <- 0
    hunted.i.m[which(is.na(hunted.i.m))] <- 0
    St.f[, t] <- St.f[, t] - (hunted.f - hunted.i.f)
    St.m[, t] <- St.m[, t] - (hunted.m - hunted.i.m)
    It.f[, t, ] <- It.f[, t, ] * (1 - hunted.i.f/Iall.f)
    It.m[, t, ] <- It.m[, t, ] * (1 - hunted.i.m/Iall.m)
  }
  f.move <- It.f[, t, ] * p
  m.move <- It.m[, t, ] * p
  It.f[, t, 1] <- It.f[, t, 1] - f.move[, 1]
  It.f[, t, 2:10] <- It.f[, t, 2:10] - f.move[, 2:10] + 
    f.move[, 1:9]
  It.m[, t, 1] <- It.m[, t, 1] - m.move[, 1]
  It.m[, t, 2:10] <- It.m[, t, 2:10] - m.move[, 2:10] + 
    m.move[, 1:9]
  CWDt.f[, t] <- f.move[, 10]
  CWDt.m[, t] <- m.move[, 10]
  Iall <- sum(It.f[, t, ] + It.m[, t, ])
  Nall <- sum(St.f[, t] + St.m[, t]) + Iall
  cases.f <- St.f[, t] * (1 - exp(-(beta.f * (Iall/Nall^theta))))
  cases.m <- St.m[, t] * (1 - exp(-(beta.m * (Iall/Nall^theta))))
  St.f[, t] <- St.f[, t] - cases.f
  St.m[, t] <- St.m[, t] - cases.m
  It.f[, t, 1] <- It.f[, t, 1] + cases.f
  It.m[, t, 1] <- It.m[, t, 1] + cases.m
  envcases.f <- St.f[, t] * env.foi
  envcases.m <- St.m[, t] * env.foi
  St.f[, t] <- St.f[, t] - envcases.f
  St.m[, t] <- St.m[, t] - envcases.m
  It.f[, t, 1] <- It.f[, t, 1] + envcases.f
  It.m[, t, 1] <- It.m[, t, 1] + envcases.m
}
counts <- list(St.f = St.f, St.m = St.m, I1t.f = It.f[, , 
                                                      1], I1t.m = It.m[, , 1], I2t.f = It.f[, , 2], I2t.m = It.m[, 
                                                                                                                 , 2], I3t.f = It.f[, , 3], I3t.m = It.m[, , 3], I4t.f = It.f[, 
                                                                                                                                                                              , 4], I4t.m = It.m[, , 4], I5t.f = It.f[, , 5], I5t.m = It.m[, 
                                                                                                                                                                                                                                           , 5], I6t.f = It.f[, , 6], I6t.m = It.m[, , 6], I7t.f = It.f[, 
                                                                                                                                                                                                                                                                                                        , 7], I7t.m = It.m[, , 7], I8t.f = It.f[, , 8], I8t.m = It.m[, 
                                                                                                                                                                                                                                                                                                                                                                     , 8], I9t.f = It.f[, , 9], I9t.m = It.m[, , 9], I10t.f = It.f[, 
                                                                                                                                                                                                                                                                                                                                                                                                                                   , 10], I10t.m = It.m[, , 10])
deaths <- list(Ht.f = Ht.f, Ht.m = Ht.m, Dt.f = Dt.f, Dt.m = Dt.m, 
               CWDt.f = CWDt.f, CWDt.m = CWDt.m)
counts.long <- melt(counts) %>% rename(age = Var1, month = Var2, 
                                       population = value, category = L1) %>% mutate(year = (month - 
                                                                                               1)/12, sex = as.factor(str_sub(category, -1)), disease = "no")
counts.long$disease[str_sub(counts.long$category, 1, 1) == 
                      "I"] <- "yes"
counts.long$disease <- as.factor(counts.long$disease)
deaths.long <- melt(deaths) %>% rename(age = Var1, month = Var2, 
                                       population = value, category = L1) %>% mutate(year = (month - 
                                                                                               1)/12, sex = as.factor(str_sub(category, -1)))
output <- list(counts = counts.long, deaths = deaths.long)