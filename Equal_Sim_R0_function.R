# Function to return Equal-R0 scenarios for betaff
equal_sims <- function (params) {
  for (v in 1:length(params)) assign(names(params)[v], params[[v]])
  
  #creating output df
  df <- matrix(NA, nrow = 6, ncol = 8)
  df[,1] <- c("Male.S", "Male.D", "Female.D", "DD-Equal", "DD-Male.D", "DD-Female.D")
  
  # this all comes from R0 in prior functions
  N <- n0 # population
  S <- rep(N/(n.age.cats*2), n.age.cats*2) #equal age distribution of disease free equilibrium 
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
  F. <- rep(F.mat, each = (2*10*n.age.cats))
  V. <- rep(V.mat, each = (2*10*n.age.cats))
  I. <- rep(I, (2*10*n.age.cats))
  partial.df <- rep(list(NA), ((2*10*n.age.cats)^2))
  partial.dv <- rep(list(NA), ((2*10*n.age.cats)^2))
  for (i in 1:((2*10*n.age.cats)^2)){
    eqf <- F.[[i]] #calling gain equations
    eqv <- V.[[i]]
    disease.class <- I.[i] #calling infectious classes
    partial.df[[i]] <- D(eqf, disease.class) # creating a partial derivivative of the "q-th" eq with respect to the "r-th" infectious class
    partial.dv[[i]] <- D(eqv, disease.class)
  }
  # So now that we have the derivatives we hold onto these - they wont change
  
  #based on these derivatives, we can decide what a basic population looks like with equal transmission
  beta.fm <- beta.ff * gamma.fm.e # contribution of fm
  beta.mm <- beta.ff * gamma.mm.e # contribution of mm
  beta.mf <- beta.ff * gamma.mf.e # contribution of mf
  fawn.sur <- fawn.an.sur^(1/12) # monthly survival
  juv.sur <- juv.an.sur^(1/12) # monthly survival
  ad.f.sur <- ad.an.f.sur^(1/12) # monthly survival
  ad.m.sur <- ad.an.m.sur^(1/12) # monthly survival
  theta <- 1 #frequency dependence
  
  # blank storage vectors
  eval.f <- rep(NA, ((2*10*n.age.cats)^2)) 
  eval.v <- rep(NA, ((2*10*n.age.cats)^2))
  
  # filling vectors
  for (i in 1:((2*10*n.age.cats)^2)){
    eval.f[i] <- eval(partial.df[[i]])
    eval.v[i] <- eval(partial.dv[[i]])
  }
  
  #creating respective matrices
  f.matrix <- matrix(eval.f, 
                     nrow = (2*10*n.age.cats), 
                     ncol = (2*10*n.age.cats),
                     byrow = T)
  v.matrix <- matrix(eval.v, 
                     nrow = (2*10*n.age.cats), 
                     ncol = (2*10*n.age.cats),
                     byrow = T)
  
  # R0 for equal transmission
  R0 <- max(Re(eigen(f.matrix %*% solve(v.matrix))$values))
  
  # filling in output matrix with proper info
  df[1,2] <- R0
  df[1,3] <- beta.ff
  df[1,4] <- gamma.fm.e
  df[1,5] <- gamma.mm.e
  df[1,6] <- gamma.mf.e
  df[1,7] <- theta
  df[1,8] <- env.foi
  
  # now we can find what values make sense for beta.ff with different transmission scenarios 
  
  #storage vectors
  R0. <- rep(NA, length.steps)
  bff <- rep(NA, length.steps)
  
  #settinging intial guess for bff
  bff[1] <- min.bff.md 
  
  #starting guessing, calculating, evaluating
  for (i in 1:length.steps) {
    beta.ff <- bff[i] #guess = beta
    
    # what each transmission type change to based on guess
    beta.fm <- beta.ff * gamma.fm.md 
    beta.mm <- beta.ff * gamma.mm.md
    beta.mf <- beta.ff * gamma.mf.md
    
    # evaluate former partial derivatives
    for (j in 1:((2*10*n.age.cats)^2)){
      eval.f[j] <- eval(partial.df[[j]])
      eval.v[j] <- eval(partial.dv[[j]])
    }
    
    # make respective matrices again
    f.matrix <- matrix(eval.f, 
                       nrow = (2*10*n.age.cats), 
                       ncol = (2*10*n.age.cats),
                       byrow = T)
    v.matrix <- matrix(eval.v, 
                       nrow = (2*10*n.age.cats), 
                       ncol = (2*10*n.age.cats),
                       byrow = T)
    #get the respective R0 per beta
    R0.[i] <- max(Re(eigen(f.matrix %*% solve(v.matrix))$values))
    
    #if the final guess is still too small, we know we need to change our inputs
    if(R0 - R0.[i] > 0 &
       i == length.steps) {
      print(expression("Male dominated step too small or initial guess too low"))
      break
    }
    
    #if the first guess is too big, we know we need to change our inputs
    if(R0 - R0.[i] < 0 &
       i == 1) {
      print(expression("Male dominated step too large or initial guess too high"))
      break
    }
    
    #if a guess is greater than R0, but not because of input parameters, 
    # we know that the former R0 or this one were pretty close
    if(R0 - R0.[i] < 0) {
      R0.m <- (R0.[i-1] + R0.[i])/2
      bff.m <- (bff[i-1] + bff[i])/2
      break
    }
    
    # if the guess is still small but not the final guess, keep trudging on
    else bff[i+1] <- bff[i] + step.size
  }
  
  #storing the average R0.m and average bff.m as well as the multipliers associated
  df[2,2] <- R0.m
  df[2,3] <- bff.m
  df[2,4] <- gamma.fm.md
  df[2,5] <- gamma.mm.md
  df[2,6] <- gamma.mf.md
  df[2,7] <- theta
  df[2,8] <- env.foi
  
  # doing this all over again but for a female dominated scenario
  R0. <- rep(NA, length.steps)
  bff <- rep(NA, length.steps)
  bff[1] <- min.bff.fd
  for (i in 1:length.steps) {
    beta.ff <- bff[i]
    beta.fm <- beta.ff * gamma.fm.fd
    beta.mm <- beta.ff * gamma.mm.fd
    beta.mf <- beta.ff * gamma.mf.fd
    for (j in 1:((2*10*n.age.cats)^2)){
      eval.f[j] <- eval(partial.df[[j]])
      eval.v[j] <- eval(partial.dv[[j]])
    }
    f.matrix <- matrix(eval.f, 
                       nrow = (2*10*n.age.cats), 
                       ncol = (2*10*n.age.cats),
                       byrow = T)
    v.matrix <- matrix(eval.v, 
                       nrow = (2*10*n.age.cats), 
                       ncol = (2*10*n.age.cats),
                       byrow = T)
    R0.[i] <- max(Re(eigen(f.matrix %*% solve(v.matrix))$values))
    if(R0 - R0.[i] > 0 &
       i == length.steps) {
      print(expression("Female dominated step too small or initial guess too low"))
      break
    }
    if(R0 - R0.[i] < 0 &
       i == 1) {
      print(expression("Female dominated step too large or initial guess too high"))
      break
    }
    if(R0 - R0.[i] < 0) {
      R0.f <- (R0.[i-1] + R0.[i])/2
      bff.f <- (bff[i-1] + bff[i])/2
      break
    }
    else bff[i+1] <- bff[i] + step.size
  }
  df[3,2] <- R0.f
  df[3,3] <- bff.f
  df[3,4] <- gamma.fm.fd
  df[3,5] <- gamma.mm.fd
  df[3,6] <- gamma.mf.fd
  df[3,7] <- theta
  df[3,8] <- env.foi
  
  # R0. <- rep(NA, length.steps)
  # bff <- rep(NA, length.steps)
  # bff[1] <- beta.ff.dens
  # theta <- 0
  # for (i in 1:length.steps) {
  #   beta.ff <- bff[i]
  #   beta.fm <- beta.ff * gamma.fm.e
  #   beta.mm <- beta.ff * gamma.mm.e
  #   beta.mf <- beta.ff * gamma.mf.e
  #   for (j in 1:((2*10*n.age.cats)^2)){
  #     eval.f[j] <- eval(partial.df[[j]])
  #     eval.v[j] <- eval(partial.dv[[j]])
  #   }
  #   f.matrix <- matrix(eval.f, 
  #                      nrow = (2*10*n.age.cats), 
  #                      ncol = (2*10*n.age.cats),
  #                      byrow = T)
  #   v.matrix <- matrix(eval.v, 
  #                      nrow = (2*10*n.age.cats), 
  #                      ncol = (2*10*n.age.cats),
  #                      byrow = T)
  #   R0.[i] <- max(Re(eigen(f.matrix %*% solve(v.matrix))$values))
  #   if(R0 - R0.[i] > 0 &
  #      i == length.steps) {
  #     print(expression("Density step too small or initial guess too low"))
  #     break
  #   }
  #   if(R0 - R0.[i] < 0 &
  #      i == 1) {
  #     print(expression("Density step too large or initial guess too high"))
  #     break
  #   }
  #   if(R0 - R0.[i] < 0) {
  #     R0.f <- (R0.[i-1] + R0.[i])/2
  #     bff.f <- (bff[i-1] + bff[i])/2
  #     break
  #   }
  #   else bff[i+1] <- bff[i] + step.size.dens
  # }
  # df[4,2] <- R0.f
  # df[4,3] <- bff.f
  # df[4,4] <- gamma.fm.e
  # df[4,5] <- gamma.mm.e
  # df[4,6] <- gamma.mf.e
  # df[4,7] <- theta
  # df[4,8] <- env.foi
  # 
  # R0. <- rep(NA, length.steps)
  # bff <- rep(NA, length.steps)
  # bff[1] <- min.bff.md.dens
  # for (i in 1:length.steps) {
  #   beta.ff <- bff[i]
  #   beta.fm <- beta.ff * gamma.fm.md 
  #   beta.mm <- beta.ff * gamma.mm.md
  #   beta.mf <- beta.ff * gamma.mf.md
  #   for (j in 1:((2*10*n.age.cats)^2)){
  #     eval.f[j] <- eval(partial.df[[j]])
  #     eval.v[j] <- eval(partial.dv[[j]])
  #   }
  #   f.matrix <- matrix(eval.f, 
  #                      nrow = (2*10*n.age.cats), 
  #                      ncol = (2*10*n.age.cats),
  #                      byrow = T)
  #   v.matrix <- matrix(eval.v, 
  #                      nrow = (2*10*n.age.cats), 
  #                      ncol = (2*10*n.age.cats),
  #                      byrow = T)
  #   R0.[i] <- max(Re(eigen(f.matrix %*% solve(v.matrix))$values))
  #   if(R0 - R0.[i] > 0 &
  #      i == length.steps) {
  #     print(expression("Male dominated (density) step too small or initial guess too low"))
  #     break
  #   }
  #   if(R0 - R0.[i] < 0 &
  #      i == 1) {
  #     print(expression("Male dominated (density) step too large or initial guess too high"))
  #     break
  #   }
  #   if(R0 - R0.[i] < 0) {
  #     R0.m <- (R0.[i-1] + R0.[i])/2
  #     bff.m <- (bff[i-1] + bff[i])/2
  #     break
  #   }
  #   else bff[i+1] <- bff[i] + step.size.dens
  # }
  # df[5,2] <- R0.m
  # df[5,3] <- bff.m
  # df[5,4] <- gamma.fm.md
  # df[5,5] <- gamma.mm.md
  # df[5,6] <- gamma.mf.md
  # df[5,7] <- theta
  # df[5,8] <- env.foi
  # 
  # # doing this all over again but for a female dominated scenario
  # R0. <- rep(NA, length.steps)
  # bff <- rep(NA, length.steps)
  # bff[1] <- min.bff.fd.dens
  # for (i in 1:length.steps) {
  #   beta.ff <- bff[i]
  #   beta.fm <- beta.ff * gamma.fm.fd
  #   beta.mm <- beta.ff * gamma.mm.fd
  #   beta.mf <- beta.ff * gamma.mf.fd
  #   for (j in 1:((2*10*n.age.cats)^2)){
  #     eval.f[j] <- eval(partial.df[[j]])
  #     eval.v[j] <- eval(partial.dv[[j]])
  #   }
  #   f.matrix <- matrix(eval.f, 
  #                      nrow = (2*10*n.age.cats), 
  #                      ncol = (2*10*n.age.cats),
  #                      byrow = T)
  #   v.matrix <- matrix(eval.v, 
  #                      nrow = (2*10*n.age.cats), 
  #                      ncol = (2*10*n.age.cats),
  #                      byrow = T)
  #   R0.[i] <- max(Re(eigen(f.matrix %*% solve(v.matrix))$values))
  #   if(R0 - R0.[i] > 0 &
  #      i == length.steps) {
  #     print(expression("Female dominated (density) step too small or initial guess too low"))
  #     break
  #   }
  #   if(R0 - R0.[i] < 0 &
  #      i == 1) {
  #     print(expression("Female dominated (density) step too large or initial guess too high"))
  #     break
  #   }
  #   if(R0 - R0.[i] < 0) {
  #     R0.f <- (R0.[i-1] + R0.[i])/2
  #     bff.f <- (bff[i-1] + bff[i])/2
  #     break
  #   }
  #   else bff[i+1] <- bff[i] + step.size.dens
  # }
  # df[6,2] <- R0.f
  # df[6,3] <- bff.f
  # df[6,4] <- gamma.fm.fd
  # df[6,5] <- gamma.mm.fd
  # df[6,6] <- gamma.mf.fd
  # df[6,7] <- theta
  # df[6,8] <- env.foi
  
  # columns named for future use
  df <- data.frame(df)
  colnames(df) <- c("Scenario", "R0", "B.ff", "G.fm", "G.mm", "G.mf", "theta", "env.foi")
  return(df)
}

params <- list(beta.ff = 0.045, # equal transmission beta one is willing to consider
               min.bff.md = 0.02961, # minimum guess for beta.ff under male dominated scenario
               min.bff.fd = 0.07572, # minimum guess for beta.ff under female dominated scenario
               # beta.ff.dens = 0.045/2000, # equal transmission beta one is willing to consider
               # min.bff.md.dens = 0.029615/2000, # minimum guess for beta.ff under male dominated scenario
               # min.bff.fd.dens = 0.075725/2000, # minimum guess for beta.ff under female dominated scenario
               #best to start with bigger steps and shorter lengths
               step.size = 0.000001, # how guesses increase (guess[prior]+step=guess[now])
               step.size.dens = 0.000000001,
               length.steps = 1000, # how far out to go with steps
               gamma.mm.e = 2, # how will gamma.mm be changed during equal scenario 
               gamma.mf.e = 1, # how will gamma.mf be changed during equal scenario
               gamma.fm.e = 2, # how will gamma.fm be changed during equal scenario
               gamma.mm.md = 4, # how will gamma.mm be changed during male dominated scenario 
               gamma.mf.md = 2, # how will gamma.mf be changed during male dominated scenario
               gamma.fm.md = 1, # how will gamma.fm be changed during male dominated scenario
               gamma.mm.fd = .5, # how will gamma.mm be changed during female dominated scenario
               gamma.mf.fd = .5, # how will gamma.mf be changed during female dominated scenario
               gamma.fm.fd = 2, # how will gamma.fm be changed during female dominated scenario
               # basic values underlying population and hunting
               fawn.an.sur = 0.354, 
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
               env.foi = 0.0005,
               n0 = 2000,
               n.years = 1)
df <- equal_sims(params)
options(digits=10) 
write.csv(df, "R0.Scenarios.10.5.csv")
# View(read_csv("R0.Scenarios.extreme.6.8.csv"))
