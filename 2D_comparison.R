require(readr)
require(CWDsims)
require(tidyverse)
require(ggpubr)

#calling the R0 validated scenarios
scenarios <- read.csv("R0.Scenarios.6.29.csv") 

#creating a dataframe with k*2 reps per scenario
k <- 50
df <- data.frame(beta.ff = rep(scenarios$B.ff, each = k*2)) 
df$gamma.mm <- rep(scenarios$G.mm, each = k*2)
df$gamma.mf <- rep(scenarios$G.mf, each = k*2)
df$gamma.fm <- rep(scenarios$G.fm, each = k*2)
df$theta <- rep(scenarios$theta, each = k*2)
df$env.foi <- rep(rep(0, 6), each = k*2)

#spreading hunting levels for four categories between the three scenarios
df$hunt.f <- rep(c(seq(0, .49, length.out= k), rep(0.05, k)), 6) 
df$hunt.m <- rep(c(rep(0.05, k), seq(0, .49, length.out= k)), 6)
# df$hunt.j.m <- rep(c(rep(.2, k*2), seq(0, .49, length.out= k),rep(.2, k)), 3)
# df$hunt.j.f <- rep(c(rep(.1, k*3), seq(0, .49, length.out= k)), 3)

#for sim - items to hold intermediates and output
df$fin.pop4.75   <- rep(NA, length(df$hunt.f))
df$fin.prev4.75   <- rep(NA, length(df$hunt.f))
df$fin.pop9.75   <- rep(NA, length(df$hunt.f))
df$fin.prev9.75  <- rep(NA, length(df$hunt.f))
df$fin.pop14.75  <- rep(NA, length(df$hunt.f))
df$fin.prev14.75 <- rep(NA, length(df$hunt.f))
df$fin.pop24.75  <- rep(NA, length(df$hunt.f))
df$fin.prev24.75 <- rep(NA, length(df$hunt.f))
df$R0 <- rep(NA, length(df$hunt.f))

# running through all scenarios
for (i in 1:(nrow(df))){
  out <- CWDsims::cwd_det_model_wiw(params = list(
    fawn.an.sur = 0.354, 
    juv.an.sur = 0.883, 
    ad.an.f.sur = 0.972, 
    ad.an.m.sur = 0.968, 
    fawn.repro = 0, 
    juv.repro = 0.198, 
    ad.repro = .928, 
    hunt.mort.fawn = 0.01,
    ini.fawn.prev = 0.02,
    ini.juv.prev = 0.03, 
    ini.ad.f.prev = 0.04,  
    ini.ad.m.prev = 0.04,
    n.age.cats = 12,
    p = 0.43,
    env.foi = 0, 
    theta = df$theta[i], 
    n0 = 2000, 
    n.years = 25, 
    rel.risk = 1.0, 
    hunt.mort.juv.f = 0.05,
    hunt.mort.juv.m = 0.05,
    hunt.mort.ad.f = df$hunt.f[i],
    hunt.mort.ad.m = df$hunt.m[i],
    beta.ff = df$beta.ff[i],
    gamma.mm = df$gamma.mm[i],
    gamma.mf = df$gamma.mf[i],
    gamma.fm = df$gamma.fm[i]))
  
  pop.sum <- out$counts %>%
    filter(month %% 12 == 10) %>%
    group_by(year) %>%
    summarize(n = sum(population)) #Total pop at time
  
  prev.sum   <- out$counts %>%
    filter(month %% 12 == 10) %>%
    group_by(year, disease) %>%
    summarize(n = sum(population)) %>%
    pivot_wider(names_from = disease, values_from = n) %>%
    mutate(prev = yes/(no + yes)) #Prev at time

  df$fin.pop4.75[i]    <- pop.sum$n[pop.sum$year == 4.75]
  df$fin.pop9.75[i]    <- pop.sum$n[pop.sum$year == 9.75]
  df$fin.pop14.75[i]   <- pop.sum$n[pop.sum$year == 14.75]
  df$fin.pop24.75[i]   <- pop.sum$n[pop.sum$year == 24.75]
  df$fin.prev4.75[i]   <- prev.sum$prev[prev.sum$year == 4.75]
  df$fin.prev9.75[i]   <- prev.sum$prev[prev.sum$year == 9.75]
  df$fin.prev14.75[i]  <- prev.sum$prev[prev.sum$year == 14.75]
  df$fin.prev24.75[i]  <- prev.sum$prev[prev.sum$year == 24.75]
}

#naming scenarios
df$scenario <- rep(rep(c("equal", "male.dominated", "female.dominated"), each = k*2), 2)

#naming the changing hunting parameter
df$rate.var <- rep(rep(c("af", "am"), each = k), 6)
df$rate.val <- NA
df$rate.val <- ifelse(df$rate.var == "af", df$hunt.f, df$rate.val)
df$rate.val <- ifelse(df$rate.var == "am", df$hunt.m, df$rate.val)
# df$rate.val <- ifelse(df$rate.var == "jm", df$hunt.j.m, df$rate.val)
# df$rate.val <- ifelse(df$rate.var == "jf", df$hunt.j.f, df$rate.val)

df$trans <- ifelse(df$theta == 1, "Freq. Dep.", "Dens. Dep.")

write.csv(df, "2D_comparison_0.5.9.2.csv")
df <- read.csv("2D_comparison_0.5.9.2.csv")

df$scenario <- factor(df$scenario, levels = c("equal", "male.dominated", "female.dominated"), ordered = T)
scenario.labs <- c("Equal", "Male-Dominated", "Female-Dominated")
names(scenario.labs) <- levels(df$scenario)

prev <- df %>% 
  filter(trans == "Freq. Dep.") %>% 
  ggplot(aes(rate.val, fin.prev14.75, color = rate.var)) + 
  geom_line(size = 1) +
  facet_grid(.~scenario,
  labeller = labeller(scenario = scenario.labs)) +
  scale_color_discrete(name = "Harvested Sex",
                       labels = c("Adult Female",
                                  "Adult Male")) +
  labs(x = "Independent Harvest Rate",
       y = "Prevalence at 15 Years") +
  theme_classic() + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) 

prev

pop <- df %>% 
  filter(trans == "Freq. Dep.") %>%  
  ggplot(aes(rate.val, fin.pop14.75, color = rate.var)) + 
  geom_line(size = 1) +
  facet_grid(.~scenario,
             labeller = labeller(scenario = scenario.labs)) +
  scale_color_discrete(name = "Harvested Sex",
                       labels = c("Adult Female",
                                  "Adult Male")) +
  labs(x = "Independent Harvest Rate",
       y = "Population at 15 Years") +
  theme_classic() + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) 

pop

ggarrange(prev, pop, ncol=1, common.legend = T, legend = "right", labels = c("A", "B"))
ggsave("Plots/2D_comparisons_pop.png", dpi=300, width = 6, height = 6)

hunts <- read.csv("2D_comparison_0.5.6.24.csv")

n0 <- 2000
n.age.cats <- 12
fawn.an.sur = 0.354
juv.an.sur = 0.883
ad.an.f.sur = 0.972 
ad.an.m.sur = 0.968
fawn.sur <- fawn.an.sur^(1/12)
juv.sur <- juv.an.sur^(1/12) 
ad.f.sur <- ad.an.f.sur^(1/12) 
ad.m.sur <- ad.an.m.sur^(1/12)
hunt.mort.fawn <- 0.01
p <- 0.43

#creating output df
df <- matrix(NA, nrow = 8, ncol = 7)
df[,1] <- c("Equal", "Male.D", "Female.D", "Envi", "DD-Equal", "DD-Male.D", "DD-Female.D", "DD-Envi")

# this all comes from R0 in prior functions
N <- n0 # population
S <- rep(N/(n.age.cats*2), n.age.cats*2) #equal age distribution of disease free equilibrium 
age <- rep(c(1:n.age.cats), 2) #age naming vector
sex <- rep(c("f","m"), each = n.age.cats) # sex naming vector
S. <- paste0("S", age, sex) # making a vector of variable names
for (i in 1:100){
  val <- S[i] #calling the value from stable DFE
  var <- S.[i] #calling the name from above
  assign(var, val) #creating a variable in the local environment
}
age <- rep(c(rep(c(1:n.age.cats), each = 10)),2) # for infectious classes
sex <- rep(c("f","m"), each = n.age.cats*10) # for infectious classes
cat <- rep(c(rep(c(1:10), n.age.cats)),2) # for infectious classes
I <- paste0("I", age, sex, cat) # for infectious classes
for (i in 1:100){
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
                                            /(N^theta))) + (env.foi.f*((A+B)/n0))), # this equation comes from the model vignette
                             list(A = A, B = B, C = str2lang(S.[1]))) #if fawn
  }
  if (i == 11){
    F.mat[[i]] <- substitute(C * (1 - exp(-(((beta.ff*(A))+
                                               (beta.mf*(B)))
                                            /(N^theta))) + (env.foi.f*((A+B)/n0))),
                             list(A = A, B = B, C = str2lang(S.[2]))) #if juvenile f
  }
  if (i %in% seq(21, 10*n.age.cats, by = 10)){
    for (j in 3:n.age.cats){
      F.mat[[i]] <- substitute(C * (1 - exp(-(((beta.ff*(A))+
                                                 (beta.mf*(B)))
                                              /(N^theta))) + (env.foi.f*((A+B)/n0))),
                               list(A = A, B = B, C = str2lang(S.[j]))) # if adult female in 3:n.age.cat
    }
  }
  if (i == 10*n.age.cats+1){
    F.mat[[i]] <- substitute(C * (1 - exp(-(((beta.mm*(B))+
                                               (beta.fm*(A)))
                                            /(N^theta))) + (env.foi.m*((A+B)/n0))),
                             list(A = A, B = B, C = str2lang(S.[1+n.age.cats]))) # for males
  }
  if (i == 10*n.age.cats+11){
    F.mat[[i]] <- substitute(C * (1 - exp(-(((beta.mm*(B))+
                                               (beta.fm*(A)))
                                            /(N^theta))) + (env.foi.m*((A+B)/n0))),
                             list(A = A, B = B, C = str2lang(S.[2+n.age.cats]))) # for males
  }
  if (i %in% seq(10*n.age.cats + 21, 2*10*n.age.cats, by = 10)){
    for (j in (n.age.cats+3):(2*n.age.cats)) {
      F.mat[[i]] <- substitute(C * (1 - exp(-(((beta.mm*(B))+
                                                 (beta.fm*(A)))
                                              /(N^theta))) + (env.foi.m*((A+B)/n0))),
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

eval.f <- rep(NA, ((2*10*n.age.cats)^2)) 
eval.v <- rep(NA, ((2*10*n.age.cats)^2))

hunts$R0 <- NA
hunts$l <- NA

for (i in 1:nrow(hunts)){
  beta.ff <- hunts$beta.ff[i]
  gamma.mm <- hunts$gamma.mm[i]
  beta.mm <- beta.ff * gamma.mm
  gamma.mf <- hunts$gamma.mf[i]
  beta.mf <- beta.ff * gamma.mf
  gamma.fm <- hunts$gamma.fm[i]
  beta.fm <- beta.ff * gamma.fm
  hunt.mort.juv.f = 0.1
  hunt.mort.juv.m = 0.1
  hunt.mort.ad.f = hunts$hunt.f[i]
  hunt.mort.ad.m = hunts$hunt.m[i]
  theta <- hunts$theta[i]
  env.foi.f <- hunts$evi.f[i]
  env.foi.m <- hunts$evi.m[i]
  beta.mm <- beta.ff * gamma.mm
  beta.mf <- beta.ff * gamma.mf
  beta.fm <- beta.ff * gamma.fm
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
  hunts$l[i]<-popbio::lambda(M)
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
  hunts$R0[i] <- max(Re(eigen(f.matrix %*% solve(v.matrix))$values))
}

write.csv(hunts, "2D_comparison_R0.6.24.csv")
hunts <- read.csv("2D_comparison_R0.6.24.csv")


trans.labs <- c("Density Dependent", "Frequency Dependent")
names(trans.labs) <- levels(hunts$trans)

R0 <- hunts %>% 
  filter(trans == "Freq. Dep.") %>%  
  ggplot(aes(rate.val, R0, color = rate.var)) + 
  geom_line() +
  facet_grid(.~scenario,
             labeller = labeller(scenario = scenario.labs)) +
  scale_color_discrete(name = "Harvest Parameter",
                       labels = c("Adult Female",
                                  "Adult Male")) +
  labs(x = "Independent Harvest Rate",
       y = expression(R[0])) +
  theme_classic() + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
R0
ggsave("Plots/2D_comparisons_R0.tiff", dpi=300, width = 6, height = 3)

R0vl <- hunts %>% 
  ggplot(aes(rate.val, color = rate.var)) + 
  geom_line(aes(y = log(1-fin.prev9.75)/(-fin.prev9.75))) +
  geom_line(aes(y = l)) +
  geom_hline(yintercept = 1) +
  scale_y_continuous(expression(R[0]), sec.axis = sec_axis(~ ., name = expression(lambda)))+
  facet_grid(trans~scenario,
             labeller = labeller(scenario = scenario.labs)) +
  scale_color_discrete(name = "Hunting Parameter",
                       labels = c("Adult Female",
                                  "Adult Male",
                                  "Juvenile Female",
                                  "Juvenile Male")) +
  labs(x = "Hunting Rate") +
  theme_classic() + 
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
R0vl
ggsave("Plots/2D_comparisons_R0.tiff", dpi=300, width = 6, height = 3)