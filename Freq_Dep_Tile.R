require(devtools)
require(CWDsims)
require(tidyr)
require(dplyr)
require(ggplot2)
require(gridExtra)
require(grid)
require(gtable)
require(RColorBrewer)
require(readr)
require(reshape2)
require(stringr)

#R0 validated scenarios
mat <- read_csv("R0.Scenarios.9.26.csv")

#basic parameters
fawn.an.sur = 0.354
juv.an.sur = 0.883 
ad.an.f.sur = 0.972 
ad.an.m.sur = 0.968 
fawn.repro = 0
juv.repro = 0.198
ad.repro = .928
hunt.mort.fawn = 0.01
hunt.mort.juv.f = 0.05
hunt.mort.juv.m = 0.05
hunt.mort.ad.f = 0.05
hunt.mort.ad.m = 0.05
ini.fawn.prev = 0.02
ini.juv.prev = 0.03
ini.ad.f.prev = 0.04  
ini.ad.m.prev = 0.04
n.age.cats = 12  
p = 0.43 
env.foi = 0
theta = 1 
n0 = 2000 
n.years = 25
rel.risk = 1.0


#year time points to slice
cut2.75 <- 2.75
cut4.75 <- 4.75
cut9.75 <- 9.75
cut14.75 <- 14.75
cut19.75 <- 19.75
cut24.75 <- 24.75

#sequences of hunting combinations to consider
hunt.mort.ad.f    <- seq(0, .5, length.out = 99) #vector of female harvest to run through
hunt.mort.ad.m    <- seq(0, .5, length.out = 99) #vector of male harvest to run through

#blank storage vectors
hunt.f <- c()
hunt.m <- c()
beta.ff.1 <- c()
gamma.mm.1 <- c()
gamma.mf.1 <- c()
gamma.fm.1 <- c()
scenario <- c()

# for each scenario, what is the hunting combination
for (l in 1:3){
  for (j in 1:length(hunt.mort.ad.m)){
    for (i in 1:length(hunt.mort.ad.f)){
      hunt.f <- c(hunt.f, hunt.mort.ad.f[i])
      hunt.m <- c(hunt.m, hunt.mort.ad.m[j])
      beta.ff.1 <- c(beta.ff.1, mat$B.ff[l])
      gamma.mm.1 <- c(gamma.mm.1, mat$G.mm[l])
      gamma.mf.1 <- c(gamma.mf.1, mat$G.mf[l])
      gamma.fm.1 <- c(gamma.fm.1, mat$G.fm[l])
      scenario <- c(scenario, mat$Scenario[l])
    }
  }
}

#storing these in a dataframe
df.freq <- data.frame(hunt.f,hunt.m, beta.ff.1, gamma.mm.1, gamma.mf.1, gamma.fm.1, scenario)

#storage vectors in df
df.freq$fin.pop2.75   <- rep(NA, length(df.freq$hunt.f))
df.freq$fin.prev2.75  <- rep(NA, length(df.freq$hunt.f))
df.freq$fin.pop4.75   <- rep(NA, length(df.freq$hunt.f))
df.freq$fin.prev4.75  <- rep(NA, length(df.freq$hunt.f))
df.freq$fin.pop9.75   <- rep(NA, length(df.freq$hunt.f))
df.freq$fin.prev9.75  <- rep(NA, length(df.freq$hunt.f))
df.freq$fin.pop14.75  <- rep(NA, length(df.freq$hunt.f))
df.freq$fin.prev14.75 <- rep(NA, length(df.freq$hunt.f))
df.freq$fin.pop19.75  <- rep(NA, length(df.freq$hunt.f))
df.freq$fin.prev19.75 <- rep(NA, length(df.freq$hunt.f))
df.freq$fin.pop24.75  <- rep(NA, length(df.freq$hunt.f))
df.freq$fin.prev24.75 <- rep(NA, length(df.freq$hunt.f))

for (i in 1:nrow(df.freq)){
  out <- CWDsims::cwd_det_model_wiw(params = list(fawn.an.sur = fawn.an.sur,
                                         juv.an.sur = juv.an.sur,
                                         ad.an.f.sur = ad.an.f.sur,
                                         ad.an.m.sur = ad.an.m.sur,
                                         fawn.repro = fawn.repro,
                                         juv.repro = juv.repro,
                                         ad.repro = ad.repro,
                                         hunt.mort.fawn = hunt.mort.fawn,
                                         hunt.mort.juv.f = hunt.mort.juv.f,
                                         hunt.mort.juv.m = hunt.mort.juv.m,
                                         hunt.mort.ad.f = df.freq$hunt.f[i], #ith female hunt
                                         hunt.mort.ad.m = df.freq$hunt.m[i], #ith male hunt
                                         ini.fawn.prev = ini.fawn.prev,
                                         ini.juv.prev = ini.juv.prev,
                                         ini.ad.f.prev = ini.ad.f.prev,
                                         ini.ad.m.prev = ini.ad.m.prev,
                                         n.age.cats = n.age.cats,
                                         p = p,
                                         env.foi = 0,
                                         beta.ff = df.freq$beta.ff.1[i], # ith beta ff
                                         gamma.mm = df.freq$gamma.mm.1[i], #""
                                         gamma.mf = df.freq$gamma.mf.1[i], #""
                                         gamma.fm = df.freq$gamma.fm.1[i], #""
                                         theta = 1, #frequency dep
                                         n0 = n0,
                                         n.years = n.years,
                                         rel.risk = rel.risk))
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
  
  df.freq$fin.pop2.75[i]    <- pop.sum$n[pop.sum$year==cut2.75] # population total at 2.75 years
  df.freq$fin.pop4.75[i]    <- pop.sum$n[pop.sum$year==cut4.75] # at 4.75
  df.freq$fin.pop9.75[i]    <- pop.sum$n[pop.sum$year==cut9.75] # at 9.75
  df.freq$fin.pop14.75[i]   <- pop.sum$n[pop.sum$year==cut14.75] # at 14.75
  df.freq$fin.pop19.75[i]   <- pop.sum$n[pop.sum$year==cut19.75] # at 19.75
  df.freq$fin.pop24.75[i]   <- pop.sum$n[pop.sum$year==cut24.75] # at 24.75
  df.freq$fin.prev2.75[i]   <- prev.sum$prev[prev.sum$year==cut2.75] # same for prevalence
  df.freq$fin.prev4.75[i]   <- prev.sum$prev[prev.sum$year==cut4.75]
  df.freq$fin.prev9.75[i]   <- prev.sum$prev[prev.sum$year==cut9.75]
  df.freq$fin.prev14.75[i]  <- prev.sum$prev[prev.sum$year==cut14.75]
  df.freq$fin.prev19.75[i]  <- prev.sum$prev[prev.sum$year==cut19.75]
  df.freq$fin.prev24.75[i]  <- prev.sum$prev[prev.sum$year==cut24.75]
}

pop.cut <- 1000 #population critical value
prev.cut <- 0.15 #prevalence critical value

# prev <- c("fin.prev2.75","fin.prev4.75","fin.prev9.75",
#           "fin.prev14.75", "fin.prev19.75", "fin.prev24.75")
# pop <- c("fin.pop2.75","fin.pop4.75","fin.pop9.75",
#          "fin.pop14.75", "fin.pop19.75", "fin.pop24.75")
# 
# df.freq <- melt(df.freq, measure.vars = prev)
# df.freq <- melt(df.freq, measure.vars = pop)
# colnames(df.freq)[c(8,10)] <- c("prevalence.t", "population.t")
# colnames(df.freq)[c(9,11)] <- c("prevalence", "population")

df.freq <- df.freq %>%  
  mutate(class = case_when(prevalence <= prev.cut &
                             population >= pop.cut ~ "Criteria Met",
                           prevalence > prev.cut &
                             population >= pop.cut ~ "High Prevalence",
                           prevalence <= prev.cut &
                             population < pop.cut ~ "Low Population",
                           prevalence > prev.cut &
                             population < pop.cut ~ "Fails Critera",
                           is.na(prevalence) |
                             is.na(population) ~ "Invalid"))

df.freq$class <- factor(df.freq$class)
order <- c("Criteria Met", "Low Population", "High Prevalence", "Fails Critera", "Invalid")
df.freq$class  <- factor(df.freq$class, levels = order)

df.freq <- df.freq %>%  
  mutate(year = case_when(str_detect(prevalence.t, "v2.75") &
                            str_detect(population.t, "p2.75") ~ "2.75",
                          str_detect(prevalence.t, "v4.75") &
                            str_detect(population.t, "p4.75") ~ "Year 5",
                          str_detect(prevalence.t, "v9.75") &
                            str_detect(population.t, "p9.75") ~ "Year 10",
                          str_detect(prevalence.t, "v14.75") &
                            str_detect(population.t, "p14.75") ~ "Year 15",
                          str_detect(prevalence.t, "v19.75") &
                            str_detect(population.t, "p19.75") ~ "Year 20",
                          str_detect(prevalence.t, "v24.75") &
                            str_detect(population.t, "p24.75") ~ "24.75"))

df.freq$year <- factor(df.freq$year)
order <- c("2.75", "Year 5", "Year 10", "Year 15", "Year 20", "24.75")
df.freq$year  <- ordered(df.freq$year, levels = order)

df.freq <- df.freq %>% 
  filter(! is.na(year))

df.freq$scenario <- NA
df.freq$scenario <- ifelse(df.freq$gamma.mm.1 == 0.5, "Female-Dominated", df.freq$scenario)
df.freq$scenario <- ifelse(df.freq$gamma.mm.1 == 4, "Male-Dominated", df.freq$scenario)
df.freq$scenario <- ifelse(df.freq$gamma.mm.1 == 2, "Equal", df.freq$scenario)

# write.csv(df.freq, "Tile_sim_FD.csv")
# write.csv(df.freq, "Tile_sim_FD.10.5.csv")
# df.freq <- read.csv("Tile_sim_FD.10.5.csv")

myColors <- c("#3288BD", "#ABDDA4", "#FFFFBF", "#D53E4F", NA)
names(myColors) <- c("Criteria Met", "Low Population", "High Prevalence", "Fails Critera", "Invalid")
colScalef <- scale_fill_manual(name = "Classification", values = myColors)

df.freq$year <- factor(df.freq$year)
order <- c("2.75", "Year 5", "Year 10", "Year 15", "Year 20", "24.75")
df.freq$year  <- ordered(df.freq$year, levels = order)

df.freq$scenario. <- factor(df.freq$scenario, ordered = T, levels = c("Equal", "Male-Dominated", "Female-Dominated"))

df.freq %>%
  filter(class != "Invalid" &
           year != "2.75" &
           year != "24.75" ) %>% 
  ggplot(aes(x = hunt.f, y = hunt.m, fill = class)) +
  geom_tile() +
  labs(x = "Adult Female Harvest Rate", 
       y = "Adult Male Harvest Rate") +
  facet_grid(scenario. ~ year) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  colScalef
beepr::beep()
# df.freq %>%
#   filter(class != "Invalid" &
#            year != "2.75" &
#            year != "24.75" ) %>% 
#   ggplot(aes(x = hunt.f, y = hunt.m, fill = class)) +
#   # geom_tile() +
#   geom_contour(aes(z = prevalence, linetype = factor(prevalence)), breaks = c(.1, .15, .2, .25)) +
#   # geom_contour(aes(z = population),
#   #              breaks = c(1000),
#   #              color = "Yellow") +
#   labs(x = "Female Hunting", 
#        y = "Male Hunting") +
#   facet_grid(scenario ~ year) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#   scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
#   colScalef

ggsave("Plots/Freq_Tile_Plots.cons.png", dpi=300, width = 8, height = 5)

invalid <- df.freq %>% 
  filter(class == "Invalid") %>% 
  group_by(year, scenario, hunt.f, hunt.m) 

breaks.prev <- c(0,0.05, .1, .15, .2, .25,.3,0.35, 0.4, 0.451)
breaks.pop <- c(250, 500, 750, 1000, 1250)

df.freq.m <- melt(df.freq, measur.vars = c("population", "prevalence"), id.vars = colnames(df.freq)[c(1:9,11,13,14)])
  
a.plot <- df.freq %>%
  filter(class != "Invalid" &
           year == "24.75" ) %>% 
  ggplot() +
  geom_contour_filled(aes(x = hunt.f, y = hunt.m, z = prevalence), breaks = c(seq(0, .35, 0.025),1)) +
  labs(x = "Adult Female Harvest Rate", 
       y = "Adult Male Harvest Rate") +
  facet_grid(. ~ scenario.) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
  scale_fill_viridis_d(name = "Prevalence", 
                       guide = guide_legend(keywidth = 0.5, keyheight = 0.25)) 

b<- df.freq %>%
  filter(class != "Invalid" &
           year == "Year 15" ) %>% 
  ggplot() +
  geom_contour_filled(aes(x = hunt.f, y = hunt.m, z = population), breaks = c(seq(0,2000, by = 250),7000)) +
  labs(x = "Adult Female Harvest Rate", 
       y = "Adult Male Harvest Rate") +
  facet_grid(year ~ scenario) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_viridis_d(name = "Population")

ggpubr::ggarrange(a, b, labels = c("A", "B"), nrow = 2)
ggsave("Plots/Freq_Contour_PlotA.conser.png", dpi=300, width = 8, height = 4.5)

require(metR)
require(tidyverse)
require(popbio)
df.freq <- read.csv("Tile_sim_FD.10.5.csv")

df.freq. <- df.freq %>% 
  filter(year == "2.75")
fawn.an.sur = 0.354
juv.an.sur = 0.883 
ad.an.f.sur = 0.972 
ad.an.m.sur = 0.968 
fawn.repro = 0
juv.repro = 0.198
ad.repro = .928
hunt.mort.fawn = 0.01
hunt.mort.juv.f = 0.05
hunt.mort.juv.m = 0.05
hunt.mort.ad.f = 0.05
hunt.mort.ad.m = 0.05
ini.fawn.prev = 0.02
ini.juv.prev = 0.03
ini.ad.f.prev = 0.04  
ini.ad.m.prev = 0.04
n.age.cats = 12  
p = 0.43 
env.foi = 0
theta = 1 
n0 = 2000 
n.years = 25
rel.risk = 1.0
fawn.sur <- fawn.an.sur^(1/12)
juv.sur <- juv.an.sur^(1/12) 
ad.f.sur <- ad.an.f.sur^(1/12) 
ad.m.sur <- ad.an.m.sur^(1/12)
p <- 0.43

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

eval.f <- rep(NA, ((2*10*n.age.cats)^2)) 
eval.v <- rep(NA, ((2*10*n.age.cats)^2))

df.freq.$R0 <- NA

for (i in 1:nrow(df.freq.)){
  beta.ff <- df.freq.$beta.ff[i]
  gamma.mm <- df.freq.$gamma.mm[i]
  beta.mm <- beta.ff * gamma.mm
  gamma.mf <- df.freq.$gamma.mf[i]
  beta.mf <- beta.ff * gamma.mf
  gamma.fm <- df.freq.$gamma.fm[i]
  beta.fm <- beta.ff * gamma.fm
  hunt.mort.juv.f = 0.05
  hunt.mort.juv.m = 0.05
  hunt.mort.ad.f = df.freq.$hunt.f[i]
  hunt.mort.ad.m = df.freq.$hunt.m[i]
  theta <- 1
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
  df.freq.$R0[i] <- max(Re(eigen(f.matrix %*% solve(v.matrix))$values))
}

write.csv(df.freq., "R0_Freq_Tile.csv")
df.freq. <- read.csv("R0_Freq_Tile.csv")

df.freq.$scenario <- factor(df.freq.$scenario, levels = c("Equal", "Male Dominated", "Female Dominated"),labels = c("Equal", "Male-Dominated", "Female-Dominated"), ordered = T)

# we need to white out the space where we know the population cannot sustain itself
df.freq.$R0 <- ifelse(df.freq.$class == "Invalid", NA, df.freq.$R0)

c.plot <- df.freq. %>%
  ggplot(aes(x = hunt.f, y = hunt.m, z = R0, fill= R0)) +
  geom_tile(alpha = 0.5) +
  stat_contour(color = "black",
               linetype = 4) +
  stat_contour(breaks = 1,
               color = "black",
               linetype = 1,
               size = 1) +
  labs(x = "Female Harvest Rate", 
       y = "Male Harvest Rate") +
  facet_grid(. ~ scenario) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_binned(type="viridis",na.value="white", n.breaks = 10)

ggsave("Plots/Freq_R0_Tile_Plots.png", dpi=300, width = 7, height = 3)

ggarrange(a.plot, c.plot, labels = c("A", "B"), nrow = 2, align = "v")
ggsave("Plots/Combined_Prev_R0.png", dpi=300, width = 7, height = 6.5)








