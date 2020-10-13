require(reshape2)
require(ggplot2)
require(dplyr)
require(CWDsims)
require(ggplot2)
require(reshape2)
require(tidyverse)
mat <- read.csv("R0.Scenarios.6.29.csv")

#female harvest styles
harvests <- c(0.00001, 0.05, 0.1, 0.15, 0.2)

df <- data.frame(B.ff = rep(rep(mat$B.ff[1:3], each = length(harvests)),2),
                 G.mm = rep(rep(mat$G.mm[1:3], each = length(harvests)),2),
                 G.mf = rep(rep(mat$G.mf[1:3], each = length(harvests)),2),
                 G.fm = rep(rep(mat$G.fm[1:3], each = length(harvests)),2))
df$harvest.m <- c(rep(harvests, 3), rep(0.05, length(harvests)*3))
df$harvest.f <- c(rep(0.05, length(harvests)*3), rep(harvests, 3))
df$harvest.sex <- c(rep("Male Harvests", length(harvests)*3), rep("Female Harvests", length(harvests)*3))

n <- 50
df <- do.call("rbind", replicate(n, df, simplify = FALSE))

df$y1 <- NA
df$y2 <- NA
df$y3 <- NA
df$y4 <- NA
df$y5 <- NA
df$y6 <- NA
df$y7 <- NA
df$y8 <- NA
df$y9 <- NA
df$y10 <- NA
df$y11 <- NA
df$y12 <- NA
df$y13 <- NA
df$y14 <- NA
df$y15 <- NA
df$y16 <- NA
df$y17 <- NA
df$y18 <- NA
df$y19 <- NA
df$y20 <- NA
df$y21 <- NA
df$y22 <- NA
df$y23 <- NA
df$y24 <- NA
df$y25 <- NA

df$yn1 <- NA
df$yn2 <- NA
df$yn3 <- NA
df$yn4 <- NA
df$yn5 <- NA
df$yn6 <- NA
df$yn7 <- NA
df$yn8 <- NA
df$yn9 <- NA
df$yn10 <- NA
df$yn11 <- NA
df$yn12 <- NA
df$yn13 <- NA
df$yn14 <- NA
df$yn15 <- NA
df$yn16 <- NA
df$yn17 <- NA
df$yn18 <- NA
df$yn19 <- NA
df$yn20 <- NA
df$yn21 <- NA
df$yn22 <- NA
df$yn23 <- NA
df$yn24 <- NA
df$yn25 <- NA

df$ybd1 <- NA
df$ybd2 <- NA
df$ybd3 <- NA
df$ybd4 <- NA
df$ybd5 <- NA
df$ybd6 <- NA
df$ybd7 <- NA
df$ybd8 <- NA
df$ybd9 <- NA
df$ybd10 <- NA
df$ybd11 <- NA
df$ybd12 <- NA
df$ybd13 <- NA
df$ybd14 <- NA
df$ybd15 <- NA
df$ybd16 <- NA
df$ybd17 <- NA
df$ybd18 <- NA
df$ybd19 <- NA
df$ybd20 <- NA
df$ybd21 <- NA
df$ybd22 <- NA
df$ybd23 <- NA
df$ybd24 <- NA
df$ybd25 <- NA

df$ycc1 <- NA
df$ycc2 <- NA
df$ycc3 <- NA
df$ycc4 <- NA
df$ycc5 <- NA
df$ycc6 <- NA
df$ycc7 <- NA
df$ycc8 <- NA
df$ycc9 <- NA
df$ycc10 <- NA
df$ycc11 <- NA
df$ycc12 <- NA
df$ycc13 <- NA
df$ycc14 <- NA
df$ycc15 <- NA
df$ycc16 <- NA
df$ycc17 <- NA
df$ycc18 <- NA
df$ycc19 <- NA
df$ycc20 <- NA
df$ycc21 <- NA
df$ycc22 <- NA
df$ycc23 <- NA
df$ycc24 <- NA
df$ycc25 <- NA

a <- list(fawn.an.sur = 0.354,
          juv.an.sur = 0.883 ,
          ad.an.f.sur = 0.972 ,
          ad.an.m.sur = 0.968, 
          fawn.repro = 0,
          juv.repro = 0.198,
          ad.repro = .928,
          hunt.mort.fawn = 0.01,
          hunt.mort.juv.f = 0.05,
          hunt.mort.juv.m = 0.05,
          hunt.mort.ad.f = 0.05,
          hunt.mort.ad.m = 0.05,
          ini.fawn.prev = 0.02,
          ini.juv.prev = 0.03,
          ini.ad.f.prev = 0.04  ,
          ini.ad.m.prev = 0.04,
          n.age.cats = 12  ,
          p = 0.43 ,
          env.foi = 0 ,
          theta = 1 ,
          n0 = 2000 ,
          n.years = 25,
          rel.risk = 1.0,
          repro.var = 0.005, # 0.00126, 
          fawn.sur.var = 0.005, # 0.03854, 
          sur.var = 0.005, # 0.03854 , 
          hunt.var = 0.005, # 0.05,  
          p = 0.43, 
          env.foi = 0, 
          beta.ff = 0.07,
          gamma.mm = 1,
          gamma.mf = 1,
          gamma.fm = 1,
          rel.risk = 1)

set.seed(123456)
for (i in 1:nrow(df)){
  tryCatch({
  params <- a
  params[["beta.ff"]] <- df$B.ff[i]
  params[["gamma.mm"]] <- df$G.mm[i]
  params[["gamma.mf"]] <- df$G.mf[i]
  params[["gamma.fm"]] <- df$G.fm[i]
  params[["hunt.mort.ad.m"]] <- df$harvest.m[i]
  params[["hunt.mort.ad.f"]] <- df$harvest.f[i]
  
  out <- cwd_stoch_model_wiw(params)
  
  prev.sum   <- out$counts %>%
    filter(month %% 12 == 10) %>%
    group_by(year, disease) %>%
    summarize(n = sum(population)) %>%
    pivot_wider(names_from = disease, values_from = n) %>%
    mutate(prev = yes/(no + yes), pop = no + yes) #Prev at time
  df[i, 8:32] <- prev.sum$prev
  df[i, 33:57] <- prev.sum$pop
  bdratio <- out$counts %>%
    filter(month %% 12 == 10 &
             age > 1) %>%
    group_by(year, sex) %>%
    summarize(n = sum(population)) %>%
    pivot_wider(names_from = sex, values_from = n) %>%
    mutate(bd = m/f)
  df[i, 58:82] <- bdratio$bd
  ccratio <- out$counts %>%
    filter(month %% 12 == 10) %>%
    mutate(cow.calf = ifelse(age==1, "calf", NA)) %>% 
    mutate(cow.calf = ifelse(age>1 &
                               sex=="f", "cow", cow.calf)) %>% 
    group_by(year, cow.calf) %>%
    summarize(n = sum(population)) %>%
    pivot_wider(names_from = cow.calf, values_from = n) %>%
    mutate(cc = calf/cow)
  df[i, 83:107] <- ccratio$cc
  }, error = function(e){print(paste(i, "crashed"))})
}

df$group <- 1:nrow(df)
df$trans <- ifelse(df$G.mm == 0.5, "Female-Dominated", NA)
df$trans <- ifelse(df$G.mm == 4, "Male-Dominated", df$trans)
df$trans <- ifelse(df$G.mm == 2, "Equal", df$trans)
melted <- melt(df, id.vars = colnames(df[,c(1:7,108:109)]), measure.vars = colnames(df[,8:107]))
melted$year <- sub("n", "", melted$variable)
melted$year <- sub("y", "", melted$year)
melted$year <- sub("bd", "", melted$year)
melted$year <- sub("cc", "", melted$year)
melted$year <- as.numeric(melted$year)
melteda <- melted[1:(nrow(df)*25),]
meltedb <- melted[((nrow(df)*25)+1):(nrow(df)*25*2),]
meltedc <- melted[((nrow(df)*25*2)+1):(nrow(df)*25*3),]
meltedd <- melted[((nrow(df)*25*3)+1):(nrow(df)*25*4),]
melted.1 <- merge(melteda, meltedb, by = colnames(melted[,c(1:9,12)]))
melted.2 <- merge(meltedc, meltedd, by = colnames(melted[,c(1:9,12)]))
melted <- merge(melted.1, melted.2, by = colnames(melted[,c(1:9,12)]))
melted$value.x <- ifelse(melted$value.y.x < 10, NA, melted$value.x.x)

melted$harvest.rate <- ifelse(melted$harvest.m == 0.05, melted$harvest.f, melted$harvest.m)
melted$harvest.rate <- as.factor(melted$harvest.rate)
melted$trans <- as.factor(melted$trans)

melted$value.x.x <- ifelse(melted$harvest.rate == 0.2 &
                             melted$harvest.sex == "Female Harvests",
                           NA, 
                           melted$value.x.x)

melted$harvest.rate <- factor(melted$harvest.rate, labels = c("0%", "5%", "10%", "15%", "20%"))
  
melted <- melted %>% 
  group_by(year, trans, harvest.sex, harvest.rate) %>% 
  mutate(mean = mean(value.x.x, na.rm = T), 
         lower = mean(value.x.x, na.rm = T) - sd(value.x.x, na.rm = T),
         upper = mean(value.x.x, na.rm = T) + sd(value.x.x, na.rm = T))

melted %>% 
  filter(harvest.rate != "20%") %>% 
  ggplot(aes(year, value.x.x, color = trans, fill = trans)) +
  # geom_line(alpha=0.05, aes(group = group)) +
  geom_line(aes(y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, color = NA) +
  facet_grid(harvest.sex~harvest.rate) +
  labs(x = "Year",
       y = "Prevalence") +
  scale_color_discrete(name = "Transmission") +
  scale_fill_discrete(name = "Transmission") +
  theme_classic() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

melted %>% 
  ggplot(aes(year, value.x.y, color = harvest.rate, fill = harvest.rate)) +
  geom_line(alpha=0.05, aes(group = group)) +
  geom_smooth(se = 0.95) +
  facet_grid(harvest.sex~trans) +
  labs(x = "Year",
       y = "Male:Female Ratio") +
  scale_color_discrete(name = "Harvest Rate",
                       labels = c("1%", "7.5%", "15%")) +
  scale_fill_discrete(name = "Harvest Rate",
                      labels = c("1%", "7.5%", "15%")) +
  theme_classic()

melted %>% 
  ggplot(aes(year, value.y.y, color = harvest.rate, fill = harvest.rate)) +
  geom_line(alpha=0.05, aes(group = group)) +
  geom_smooth(se = 0.95) +
  facet_grid(harvest.sex~trans) +
  labs(x = "Year",
       y = "Young of Year:Adult Female Ratio") +
  scale_color_discrete(name = "Harvest Rate",
                       labels = c("1%", "7.5%", "15%")) +
  scale_fill_discrete(name = "Harvest Rate",
                      labels = c("1%", "7.5%", "15%")) +
  ylim (0,1) +
  theme_classic()
  
ggsave("/Users/willrogers/Downloads/cwdsims-master/Plots/Hunting_Plots_sd.png", dpi=300, width = 7, height = 4)



