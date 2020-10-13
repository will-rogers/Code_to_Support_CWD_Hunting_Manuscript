params <- list(fawn.an.sur = 0.6, juv.an.sur = 0.8, ad.an.f.sur = 0.95,
ad.an.m.sur = 0.9, fawn.repro = 0, juv.repro = 0.6, ad.repro = 1,
hunt.mort.fawn = 0.01, hunt.mort.juv.f = 0.1, hunt.mort.juv.m = 0.1,
hunt.mort.ad.f = 0.1, hunt.mort.ad.m = 0.2, ini.fawn.prev = 0.02,
ini.juv.prev = 0.03, ini.ad.f.prev = 0.04,  ini.ad.m.prev = 0.04,
n.age.cats = 12,  p = 0.43, env.foi = 0.0,  beta.ff = 0.045,
gamma.mm = 2.3, gamma.mf = 2.3, gamma.fm = 1,
theta = 1, n0 = 10000, n.years = 50, rel.risk = 1.0)
out <- cwd_wiw_track_R0(params = params)

a <- plot_track_sex_trans(out$tracking.inf, T)
b <- plot_R0(out$R0)

ggarrange(a,b, ncol=2, nrow = 1)

df <- matrix(NA, ncol = 8, nrow = 3)
colnames(df) <- c("Next.Gen", "Fin.Prev", "Min.r", "Max.r", "B.ff", "G.fm", "G.mf", "G.mm")
rownames(df) <- c("Equal", "Male.Dom", "Female.Dom")
df[1,1] <- as.numeric(out$R0[[1]][1,1])
df[1,2] <- as.numeric(out$R0[[1]][2,1])
df[1,3] <- as.numeric(out$R0[[2]][1,1])
df[1,4] <- as.numeric(out$R0[[2]][1,2])
df[1,5] <- params$beta.ff
df[1,6] <- params$gamma.fm
df[1,7] <- params$gamma.mf
df[1,8] <- params$gamma.mm
write.csv(df, file = "R0_validated.csv")

 
