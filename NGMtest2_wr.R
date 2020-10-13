fawn.an.sur = 0.6
juv.an.sur = 0.8
ad.an.f.sur = 0.95 
ad.an.m.sur = 0.9
fawn.repro = 0
juv.repro = 0.6
ad.repro = 1 

hunt.mort.fawn = 0.01/12
hunt.mort.juv.f = 0.1/12
hunt.mort.juv.m = 0.1/12
hunt.mort.ad.f = 0.1/12
hunt.mort.ad.m = 0.2/12

p = 0.075
env.foi = 0
beta.ff = 0.06^(1/12)
gamma.mm = 1
gamma.mf = 1
gamma.fm = 1
theta = 1
rel.risk = 1.0

beta.mm <- beta.ff * gamma.mm
beta.mf <- beta.ff * gamma.mf
beta.fm <- beta.ff * gamma.fm

fawn.sur <- fawn.an.sur^(1/12)
juv.sur <- juv.an.sur^(1/12)
ad.f.sur <- ad.an.f.sur^(1/12)
ad.m.sur <- ad.an.m.sur^(1/12)

# PCC: note that the starting conditions should be everyone is susceptible
N <- 2000

I1f1 <- 0
I2f1 <- 0
I3f1 <- 0
I1m1 <- 0
I2m1 <- 0
I3m1 <- 0

I1f2 <- 0
I2f2 <- 0
I3f2 <- 0
I1m2 <- 0
I2m2 <- 0
I3m2 <- 0

S1f <- N/4
S2f <- N/8
S3f <- N/8

S1m <- N/4
S2m <- N/8
S3m <- N/8

#scaling transmission by the total population size not the subcategory...need to check this in the main code. 
# need a negative beta inserted, changed parantheses
# env.foi goes into the exp? I think this disappears in the derivative as not a function of I. 

F1f1 <- quote(S1f * (1 - exp(-(((beta.ff*(I1f1+I1f2+I2f1+I2f2+I3f1+I3f2))+
                                 (beta.mf*(I1m1+I1m2+I2m1+I2m2+I3m1+I3m2)))/(N^theta)))))
F2f1 <- quote(S2f * (1 - exp(-(((beta.ff*(I1f1+I1f2+I2f1+I2f2+I3f1+I3f2))+
                                 (beta.mf*(I1m1+I1m2+I2m1+I2m2+I3m1+I3m2)))/(N^theta)))))
F3f1 <- quote(S3f * (1 - exp(-(((beta.ff*(I1f1+I1f2+I2f1+I2f2+I3f1+I3f2))+
                                 (beta.mf*(I1m1+I1m2+I2m1+I2m2+I3m1+I3m2)))/(N^theta)))))
F1m1 <- quote(S1m * (1 - exp(-(((beta.mm*(I1f1+I1f2+I2f1+I2f2+I3f1+I3f2))+
                                 (beta.fm*(I1m1+I1m2+I2m1+I2m2+I3m1+I3m2)))/(N^theta)))))
F2m1 <- quote(S2m * (1 - exp(-(((beta.mm*(I1f1+I1f2+I2f1+I2f2+I3f1+I3f2))+
                                 (beta.fm*(I1m1+I1m2+I2m1+I2m2+I3m1+I3m2)))/(N^theta)))))
F3m1 <- quote(S3m * (1 - exp(-(((beta.mm*(I1f1+I1f2+I2f1+I2f2+I3f1+I3f2))+
                                 (beta.fm*(I1m1+I1m2+I2m1+I2m2+I3m1+I3m2)))/(N^theta)))))

F1f2 <- 0
F2f2 <- 0
F3f2 <- 0
F1m2 <- 0
F2m2 <- 0
F3m2 <- 0

#all loses
#  PCC: changed survival rates to account for the different ages. 
# NOTE: survival rates changed to a monthly timescale
# NOTE fixed some m-f errors here

V1f1 <- quote(I1f1*(hunt.mort.fawn+(1-fawn.sur) + p))
V2f1 <- quote(I2f1*(hunt.mort.juv.f+(1-juv.sur) + p))
V3f1 <- quote(I3f1*(hunt.mort.ad.f+(1-ad.f.sur) + p))

V1m1 <- quote(I1m1*(hunt.mort.fawn+(1-fawn.sur) + p))
V2m1 <- quote(I2m1*(hunt.mort.juv.m+(1-juv.sur) + p))
V3m1 <- quote(I3m1*(hunt.mort.ad.m+(1-ad.m.sur) + p))

V1f2 <- quote(I1f2*(hunt.mort.fawn+(1-fawn.sur) + p))
V2f2 <- quote(I2f2*(hunt.mort.juv.f+(1-juv.sur) + p))
V3f2 <- quote(I3f2*(hunt.mort.ad.f+(1-ad.f.sur) + p))

V1m2 <- quote(I1m2*(hunt.mort.fawn+(1-fawn.sur) + p))
V2m2 <- quote(I2m2*(hunt.mort.juv.m+(1-juv.sur) + p))
V3m2 <- quote(I3m2*(hunt.mort.ad.m+(1-ad.m.sur) + p))

# gained transfers
V.1f1 <- 0
V.2f1 <- 0
V.3f1 <- 0
V.1m1 <- 0
V.2m1 <- 0
V.3m1 <- 0
V.1f2 <- quote(I1f1*p)
V.2f2 <- quote(I2f1*p)
V.3f2 <- quote(I3f1*p)
V.1m2 <- quote(I1m1*p)
V.2m2 <- quote(I2m1*p)
V.3m2 <- quote(I3m1*p) #fixed error here

# matrix is ordered as: (f11,f21,f31, m11,m21,31, f12,f22,f32, m12,m22, m32)
V1 = substitute(a - b, list(a = V1f1, b = V.1f1))
V2 = substitute(a - b, list(a = V2f1, b = V.2f1))
V3 = substitute(a - b, list(a = V3f1, b = V.3f1))
V4 = substitute(a - b, list(a = V1m1, b = V.1m1))
V5 = substitute(a - b, list(a = V2m1, b = V.2m1))
V6 = substitute(a - b, list(a = V3m1, b = V.3m1))

V7 = substitute(a - b, list(a = V1f2, b = V.1f2))
V8 = substitute(a - b, list(a = V2f2, b = V.2f2))
V9 = substitute(a - b, list(a = V3f2, b = V.3f2))
V10 = substitute(a - b, list(a = V1m2, b = V.1m2))
V11 = substitute(a - b, list(a = V2m2, b = V.2m2))
V12 = substitute(a - b, list(a = V3m2, b = V.3m2))

f11 = D(F1f1, "I1f1");
f12 = D(F1f1, "I2f1");
f13 = D(F1f1, "I3f1");
f14 = D(F1f1, "I1f2");
f15 = D(F1f1, "I2f2");
f16 = D(F1f1, "I3f2");

f17 = D(F1f1, "I1m1");
f18 = D(F1f1, "I1m2");
f19 = D(F1f1, "I2m1");
f110 = D(F1f1, "I2m2");
f111 = D(F1f1, "I3m1");
f112 = D(F1f1, "I3m2")

# simplifying (only fecundities on the top row)
f <- matrix(0, nrow = 12, ncol = 12)
f[1,] <-c(eval(f11),eval(f12),eval(f13),eval(f14),eval(f15),eval(f16),
              eval(f17),eval(f18),eval(f19),eval(f110),eval(f111),eval(f112))

# check that f entries are positive. 
f

###########
# PCC: changed the ordering to match f 

v11 = D(V1, "I1f1");v12 = D(V1, "I2f1");v13 = D(V1, "I3f1")
v14 = D(V1, "I1m1");v15 = D(V1, "I2m1");v16 = D(V1, "I3m1")
v17 = D(V1, "I1f2");v18 = D(V1, "I2f2");v19 = D(V1, "I3f2")
v110 = D(V1, "I1m2");v111 = D(V1, "I2m2");v112 = D(V1, "I3m2")

v21 = D(V2, "I1f1");v22 = D(V2, "I2f1");v23 = D(V2, "I3f1")
v24 = D(V2, "I1m1");v25 = D(V2, "I2m1");v26 = D(V2, "I3m1")
v27 = D(V2, "I1f2");v28 = D(V2, "I2f2");v29 = D(V2, "I3f2")
v210 = D(V2, "I1m2");v211 = D(V2, "I2m2");v212 = D(V2, "I3m2")

v31 = D(V3, "I1f1");v32 = D(V3, "I2f1");v33 = D(V3, "I3f1")
v34 = D(V3, "I1m1");v35 = D(V3, "I2m1");v36 = D(V3, "I3m1")
v37 = D(V3, "I1f2");v38 = D(V3, "I2f2");v39 = D(V3, "I3f2")
v310 = D(V3, "I1m2");v311 = D(V3, "I2m2");v312 = D(V3, "I3m2")

v41 = D(V4, "I1f1");v42 = D(V4, "I2f1");v43 = D(V4, "I3f1")
v44 = D(V4, "I1m1");v45 = D(V4, "I2m1");v46 = D(V4, "I3m1")
v47 = D(V4, "I1f2");v48 = D(V4, "I2f2");v49 = D(V4, "I3f2")
v410 = D(V4, "I1m2");v411 = D(V4, "I2m2");v412 = D(V4, "I3m2")

v51 = D(V5, "I1f1");v52 = D(V5, "I2f1");v53 = D(V5, "I3f1")
v54 = D(V5, "I1m1");v55 = D(V5, "I2m1");v56 = D(V5, "I3m1")
v57 = D(V5, "I1f2");v58 = D(V5, "I2f2");v59 = D(V5, "I3f2")
v510 = D(V5, "I1m2");v511 = D(V5, "I2m2");v512 = D(V5, "I3m2")

v61 = D(V6, "I1f1");v62 = D(V6, "I2f1");v63 = D(V6, "I3f1")
v64 = D(V6, "I1m1");v65 = D(V6, "I2m1");v66 = D(V6, "I3m1")
v67 = D(V6, "I1f2");v68 = D(V6, "I2f2");v69 = D(V6, "I3f2")
v610 = D(V6, "I1m2");v611 = D(V6, "I2m2");v612 = D(V6, "I3m2")

v71 = D(V7, "I1f1");v72 = D(V7, "I2f1");v73 = D(V7, "I3f1")
v74 = D(V7, "I1m1");v75 = D(V7, "I2m1");v76 = D(V7, "I3m1")
v77 = D(V7, "I1f2");v78 = D(V7, "I2f2");v79 = D(V7, "I3f2")
v710 = D(V7, "I1m2");v711 = D(V7, "I2m2");v712 = D(V7, "I3m2")

v81 = D(V8, "I1f1");v82 = D(V8, "I2f1");v83 = D(V8, "I3f1")
v84 = D(V8, "I1m1");v85 = D(V8, "I2m1");v86 = D(V8, "I3m1")
v87 = D(V8, "I1f2");v88 = D(V8, "I2f2");v89 = D(V8, "I3f2")
v810 = D(V8, "I1m2");v811 = D(V8, "I2m2");v812 = D(V8, "I3m2")

v91 = D(V9, "I1f1");v92 = D(V9, "I2f1");v93 = D(V9, "I3f1")
v94 = D(V9, "I1m1");v95 = D(V9, "I2m1");v96 = D(V9, "I3m1")
v97 = D(V9, "I1f2");v98 = D(V9, "I2f2");v99 = D(V9, "I3f2")
v910 = D(V9, "I1m2");v911 = D(V9, "I2m2");v912 = D(V9, "I3m2")

v101 = D(V10, "I1f1");v102 = D(V10, "I2f1");v103 = D(V10, "I3f1")
v104 = D(V10, "I1m1");v105 = D(V10, "I2m1");v106 = D(V10, "I3m1")
v107 = D(V10, "I1f2");v108 = D(V10, "I2f2");v109 = D(V10, "I3f2")
v1010 = D(V10, "I1m2");v1011 = D(V10, "I2m2");v1012 = D(V10, "I3m2")

v111 = D(V11, "I1f1");v112 = D(V11, "I2f1");v113 = D(V11, "I3f1")
v114 = D(V11, "I1m1");v115 = D(V11, "I2m1");v116 = D(V11, "I3m1")
v117 = D(V11, "I1f2");v118 = D(V11, "I2f2");v119 = D(V11, "I3f2")
v1110 = D(V11, "I1m2");v1111 = D(V11, "I2m2");v1112 = D(V11, "I3m2")

v121 = D(V12, "I1f1");v122 = D(V12, "I2f1");v123 = D(V12, "I3f1")
v124 = D(V12, "I1m1");v125 = D(V12, "I2m1");v126 = D(V12, "I3m1")
v127 = D(V12, "I1f2");v128 = D(V12, "I2f2");v129 = D(V12, "I3f2")
v1210 = D(V12, "I1m2");v1211 = D(V12, "I2m2");v1212 = D(V12, "I3m2")

v <- matrix(c(eval(v11),eval(v12),eval(v13),eval(v14),eval(v15),eval(v16),
              eval(v17),eval(v18),eval(v19),eval(v110),eval(v111),eval(v112),
              eval(v21),eval(v22),eval(v23),eval(v24),eval(v25),eval(v26),
              eval(v27),eval(v28),eval(v29),eval(v210),eval(v211),eval(v212),
              eval(v31),eval(v32),eval(v33),eval(v34),eval(v35),eval(v36),
              eval(v37),eval(v38),eval(v39),eval(v310),eval(v311),eval(v312),
              eval(v41),eval(v42),eval(v43),eval(v44),eval(v45),eval(v46),
              eval(v47),eval(v48),eval(v49),eval(v410),eval(v411),eval(v412),
              eval(v51),eval(v52),eval(v53),eval(v54),eval(v55),eval(v56),
              eval(v57),eval(v58),eval(v59),eval(v510),eval(v511),eval(v512),
              eval(v61),eval(v62),eval(v63),eval(v64),eval(v65),eval(v66),
              eval(v67),eval(v68),eval(v69),eval(v610),eval(v611),eval(v612),
              eval(v71),eval(v72),eval(v73),eval(v74),eval(v75),eval(v76),
              eval(v77),eval(v78),eval(v79),eval(v710),eval(v711),eval(v712),
              eval(v81),eval(v82),eval(v83),eval(v84),eval(v85),eval(v86),
              eval(v87),eval(v88),eval(v89),eval(v810),eval(v811),eval(v812),
              eval(v91),eval(v92),eval(v93),eval(v94),eval(v95),eval(v96),
              eval(v97),eval(v98),eval(v99),eval(v910),eval(v911),eval(v912),
              eval(v101),eval(v102),eval(v103),eval(v104),eval(v105),eval(v106),
              eval(v107),eval(v108),eval(v109),eval(v1010),eval(v1011),eval(v1012),
              eval(v111),eval(v112),eval(v113),eval(v114),eval(v115),eval(v116),
              eval(v117),eval(v118),eval(v119),eval(v1110),eval(v1111),eval(v1112),
              eval(v121),eval(v122),eval(v123),eval(v124),eval(v125),eval(v126),
              eval(v127),eval(v128),eval(v129),eval(v1210),eval(v1211),eval(v1212)), nrow = 12, byrow=T)

v
det(v) #should not be 0 
max(eigen(f %*% solve(v))$values)

