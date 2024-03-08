#### Load in libraries ####

library(nleqslv)
library(tidyverse)
library(cowplot)
library(ggpmisc)
library(ggpubr)
devtools::load_all()

v.c = df.conc$Concentration
names(v.c) = df.conc$Metabolites

####Write a solving function####

fn.solving = function(v.c = df.conc$Concentration){
  names(v.c) = df.conc$Metabolites

  fn = function(z) {
    x = z^2
    y = c()

    #Kd Mg
    y[1] = log10(x[1]*x[5]/(Kd.app[1]*x[13]))
    y[2] = log10(x[1]*x[6]/(Kd.app[2]*x[14]))
    y[3] = log10(x[1]*x[7]/(Kd.app[3]*x[15]))
    y[4] = log10(x[1]*x[8]/(Kd.app[4]*x[16]))
    y[5] = log10(x[1]*x[9]/(Kd.app[5]*x[17]))
    y[6] = log10(x[1]*x[10]/(Kd.app[6]*x[18]))
    y[7] = log10(x[1]*x[11]/(Kd.app[7]*x[19]))
    y[8] = log10(x[1]*x[12]/(Kd.app[8]*x[20]))

    #Kd Mn
    y[9]  = log10(x[2]*x[5]/(Kd.app[9]*x[21]))
    y[10] = log10(x[2]*x[6]/(Kd.app[10]*x[22]))
    y[11] = log10(x[2]*x[7]/(Kd.app[11]*x[23]))
    y[12] = log10(x[2]*x[8]/(Kd.app[12]*x[24]))
    y[13] = log10(x[2]*x[9]/(Kd.app[13]*x[25]))
    y[14] = log10(x[2]*x[10]/(Kd.app[14]*x[26]))
    y[15] = log10(x[2]*x[11]/(Kd.app[15]*x[27]))
    y[16] = log10(x[2]*x[12]/(Kd.app[16]*x[28]))


    #Kd Zn
    y[17] = log10(x[3]*x[5]/(Kd.app[17]*x[29]))
    y[18] = log10(x[3]*x[6]/(Kd.app[18]*x[30]))
    y[19] = log10(x[3]*x[7]/(Kd.app[19]*x[31]))
    y[20] = log10(x[3]*x[8]/(Kd.app[20]*x[32]))
    y[21] = log10(x[3]*x[9]/(Kd.app[21]*x[33]))
    y[22] = log10(x[3]*x[10]/(Kd.app[22]*x[34]))
    y[23] = log10(x[3]*x[11]/(Kd.app[23]*x[35]))
    y[24] = log10(x[3]*x[12]/(Kd.app[24]*x[36]))

    #Kd Ca
    y[25] = log10(x[4]*x[5]/(Kd.app[25]*x[37]))
    y[26] = log10(x[4]*x[6]/(Kd.app[26]*x[38]))
    y[27] = log10(x[4]*x[7]/(Kd.app[27]*x[39]))
    y[28] = log10(x[4]*x[8]/(Kd.app[28]*x[40]))
    y[29] = log10(x[4]*x[9]/(Kd.app[29]*x[41]))
    y[30] = log10(x[4]*x[10]/(Kd.app[30]*x[42]))
    y[31] = log10(x[4]*x[11]/(Kd.app[31]*x[43]))
    y[32] = log10(x[4]*x[12]/(Kd.app[32]*x[44]))

    #Total M2+
    y[33] = (x[1] + x[13] + x[14] + x[15] + x[16] + x[17] + x[18] + x[19] + x[20] - v.c[9])        #Mg
    y[34] = 10*(x[2] + x[21] + x[22] + x[23] + x[24] + x[25] + x[26] + x[27] + x[28] - v.c[10])    #Mn
    y[35] = 100*(x[3] + x[29] + x[30] + x[31] + x[32] + x[33] + x[34] + x[35] + x[36] - v.c[11])  #Zn
    y[36] = 10000*(x[4] + x[37] + x[38] + x[39] + x[40] + x[41] + x[42] + x[43] + x[44] - v.c[12])  #Ca


    #Total metabolites
    y[37] = x[5]  + x[13] + x[21] + x[29] + x[37] - v.c[1]
    y[38] = x[6]  + x[14] + x[22] + x[30] + x[38] - v.c[2]
    y[39] = x[7]  + x[15] + x[23] + x[31] + x[39] - v.c[3]
    y[40] = x[8]  + x[16] + x[24] + x[32] + x[40] - v.c[4]
    y[41] = x[9]  + x[17] + x[25] + x[33] + x[41] - v.c[5]
    y[42] = x[10] + x[18] + x[26] + x[34] + x[42] - v.c[6]
    y[43] = x[11] + x[19] + x[27] + x[35] + x[43] - v.c[7]
    y[44] = x[12] + x[20] + x[28] + x[36] + x[44] - v.c[8]

    return(y)
  }

  ####Names of each species in x####
  v.names = c("Mg",         #1
              "Mn",                   #2
              "Zn",                   #3
              "Ca",                   #4
              "L-Aspartic acid",      #5
              "Glutathione",          #6
              "Glucose 1-P",          #7
              "ATP",                  #8
              "AMP",                  #9
              "L-Alanine",            #10
              "L-Asparagine",         #11
              "Pyruvic acid",         #12
              "L-Aspartic acid_Mg",   #13
              "Glutathione_Mg",       #14
              "Glucose 1-P_Mg",       #15
              "ATP_Mg",               #16
              "AMP_Mg",               #17
              "L-Alanine_Mg",         #18
              "L-Asparagine_Mg",      #19
              "Pyruvic acid_Mg",      #20
              "L-Aspartic acid_Mn",   #21
              "Glutathione_Mn",       #22
              "Glucose 1-P_Mn",       #23
              "ATP_Mn",               #24
              "AMP_Mn",               #25
              "L-Alanine_Mn",         #26
              "L-Asparagine_Mn",      #27
              "Pyruvic acid_Mn",      #28
              "L-Aspartic acid_Zn",   #29
              "Glutathione_Zn",       #30
              "Glucose 1-P_Zn",       #31
              "ATP_Zn",               #32
              "AMP_Zn",               #33
              "L-Alanine_Zn",         #34
              "L-Asparagine_Zn",      #35
              "Pyruvic acid_Zn",      #36
              "L-Aspartic acid_Ca",   #37
              "Glutathione_Ca",       #38
              "Glucose 1-P_Ca",       #39
              "ATP_Ca",               #40
              "AMP_Ca",               #41
              "L-Alanine_Ca",         #42
              "L-Asparagine_Ca",      #43
              "Pyruvic acid_Ca")      #44

  #### Starting variables for the solution ####

  v.start = c(2.7,                    #1  "Mg"
              0.07,                   #2  "Mn"
              0.00071,                #3  "Zn"
              0.000078,               #4  "Ca"
              100,                    #5  "L-Aspartic acid"
              17,                   #6  "Glutathione"
              18,                     #7  "Glucose 1-P"
              0.64,                      #8  "ATP"
              7.3,                      #9  "AMP"
              4,                      #10 "L-Alanine"
              3.8,                    #11 "L-Asparagine"
              3.5,                      #12 "Pyruvic acid"
              0.15,                      #13 "L-Asp_Mg"
              0.00043,                  #14 "Glutathione_Mg"
              12,                      #15 "Glucose 1-P_Mg"
              23,                     #16 "ATP_Mg"
              1.8,                      #17 "AMP_Mg"
              0.0014,                  #18 "L-Alanine_Mg"
              0.00001,                  #19 "L-Asparagine_Mg"
              0.14,                      #20 "Pyruvic acid_Mg"
              0.072,                   #21 "L-Aspartic acid_Mn"
              0.0000011,                   #22 "Glutathione_Mn"
              0.16,                    #23 "Glucose 1-P_Mn"
              3.6,                    #24 "ATP_Mn"
              0.12,                   #25 "AMP_Mn"
              0.00014,                   #26 "L-Alanine_Mn"
              0.000015,                   #27 "L-Asparagine_Mn"
              0.0032,                   #28 "Pyruvic acid_Mn"
              0.11,                 #29 "L-Aspartic acid_Zn"
              0.044,                 #30 "Glutathione_Zn"
              0.0023,                 #31 "Glucose 1-P_Zn"
              0.085,                 #32 "ATP_Zn"
              0.0022,                 #33 "AMP_Zn"
              0.0000036,                 #34 "L-Alanine_Zn"
              0.000023,                 #35 "L-Asparagine_Zn"
              0.000038,                 #36 "Pyruvic acid_Zn"
              0.00000088,                #37 "L-Aspartic acid_Ca"
              0.00000034,                #38 "Glutathione_Ca"
              0.00011,                #39 "Glucose 1-P_Ca"
              0.00059,                #40 "ATP_Ca"
              0.000051,                #41 "AMP_Ca"
              0.00000011,                #42 "L-Alanine_Ca"
              0.00000012,                #43 "L-Asparagine_Ca"
              0.0000026)                #44 "Pyruvic acid_Ca"

  #### Solve the equations ####

  p = nleqslv(sqrt(v.start),
              fn,
              control = list(xtol = 10^-20,
                             ftol = 10^-20,
                             maxit = 1000000,
                             allowSingular = T))
  ####Format and compile the results####

  v.s = p$x^2
  names(v.s) = v.names

  M = c()
  L = c()
  Conc = c()

  for (i in 1:length(v.s)){
    #print(i)
    v.name = strsplit(names(v.s[i]), split = "_")[[1]]
    if (length(v.name) > 1){
      M[i] = v.name[2]
      L[i] = v.name[1]
      Conc[i] = v.s[i]
    }else{
      if (v.name %in% c("Mg", "Mn", "Zn", "Ca")){
        M[i] = v.name[1]
        L[i] = "free"
        Conc[i] = v.s[i]
      }else{
        M[i] = "free"
        L[i] = v.name
        Conc[i] = v.s[i]
      }
    }
  }

  df = data.frame(M, L, Conc)

  M = c("free", "Mg", "Mn", "Zn", "Ca")
  L = c("free",
        "L-Aspartic acid",      #4
        "Glutathione",          #5
        "Glucose 1-P",          #6
        "ATP",                  #7
        "AMP",                  #8
        "L-Alanine",            #9
        "L-Asparagine",         #10
        "Pyruvic acid")

  df$L = factor(df$L, levels = L)
  df$M = factor(df$M, levels = M)

  df$Mg.T = v.c[9]
  df$Mn.T = v.c[10]
  df$Zn.T = v.c[11]
  df$Ca.T = v.c[12]

  return(df)
}

####Test solving function####

df.test = fn.solving()
print(df.test)

####Define a list of perturbed systems####

#Mg

log(v.c[9], 2)

v.Mg.T = 2^seq(log(v.c[9], 2)-0.5, log(v.c[9], 2)+0.5, length.out = 10)
v.Mn.T = 2^seq(log(v.c[10], 2)-0.5, log(v.c[10], 2)+0.5, length.out = 10)
v.Zn.T = 2^seq(log(v.c[11], 2)-0.5, log(v.c[11], 2)+0.5, length.out = 10)
v.Ca.T = 2^seq(log(v.c[12], 2)-0.5, log(v.c[12], 2)+0.5, length.out = 10)

l.v.c = {}

for (i in 1:10){
  l.v.c[[i]] = v.c
  l.v.c[[i]][9] = v.Mg.T[i]
}
for (i in 11:20){
  l.v.c[[i]] = v.c
  l.v.c[[i]][10] = v.Mn.T[(i-10)]
}
for (i in 21:30){
  l.v.c[[i]] = v.c
  l.v.c[[i]][11] = v.Zn.T[(i-20)]
}
for (i in 31:40){
  l.v.c[[i]] = v.c
  l.v.c[[i]][12] = v.Ca.T[(i-30)]
}

####Solve the perterbed systems####

l.df = {}
for (i in 1:40){
  l.df[[i]] = fn.solving(l.v.c[[i]])
}

df = bind_rows(l.df)

head(df)

P.Mg = ggplot(df %>%
         filter(M == "Mg") %>%
         filter(L == "free") %>%
         filter(Mn.T == 4) %>%
         filter(Zn.T == 0.24) %>%
         filter(Ca.T == 0.0005),
       aes(x = Mg.T, y = Conc)) +
  stat_poly_line() +
  geom_point() +
  stat_regline_equation(aes(label = ..eq.label..)) +
  theme_classic() +
  xlab("Total [Mg2+] (mM)") +
  ylab("Free [Mg2+] (mM)") +
  scale_x_continuous(limits = c(27, 62), breaks = c(30, 40, 50, 60))

P.Mn = ggplot(df %>%
         filter(M == "Mn") %>%
         filter(L == "free") %>%
         filter(Mg.T == 40) %>%
         filter(Zn.T == 0.24) %>%
         filter(Ca.T == 0.0005),
       aes(x = Mn.T, y = Conc)) +
  stat_poly_line() +
  geom_point() +
  stat_regline_equation(aes(label = ..eq.label..)) +
  theme_classic() +
  xlab("Total [Mn2+] (mM)") +
  ylab("Free [Mn2+] (mM)") +
  scale_x_continuous(limits = c(2.75, 6.25), breaks = c(3, 4, 5, 6))

P.Zn = ggplot(df %>%
         filter(M == "Zn") %>%
         filter(L == "free") %>%
         filter(Mg.T == 40) %>%
         filter(Mn.T == 4) %>%
         filter(Ca.T == 0.0005),
       aes(x = Zn.T, y = Conc)) +
  stat_poly_line() +
  geom_point() +
  stat_regline_equation(aes(label = ..eq.label..)) +
  theme_classic() +
  xlab("Total [Zn2+] (mM)") +
  ylab("Free [Zn2+] (mM)") 
  #scale_x_continuous(limits = c(0.00275, 0.00625), breaks = c(0.003, 0.004, 0.005, 0.006))

P.Ca = ggplot(df %>%
         filter(M == "Ca") %>%
         filter(L == "free") %>%
         filter(Mg.T == 40) %>%
         filter(Mn.T == 4) %>%
         filter(Zn.T == 0.24),
       aes(x = Ca.T, y = Conc)) +
  stat_poly_line() +
  geom_point() +
  stat_regline_equation(aes(label = ..eq.label..)) +
  theme_classic() +
  xlab("Total [Ca2+] (mM)") +
  ylab("Free [Ca2+] (mM)")
  #scale_x_continuous(limits = c(0.00275, 0.006), breaks = c(0.003, 0.004, 0.005, 0.006))


P = plot_grid(P.Mg, P.Mn, P.Zn, P.Ca, labels = c("A", "B", "C", "D"), align = "v")

ggsave("Figures_Tables/Figure_3/Figure_3.svg",
       P, width = 5, height = 3.5, scale = 2.0)

print(df %>%
         filter(M == "Ca") %>%
         filter(L == "free") %>%
         filter(Mg.T == 40) %>%
         filter(Mn.T == 4) %>%
         filter(Zn.T == 0.24))

print(v.Ca.T)