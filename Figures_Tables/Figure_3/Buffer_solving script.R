#### Load in libraries ####

library(nleqslv)
library(tidyverse)
library(viridis)
library(cowplot)
library(ggpmisc)
library(ggpubr)
devtools::load_all()

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
    y[35] = 1000*(x[3] + x[29] + x[30] + x[31] + x[32] + x[33] + x[34] + x[35] + x[36] - v.c[11])  #Zn
    y[36] = 1000*(x[4] + x[37] + x[38] + x[39] + x[40] + x[41] + x[42] + x[43] + x[44] - v.c[12])  #Ca


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

  v.start = c(2,                      #1
              0.2,                    #2
              0.0004,                 #3
              0.000004,               #4
              100,                    #5
              16.6,                   #6
              20,                     #7
              1,                      #8
              4,                      #9
              4,                      #10
              3.8,                    #11
              1,                      #12
              1,                      #13
              10^-4,                  #14
              2,                      #15
              22,                     #16
              2,                      #17
              10^-4,                  #18
              10^-4,                  #19
              1,                      #20
              0.01,                   #21
              0.01,                   #22
              0.1,                    #23
              3.0,                    #24
              0.01,                   #25
              0.01,                   #26
              0.01,                   #27
              0.01,                   #28
              0.0001,                 #29
              0.0001,                 #30
              0.0001,                 #31
              0.0001,                 #32
              0.0001,                 #33
              0.0001,                 #34
              0.0001,                 #35
              0.0001,                 #36
              0.00001,                #37
              0.00001,                #38
              0.00001,                #39
              0.00001,                #40
              0.00001,                #41
              0.00001,                #42
              0.00001,                #43
              0.00001)                #44

  #### Solve the equations ####

  p = nleqslv(sqrt(v.start),
              fn,
              control = list(xtol = 10^-15,
                             ftol = 10^-15,
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
df.test

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
  l.v.c[[i]][12] = v.Zn.T[(i-30)]
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
         filter(Zn.T == 0.004) %>%
         filter(Ca.T == 0.0001),
       aes(x = Mg.T, y = Conc)) +
  facet_wrap(~M, ncol = 1, scales = "free") +
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
         filter(Zn.T == 0.004) %>%
         filter(Ca.T == 0.0001),
       aes(x = Mn.T, y = Conc)) +
  facet_wrap(~M, ncol = 1, scales = "free") +
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
         filter(Ca.T == 0.0001),
       aes(x = Zn.T, y = Conc)) +
  facet_wrap(~M, ncol = 1, scales = "free") +
  stat_poly_line() +
  geom_point() +
  stat_regline_equation(aes(label = ..eq.label..)) +
  theme_classic() +
  xlab("Total [Zn2+] (mM)") +
  ylab("Free [Zn2+] (mM)") +
  scale_x_continuous(limits = c(0.00275, 0.00625), breaks = c(0.003, 0.004, 0.005, 0.006))

P.Ca = ggplot(df %>%
         filter(M == "Ca") %>%
         filter(L == "free") %>%
         filter(Mg.T == 40) %>%
         filter(Mn.T == 4) %>%
         filter(Zn.T == 0.004),
       aes(x = Ca.T, y = Conc)) +
  facet_wrap(~M, ncol = 1, scales = "free") +
  stat_poly_line() +
  geom_point() +
  stat_regline_equation(aes(label = ..eq.label..)) +
  theme_classic() +
  xlab("Total [Ca2+] (mM)") +
  ylab("Free [Ca2+] (mM)")+
  scale_x_continuous(limits = c(0.00275, 0.006), breaks = c(0.003, 0.004, 0.005, 0.006))


P = plot_grid(P.Mg, P.Mn, P.Zn, P.Ca, labels = c("A", "B", "C", "D"), align = "v")

ggsave("Figures_Tables/Figure_3/Figure_3.svg",
       P, width = 5, height = 3.5, scale = 2.0)
