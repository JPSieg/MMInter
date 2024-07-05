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
    y[1] = log10(x[1]*x[5]/(Kd.app[3]*x[7]))
    y[2] = log10(x[1]*x[6]/(Kd.app[5]*x[8]))

    #Kd Mn
    y[3] = log10(x[2]*x[5]/(Kd.app[11]*x[9]))
    y[4] = log10(x[2]*x[6]/(Kd.app[13]*x[10]))


    #Kd Zn
    y[5] = log10(x[3]*x[5]/(Kd.app[19]*x[11]))
    y[6] = log10(x[3]*x[6]/(Kd.app[21]*x[12]))

    #Kd Ca
    y[7] = log10(x[4]*x[5]/(Kd.app[27]*x[13]))
    y[8] = log10(x[4]*x[6]/(Kd.app[29]*x[14]))

    #Total M2+
    y[9] = (x[1] + x[7] + x[8] - v.c[9])          #Mg
    y[10] = 10*(x[2] + x[9] + x[10] - v.c[10])      #Mn
    y[11] = 1000*(x[3] + x[11] + x[12]  - v.c[11])    #Zn
    y[12] = 100000*(x[4] + x[13] + x[14] - v.c[12])  #Ca


    #Total metabolites
    y[13] = (v.c[1]/v.c[3])*(x[5] + x[7] + x[9]  + x[11] + x[13] - v.c[3])
    y[14] = (v.c[1]/v.c[5])*(x[6] + x[8] + x[10] + x[12] + x[14] - v.c[5])

    return(y)
  }

  ####Names of each species in x####
  v.names = c("Mg",                   #1
              "Mn",                   #2
              "Zn",                   #3
              "Ca",                   #4
              "Glucose 1-P",          #5
              "AMP",                  #6
              "Glucose 1-P_Mg",       #7
              "AMP_Mg",               #8
              "Glucose 1-P_Mn",       #9
              "AMP_Mn",               #10
              "Glucose 1-P_Zn",       #11
              "AMP_Zn",               #12
              "Glucose 1-P_Ca",       #13
              "AMP_Ca")               #14

  #### Starting variables for the solution ####

  v.start = c(2.6,                    #1  "Mg"
              0.066,                   #2  "Mn"
              0.00064,                #3  "Zn"
              0.000011,               #4  "Ca"
              17,                     #5  "Glucose 1-P"
              7.4,                      #6  "AMP"
              12,                      #7 "Glucose 1-P_Mg"
              1.7,                      #8 "AMP_Mg"
              0.16,                    #9 "Glucose 1-P_Mn"
              0.13,                   #10 "AMP_Mn"
              0.0022,                 #11 "Glucose 1-P_Zn"
              0.0021,                 #12 "AMP_Zn"
              0.00013,                #13 "Glucose 1-P_Ca"
              0.0000065)              #14 "AMP_Ca"

  #### Solve the equations ####

  p = nleqslv(sqrt(v.start),
              fn,
              control = list(xtol = 10^-30,
                             ftol = 10^-30,
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

  C = as.character(df$L)

  C[df$L %in% c("L-Aspartic acid", "L-Alanine", "L-Asparagine", "Pyruvic acid")] = "Amino acid/other\ncarboxylate ligands"
  C[df$L %in% c("Glutathione")] = "Glutathione"
  C[df$L %in% c("Glucose 1-P", "AMP")] = "Mono/Di\nphosphate ligands"
  C[df$L %in% c("ATP")] = "NTP ligands"

  df$C = factor(C, levels = c("free", "Amino acid/other\ncarboxylate ligands", "Glutathione", "Mono/Di\nphosphate ligands", "NTP ligands"))

  return(df)
}

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
                filter(Ca.T == 0.0001),
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
                filter(Ca.T == 0.0001),
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
                filter(Ca.T == 0.0001),
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

P.Mg
P.Mn
P.Zn
P.Ca

####Determine slope####

df.Mg = df %>%
  filter(M == "Mg") %>%
  filter(L == "free") %>%
  filter(Mn.T == 4) %>%
  filter(Zn.T == 0.24) %>%
  filter(Ca.T == 0.0001)
df.Mn = df %>%
  filter(M == "Mn") %>%
  filter(L == "free") %>%
  filter(Mg.T == 40) %>%
  filter(Zn.T == 0.24) %>%
  filter(Ca.T == 0.0001)
df.Zn = df %>%
  filter(M == "Zn") %>%
  filter(L == "free") %>%
  filter(Mg.T == 40) %>%
  filter(Mn.T == 4) %>%
  filter(Ca.T == 0.0001)
df.Ca = df %>%
  filter(M == "Ca") %>%
  filter(L == "free") %>%
  filter(Mg.T == 40) %>%
  filter(Mn.T == 4) %>%
  filter(Zn.T == 0.24)

df.Mg

M = c("Mg", "Mn", "Zn", "Ca")
slope = c(lm(Conc ~ Mg.T, df.Mg)[[1]][2],
          lm(Conc ~ Mn.T, df.Mn)[[1]][2],
          lm(Conc ~ Zn.T, df.Zn)[[1]][2],
          lm(Conc ~ Ca.T, df.Ca)[[1]][2])
System = "MonoDiphosphate_ligands"

df = data.frame(M, slope, System)

list.files("Figures_Tables/Figure_S3_simplified_buffer")

write.csv(df, "Figures_Tables/Figure_S3_simplified_buffer/MonoDiphosphate_ligands.csv", row.names = F, quote = F)
