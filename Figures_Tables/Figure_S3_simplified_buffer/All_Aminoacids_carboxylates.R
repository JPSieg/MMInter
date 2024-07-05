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
    y[1] = log10(x[1]*x[5]/(Kd.app[1]*x[9]))
    y[2] = log10(x[1]*x[6]/(Kd.app[6]*x[10]))
    y[3] = log10(x[1]*x[7]/(Kd.app[7]*x[11]))
    y[4] = log10(x[1]*x[8]/(Kd.app[8]*x[12]))

    #Kd Mn
    y[5]  = log10(x[2]*x[5]/(Kd.app[9]*x[13]))
    y[6] = log10(x[2]*x[6]/(Kd.app[14]*x[14]))
    y[7] = log10(x[2]*x[7]/(Kd.app[15]*x[15]))
    y[8] = log10(x[2]*x[8]/(Kd.app[16]*x[16]))

    #Kd Zn
    y[9] = log10(x[3]*x[5]/(Kd.app[17]*x[17]))
    y[10] = log10(x[3]*x[6]/(Kd.app[22]*x[18]))
    y[11] = log10(x[3]*x[7]/(Kd.app[23]*x[19]))
    y[12] = log10(x[3]*x[8]/(Kd.app[24]*x[20]))

    #Kd Ca
    y[13] = log10(x[4]*x[5]/(Kd.app[25]*x[21]))
    y[14] = log10(x[4]*x[6]/(Kd.app[30]*x[22]))
    y[15] = log10(x[4]*x[7]/(Kd.app[31]*x[23]))
    y[16] = log10(x[4]*x[8]/(Kd.app[32]*x[24]))

    #Total M2+
    y[17] = (x[1] + x[9] + x[10] + x[11] + x[12] - v.c[9])        #Mg
    y[18] = 10*(x[2] + x[13] + x[14] + x[15] + x[16] - v.c[10])    #Mn
    y[19] = 1000*(x[3]  + x[17] + x[18] + x[19] + x[20] - v.c[11])  #Zn
    y[20] = 100000*(x[4]  + x[21] + x[22] + x[23] + x[24] - v.c[12])  #Ca

    #Total metabolites
    y[21] = (v.c[1]/v.c[1])*(x[5]  + x[9] + x[13] + x[17] + x[21] - v.c[1])
    y[22] = (v.c[1]/v.c[6])*(x[6] + x[10] + x[14] + x[18] + x[22] - v.c[6])
    y[23] = (v.c[1]/v.c[7])*(x[7] + x[11] + x[15] + x[19] + x[23] - v.c[7])
    y[24] = (v.c[1]/v.c[8])*(x[8] + x[12] + x[16] + x[20] + x[24] - v.c[8])

    return(y)
}

  ####Names of each species in x####
  v.names = c("Mg",                   #1
              "Mn",                   #2
              "Zn",                   #3
              "Ca",                   #4
              "L-Aspartic acid",      #5
              "L-Alanine",            #6
              "L-Asparagine",         #7
              "Pyruvic acid",         #8
              "L-Aspartic acid_Mg",   #9
              "L-Alanine_Mg",         #10
              "L-Asparagine_Mg",      #11
              "Pyruvic acid_Mg",      #12
              "L-Aspartic acid_Mn",   #13
              "L-Alanine_Mn",         #14
              "L-Asparagine_Mn",      #15
              "Pyruvic acid_Mn",      #16
              "L-Aspartic acid_Zn",   #17
              "L-Alanine_Zn",         #18
              "L-Asparagine_Zn",      #19
              "Pyruvic acid_Zn",      #20
              "L-Aspartic acid_Ca",   #21
              "L-Alanine_Ca",         #22
              "L-Asparagine_Ca",      #23
              "Pyruvic acid_Ca")      #24

  #### Starting variables for the solution ####

  v.start = c(2.6,                    #1  "Mg"
              0.066,                  #2  "Mn"
              0.00064,                #3  "Zn"
              0.000011,               #4  "Ca"
              100,                    #5  "L-Aspartic acid"
              4,                      #6 "L-Alanine"
              3.8,                    #7 "L-Asparagine"
              3.5,                      #8 "Pyruvic acid"
              0.15,                      #9 "L-Asp_Mg"
              0.0015,                  #10 "L-Alanine_Mg"
              0.00001,                  #11 "L-Asparagine_Mg"
              0.13,                      #12 "Pyruvic acid_Mg"
              0.074,                   #13 "L-Aspartic acid_Mn"
              0.00011,                   #14 "L-Alanine_Mn"
              0.0014,                   #15 "L-Asparagine_Mn"
              0.0032,                   #16 "Pyruvic acid_Mn"
              0.11,                     #17 "L-Aspartic acid_Zn"
              0.000014,                 #18 "L-Alanine_Zn"
              0.000022,                 #19 "L-Asparagine_Zn"
              0.000038,                 #20 "Pyruvic acid_Zn"
              0.00000012,               #21 "L-Aspartic acid_Ca"
              0.000000015,              #22 "L-Alanine_Ca"
              0.000000000042,           #23 "L-Asparagine_Ca"
              0.00000036)               #24 "Pyruvic acid_Ca"
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
System = "Aminoacids_carboxylates"

df = data.frame(M, slope, System)

list.files("Figures_Tables/Figure_S3_simplified_buffer")

write.csv(df, "Figures_Tables/Figure_S3_simplified_buffer/Aminoacids_carboxylates.csv", row.names = F, quote = F)
