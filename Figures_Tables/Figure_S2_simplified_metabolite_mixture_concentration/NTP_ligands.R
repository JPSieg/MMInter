#### Load in libraries ####

library(nleqslv)
library(tidyverse)
library(viridis)
library(cowplot)
library(ggpmisc)
library(ggpubr)
devtools::load_all()

####Write a solving function####

fn.solving = function(v.c = df.conc$Concentration[c(4,9:12)]){
  names(v.c) = df.conc$Metabolites[c(4,9:12)]

  fn = function(z) {
    x = z^2
    y = c()

    #Kd Mg
    y[1] = log10(x[1]*x[5]/(Kd.app[4]*x[6])) #ATP

    #Kd Mn
    y[2]  = log10(x[2]*x[5]/(Kd.app[12]*x[7])) #ATP

    #Kd Zn
    y[3] = log10(x[3]*x[5]/(Kd.app[20]*x[8])) #ATP

    #Kd Ca
    y[4] = log10(x[4]*x[5]/(Kd.app[28]*x[9])) #ATP

    #Total M2+
    y[5] = (x[1] + x[6]  - v.c[2])        #Mg
    y[6] = 10*(x[2] + x[7] - v.c[3])    #Mn
    y[7] = 1000*(x[3] + x[8] - v.c[4])  #Zn
    y[8] = 100000*(x[4] + x[9] - v.c[5])  #Ca


    #Total metabolites
    y[9] = (x[5]  + x[6] + x[7] + x[8] + x[9] - v.c[1]) #ATP


    return(y)
  }

  ####Names of each species in x####
  v.names = c("Mg",                   #1
              "Mn",                   #2
              "Zn",                   #3
              "Ca",                   #4
              "ATP",                  #5
              "ATP_Mg",               #6
              "ATP_Mn",               #7
              "ATP_Zn",               #8
              "ATP_Ca")               #9

  #### Starting variables for the solution ####

  v.start = c(3,                    #1
              4,                   #2
              0.24,                #3
              0.0001,               #4
              0.01,                    #5
              27,                   #6
              0.1,                     #7
              0.1,                      #8
              0.00001)
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

  df$Mg.T = v.c[2]
  df$Mn.T = v.c[3]
  df$Zn.T = v.c[4]
  df$Ca.T = v.c[5]

  C = as.character(df$L)

  C[df$L %in% c("L-Aspartic acid", "L-Alanine", "L-Asparagine", "Pyruvic acid")] = "Amino acid/other\ncarboxylate ligands"
  C[df$L %in% c("Glutathione")] = "Glutathione"
  C[df$L %in% c("Glucose 1-P", "AMP")] = "Mono/Di\nphosphate ligands"
  C[df$L %in% c("ATP")] = "NTP ligands"

  df$C = factor(C, levels = c("free", "Amino acid/other\ncarboxylate ligands", "Glutathione", "Mono/Di\nphosphate ligands", "NTP ligands"))

  return(df)
}

####Test solving function####

df = fn.solving()
df$System = "NTPs only"

write.csv(df, "Figures_Tables/Figure_S2_simplified_metabolite_mixture_concentration/NTPs.csv", row.names = F, quote = F)
