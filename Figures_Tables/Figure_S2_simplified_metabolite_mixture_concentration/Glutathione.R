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
    y[1] = log10(x[1]*x[6]/(Kd.app[2]*x[6]))

    #Kd Mn
    y[2] = log10(x[2]*x[6]/(Kd.app[10]*x[7]))

    #Kd Zn
    y[3] = log10(x[3]*x[6]/(Kd.app[18]*x[8]))

    #Kd Ca
    y[4] = log10(x[4]*x[6]/(Kd.app[26]*x[9]))

    #Total M2+
    y[5] = (x[1] + x[6] - v.c[9])        #Mg
    y[6] = 10*(x[2] + x[7] - v.c[10])    #Mn
    y[7] = 1000*(x[3] + x[8] - v.c[11])  #Zn
    y[8] = 100000*(x[4] + x[9] - v.c[12])  #Ca


    #Total metabolites
    y[9] = (v.c[1]/v.c[2])*(x[5]  + x[6] + x[7] + x[8] + x[9] - v.c[2])

    return(y)
  }

  ####Names of each species in x####
  v.names = c("Mg",                   #1
              "Mn",                   #2
              "Zn",                   #3
              "Ca",                   #4
              "Glutathione",          #5
              "Glutathione_Mg",       #6
              "Glutathione_Mn",       #7
              "Glutathione_Zn",       #8
              "Glutathione_Ca")       #9

  #### Starting variables for the solution ####

  v.start = c(40,                    #1  "Mg"
              4,                  #2  "Mn"
              0.01,                #3  "Zn"
              0.0001,               #4  "Ca"
              16,                     #5  "Glutathione"
              0.0001,                #6  "Glutathione_Mg"
              0.0001,              #7  "Glutathione_Mn"
              0.23,                  #8  "Glutathione_Zn"
              0.000000001)            #9  "Glutathione_Ca"

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

####Test solving function####

df = fn.solving()
df$System = "Glutathione"

write.csv(df, "Figures_Tables/Figure_S2_simplified_metabolite_mixture_concentration/Glutathione.csv", row.names = F, quote = F)

