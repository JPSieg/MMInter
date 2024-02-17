#### Load in libraries ####

library(tidyverse)
devtools::load_all()

#### Define constants ####

v.c = df.conc$Concentration
names(v.c) = df.conc$Metabolites
Kd.app

df.conc

#### Name variables####

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

#### Write system of equations ####

fn = function(z, pH = 7) {
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

#### Check the results ####

p
v.s = p$x^2
names(v.s) = v.names
v.s
p

#### Parse the results into a graphing data frame ####

M = c()
L = c()
Conc = c()

for (i in 1:length(v.s)){
  print(i)
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

#### Determine how precisely M was solved for ####

df.M = df %>% group_by(M) %>% summarise(Sum = sum(Conc))
df.M$Sum[1] = NA
df.M$M = factor(df.M$M,
                levels = c("free", "Mg", "Mn", "Zn"))
df.M$Known = c(NA, df.conc$Concentration[9:12])
df.M$Norm.resid = (df.M$Sum - df.M$Known)/df.M$Known

#### Determine how precisely L was solved for ####

df.L = df %>% group_by(L) %>% summarise(Sum = sum(Conc))
df.L$Sum[1] = NA
df.L$L = factor(df.L$L,
                levels = c("free",
                           "L-Aspartic acid",      #4
                           "Glutathione",          #5
                           "Glucose 1-P",          #6
                           "ATP",                  #7
                           "AMP",                  #8
                           "L-Alanine",            #9
                           "L-Asparagine",         #10
                           "Pyruvic acid"))
df.L$Known = c(NA, df.conc$Concentration[1:8])
df.L$Norm.resid = (df.L$Sum - df.L$Known)/df.L$Known

#### Determine how precicely Kds were solved for ####

Kd.s = c()
free.L = c()
free.M = c()
for (i in 1:nrow(df)){
  Kd.pos = Kd.app[which(grepl(df$M[i], names(Kd.app)))]
  Kd.pos = Kd.pos[which(grepl(df$L[i], names(Kd.pos)))]
  Kd.s[i] = Kd.pos
  df.i = df %>% filter(L == df$L[i]) %>% filter(M == "free")
  free.L[i] = df.i$Conc[1]
  df.i = df %>% filter(M == df$M[i]) %>% filter(L == "free")
  free.M[i] = df.i$Conc[1]
}

df$Kd.s  = Kd.s

df$Kd.obs = free.M*free.L/df$Conc

df$norm.resid.Kd = (df$Kd.obs - df$Kd.s)/df$Kd.s

#### Make the plots ####

P.mid = ggplot(df, aes(x = M, y = L, fill = Conc, label = formatC(Conc, format = "e", digits = 1))) +
  geom_tile() +
  geom_text() +
  theme_classic() +
  scale_fill_viridis(trans = "log10", direction = -1, name = "Conc.\n(mM)",
                     breaks = c(10^-12, 10^-10, 10^-8, 10^-6, 10^-4, 10^-2, 1, 10^2),
                     limits = c(10^-12, 1.1*10^2)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none",
        legend.key.size = unit(5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(3, 'cm')) +
  ylab("Metabolite\nconcentration (mM)") +
  xlab("Divalent metal ion\n concentration (mM)")
P.mid

P.bottom = ggplot(df %>% filter(M != "free") %>% filter(L != "free"),
       aes(x = L, y = norm.resid.Kd, color = M, shape = M)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  theme_classic() +
  scale_y_continuous(breaks = c(-10^-15, 0, 10^-15), limits = c(-10^-15, 10^-15)) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "right") +
  scale_color_manual(values = c("grey", "red", "black", "blue"), name = "") +
  scale_shape(name = "") +
  ylab("Normalized\nresiduals KD")
P.bottom

P.Top = ggplot(df.M, aes(x = M, y = Norm.resid)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  theme_classic() +
  scale_y_continuous(limits = c(-10^-15, 10^-15), breaks = c(-10^-15, 0, 10^-15)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.line.x = element_blank()) +
  ylab("Normalized\nresiduals")
P.Top

p.side = ggplot(df.L, aes(x = L, y = Norm.resid)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  geom_point() +
  geom_point() +
  theme_classic() +
  scale_y_continuous(limits = c(-10^-15, 10^-15), breaks = c(-10^-15, 0, 10^-15)) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none") +
  ylab("Normalized\nresiduals")
p.side

p.blank = ggplot() + theme_classic() + theme(axis.line = element_blank())

P.left = plot_grid(P.Top, P.mid, ncol = 1, align = "v", rel_heights = c(1, 6))
P.right = plot_grid(p.blank, p.side, ncol = 1, align = "h", rel_heights = c(1, 6))

PA = plot_grid(P.left, P.right, nrow = 1, rel_widths = c(4.5, 1), align = "h")
PA = plot_grid(PA, P.bottom, ncol = 1, rel_heights = c(4, 1), labels = c("A", "B"))

C = as.character(df$L)

C[df$L %in% c("L-Aspartic acid", "Glutathione", "L-Alanine", "L-Asparagine", "Pyruvic acid")] = "Amino acid/other\ncarboxylate ligands"
C[df$L %in% c("Glucose 1-P", "AMP")] = "Mono/Di\nphosphate ligands"
C[df$L %in% c("ATP")] = "NTP ligands"

df$C = factor(C, levels = c("free", "Amino acid/other\ncarboxylate ligands", "Mono/Di\nphosphate ligands", "NTP ligands"))

PB = ggplot(df %>% filter(M != "free"),
       aes(x = C, y = Conc, fill = C)) +
  facet_wrap(~M, ncol = 1, scales = "free") +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 25, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Concentration (mM)")

#### Final plot ####

P = plot_grid(P.mid, PB, nrow = 1, labels = c("A", "B"), rel_widths = c(2,1))

ggsave("Figures_Tables/Figure_2/Figure_2.svg",
       P, width = 5, height = 3.5, scale = 2.0)
