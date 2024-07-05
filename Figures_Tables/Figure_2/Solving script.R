#### Load in libraries ####

library(nleqslv)
library(tidyverse)
library(viridis)
library(cowplot)
devtools::load_all()

#### Define constants ####

v.c = df.conc$Concentration
names(v.c) = df.conc$Metabolites

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
            0.000043,                  #14 "Glutathione_Mg"
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
            0.0000026)

#### Solve the equations ####

p = nleqslv(sqrt(v.start),
            omega,
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
        legend.position = "top",
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
  theme() +
  scale_color_manual(values = c("grey", "red", "black", "blue"), name = "") +
  scale_shape(name = "") +
  ylab("Normalized\nresiduals KD")
P.bottom

P.Top = ggplot(df.M, aes(x = M, y = Norm.resid)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  theme_classic() +
  scale_y_continuous(limits = c(-10^-15, 10^-15), breaks = c(-10^-15, 0, 10^-15)) +
  theme() +
  ylab("Normalized\nresiduals")
P.Top

p.side = ggplot(df.L, aes(x = L, y = Norm.resid)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  geom_point() +
  geom_point() +
  theme_classic() +
  #scale_y_continuous(limits = c(-10^-15, 10^-15), breaks = c(-10^-15, 0, 10^-15)) +
  theme(legend.position = "none") +
  ylab("Normalized\nresiduals")
p.side

p.blank = ggplot() + theme_classic() + theme(axis.line = element_blank())

P.left = plot_grid(P.Top, P.mid, ncol = 1, align = "v", rel_heights = c(1, 6))
P.right = plot_grid(p.blank, p.side, ncol = 1, align = "h", rel_heights = c(1, 6))

PA = plot_grid(P.left, P.right, nrow = 1, rel_widths = c(4.5, 1), align = "h")
PA = plot_grid(PA, P.bottom, ncol = 1, rel_heights = c(4, 1), labels = c("A", "B"))

C = as.character(df$L)

C[df$L %in% c("L-Aspartic acid", "Glutathione", "L-Alanine", "L-Asparagine", "Pyruvic acid")] = "Amino acid/other\ncarboxylate ligands"
C[df$L %in% c("Glutathione")] = "Glutathione"
C[df$L %in% c("Glucose 1-P", "AMP")] = "Mono/Di\nphosphate ligands"
C[df$L %in% c("ATP")] = "NTP ligands"

df$C = factor(C, levels = c("free", "Amino acid/other\ncarboxylate ligands", "Glutathione", "Mono/Di\nphosphate ligands", "NTP ligands"))

PB = ggplot(df %>% filter(M != "free"),
       aes(x = C, y = Conc, fill = C)) +
  facet_wrap(~M, ncol = 1, axis.labels = "margins", scales = "free_y") +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 25, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Concentration (mM)")

PB

#### Final plot ####

P = plot_grid(P.mid, PB, nrow = 1, labels = c("A", "B"), rel_widths = c(2,1.2))

ggsave("Figures_Tables/Figure_2/Figure_2.svg",
       P, width = 5.2, height = 3.5, scale = 1.6)

write.csv(df, "Figures_Tables/Worked_example/Solution.csv", row.names = F)
