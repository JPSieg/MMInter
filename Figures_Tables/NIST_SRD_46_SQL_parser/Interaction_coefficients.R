df.is = data.frame("IS" = c(0.05, 0.10, 0.15, 0.2, 0.5, 1.0, 2.0, 3.0),
                   "Correction" = c(0.09, 0.11, 0.12, 0.13, 0.15, 0.14, 0.11, 0.07))

Kd.app = c()

####Aspartate Mg####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Aspartic", df$L))]
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "Mg<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     0)
}

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[1] = 1000/(a.L*(10^mean(log10K.0[1])))

####Glutathione Mg####

Kd.app[2] = 100000

####Glucose 1-P Mg####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Glucose", df$L))]
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "Mg<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = as.numeric(df$constant[1])

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(gsub(")", "", gsub("(", "", df$constant[i], fixed = T), fixed = T)) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[3] = 1000/(a.L*(10^mean(log10K.0[1])))

####ATP Mg####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("ATP", df$L))]
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "Mg<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -4,
                                     -2)
}

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[4] = 1000/(a.L*(10^mean(log10K.0[1])))

####AMP Mg####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("AMP", df$L))]
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "Mg<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     -0)
}

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[5] = 1000/(a.L*(10^mean(log10K.0[1])))

####Alanine Mg####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Alanine", df$L))]
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "Mg<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = as.numeric(df$constant)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[6] = 1000/(a.L*(10^mean(log10K.0[1])))

####Glutamine Mg####

Kd.app[7] = 100000

####Pyruvate Mg####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Pyruvic", df$L))]
df = subset(df, df$L == "2-Oxopropanoic acid (Pyruvic acid)")
df = subset(df, df$M == "Mg<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -1,
                                     1)
}
log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "2-Oxopropanoic acid (Pyruvic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))


#Correct for pH 7.0
pH = 7
a.L = 1/(1 + (10^(pKa.HL - pH)))
a.HL = 1 - a.L
Kd.app[8] = 1000/(a.L*(10^mean(log10K.0[1])))

####Aspartate Mn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Aspartic", df$L))]
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "Mn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")
df$constant = gsub(")", "", gsub("(", "", df$constant, fixed = T), fixed = T)
#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     0)
}

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[9] = 1000/(a.L*(10^mean(log10K.0[1])))

####Glutathione Mn####
#Pull out the relevant data

Kd.app[10] = 500000

####Glucose 1-P Mn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Glucose", df$L))]
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "Mn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = as.numeric(df$constant[1])

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(gsub(")", "", gsub("(", "", df$constant[i], fixed = T), fixed = T)) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[11] = 1000/(a.L*(10^mean(log10K.0[1])))

####ATP Mn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("ATP", df$L))]
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "Mn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -4,
                                     -2)
}

log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[12] = 1000/(a.L*(10^mean(log10K.0[1])))

####AMP Mn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("AMP", df$L))]
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "Mn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     -0)
}

log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[13] = 1000/(a.L*(10^mean(log10K.0[1])))

####Alanine Mn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Alanine", df$L))]
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "Mn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$temperature == "25")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     -0)
}

log10K.0 = mean(log10K.0)

log10K.0 = as.numeric(df$constant)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[14] = 1000/(a.L*(10^mean(log10K.0[1])))

####Glutamine Mn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Glutamine", df$L))]
df = subset(df, df$L == "L-2-Aminopentanedioic acid 5-amide (Glutamine)")
df = subset(df, df$M == "Mn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$temperature == "25")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     -0)
}

log10K.0 = mean(log10K.0)

log10K.0 = as.numeric(df$constant)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopentanedioic acid 5-amide (Glutamine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopentanedioic acid 5-amide (Glutamine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[15] = 1000/(a.L*(10^mean(log10K.0[1])))


####Pyruvate Mn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Pyruvic", df$L))]
df = subset(df, df$L == "2-Oxopropanoic acid (Pyruvic acid)")
df = subset(df, df$M == "Mn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -1,
                                     1)
}
log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "2-Oxopropanoic acid (Pyruvic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))


#Correct for pH 7.0
pH = 7
a.L = 1/(1 + (10^(pKa.HL - pH)))
a.HL = 1 - a.L
Kd.app[16] = 1000/(a.L*(10^mean(log10K.0[1])))

####Aspartate Zn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Aspartic", df$L))]
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "Zn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")
#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     0)
}

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[17] = 1000/(a.L*(10^mean(log10K.0[1])))

####Glutathione Zn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Glutathione", df$L))]
df = subset(df, df$L == "L-5-Glutamyl-L-cysteinylglycine (Glutathione)")
df = subset(df, df$M == "Zn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.15")
#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     0)
}

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-5-Glutamyl-L-cysteinylglycine (Glutathione)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-5-Glutamyl-L-cysteinylglycine (Glutathione)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[18] = 1000/(a.L*(10^mean(log10K.0[1])))

####Glucose 1-P Zn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Glucose", df$L))]
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "Zn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     0)
}

log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(gsub(")", "", gsub("(", "", df$constant[i], fixed = T), fixed = T)) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[19] = 1000/(a.L*(10^mean(log10K.0[1])))

####ATP Zn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("ATP", df$L))]
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "Zn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -4,
                                     -2)
}

log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[20] = 1000/(a.L*(10^mean(log10K.0[1])))

####AMP Zn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("AMP", df$L))]
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "Zn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")
df$constant = gsub(")", "", gsub("(", "", df$constant, fixed = T), fixed = T)

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     -0)
}

log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[21] = 1000/(a.L*(10^mean(log10K.0[1])))

####Alanine Zn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Alanine", df$L))]
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "Mn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$temperature == "25")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     -0)
}

log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[22] = 1000/(a.L*(10^mean(log10K.0[1])))

####Glutamine Zn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Glutamine", df$L))]
df = subset(df, df$L == "L-2-Aminopentanedioic acid 5-amide (Glutamine)")
df = subset(df, df$M == "Mn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$temperature == "25")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     -0)
}

log10K.0 = mean(log10K.0)

log10K.0 = as.numeric(df$constant)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopentanedioic acid 5-amide (Glutamine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopentanedioic acid 5-amide (Glutamine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[23] = 1000/(a.L*(10^mean(log10K.0[1])))

####Pyruvate Zn####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Pyruvic", df$L))]
df = subset(df, df$L == "2-Oxopropanoic acid (Pyruvic acid)")
df = subset(df, df$M == "Zn<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -1,
                                     1)
}
log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "2-Oxopropanoic acid (Pyruvic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))


#Correct for pH 7.0
pH = 7
a.L = 1/(1 + (10^(pKa.HL - pH)))
a.HL = 1 - a.L
Kd.app[24] = 1000/(a.L*(10^mean(log10K.0[1])))

####Aspartate Ca####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Aspartic", df$L))]
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "Ca<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")
#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     0)
}
log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-Aminobutanedioic acid (Aspartic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[25] = 1000/(a.L*(10^mean(log10K.0[1])))

####Glutathione Ca####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Glutathione", df$L))]
df = subset(df, df$L == "L-5-Glutamyl-L-cysteinylglycine (Glutathione)")
df = subset(df, df$M == "Ca<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     0)
}

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-5-Glutamyl-L-cysteinylglycine (Glutathione)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-5-Glutamyl-L-cysteinylglycine (Glutathione)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[26] = 1000/(a.L*(10^mean(log10K.0[1])))

####Glucose 1-P Ca####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Glucose", df$L))]
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "Ca<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     0)
}

log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "alpha-D(-)-Glucose-1-dihydrogenphosphate (Glucose-1-phosphate)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(gsub(")", "", gsub("(", "", df$constant[i], fixed = T), fixed = T)) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[27] = 1000/(a.L*(10^mean(log10K.0[1])))

####ATP Ca####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("ATP", df$L))]
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "Ca<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -4,
                                     -2)
}

log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(tetrahydrogentriphosphate) (ATP)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[28] = 1000/(a.L*(10^mean(log10K.0[1])))

####AMP Ca####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("AMP", df$L))]
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "Ca<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")
df = subset(df, df$ionicstrength == "0.1")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     -0)
}

log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "Adenosine-5'-(dihydrogenphosphate) (AMP-5)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[29] = 1000/(a.L*(10^mean(log10K.0[1])))

####Alanine Ca####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Alanine", df$L))]
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "Ca<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$temperature == "25")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -2,
                                     -0)
}

log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Determine pKa.H2L

df = df.NISTSRD46
df = subset(df, df$L == "L-2-Aminopropanoic acid (Alanine)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[H<sub>2</sub>L]/[HL][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.H2L = as.numeric(df$constant)

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
Kd.app[30] = 1000/(a.L*(10^mean(log10K.0[1])))

####Glutamine Ca####

Kd.app[31] = 500000

####Pyruvate Ca####

#Pull out the relevant data

df = df.NISTSRD46
df$L[which(grepl("Pyruvic", df$L))]
df = subset(df, df$L == "2-Oxopropanoic acid (Pyruvic acid)")
df = subset(df, df$M == "Ca<sup>2+</sup>")
df = subset(df, df$Beta == "[ML]/[M][L]")
df = subset(df, df$Type == "K")

#Calculate the log10K at IS = 0.15

log10K.0 = c()
for (i in 1:nrow(df)){
  log10K.0[i] = log10K.IS.calculator(as.numeric(df$constant[i]),
                                     as.numeric(df$ionicstrength[i]),
                                     0.15,
                                     2,
                                     -1,
                                     1)
}
log10K.0 = mean(log10K.0)

#Determine pKa.HL

df = df.NISTSRD46
df = subset(df, df$L == "2-Oxopropanoic acid (Pyruvic acid)")
df = subset(df, df$M == "H<sup>+</sup>")
df = subset(df, df$Beta == "[HL]/[L][H]")
df = subset(df, df$Type == "K")
df = subset(df, df$temperature == "25")
df = subset(df, df$ionicstrength == "0.1")

for (i in 1:nrow(df)){
  df$constant[i] = as.numeric(df$constant[i]) + df.is$Correction[which.min(abs(df.is$IS - as.numeric(df$ionicstrength[i])))]
}

pKa.HL = mean(as.numeric(df$constant))

#Correct for pH 7.0
pH = 7
a.L = 1/(1 + (10^(pKa.HL - pH)))
a.HL = 1 - a.L
Kd.app[32] = 1000/(a.L*(10^mean(log10K.0[1])))

names(Kd.app) = c(paste(df.conc$Metabolites[1:8], "Mg", sep = "-"),
                  paste(df.conc$Metabolites[1:8], "Mn", sep = "-"),
                  paste(df.conc$Metabolites[1:8], "Zn", sep = "-"),
                  paste(df.conc$Metabolites[1:8], "Ca", sep = "-"))

write.csv(data.frame(Kd.app), "data/Kd_app.csv")
