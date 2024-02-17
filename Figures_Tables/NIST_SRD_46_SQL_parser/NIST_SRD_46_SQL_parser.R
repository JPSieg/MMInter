library(tidyverse)

list.files("data/SRD_46_SQL")

####Read in constant data####

df = read.delim("data/SRD_46_SQL/verkn_ligand_metal.txt", header = F)
colnames(df) = c("verkn_ligand_metalID", "ligandenNr", "metalNr", "beta_definitionNr", "constanttypNr", "temperature", "ionicstrength",
                 "constant", "constant_sic", "error", "footnoteNr", "solventNr", "electrolyte", "comment", "Importhilfe", "Acc_feld")

####Read in ligands data####

df.L = read.delim("data/SRD_46_SQL/liganden.txt", header = F)
colnames(df.L) = c("name_ligand", "formula", "link_pictures", "pictures_description", "figure_definition", "ligand_classNr", "comment",
                 "Importhilfe", "Acc_feld")
head(df.L)

####Read in the metals data####

df.M = read.delim("data/SRD_46_SQL/metal.txt", header = F)
colnames(df.M) = c("metalID", "name_metal", "name_metal_pur", "type", "order", "part_of", "preselection",
                   "comment", "Acc_feld")
head(df.M)

####Read in the Literature data####

df.Lit = read.delim("data/SRD_46_SQL/literature_alt.txt", header = F)
colnames(df.Lit) = c("literature_altID", "literature_alt", "literature_shortcut", "comment", "Acc_feld")
head(df.Lit)
####Read in the Constant type data####

df.Type = read.delim("data/SRD_46_SQL/constanttyp.txt", header = F)
colnames(df.Type) = c("constanttypID", "name_constanttyp", "comment", "Acc_feld")

####Read in equilibrium####

df.beta = read.delim("data/SRD_46_SQL/beta_definition.txt", header = F)
colnames(df.beta) = c("beta_definitionID", "name_beta_definition", "comment", "Acc_feld")

####Extract the information I need for each constant####

L = c()
M = c()
Lit = c()
Type = c()
Beta = c()

for (i in 1:nrow(df)){
  print(i)
  if (length(which(df.L$name_ligand == df$ligandenNr[i])) > 0){
    L[i] = df.L$formula[which(df.L$name_ligand == df$ligandenNr[i])]
  }else{
    L[i] = NA
  }
  M[i] = df.M$name_metal[which(df.M$metalID == df$metalNr[i])]
  Lit[i] = df.Lit$literature_altID[which(df.Lit$literature_altID == df$beta_definitionNr[i])]
  Type[i] = df.Type$name_constanttyp[which(df.Type$constanttypID == df$constanttypNr[i])]
  if (length(which(df.beta$beta_definitionID == df$beta_definitionNr[i])) > 0){
    Beta[i] = df.beta$name_beta_definition[which(df.beta$beta_definitionID == df$beta_definitionNr[i])]
  }else{
    Beta[i] = NA
  }
}

df$L = L
df$M = M
df$Lit = Lit
df$Type = Type
df$Beta = Beta


list.files("data")
write.csv(df %>% select(M, L, Beta, Type, constant, temperature, ionicstrength, ligandenNr, metalNr, Lit),
          "data/SRD_46_SQL_JPS_parsed.csv", row.names = F)


