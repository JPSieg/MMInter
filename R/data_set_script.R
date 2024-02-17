df.conc = read.csv("data/Concentrations.csv")
usethis::use_data(df.conc, overwrite = T)
df.NISTSRD46 = read.csv("data/SRD_46_SQL_JPS_parsed.csv")
usethis::use_data(df.NISTSRD46, overwrite = T)
df = read.csv("data/Kd_app.csv")
Kd.app = df$Kd.app
names(Kd.app) = df$X
usethis::use_data(Kd.app, overwrite = T)
