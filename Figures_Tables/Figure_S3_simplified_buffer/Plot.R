list.files("Figures_Tables/Figure_S3_simplified_buffer")
v.files = list.files("Figures_Tables/Figure_S3_simplified_buffer", pattern = ".csv")
v.files = paste("Figures_Tables/Figure_S3_simplified_buffer", v.files, sep = "/")

l.df = lapply(v.files, read.csv)

df = bind_rows(l.df) %>% filter(System != "")

unique(df$System)

df$M = factor(df$M,
              levels = c("Mg", "Mn", "Zn", "Ca"))

df$System = factor(df$System,
                   levels = c("All metabolites", "Aminoacids_carboxylates", "Glutathione", "MonoDiphosphate_ligands", "NTPs only"),
                   labels = c("All metabolites", "Amino acid/other\ncarboxylate ligands", "Glutathione", "Mono/Di\nphosphate ligands", "NTP ligands"))

P = ggplot(df, aes(x = System, y = slope)) +
  facet_wrap(~M, nrow = 1) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("[Free M2+]/[Total M2+] (mM/mM)")

P

list.files()

ggsave("Figures_Tables/Figure_S3_simplified_buffer/Figure_S3.svg", P, width = 8, height = 3, scale = 1.2)
