v.files = list.files("Figures_Tables/Figure_S2_simplified_metabolite_mixture_concentration", pattern = ".csv")
v.files = paste("Figures_Tables/Figure_S2_simplified_metabolite_mixture_concentration", v.files, sep = "/")

l.df = lapply(v.files, read.csv)

df = bind_rows(l.df) %>% filter(System != "")

unique(df$System)

df$M = factor(df$M,
              levels = c("Mg", "Mn", "Zn", "Ca"))

df$System = factor(df$System,
                   levels = c("All metabolites", "Aminoacids_carboxylates", "Glutathione", "MonoDiphosphate_ligands", "NTPs only"),
                   labels = c("All metabolites", "Amino acid/other\ncarboxylate ligands", "Glutathione", "Mono/Di\nphosphate ligands", "NTP ligands"))

P = ggplot(df %>% filter(L == "free"), aes(x = System, y = Conc)) +
  facet_wrap(~M, scales = "free", nrow = 1) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("[Free M2+] (mM)")

P

list.files()

ggsave("Figures_Tables/Figure_S2_simplified_metabolite_mixture_concentration/Figure_S2.svg", P, width = 8, height = 3, scale = 1.5)
