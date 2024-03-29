library(rstatix)   #basic statistics
library(multcomp)  #multiple comparison statistics
library(ggplot2)   #basic ggplot
library(ggpubr)    #provides plots and some statistics
library(ggtext)    #allows HTML-tags in text (element_markdown)
library(broom)     #simplifies statistical data (tidy())
library(svglite)   #for saving images in svg format

setwd("/media/sf_Genomik-data/user_data/Zemlianski Viacheslav/Cut project/Manuscript/Fig1")

theme_set(theme_classic() + theme(text = element_text(family = "Arial")))

################################################################################
# Fig1B, chloride vs sulfate

Fig1B_data <- read.csv(file = "Fig1B_data.csv", header = TRUE)                  # loading raw data
Fig1B_data$genotype <- factor(Fig1B_data$genotype, 
                               levels = c("WT", "Dcbf11"))
Fig1B_data$treatment <- factor(Fig1B_data$treatment, 
                               levels = c("YES","EMM-NH4Cl","EMM-NH4SO4"))
Fig1B_data$cat_freq <- Fig1B_data$cat_freq*100                                  # converting into percentage

Fig1B_stats <- aov(cat_freq ~ treatment,                                        # statistical tests
                 data = Fig1B_data[Fig1B_data$genotype == "Dcbf11",]) %>%
  glht(mcp(treatment = "Dunnet"), alternative = "less") %>%
  tidy() %>%
  add_significance("adj.p.value")                                               # adding stars

ggplot(Fig1B_data, aes(x = treatment, fill = genotype, y = cat_freq)) +
  geom_bar(position="dodge", stat="summary", fun = "mean", 
           na.rm = T, colour="black", width = 0.5) +
  scale_fill_manual(values = c("white", "gray50"),
                    labels = c("WT", expression(italic("\u0394cbf11")))) +      # colors of the bars
  scale_x_discrete(labels = c("YES", 
                              "EMM+<br>NH<sub>4</sub>Cl", 
                              "EMM+<br>NH<sub>4</sub>SO<sub>4</sub>")) +        # labels of the bars
  scale_y_continuous(expand=c(0, 0), limits=c(0,12), breaks=c(0,3,6,9,12)) +
  ylab("Frequency of catastrophic mitosis (%)") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black"),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = c(0.5, 1.07),                                         # legend design
        legend.direction= "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.text.align = 0,
        legend.key.size = unit(0.3, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(0.7, 0.2, 0, 0.7), "cm")
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(cat_freq~genotype+treatment,                   # adding errorbars
                                 data = Fig1B_data, 
                                 FUN = function(Fig1B_data) 
                                   c(AVG = mean(Fig1B_data, na.rm=T),
                                     SD = sd(Fig1B_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                position = position_dodge(0.5), width=0.3,
                colour="black", linewidth=0.75
                ) +
  geom_signif(xmin = c(1, seq(2, nrow(Fig1B_stats))),                           # adding statistics
              xmax = c(mean(c(2, nrow(Fig1B_stats)+1)), 
                       seq(3, nrow(Fig1B_stats)+1)),
              annotations = "",
              position = position_nudge(0.12, 0),
              y_position = c(max(Fig1B_data$cat_freq) + 0.6, 
                             max(Fig1B_data$cat_freq) + 0.21)) +
  geom_text(data = Fig1B_stats, 
            aes(x = seq(2, nrow(Fig1B_stats)+1), 
                y = max(Fig1B_data$cat_freq) - 0.5, 
                label = paste0("p=", format(adj.p.value, digits=2)), 
                fill = NULL), 
            position = position_nudge(0.12, 0),
            colour = symnum(Fig1B_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt)

ggsave("Fig1B.png", width = 8, height = 7, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig1B.svg", width=8, height=7, units="cm", dpi=600, fix_text_size=F)



################################################################################
# Fig1C, ammonium concentrations panel

Fig1C_data <- read.csv(file = "Fig1C_data.csv", header = TRUE)                  # loading raw data
Fig1C_data$genotype <- factor(Fig1C_data$genotype, 
                              levels = c("WT", "Dcbf11"))
Fig1C_data$treatment <- factor(Fig1C_data$treatment, 
                               levels = c("YES","93mM","50mM","20mM"))
Fig1C_data$cat_freq <- Fig1C_data$cat_freq*100                                  # converting into percentage

Fig1C_stats <- aov(cat_freq ~ treatment,                                        # statistical tests
                   data = Fig1C_data[Fig1C_data$genotype == "Dcbf11",]) %>%
  glht(mcp(treatment = "Dunnet"), alternative = "less") %>%
  tidy() %>%
  add_significance("adj.p.value")                                               # adding stars

ggplot(Fig1C_data, aes(x = treatment, fill = genotype, y = cat_freq)) +
  geom_bar(position="dodge", stat="summary", fun = "mean", 
           na.rm = T, colour="black", width = 0.5) +
  scale_fill_manual(values = c("white", "gray50"),
                    labels = c("WT", expression(italic("\u0394cbf11")))) +      # colors of the bars
  scale_x_discrete(labels = c("YES", 
                              "EMM+<br>NH<sub>4</sub>Cl 93mM",                  # labels of the bars
                              "EMM+<br>NH<sub>4</sub>Cl 50mM",
                              "EMM+<br>NH<sub>4</sub>Cl 20mM")) +
  scale_y_continuous(expand=c(0,0), limits=c(0,12), breaks=c(0,3,6,9,12)) +
  ylab("Frequency of catastrophic mitosis (%)") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black", 
                                       angle = 30, hjust = 0.5, vjust = 0.9),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = c(0.5, 1.07),                                         # legend design
        legend.direction= "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 9, colour = "black"),
        legend.text.align = 0,
        legend.key.size = unit(0.3, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(0.7, 0.2, 0, 0.2), "cm")
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(cat_freq~genotype+treatment,                   # adding errorbars
                                 data = Fig1C_data, 
                                 FUN = function(Fig1C_data) 
                                   c(AVG = mean(Fig1C_data, na.rm=T),
                                     SD = sd(Fig1C_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                position = position_dodge(0.5), width=0.3,
                colour="black", linewidth=0.75
  ) +
  geom_signif(xmin = c(1, seq(2, nrow(Fig1C_stats))),                           # adding statistics
              xmax = c(mean(c(2, nrow(Fig1C_stats)+1)), 
                       seq(3, nrow(Fig1C_stats)+1)),
              annotations = "",
              position = position_nudge(0.12, 0),
              y_position = c(11.9, rep(11.6, nrow(Fig1C_stats)-1))) +
  geom_text(data = Fig1C_stats, 
            aes(x = seq(2, nrow(Fig1C_stats)+1), 
                y = 9.3, 
                label = paste0("p=", format(adj.p.value, 
                                            digits=2, 
                                            scientific=T)), 
                fill = NULL), 
            position = position_nudge(0.12, 0),
            angle = 90,
            colour = symnum(Fig1C_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt)

ggsave("Fig1C.png", width = 8, height = 8, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig1C.svg", width=8, height=8, units="cm", dpi=600, fix_text_size=F)



################################################################################
# Fig1D, glutamate concentrations panel

Fig1D_data <- read.csv(file = "Fig1D_data.csv", header = TRUE)                  # loading raw data
Fig1D_data$genotype <- factor(Fig1D_data$genotype, 
                              levels = c("WT", "Dcbf11"))
Fig1D_data$treatment <- factor(Fig1D_data$treatment, 
                               levels = c("YES", 
                                          "Ammonium", 
                                          "Glu_20mM", 
                                          "Glu_100mM"))
Fig1D_data$cat_freq <- Fig1D_data$cat_freq*100                                  # converting into percentage

Fig1D_stats <- aov(cat_freq ~ treatment,                                        # statistical tests
                   data = Fig1D_data[Fig1D_data$genotype == "Dcbf11",]) %>%
  glht(mcp(treatment = "Dunnet"), alternative = "less") %>%
  tidy() %>%
  add_significance("adj.p.value")                                               # adding stars

ggplot(Fig1D_data, aes(x = treatment, fill = genotype, y = cat_freq)) +
  geom_bar(position="dodge", stat="summary", fun = "mean", 
           na.rm = T, colour="black", width = 0.5) +
  scale_fill_manual(values = c("white", "gray50"),
                    labels = c("WT", expression(italic("\u0394cbf11")))) +      # colors of the bars
  scale_x_discrete(labels = c("YES",                                            # labels of the bars
                              "EMM+<br>NH<sub>4</sub>Cl 93mM", 
                              "EMM+<br>Glu 20mM",
                              "EMM+<br>Glu 100mM")) +
  scale_y_continuous(expand=c(0,0), limits=c(0,12), breaks=c(0,3,6,9,12)) +
  #coord_cartesian(ylim = c(0, 12), expand = T, clip = "off") +
  ylab("Frequency of catastrophic mitosis (%)") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black", 
                                       angle = 30, hjust = 0.5, vjust = 0.9),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = c(0.5, 1.07),                                         # legend design
        legend.direction= "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 9, colour = "black"),
        legend.text.align = 0,
        legend.key.size = unit(0.3, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(0.7, 0.2, 0, 0.2), "cm")
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(cat_freq~genotype+treatment,                   # adding errorbars
                                 data = Fig1D_data, 
                                 FUN = function(Fig1D_data) 
                                   c(AVG = mean(Fig1D_data, na.rm=T),
                                     SD = sd(Fig1D_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                position = position_dodge(0.5), width=0.3,
                colour="black", linewidth=0.75
  ) + 
  
  geom_signif(xmin = c(1, seq(2, nrow(Fig1D_stats))),                           # adding statistics
              xmax = c(mean(c(2, nrow(Fig1D_stats)+1)), 
                       seq(3, nrow(Fig1D_stats)+1)),
              annotations = "",
              position = position_nudge(0.12, 0),
              y_position = c(11.9, rep(11.6, nrow(Fig1D_stats)-1))) +
  geom_text(data = Fig1D_stats, 
            aes(x = seq(2, nrow(Fig1D_stats)+1), 
                y = 9.6, 
                label = paste0("p=", format(adj.p.value, 
                                            digits=2, 
                                            scientific=T)), 
                fill = NULL), 
            position = position_nudge(0.12, 0),
            angle = 90,
            colour = symnum(Fig1D_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt)

ggsave("Fig1D.png", width = 8, height = 8, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig1D.svg", width=8, height=8, units="cm", dpi=600, fix_text_size=F)



################################################################################
#Fig1E, poor nitrogen sources panel

Fig1E_data <- read.csv(file = "Fig1E_data.csv", header = TRUE)                  # loading raw data
Fig1E_data$genotype <- factor(Fig1E_data$genotype, 
                              levels = c("WT", "Dcbf11"))
Fig1E_data$treatment <- factor(Fig1E_data$treatment, 
                               levels = c("YES","Amm","Glu","Pro", "Ura"))
Fig1E_data$cat_freq <- Fig1E_data$cat_freq*100                                  # converting into percentage

Fig1E_stats <- aov(cat_freq ~ treatment,                                        # statistical tests
                   data = Fig1E_data[Fig1E_data$genotype == "Dcbf11",]) %>%
  glht(mcp(treatment = "Dunnet"), alternative = "less") %>%
  tidy() %>%
  add_significance("adj.p.value")

ggplot(Fig1E_data, aes(x = treatment, fill = genotype, y = cat_freq)) +
  geom_bar(position="dodge", stat="summary", fun = "mean", 
           na.rm = T, colour="black", width = 0.5) +
  scale_fill_manual(values = c("white", "gray50"),
                    labels = c("WT", expression(italic("\u0394cbf11")))) +      # colors of the bars
  scale_x_discrete(labels = c("YES",                                            # labels of the bars
                              "YES+<br>NH<sub>4</sub>Cl 50mM",
                              "YES+<br>Glu 50mM",
                              "YES+<br>Pro 50mM",
                              "YES+<br>Ura 25mM")) +
  scale_y_continuous(expand=c(0,0), limits=c(0,12), breaks=c(0,3,6,9,12)) +
  ylab("Frequency of catastrophic mitosis (%)") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black", 
                                       angle = 30, hjust = 0.5, vjust = 0.9),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = c(0.5, 1.07),                                         # legend design
        legend.direction= "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 9, colour = "black"),
        legend.text.align = 0,
        legend.key.size = unit(0.3, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(0.7, 0.2, 0, 0.2), "cm")
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(cat_freq~genotype+treatment,                   # adding errorbars
                                 data = Fig1E_data, 
                                 FUN = function(Fig1E_data) 
                                   c(AVG = mean(Fig1E_data, na.rm=T),
                                     SD = sd(Fig1E_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                position = position_dodge(0.5), width=0.3,
                colour="black", linewidth=0.75
  ) + 
  
  geom_signif(xmin = c(1, seq(2, nrow(Fig1E_stats))),                           # adding statistics
              xmax = c(mean(c(2, nrow(Fig1E_stats)+1)), 
                       seq(3, nrow(Fig1E_stats)+1)),
              annotations = "",
              position = position_nudge(0.12, 0),
              y_position = c(11.9, rep(11.6, nrow(Fig1E_stats)-1))) +
  geom_text(data = Fig1E_stats, 
            aes(x = seq(2, nrow(Fig1E_stats)+1), 
                y = 9.8, 
                label = paste0("p=", format(adj.p.value, 
                                            digits=2, 
                                            scientific=F)), 
                fill = NULL), 
            position = position_nudge(0.12, 0),
            angle = 90,
            colour = symnum(Fig1E_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 3)

ggsave("Fig1E.png", width = 10, height = 8, units = "cm", dpi = 600)            # saving the plot
ggsave("Fig1E.svg", width=10, height=8, units="cm", dpi=600, fix_text_size=F)



################################################################################
# Fig1F, cerulenin+ammonium treatment

Fig1F_data <- read.csv(file = "Fig1F_data.csv", header = TRUE)                  # loading raw data
Fig1F_data$treatment <- factor(Fig1F_data$treatment, 
                               levels = c("Cerulenin","Cerulenin+Ammonium"))
Fig1F_data$cat_freq <- Fig1F_data$cat_freq*100                                  # converting into percentage

ggplot(Fig1F_data, aes(x = treatment, y = cat_freq)) +
  geom_bar(stat="summary", fun = "mean", na.rm = T, 
           colour="black", width = 0.3, fill = "white") +
  scale_x_discrete(labels = c("Cerulenin",
                              "Cerulenin+<br>NH<sub>4</sub>Cl")) +
  #scale_fill_manual(values = c("grey90")) +
  scale_y_continuous(expand=c(0,0), limits=c(0,12), breaks=c(0,3,6,9,12)) +
  ylab("Frequency of catastrophic mitosis (%)") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black", 
                                       angle = 0, hjust = 0.5, vjust = 0.7),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "none",                                               # legend design
        plot.margin=unit(c(0.7, 0.2, 0, 0.2), "cm")
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(cat_freq~treatment,                            # adding errorbars
                                 data = Fig1F_data,
                                 FUN = function(Fig1F_data)
                                   c(AVG = mean(Fig1F_data, na.rm=T),
                                     SD = sd(Fig1F_data, na.rm=T))),
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                width=0.1,
                colour="black", linewidth=0.75
                ) +
  stat_compare_means(method = "t.test", 
                     method.args = list(alternative = "less",
                                        var.equal = T),
                     label.x = 2, label.y = 11.4,
                     aes(label = paste0("p=", as.numeric(..p.format..))),
                     hjust = 0.5,
                     size = 9/.pt)

ggsave("Fig1F.png", width = 5, height = 7, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig1F.svg", width=5, height=7, units="cm", dpi=600, fix_text_size=F)



################################################################################
# Fig1X, running Fig1E data against combined YES controls from all experiments
# For response to the review only

Fig1X_data <- read.csv(file = "Fig1X_data.csv", header = TRUE)                  # loading raw data
Fig1X_data$genotype <- factor(Fig1X_data$genotype, 
                              levels = c("WT", "Dcbf11"))
Fig1X_data$treatment <- factor(Fig1X_data$treatment, 
                               levels = c("YES","Amm","Glu","Pro", "Ura"))
Fig1X_data$cat_freq <- Fig1X_data$cat_freq*100                                  # converting into percentage

Fig1X_stats <- aov(cat_freq ~ treatment,                                        # statistical tests
                   data = Fig1X_data[Fig1X_data$genotype == "Dcbf11",]) %>%
  glht(mcp(treatment = "Dunnet"), alternative = "less") %>%
  tidy()

ggplot(Fig1X_data, aes(x = treatment, fill = genotype, y = cat_freq)) +
  geom_bar(position="dodge", stat="summary", fun = "mean", 
           na.rm = T, colour="black", width = 0.5) +
  scale_fill_manual(values = c("white", "gray50"),
                    labels = c("WT", expression(italic("\u0394cbf11")))) +      # colors of the bars
  scale_x_discrete(labels = c("YES",                                            # labels of the bars
                              "YES+<br>NH<sub>4</sub>Cl 50mM",
                              "YES+<br>Glu 50mM",
                              "YES+<br>Pro 50mM",
                              "YES+<br>Ura 25mM")) +
  scale_y_continuous(expand=c(0,0), limits=c(0,12), breaks=c(0,3,6,9,12)) +
  ylab("Frequency of catastrophic mitosis (%)") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black", 
                                       angle = 30, hjust = 0.5, vjust = 0.9),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = c(0.5, 1.07),                                         # legend design
        legend.direction= "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 9, colour = "black"),
        legend.text.align = 0,
        legend.key.size = unit(0.3, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(0.7, 0.2, 0, 0.2), "cm")
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(cat_freq~genotype+treatment,                   # adding errorbars
                                 data = Fig1X_data, 
                                 FUN = function(Fig1X_data) 
                                   c(AVG = mean(Fig1X_data, na.rm=T),
                                     SD = sd(Fig1X_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                position = position_dodge(0.5), width=0.3,
                colour="black", linewidth=0.75
  ) + 
  
  geom_signif(xmin = c(1, seq(2, nrow(Fig1X_stats))),                           # adding statistics
              xmax = c(mean(c(2, nrow(Fig1X_stats)+1)), 
                       seq(3, nrow(Fig1X_stats)+1)),
              annotations = "",
              position = position_nudge(0.12, 0),
              y_position = c(11.9, rep(11.6, nrow(Fig1X_stats)-1))) +
  geom_text(data = Fig1X_stats, 
            aes(x = seq(2, nrow(Fig1X_stats)+1), 
                y = 9.8, 
                label = paste0("p=", format(adj.p.value, 
                                            digits=2, 
                                            scientific=F)), 
                fill = NULL), 
            position = position_nudge(0.12, 0),
            angle = 90,
            colour = symnum(Fig1X_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 3)

ggsave("Fig1X.png", width = 10, height = 8, units = "cm", dpi = 600)            # saving the plot
