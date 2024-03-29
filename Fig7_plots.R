library(multcomp)  #multiple comparison statistics
library(ggplot2)   #basic ggplot
library(ggpubr)    #provides plots and some statistics
library(ggtext)    #allows HTML-tags in text (element_markdown)
library(broom)     #simplifies statistical data (tidy())
library(svglite)   #for saving images in svg format

setwd("/media/sf_Genomik-data/user_data/Zemlianski Viacheslav/Cut project/Manuscript/Fig7")

theme_set(theme_classic() + theme(text = element_text(family = "Arial")))

################################################################################
#Fig7A, rapamycin treatment

Fig7A_data <- read.csv(file = "Fig7A_data.csv", header = TRUE) %>%              # loading raw data
  mutate(treatment = factor(treatment,levels = c("DMSO", "Amm",
                                                 "Rapa", "Rapa+Amm"))) %>%
  mutate(cat_freq = cat_freq * 100)                                             # converting into percentage

Fig7A_stats <- aov(cat_freq ~ treatment, Fig7A_data) %>%                        # statistical tests
  glht(mcp(treatment = "Dunnet"), alternative = "less") %>%
  tidy()

ggplot(Fig7A_data, aes(x = treatment, y = cat_freq, fill = treatment)) +        # plot initialization
  geom_bar(stat="summary", fun = "mean",                                        # setting bars
           na.rm = T, colour="black", width = 0.5) +
  scale_x_discrete(labels = c("DMSO",                                           # x-axis labels
                              "NH<sub>4</sub>Cl",
                              "Rapa",
                              "NH<sub>4</sub>Cl+<br>Rapa")) +
  scale_fill_manual(values = rep("grey90", nrow(Fig7A_stats) + 1)) +            # bar colors
  scale_y_continuous(expand=c(0,0), limits=c(0,13.5), breaks=c(0,3,6,9,12)) +   # y-axis settings
  ylab("Frequency of catastrophic mitosis (%)") +
  theme_classic() +                                                             # removing background
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black"),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "none",                                               # removing legend
        plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")                           # margins
  ) +
  geom_point(na.rm=T, show.legend=F, size=2) +                                  # adding datapoints
  geom_errorbar(data = aggregate(cat_freq ~ treatment,                          # adding errorbars
                                 data = Fig7A_data, 
                                 FUN = function(Fig7A_data) 
                                   c(AVG = mean(Fig7A_data, na.rm=T),
                                     SD = sd(Fig7A_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                width=0.3, colour="black", linewidth=0.75
  ) +
  geom_signif(xmin = c(1, seq(2, nrow(Fig7A_stats))),                           # adding statistics
              xmax = c(mean(c(2, nrow(Fig7A_stats)+1)), 
                       seq(3, nrow(Fig7A_stats)+1)),
              annotations = "",
              y_position = c(13, rep(12.75, nrow(Fig7A_stats)-1))
  ) +
  geom_text(data = Fig7A_stats, 
            aes(x = seq(2, nrow(Fig7A_stats)+1), 
                y = 10.7, 
                label = paste0("p=", format(adj.p.value, 
                                            digits=2, 
                                            scientific=T)), 
                fill = NULL), 
            angle = 90, colour = "black", size = 9/.pt)

ggsave("Fig7A.png", width = 8, height = 8, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig7A.svg", width=8, height=8, units="cm", dpi=600, fix_text_size=F)



################################################################################
#Fig7B, tor1 mutants

Fig7B_data <- read.csv("Fig7B_data.csv", header = T) %>%                        # loading raw data
  mutate(sample = paste(genotype, temperature, sep = "_"),
         sample = factor(sample, levels = c("WT_32", "WT_34",
                                            "Dcbf11_32", "Dcbf11_34",
                                            "tor1ts_32", "tor1ts_34",
                                            "tor1KO_32", "tor1KO_34")),
         genotype = factor(genotype, levels = c("WT", "Dcbf11",
                                                "tor1ts", "tor1KO")),
         temperature = factor(temperature, levels = c("32", "34")),
         cat_freq = cat_freq * 100
         )

# Fig7B_t_1 <- t.test(c(Fig7B_data$cat_freq[Fig7B_data$sample == "Dcbf11_34"]),   #statistical tests
#                     y = c(Fig7B_data$cat_freq[Fig7B_data$sample == "tor1ts_34"]),
#                     alternative = "greater",
#                     var.equal = T)
# Fig7B_t_2 <- t.test(c(Fig7B_data$cat_freq[Fig7B_data$sample == "Dcbf11_32"]),
#                     y = c(Fig7B_data$cat_freq[Fig7B_data$sample == "tor1KO_32"]),
#                     alternative = "greater",
#                     var.equal = T)
# Fig7B_stats <- rbind(tidy(Fig7B_t_1), tidy(Fig7B_t_2)) %>%
#   mutate(p.adj = p.adjust(p.value, method = "holm"))
# rm(Fig7B_t_1, Fig7B_t_2)

Fig7B_stats <- aov(cat_freq ~ sample,                                           # statistical tests
                   data = Fig7B_data) %>%
  glht(mcp(sample = c("Dcbf11_34 - tor1ts_34 <= 0", 
                      "Dcbf11_32 - tor1KO_32 <= 0"
  ))) %>%
  tidy()

ggplot(Fig7B_data, aes(x = genotype, fill = temperature, y = cat_freq)) +       # plot initialization
  geom_bar(position = position_dodge(preserve = "single"),                      # adding bars
           stat="summary", fun = "mean", na.rm = F, 
           colour="black", width = 0.5) +
  scale_fill_manual(values = c("white", "gray50"),                              # bar colors
                    labels = c("32\u00B0C", "34\u00B0C")) +                     # legend labels
  scale_x_discrete(labels = c("WT",                                             # x-axis labels
                              "<i>\u0394cbf11</i>",
                              "<i>\u0394cbf11<br>tor1-D</i>",
                              "<i>\u0394cbf11<br>\u0394tor1</i>"
  )) +
  scale_y_continuous(expand=c(0,0), limits=c(0,28.7), breaks=seq(0,24,3)) +     # y-axis setings
  coord_cartesian(ylim = c(0, 24), expand = T, clip = "off") +
  ylab("Frequency of catastrophic mitosis (%)") +
  theme_classic() +                                                             # removing background
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black"),
                                       # angle = 30, hjust = 0.7, vjust = 0.9),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = c(0.15, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 9, colour = "black"),
        legend.text.align = 0,
        legend.key.size = unit(0.3, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(1.05, 0.2, 0.2, 0.2), "cm")                          # plot margins
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(cat_freq~genotype+temperature,                 # adding errorbars
                                 data = Fig7B_data, 
                                 FUN = function(Fig7B_data) 
                                   c(AVG = mean(Fig7B_data, na.rm=T),
                                     SD = sd(Fig7B_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                position = position_dodge(0.5, preserve = "single"), 
                width=0.3, colour="black", linewidth=0.75
  ) +
  geom_signif(xmin = c(2.125, 1.875),                                           # adding statistics brackets
              xmax = c(3.125, 3.875),
              y_position = c(24.8, 26.5),
              annotations = "") +
  geom_text(data = Fig7B_stats, 
            aes(x = c(2.7, 3), 
                y = c(24.8, 26.5), 
                vjust = -0.5, 
                hjust = 0.6,
                label = paste0("p=", format(adj.p.value, scientific=T, digits=2)), 
                fill = NULL),
            colour = symnum(Fig7B_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt) +
  geom_text(aes(x = 4.2, y = 0.8, label = "NA"), size = 9/.pt)

ggsave("Fig7B.png", width = 8, height = 8, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig7B.svg", width=8, height=8, units="cm", dpi=600, fix_text_size=F)



################################################################################
#Fig7C, ssp2KO

Fig7C_data <- read.csv("Fig7C_data.csv", header = T) %>%                        # loading raw data
  mutate(genotype = factor(genotype, levels = c("WT", "Dcbf11", "Dssp2")),
         cat_freq = cat_freq * 100)

Fig7C_stats <- t.test(c(Fig7C_data$cat_freq[Fig7C_data$genotype == "Dcbf11"]),  # statistical tests
                    y = c(Fig7C_data$cat_freq[Fig7C_data$genotype == "Dssp2"]),
                    alternative = "greater",
                    var.equal = T) %>%
  tidy()

ggplot(Fig7C_data, aes(x = genotype, y = cat_freq, fill = genotype)) +          # plot initialization
  geom_bar(stat="summary", fun = "mean",                                        # setting bars
           na.rm = T, colour="black", width = 0.5) + 
  scale_x_discrete(labels = c("WT", 
                              "<i>\u0394cbf11</i>",
                              "<i>\u0394cbf11<br>\u0394ssp2</i>"
  )) +
  scale_fill_manual(values = rep("grey90", 3)) +                                # bar colors
  scale_y_continuous(expand=c(0,0), limits=c(0,13.5), breaks=c(0,3,6,9,12)) +   # y-axis settings
  ylab("Frequency of catastrophic mitosis (%)") +
  theme_classic() +                                                             # removing background
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black", vjust = 1),
                                   # angle = 30, hjust = 1),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "none",                                               # removing legend
        plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")                           # plot margins
  ) +
  geom_point(na.rm=T, show.legend=F, size=2) +                                  # adding datapoints
  geom_errorbar(data = aggregate(cat_freq ~ genotype,                           # adding errorbars
                                 data = Fig7C_data, 
                                 FUN = function(Fig7C_data) 
                                   c(AVG = mean(Fig7C_data, na.rm=T),
                                     SD = sd(Fig7C_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                width=0.3, colour="black", linewidth=0.75
  ) +
  geom_signif(xmin = 2, xmax = 3, y_position = 12.7,                            # adding statistics brackets
              annotations = paste0("p=", format(Fig7C_stats$p.value, 
                                                scientific=T, 
                                                digits=2)),
              vjust = -0.5, textsize = 9/.pt)

ggsave("Fig7C.png", width = 6.5, height = 8, units = "cm", dpi = 600)           # saving the plot
ggsave("Fig7C.svg", width=6.5, height=8, units="cm", dpi=600, fix_text_size=F)



################################################################################
#Fig7X, effect of autophagy suppression
# For response to the review only

Fig7X_data <- read.csv("Fig7X_data.csv", header = T) %>%                        # loading raw data
  mutate(sample = paste(genotype, treatment, sep = "_"),
         sample = factor(sample, levels = c("WT_DMSO", "WT_BAF",
                                            "Dcbf11_DMSO", "Dcbf11_BAF")),
         genotype = factor(genotype, levels = c("WT", "Dcbf11")),
         treatment = factor(treatment, levels = c("DMSO", "BAF")),
         cat_freq = cat_freq * 100
  )

Fig7X_stats <- t.test(c(Fig7X_data$cat_freq[Fig7X_data$sample == "Dcbf11_DMSO"]),  # statistical tests
                      y = c(Fig7X_data$cat_freq[Fig7X_data$sample == "Dcbf11_BAF"]),
                      alternative = "two.sided",
                      var.equal = T) %>%
  tidy()

ggplot(Fig7X_data, aes(x = genotype, fill = treatment, y = cat_freq)) +         # plot initialization
  geom_bar(position = position_dodge(preserve = "single"),                      # adding bars
           stat="summary", fun = "mean", na.rm = F, 
           colour="black", width = 0.5) +
  scale_fill_manual(values = c("white", "gray50"),                              # bar colors
                    labels = c("DMSO", "bafilomycin")) +                        # legend labels
  scale_x_discrete(labels = c("WT",                                             # x-axis labels
                              "<i>\u0394cbf11</i>"
  )) +
  scale_y_continuous(expand=c(0,0), limits=c(0,14), breaks=c(0,3,6,9,12)) +     # y-axis setings
  # coord_cartesian(ylim = c(0, 24), expand = T, clip = "off") +
  ylab("Frequency of catastrophic mitosis (%)") +
  theme_classic() +                                                             # removing background
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black"),
        # angle = 30, hjust = 0.7, vjust = 0.9),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = c(0.3, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 9, colour = "black"),
        legend.text.align = 0,
        legend.key.size = unit(0.3, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(1.05, 0.2, 0.2, 0.2), "cm")                          # plot margins
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(cat_freq~genotype+treatment,                   # adding errorbars
                                 data = Fig7X_data, 
                                 FUN = function(Fig7X_data) 
                                   c(AVG = mean(Fig7X_data, na.rm=T),
                                     SD = sd(Fig7X_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                position = position_dodge(0.5, preserve = "single"), 
                width=0.3, colour="black", linewidth=0.75
  ) +
  geom_signif(xmin = 1.865, xmax = 2.125, y_position = 12.5,                            # adding statistics brackets
              annotations = paste0("p=", format(Fig7X_stats$p.value, 
                                                scientific=F, 
                                                digits=3)),
              vjust = -0.5, textsize = 9/.pt)


ggsave("Fig7X.png", width = 6.5, height = 8, units = "cm", dpi = 600)           # saving the plot
