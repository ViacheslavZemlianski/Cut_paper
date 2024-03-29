library(ggplot2)    #basic ggplot
library(ggpubr)     #provides plots and some statistics
library(ggtext)     #allows HTML-tags in text (element_markdown)
library(svglite)    #for saving images in svg format
library(broom)      #simplifies statistical data (tidy())
library(dplyr)      #df rows reordering (arrange)

setwd("/media/sf_Genomik-data/user_data/Zemlianski Viacheslav/Cut project/Manuscript/Fig6")

theme_set(theme_classic() + theme(text = element_text(family = "Arial")))

################################################################################
#Fig6A, 34C mutants rescue

Fig6A_data <- read.csv(file = "Fig6A_data.csv", header = TRUE) %>%              # loading raw data
  mutate(genotype = factor(genotype, levels = c("WT",
                                              "cut2-447",
                                              "cut3-477",
                                              "psm3-304"))) %>%
  arrange(genotype) %>%
  mutate(treatment = factor(treatment, levels = c("YES","YES-N"))) %>%
  mutate(cat_freq = cat_freq*100)                                               # converting into percentage

Fig6A_stats <- compare_means(cat_freq ~ treatment,                              # statistical tests
                             group.by = "genotype",
                             data = Fig6A_data,
                             method = "t.test",
                             alternative = "less",
                             var.equal = T,
                             p.adjust.method = "holm"
                             )

ggplot(data = Fig6A_data, aes(x = genotype, fill = treatment, y = cat_freq)) +  # plot initialization
  ggtitle("34\u00B0C") +                                                        # plot title
  geom_bar(position="dodge", stat="summary", fun = "mean", 
           na.rm = T, colour="black", width = 0.5) +
  scale_x_discrete(labels = c("WT", 
                              expression(italic("cut2-447")),
                              expression(italic("cut3-477")),
                              expression(italic("psm3-304"))
                              )) +
  scale_fill_manual(values = c("white", "gray50"),
                    labels = c("YES", 
                               expression("YES+NH"[4]*"Cl"))) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 45)) +
  ylab("Frequency of catastrophic mitosis (%)") +
  theme(plot.title = element_text(hjust = 0.5),                                 # plot title formatting
        axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9, colour = "black", 
                                   angle = 30, hjust = 1),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = c(0.2, 0.98),
        legend.title = element_blank(),
        legend.text = element_text(size = 9, colour = "black"),
        legend.text.align = 0,
        legend.key.height = unit(0.35, 'cm'),
        legend.key.width = unit(0.35, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(0.2, 0.2, 0.32, 0.2), "cm")                           # plot margins
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(cat_freq~genotype+treatment,                   # adding errorbars
                                 data = Fig6A_data, 
                                 FUN = function(Fig6A_data) 
                                   c(AVG = mean(Fig6A_data, na.rm=T),
                                     SD = sd(Fig6A_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                position = position_dodge(0.5), width=0.3,
                colour="black", linewidth=0.75
  ) +
  geom_signif(xmin = seq(0.875, nrow(Fig6A_stats) - 0.125, 1),                  # adding statistics
              xmax = seq(1.125, nrow(Fig6A_stats) + 0.125, 1),
              y_position = c(3, 36, 38, 15.5),
              annotations = ""
  ) +
  geom_text(data = Fig6A_stats, 
            aes(x = seq(1, nrow(Fig6A_stats)), 
                y = c(3, 36, 38, 15.5), 
                vjust = -0.5, 
                hjust = 0.6,
                label = paste0("p=", signif(p.adj)), 
                fill = NULL), 
            position = position_nudge(0.125, 0),
            colour = symnum(Fig6A_stats$p.adj,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt)

ggsave("Fig6A.png", width = 8, height = 8, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig6A.svg", width=8, height=8, units="cm", dpi=600, fix_text_size=F)



################################################################################
#Fig6B, 32C mutants rescue

Fig6B_data <- read.csv(file = "Fig6B_data.csv", header = TRUE) %>%              # loading raw data
  mutate(genotype = factor(genotype, levels = c("WT",
                                                "cut8-563",
                                                "smc6-x",
                                                "nse3-R254E"))) %>%
  arrange(genotype) %>%
  mutate(treatment = factor(treatment, levels = c("YES","YES-N"))) %>%
  mutate(cat_freq = cat_freq*100)                                               # converting into percentage

Fig6B_stats <- compare_means(cat_freq ~ treatment,                              # statistical tests
                             group.by = "genotype",
                             data = Fig6B_data,
                             method = "t.test",
                             alternative = "less",
                             var.equal = T,
                             p.adjust.method = "holm"
)

ggplot(data = Fig6B_data, aes(x = genotype, fill = treatment, y = cat_freq)) +  # plot initialization
  ggtitle("32\u00B0C") +                                                        # plot title
  geom_bar(position="dodge", stat="summary", fun = "mean", 
           na.rm = T, colour="black", width = 0.5) +
  scale_x_discrete(labels = c("WT", 
                              expression(italic("cut8-563")),
                              expression(italic("smc6-x")),
                              expression(italic("nse3-R254E"))
  )) +
  scale_fill_manual(values = c("white", "gray50"),
                    labels = c("YES", 
                               expression("YES+NH"[4]*"Cl"))) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 45)) +
  ylab("Frequency of catastrophic mitosis (%)") +
  theme(plot.title = element_text(hjust = 0.5),                                 # plot title formatting
        axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9, colour = "black", 
                                   angle = 30, hjust = 1),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "none",
        plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")                           # plot margins
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(cat_freq~genotype+treatment,                   # adding errorbars
                                 data = Fig6B_data, 
                                 FUN = function(Fig6B_data) 
                                   c(AVG = mean(Fig6B_data, na.rm=T),
                                     SD = sd(Fig6B_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                position = position_dodge(0.5), width=0.3,
                colour="black", linewidth=0.75
  ) +
  geom_signif(xmin = seq(0.875, nrow(Fig6B_stats) - 0.125, 1),                  # adding statistics
              xmax = seq(1.125, nrow(Fig6B_stats) + 0.125, 1),
              y_position = c(3, 39, 13, 6),
              annotations = ""
  ) +
  geom_text(data = Fig6B_stats, 
            aes(x = seq(1, nrow(Fig6B_stats)), 
                y = c(3, 39, 13, 6), 
                vjust = -0.5, 
                hjust = 0.6,
                label = paste0("p=", signif(p.adj)), 
                fill = NULL), 
            position = position_nudge(0.125, 0),
            colour = symnum(Fig6B_stats$p.adj,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt)

ggsave("Fig6B.png", width = 8, height = 8, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig6B.svg", width=8, height=8, units="cm", dpi=600, fix_text_size=F)



################################################################################
#Fig6C, 32C mutants rescue

Fig6C_data <- read.csv(file = "Fig6C_data.csv", header = TRUE) %>%              # loading raw data
  mutate(genotype = factor(genotype, levels = c("WT",
                                                "cut1-206",
                                                "cut4-533",
                                                "cut9-665",
                                                "cut15-85"))) %>%
  arrange(genotype) %>%
  mutate(treatment = factor(treatment, levels = c("YES","YES-N"))) %>%
  mutate(cat_freq = cat_freq*100)                                               # converting into percentage

Fig6C_stats <- compare_means(cat_freq ~ treatment,                              # statistical tests
                             group.by = "genotype",
                             data = Fig6C_data,
                             method = "t.test",
                             alternative = "less",
                             var.equal = T,
                             p.adjust.method = "holm"
)

ggplot(data = Fig6C_data, aes(x = genotype, fill = treatment, y = cat_freq)) +  # plot initialization
  ggtitle("28\u00B0C") +                                                        # plot title
  geom_bar(position="dodge", stat="summary", fun = "mean", 
           na.rm = T, colour="black", width = 0.5) +
  scale_x_discrete(labels = c("WT", 
                              expression(italic("cut1-206")),
                              expression(italic("cut4-533")),
                              expression(italic("cut9-665")),
                              expression(italic("cut15-85"))
  )) +
  scale_fill_manual(values = c("white", "gray50"),
                    labels = c("YES", 
                               expression("YES+NH"[4]*"Cl"))) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 70)) +
  ylab("Frequency of catastrophic mitosis (%)") +
  theme_classic() +                                                             # removing background
  theme(plot.title = element_text(hjust = 0.5),                                 # plot title formatting
        axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9, colour = "black", 
                                   angle = 30, hjust = 1),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "none",
        plot.margin=unit(c(0.2, 0.2, 0.45, 0.2), "cm")                          # plot margins
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(cat_freq~genotype+treatment,                   # adding errorbars
                                 data = Fig6C_data, 
                                 FUN = function(Fig6C_data) 
                                   c(AVG = mean(Fig6C_data, na.rm=T),
                                     SD = sd(Fig6C_data, na.rm=T))), 
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                position = position_dodge(0.5), width=0.3,
                colour="black", linewidth=0.75
  ) +
  geom_signif(xmin = seq(0.875, nrow(Fig6C_stats) - 0.125, 1),                  # adding statistics
              xmax = seq(1.125, nrow(Fig6C_stats) + 0.125, 1),
              y_position = c(5, 18, 12, 65, 33),
              annotations = ""
  ) +
  geom_text(data = Fig6C_stats, 
            aes(x = seq(1, nrow(Fig6C_stats)), 
                y = c(5, 18, 12, 65, 33), 
                vjust = -0.5, 
                hjust = 0.6,
                label = paste0("p=", signif(p.adj)), 
                fill = NULL), 
            position = position_nudge(0.125, 0),
            colour = symnum(Fig6C_stats$p.adj,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt)

ggsave("Fig6C.png", width = 9.5, height = 8, units = "cm", dpi = 600)           # saving the plot
ggsave("Fig6C.svg", width=9.5, height=8, units="cm", dpi=600, fix_text_size=F)
