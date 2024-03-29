library(multcomp)  #multiple comparison statistics
library(ggplot2)   #basic ggplot
library(ggpubr)    #provides plots and some statistics
library(ggtext)    #allows HTML-tags in text (element_markdown)
library(broom)     #simplifies statistical data (tidy())
library(svglite)   #for saving images in svg format

setwd("/media/sf_Genomik-data/user_data/Zemlianski Viacheslav/Cut project/Manuscript/Fig4")

theme_set(theme_classic() + theme(text = element_text(family = "Arial")))

################################################################################
#Fig4A, total FA

Fig4A_data <- read.csv(file = "Fig4A_data.csv", header = TRUE) %>%              #loading raw data
  mutate(sample = paste(genotype, treatment, sep = "_"),
         sample = factor(sample, levels = c("WT_YES", "WT_YES_N",
                                            "Dcbf11_YES", "Dcbf11_YES_N")),
         genotype = factor(genotype, levels = c("WT", "Dcbf11")),
         treatment = factor(treatment, levels = c("YES", "YES_N"))
         )

# Fig4A_t_1 <- t.test(x = c(Fig4A_data$amount[Fig4A_data$genotype == "WT" &
#                                             Fig4A_data$treatment == "YES"]),
#                     y = c(Fig4A_data$amount[Fig4A_data$genotype == "Dcbf11" &
#                                             Fig4A_data$treatment == "YES"]),
#                     alternative = "greater", var.equal = T)
# Fig4A_t_2 <- t.test(x = c(Fig4A_data$amount[Fig4A_data$genotype == "Dcbf11" &
#                                               Fig4A_data$treatment == "YES-N"]),
#                     y = c(Fig4A_data$amount[Fig4A_data$genotype == "Dcbf11" &
#                                               Fig4A_data$treatment == "YES"]),
#                     alternative = "greater", var.equal = T)
# Fig4A_stats <- rbind(tidy(Fig4A_t_1), tidy(Fig4A_t_2)) %>%
#   mutate(p.adj = p.adjust(p.value, method = "holm"))
# rm(Fig4A_t_1, Fig4A_t_2)

Fig4A_stats <- aov(amount ~ sample,                                             # statistical test
                   data = Fig4A_data) %>%
  glht(mcp(sample = c("WT_YES - Dcbf11_YES <= 0", 
                      "Dcbf11_YES_N - Dcbf11_YES <= 0"
  ))) %>%
  tidy()

ggplot(Fig4A_data, aes(x = genotype, y = amount, fill = treatment)) +           # plot initialization
  geom_bar(position = "dodge", stat="summary", fun = "mean",                    # adding bars
           na.rm = T, colour="black", width = 0.5) +
  scale_fill_manual(values = c("white", "gray50"),                              # bars colors and legend labels
                    labels = c("YES", 
                               expression("YES+NH"[4]*"Cl"))) +
  scale_x_discrete(labels = c("WT", expression(italic("\u0394cbf11")))) +       # x-labels          
  scale_y_continuous(expand=c(0, 0), limits=c(0, 40), breaks = seq(0, 35, 5)) + # y-axis settings
  ylab("Total FA (\u03BCg/mg DCW)") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "top",                                                # legend design
        legend.title = element_blank(),
        legend.text = element_text(size = 9, colour = "black"),
        legend.text.align = 0,
        legend.key.height = unit(0.35, 'cm'),
        legend.key.width = unit(0.35, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")                           # plot margins
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(amount~genotype+treatment,                     # adding errorbars
                                 data = Fig4A_data, 
                                 FUN = function(Fig4A_data) 
                                   c(AVG = mean(Fig4A_data, na.rm=T),
                                     SD = sd(Fig4A_data, na.rm=T))), 
                aes(y = amount[,"AVG"],
                    ymin = pmax(amount[,"AVG"]-amount[,"SD"], 0),
                    ymax = amount[,"AVG"]+amount[,"SD"]),
                position = position_dodge(0.5), width=0.3,
                colour="black", linewidth=0.75
  ) +
  geom_signif(xmin = c(1 - 0.125, 2 - 0.125),                                   # adding statistics
              xmax = c(2 - 0.125, 2 + 0.125),
              y_position = c(37, 30),
              tip_length = 0.07,
              annotations = ""
  ) +
  geom_text(data = Fig4A_stats, 
            aes(x = c(1.5 - 0.125, 2), 
                y = c(37, 30), 
                vjust = -0.5,
                label = paste0("p=", format(adj.p.value, scientific=T, digits=3)), 
                fill = NULL),
            colour = symnum(Fig4A_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt)

ggsave("Fig4A.png", width = 5, height = 7, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig4A.svg", width=5, height=7, units="cm", dpi=600, fix_text_size=F)



################################################################################
#Fig4B, total FA

Fig4B_data <- read.csv(file = "Fig4B_data.csv", header = TRUE) %>%              # loading raw data
  mutate(sample = paste(genotype, treatment, sep = "_"),
         sample = factor(sample, levels = c("WT_YES", "WT_YES_N",
                                            "Dcbf11_YES", "Dcbf11_YES_N")),
         genotype = factor(genotype, levels = c("WT", "Dcbf11")),
         treatment = factor(treatment, levels = c("YES", "YES_N"))
  )

# Fig4B_t_1 <- t.test(x = c(Fig4B_data$amount[Fig4B_data$genotype == "WT" &
#                                               Fig4B_data$treatment == "YES"]),
#                     y = c(Fig4B_data$amount[Fig4B_data$genotype == "Dcbf11" &
#                                               Fig4B_data$treatment == "YES"]),
#                     alternative = "two.sided", var.equal = T)
# Fig4B_t_2 <- t.test(x = c(Fig4B_data$amount[Fig4B_data$genotype == "Dcbf11" &
#                                               Fig4B_data$treatment == "YES-N"]),
#                     y = c(Fig4B_data$amount[Fig4B_data$genotype == "Dcbf11" &
#                                               Fig4B_data$treatment == "YES"]),
#                     alternative = "two.sided", var.equal = T)
# Fig4B_stats <- rbind(tidy(Fig4B_t_1), tidy(Fig4B_t_2)) %>%
#   mutate(p.adj = p.adjust(p.value, method = "holm"))
# rm(Fig4B_t_1, Fig4B_t_2)

Fig4B_stats <- aov(amount ~ sample,                                             # statistical test
                   data = Fig4B_data) %>%
  glht(mcp(sample = c("WT_YES - Dcbf11_YES >= 0", 
                      "Dcbf11_YES_N - Dcbf11_YES >= 0"
  ))) %>%
  tidy()

ggplot(Fig4B_data, aes(x = genotype, y = amount, fill = treatment)) +           # plot initialization
  geom_bar(position = "dodge", stat="summary", fun = "mean",                    # adding bars
           na.rm = T, colour="black", width = 0.5) +
  scale_fill_manual(values = c("white", "gray50"),                              # bars colors and legend labels
                    labels = c("YES", 
                               expression("YES+NH"[4]*"Cl"))) +
  scale_x_discrete(labels = c("WT", expression(italic("\u0394cbf11")))) +       # x-labels          
  scale_y_continuous(expand=c(0, 0), limits=c(0, 55)) +                         # y-axis settings
  ylab("Saturated FA (% of total)") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "none",                                               # legend design
        plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")                           # plot margins
  ) +
  geom_point(position=position_dodge(0.5), na.rm=T, show.legend=F, size=2) +    # adding individual points
  geom_errorbar(data = aggregate(amount~genotype+treatment,                     # adding errorbars
                                 data = Fig4B_data, 
                                 FUN = function(Fig4B_data) 
                                   c(AVG = mean(Fig4B_data, na.rm=T),
                                     SD = sd(Fig4B_data, na.rm=T))), 
                aes(y = amount[,"AVG"],
                    ymin = pmax(amount[,"AVG"]-amount[,"SD"], 0),
                    ymax = amount[,"AVG"]+amount[,"SD"]),
                position = position_dodge(0.5), width=0.3,
                colour="black", linewidth=0.75
  ) +
  geom_signif(xmin = c(1 - 0.125, 2 - 0.125),                                   # adding statistics
              xmax = c(2 - 0.125, 2 + 0.125),
              y_position = c(43, 49),
              tip_length = 0.04,
              annotations = ""
  ) +
  geom_text(data = Fig4B_stats, 
            aes(x = c(1.5 - 0.125, 2), 
                y = c(43, 49), 
                vjust = -0.5,
                label = paste0("p=", format(adj.p.value, scientific=T, digits=3)), 
                fill = NULL),
            colour = symnum(Fig4B_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt)

ggsave("Fig4B.png", width = 5, height = 6, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig4B.svg", width=5, height=6, units="cm", dpi=600, fix_text_size=F)



################################################################################
# Fig4C, individual FA content

Fig4C_data <- read.csv(file = "Fig4C_data.csv", header = TRUE) %>%              # loading raw data
  mutate(sample = paste(genotype, treatment, sep = "_"),
         sample = factor(sample, levels = c("WT_YES", "WT_YES-N",
                                            "Dcbf11_YES", "Dcbf11_YES-N")),
         genotype = factor(genotype, levels = c("WT", "Dcbf11")),
         treatment = factor(treatment, levels = c("YES", "YES-N")),
         FAx = as.integer(as.factor(FA))
  )

Fig4C_stats <- compare_means(amount ~ genotype,                                 # statistical tests
                             group.by = "FA",
                             data = Fig4C_data[Fig4C_data$treatment == "YES",],
                             method = "t.test",
                             alternative = "two.sided",
                             var.equal = T,
                             p.adjust.method = "holm"
)

Fig4C_plot_base <- ggplot(data = Fig4C_data,                                    # basic plot settings for Fig4C
                          aes(x = FAx, fill = sample, y = amount)) +            # plot initialization
  geom_bar(position="dodge", stat="summary", fun = "mean", 
           na.rm = T, colour="black", width = 0.8) +
  scale_fill_manual(
    values = c("red", "pink", "blue", "lightblue"),                             # bar colors
    labels = c("WT <span style='color:#ffffff'>_____</span> YES",               # legend labels
    "WT <span style='color:#ffffff'>_____</span> YES+NH<sub>4</sub>Cl",
    "<i>\u0394cbf11</i> <span stile='color:#ffffff'></span> YES",
    "<i>\u0394cbf11</i> <span stile='color:#ffffff'></span> YES+NH<sub>4</sub>Cl"
  )) +
  ylab("FA content (% of total)") +                                             # y-axis label
  scale_x_continuous(breaks = seq(1, 10, 0.5), expand=c(0.02, 0.02),
                     labels = c("C12:0", "", "C14:0", "", "C16:0", "",
                                "C16:1", "", "C18:0", "", "C18:1", "",
                                "C20:0", "", "C22:0", "", "C24:0", "",
                                "C26:0")) +                                     # removing background
  theme(legend.text = element_markdown(size = 9),
        axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_line(colour = c(rep(c(NA, "black"), 9), NA), 
                                    linewidth = 0.75),
        # axis.ticks.length.x = unit(0.4, "cm"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        plot.margin=unit(c(0, 0.2, 0.5, 0.2), "cm")                             # plot margins
  ) +
  geom_point(position=position_dodge(0.8), na.rm=T,                             # adding individual points
             show.legend=F, size=0.5, alpha = 0.3) +
  geom_errorbar(data = aggregate(amount~FAx+sample,                             # adding errorbars
                                 data = Fig4C_data,
                                 FUN = function(Fig4C_data)
                                   c(AVG = mean(Fig4C_data, na.rm=T),
                                     SD = sd(Fig4C_data, na.rm=T))),
                aes(y = amount[,"AVG"],
                    ymin = pmax(amount[,"AVG"]-amount[,"SD"], 0),
                    ymax = amount[,"AVG"]+amount[,"SD"]),
                position = position_dodge(0.8), width=0.3,
                colour="black", linewidth=0.5)

Fig4C_subset1 = c("C16:0", "C18:0", "C18:1")                                    # 1st subset of FA

Fig4C_plot_1 <- Fig4C_plot_base +                                               # top plot subpanel
  scale_y_continuous(expand=c(0, 0), limits=c(-20, 87)) +                       # y-axis scale
  coord_cartesian(ylim = c(0, 87), expand = T, clip = "on") +
  labs(tag = "\u2193ZOOM\u2193", textsize = 9) +                                # adding zoom pointer
  theme(legend.position = c(0.8, 0.5),
        legend.title = element_blank(),
        legend.key.height = unit(0.35, 'cm'),
        legend.key.width = unit(0.35, 'cm'),
        plot.tag.position = c(0.53, -0.04),                                     # zoom pointer setings
        plot.tag = element_text(size = 9)
        ) +
  geom_signif(xmin = which(Fig4C_stats$FA %in% Fig4C_subset1) - 0.3,            # adding statistics
              xmax = which(Fig4C_stats$FA %in% Fig4C_subset1) + 0.1,
              y_position = c(21, 21, 80),
              tip_length = 0.01,
              annotation = "") +
  geom_text(data = Fig4C_stats[Fig4C_stats$FA %in% Fig4C_subset1,], 
            aes(x = which(Fig4C_stats$FA %in% Fig4C_subset1), 
                y = c(21, 21, 80), 
                vjust = c(0, 0, -0.5), 
                hjust = c(-0.1, -0.1, 0.6),
                label = paste0("p=", format(p.adj, scientific=T, digits=2)), 
                fill = NULL),
            colour = symnum(Fig4C_stats$p.adj[Fig4C_stats$FA
                                              %in% Fig4C_subset1],
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt, angle = c(90, 90, 0))

Fig4C_subset2 <- c("C12:0", "C14:0", "C16:1", "C24:0")                          # 2nd subset of FA

Fig4C_plot_2 <- Fig4C_plot_base +                                               # middle plot subpanel
  scale_y_continuous(expand=c(0, 0)) +                                          # y-axis settings
  coord_cartesian(ylim = c(0, 4), clip = "on") +
  labs(tag = "\u2193ZOOM\u2193", textsize = 9) +                                # adding zoom pointer
  theme(legend.position = "none",                                               # removing legend
        plot.tag.position = c(0.53, -0.04),                                     # zoom pointer settings
        plot.tag = element_text(size = 9)) +
  geom_signif(xmin = which(Fig4C_stats$FA %in% Fig4C_subset2) - 0.3,            # adding statistics
              xmax = which(Fig4C_stats$FA %in% Fig4C_subset2) + 0.1,
              y_position = c(1.7, 1.1, 1.9, 2.4),
              tip_length = 0.0005,
              annotation = "") +
  geom_text(data = Fig4C_stats[Fig4C_stats$FA %in% Fig4C_subset2,], 
            aes(x = which(Fig4C_stats$FA %in% Fig4C_subset2), 
                y = c(1.7, 1.1, 1.9, 2.4), 
                vjust = 0, 
                hjust = -0.1,
                label = paste0("p=", format(p.adj, scientific=T, digits=2)), 
                fill = NULL),
            colour = symnum(Fig4C_stats$p.adj[Fig4C_stats$FA %in% 
                                                Fig4C_subset2],
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt, angle = 90)

Fig4C_subset3 <- c("C20:0", "C22:0", "C26:0")                                   # 3rd subset of FA

Fig4C_plot_3 <- Fig4C_plot_base +                                               # bottom plot subpanel
  scale_y_continuous(expand=c(0, 0)) +                                          # y-axis settings
  coord_cartesian(ylim = c(0, 0.65), clip = "on") +
  theme(legend.position = "none") +                                             # removing legend
  geom_signif(xmin = which(Fig4C_stats$FA %in% Fig4C_subset3) - 0.3,            # adding statistics
              xmax = which(Fig4C_stats$FA %in% Fig4C_subset3) + 0.1,
              y_position = c(0.4, 0.2, 0.2),
              tip_length = 0.0001,
              annotation = "") +
  geom_text(data = Fig4C_stats[Fig4C_stats$FA %in% Fig4C_subset3,],
            aes(x = which(Fig4C_stats$FA %in% Fig4C_subset3),
                y = c(0.4, 0.2, 0.2),
                vjust = 0,
                hjust = -0.1,
                label = paste0("p=", format(p.adj, scientific=T, digits=2)),
                fill = NULL),
            colour = symnum(Fig4C_stats$p.adj[Fig4C_stats$FA %in%
                                                Fig4C_subset3],
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt, angle = c(90, 90, 90))

ggarrange(Fig4C_plot_1, Fig4C_plot_2, Fig4C_plot_3, ncol = 1, align = "v")      # merging subpanels

ggsave("Fig4C.png", width = 11, height = 17, units = "cm", dpi = 600)           # saving the plot
ggsave("Fig4C.svg", width=11, height=17, units="cm", dpi=600, fix_text_size=F)



################################################################################
# Fig4E, neutral lipid content

Fig4E_data <- read.csv(file = "Fig4E_data.csv", header = TRUE) %>%              # loading raw data
  mutate(sample = paste(genotype, treatment, sep = "_"),
         sample = factor(sample, levels = c("WT_YES", "WT_YES_N",
                                            "Dcbf11_YES", "Dcbf11_YES_N")),
         genotype = factor(genotype, levels = c("WT", "Dcbf11")),
         treatment = factor(treatment, levels = c("YES", "YES_N"))
  )

Fig4E_stats <- compare_means(amount ~ genotype,                                 # statistical tests
                             group.by = "lipid",
                             data = Fig4E_data[Fig4E_data$treatment == "YES",],
                             method = "t.test",
                             alternative = "two.sided",
                             var.equal = T,
                             p.adjust.method = "holm"
)

ggplot(data = Fig4E_data,                                                       # basic plot settings for Fig4C
       aes(x = lipid, fill = sample, y = amount)) +                             # plot initialization
  geom_bar(position="dodge", stat="summary", fun = "mean", 
           na.rm = T, colour="black", width = 0.8) +
  scale_fill_manual(
    values = c("red", "pink", "blue", "lightblue"),                             # bar colors
    labels = c("WT <span style='color:#ffffff'>_____</span> YES",               # legend labels
    "WT <span style='color:#ffffff'>_____</span> YES+NH<sub>4</sub>Cl",
    "<i>\u0394cbf11</i> <span stile='color:#ffffff'></span> YES",
    "<i>\u0394cbf11</i> <span stile='color:#ffffff'></span> YES+NH<sub>4</sub>Cl"
    )) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 50000)) +                      # y-axis settings
  ylab("Neutral lipids content (AU)") +                                         # y-axis label
  scale_x_discrete(labels = c("SE-1", "SE-2", "SQ")) +                          # removing background
  theme(legend.text = element_markdown(size = 9),
        axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_line(colour = "black", linewidth = 0.75),
        # axis.ticks.length.x = unit(0.4, "cm"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "none",
        plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")                           # plot margins
  ) +
  geom_point(position=position_dodge(0.8), na.rm=T,                             # adding individual points
             show.legend=F, size=2, alpha = 1) +
  geom_errorbar(data = aggregate(amount~lipid+sample,                           # adding errorbars
                                 data = Fig4E_data,
                                 FUN = function(Fig4E_data)
                                   c(AVG = mean(Fig4E_data, na.rm=T),
                                     SD = sd(Fig4E_data, na.rm=T))),
                aes(y = amount[,"AVG"],
                    ymin = pmax(amount[,"AVG"]-amount[,"SD"], 0),
                    ymax = amount[,"AVG"]+amount[,"SD"]),
                position = position_dodge(0.8), width=0.3,
                colour="black", linewidth=0.5
  ) +
  geom_signif(xmin = c(1 - 0.3, 2 - 0.3, 3 - 0.3),
              xmax = c(1 + 0.1, 2 + 0.1, 3 + 0.1),
              y_position = c(44000, 39000, 20000),
              tip_length = 0.03,
              annotation = ""
  ) +
  geom_text(data = Fig4E_stats, 
            aes(x = c(1 - 0.1, 2 - 0.1, 3 - 0.1), 
                y = c(44000, 39000, 20000), 
                vjust = -0.5,
                label = paste0("p=", format(p.adj, scientific=T, digits=2)), 
                fill = NULL),
            colour = symnum(Fig4E_stats$p.adj,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray")),
            size = 9/.pt)

ggsave("Fig4E.png", width = 9, height = 6, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig4E.svg", width=9, height=6, units="cm", dpi=600, fix_text_size=F)
