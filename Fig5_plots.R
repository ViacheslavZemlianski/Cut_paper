library(ggplot2)    #basic ggplot
library(ggpubr)     #provides plots and some statistics
library(ggtext)     #allows HTML-tags in text (element_markdown)
library(svglite)    #for saving images in svg format
library(broom)      #simplifies statistical data (tidy())
library(ggbeeswarm) #make symmetric datapoints distribution (geom_beeswarm)

setwd("/media/sf_Genomik-data/user_data/Zemlianski Viacheslav/Cut project/Manuscript/Fig5")

theme_set(theme_classic() + theme(text = element_text(family = "Arial")))

Fig5_data <- read.csv(file = "Fig5_data.csv", header = TRUE) %>%                # loading raw data
  mutate(sample = paste(genotype, treatment, sep = "_"),
         sample = factor(sample, levels = c("WT_YES", "WT_YES_N",
                                            "Dcbf11_YES", "Dcbf11_YES_N")),
         genotype = factor(genotype, levels = c("WT","Dcbf11")),
         treatment = factor(treatment, levels = c("YES","YES_N")),
         phase = factor(phase, levels = c("A", "PM")))

################################################################################
#Fig5B, length of the whole mitosis

Fig5_data_PM <- Fig5_data[Fig5_data$phase == "PM",]                             # generating data subset
Fig5_data_A <- Fig5_data[Fig5_data$phase == "A",]
Fig5B_data <- Fig5_data_PM[,c("genotype","treatment","sample")] %>%
  mutate(length = Fig5_data_PM$length + Fig5_data_A$length)
rm(Fig5_data_PM, Fig5_data_A)

# Fig5B_t_1 <- t.test(c(Fig5B_data$length[Fig5_data$sample == "WT_YES"]),         #statistical tests
#                     y = c(Fig5B_data$length[Fig5_data$sample == "Dcbf11_YES"]),
#                     alternative = "less",
#                     var.equal = T)
# Fig5B_t_2 <- t.test(c(Fig5B_data$length[Fig5_data$sample == "Dcbf11_YES_N"]),
#                     y = c(Fig5B_data$length[Fig5_data$sample == "Dcbf11_YES"]),
#                     alternative = "less",
#                     var.equal = T)
# Fig5B_stats <- rbind(tidy(Fig5B_t_1), tidy(Fig5B_t_2)) %>%
#   mutate(p.adj = p.adjust(p.value, method = "holm"))
# rm(Fig5B_t_1, Fig5B_t_2)

Fig5B_stats <- aov(length ~ sample,                                             # statistical test
                   data = Fig5B_data) %>%
  glht(mcp(sample = c("WT_YES - Dcbf11_YES >= 0", 
                      "Dcbf11_YES_N - Dcbf11_YES >= 0"
  ))) %>%
  tidy()

ggplot(data = Fig5_data, aes(x = sample, fill = phase, y = length)) +           # plot initialization
  geom_bar(stat="summary", fun = "mean", colour="black") +                      # adding bars 
  scale_fill_manual(values = c("gray90", "gray40"),                             # bar colors
                    labels = c("Anaphase", "Prophase + metaphase")) +
  scale_x_discrete(labels = c("-","+","-","+")) +                               # x-axis labels
  annotate(geom = "text", x = c(1.5, 3.5), y = -7.5, 
           label = c("WT", expression(italic("\u0394cbf11"))), size = 3) +
  labs(tag = expression("NH"[4]*"Cl"), textsize = 9) +
  geom_hline(yintercept = -5.5) +
  geom_beeswarm(data = Fig5B_data,                                              # adding datapoints
                aes(x = sample, y = length, fill = NULL),
                cex = 2,
                alpha = 0.3,
                key_glyph = draw_key_blank) +
  ylab("Mitosis (P+M+A) length (min)") +                                        # y-axis title
  scale_y_continuous(expand=c(0, 0), limits=c(-10, 60),                         # y-axis scale and labels
                     breaks=c(0, 10, 20, 30, 40, 50, 60)) +
  coord_cartesian(ylim = c(0, 60), expand = T, clip = "off") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 9, vjust = 0),
        plot.tag.position = c(0, 0),
        plot.tag = element_text(size = 9, hjust = 0, vjust = 0.3),
        axis.title.y = element_markdown(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "top",
        legend.title = element_blank(),                                         # legend design
        legend.text = element_text(size = 9, colour = "black"),
        legend.key.height = unit(0.35, 'cm'),
        legend.key.width = unit(0.35, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(0.2, 0.2, 0.7, 0.2), "cm")                           # plot margins
  ) +
  geom_errorbar(data = aggregate(length~sample,                                 # adding errorbars
                                 data = Fig5B_data, 
                                 FUN = function(Fig5B_data) 
                                   c(AVG = mean(Fig5B_data, na.rm=T),
                                     SD = sd(Fig5B_data, na.rm=T))), 
                aes(y = length[,"AVG"],
                    ymin = pmax(length[,"AVG"]-length[,"SD"], 0),
                    ymax = length[,"AVG"]+length[,"SD"],
                    fill = NULL),
                width=0.3,
                colour="black", linewidth=0.75
  ) +
  geom_signif(xmin = c(1, 3),                                                   # adding statistics
              xmax = c(3, 4),
              y_position = c(59.5, 52),
              annotations = ""
  ) +
  geom_text(data = Fig5B_stats,
            aes(x = c(2, 3.5),
                y = c(59.5, 52),
                vjust = -0.5,
                fill = NULL,
                label = paste0("p=", format(adj.p.value, scientific=T, digits=2))),
            position = position_nudge(0, 0.4),
            size = 9/.pt,
            colour = symnum(Fig5B_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray"))
  )

ggsave("Fig5B.png", width = 8, height = 10, units = "cm", dpi = 600)            # saving the plot
ggsave("Fig5B.svg", width=8, height=10, units="cm", dpi=600, fix_text_size=F)



################################################################################
#Fig5C, prophase + metaphase length

Fig5C_data <- Fig5_data[Fig5_data$phase == "PM",]                               # generating data subset

Fig5C_stats <- aov(length ~ sample,                                             # statistical tests
                    data = Fig5C_data) %>%
  glht(mcp(sample = c("WT_YES - Dcbf11_YES >= 0", 
                      "Dcbf11_YES_N - Dcbf11_YES >= 0"
                      ))) %>%
  tidy()

# Fig5C_t_1 <- t.test(c(Fig5C_data$length[Fig5C_data$sample == "WT_YES"]),        #statistical tests
#                     y = c(Fig5C_data$length[Fig5C_data$sample == "Dcbf11_YES"]),
#                     alternative = "less",
#                     var.equal = T)
# Fig5C_t_2 <- t.test(c(Fig5C_data$length[Fig5C_data$sample == "Dcbf11_YES_N"]),
#                     y = c(Fig5C_data$length[Fig5C_data$sample == "Dcbf11_YES"]),
#                     alternative = "less",
#                     var.equal = T)
# Fig5C_stats <- rbind(tidy(Fig5C_t_1), tidy(Fig5C_t_2)) %>%
#   mutate(p.adj = p.adjust(p.value, method = "holm"))
# rm(Fig5C_t_1, Fig5C_t_2)

ggplot(data = Fig5C_data, aes(x = genotype, fill = treatment, y = length)) +    # plot initialization
  stat_boxplot(geom ='errorbar',                                                # adding boxes
               position = position_dodge(0.8),
               width = 0.3) + 
  geom_boxplot(position = position_dodge(0.8),
               outlier.shape=NA,
               key_glyph = draw_key_polygon,
               colour="black", width = 0.7) +
  geom_beeswarm(dodge.width = 0.8,                                              # adding datapoints
                cex = 1.7,
                alpha = 0.3,
                key_glyph = draw_key_blank) +
  scale_fill_manual(values = c("white", "gray50"),                              # boxes colors and legend labels
                    labels = c("YES",
                               expression("YES+NH"[4]*"Cl"))) +
  scale_x_discrete(labels = c("WT", "<i>\u0394cbf11</i>")) +                    # x-labels
  scale_y_continuous(expand=c(0, 0), limits=c(5, 46.5),                         # y-axis scale and labels
                     breaks=c(5, 15, 25, 35, 45)) +
  ylab("Prophase+metaphase length (min)") +                                     # y-axis title
  theme(plot.title = element_text(hjust = 0.5),                                 # plot title formatting
        panel.grid.major = element_blank(),                                     # removing background
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black"),
        axis.title.y = element_markdown(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = c(0.2, 0.9),                                          # legend design
        legend.title = element_blank(),
        legend.text = element_text(size = 9, colour = "black",
                                   vjust = 0.5, hjust = 0),
        legend.key.height = unit(0.35, 'cm'),
        legend.key.width = unit(0.35, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(0.2, 0.2, 0, 0.2), "cm")                             # plot margins
  ) +
  geom_signif(xmin = c(0.8, 1.8),                                               # adding errorbars
              xmax = c(1.8, 2.2),
              y_position = c(36, 38),
              annotations = ""
  ) +
  geom_text(data = Fig5C_stats,
            aes(x = c(1.3, 2),
                y = c(36, 38),
                vjust = -0.4,
                fill = NULL,
                label = paste0("p=", format(adj.p.value, scientific=F, digits=2))),
            position = position_nudge(0, 0.4),
            size = 9/.pt,
            colour = symnum(Fig5C_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray"))
  )

ggsave("Fig5C.png", width = 8, height = 8, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig5C.svg", width=8, height=8, units="cm", dpi=600, fix_text_size=F)



################################################################################
#Fig5D, anaphase length

Fig5D_data <- Fig5_data[Fig5_data$phase == "A",]                                # generating data subset

Fig5D_stats <- aov(length ~ sample,                                             # statistical tests
                   data = Fig5D_data) %>%
  glht(mcp(sample = c("WT_YES - Dcbf11_YES >= 0", 
                      "Dcbf11_YES_N - Dcbf11_YES >= 0"
  ))) %>%
  tidy()

# Fig5D_t_1 <- t.test(c(Fig5D_data$length[Fig5D_data$sample == "WT_YES"]),        #statistical tests
#                     y = c(Fig5D_data$length[Fig5D_data$sample == "Dcbf11_YES"]),
#                     alternative = "less",
#                     var.equal = T)
# Fig5D_t_2 <- t.test(c(Fig5D_data$length[Fig5D_data$sample == "Dcbf11_YES_N"]),
#                     y = c(Fig5D_data$length[Fig5D_data$sample == "Dcbf11_YES"]),
#                     alternative = "less",
#                     var.equal = T)
# Fig5D_stats <- rbind(tidy(Fig5D_t_1), tidy(Fig5D_t_2)) %>%
#   mutate(p.adj = p.adjust(p.value, method = "holm"))
# rm(Fig5D_t_1, Fig5D_t_2)

ggplot(data = Fig5D_data, aes(x = genotype, fill = treatment, y = length)) +    # plot initialization
  stat_boxplot(geom ='errorbar',                                                # adding boxes
               position = position_dodge(0.9),
               width = 0.3) + 
  geom_boxplot(position = position_dodge(0.9),
               outlier.shape=NA,
               key_glyph = draw_key_polygon,
               colour="black", width = 0.7) +
  geom_beeswarm(dodge.width = 0.9,                                              # adding datapoints
                cex = 1.3,
                alpha = 0.3,
                key_glyph = draw_key_blank) +
  scale_fill_manual(values = c("white", "gray50"),                              # boxes colors and legend labels
                    labels = c("YES",
                               expression("YES+NH"[4]*"Cl"))) +
  scale_x_discrete(labels = c("WT", "<i>\u0394cbf11</i>")) +                    # x-labels
  scale_y_continuous(expand=c(0, 0), limits=c(5, 35),                           # y-axis scale and labels
                     breaks=c(5, 10, 15, 20, 25, 30, 35)) +
  ylab("Anaphase length (min)") +                                               # y-axis title
  theme(plot.title = element_text(hjust = 0.5),                                 # plot title formatting
        panel.grid.major = element_blank(),                                     # removing background
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black"),
        axis.title.y = element_markdown(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "none",                                               # legend design
        plot.margin=unit(c(0.2, 0.2, 0, 0.2), "cm")                             # plot margins
  ) +
  geom_signif(xmin = c(0.77, 1.77),                                             # adding errorbars
              xmax = c(1.77, 2.27),
              y_position = c(30, 25),
              annotations = ""
  ) +
  geom_text(data = Fig5D_stats,
            aes(x = c(1.27, 2.03),
                y = c(30, 25),
                vjust = -0.3,
                fill = NULL,
                label = paste0("p=", format(adj.p.value, scientific=T, digits=2))),
            position = position_nudge(0, 0.4),
            size = 9/.pt,
            colour = symnum(Fig5D_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray"))
  )

ggsave("Fig5D.png", width = 8, height = 8, units = "cm", dpi = 600)             # saving the plot
ggsave("Fig5D.svg", width=8, height=8, units="cm", dpi=600, fix_text_size=F)
