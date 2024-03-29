library(multcomp)  #multiple comparison statistics
library(ggplot2)   #basic ggplot
library(ggpubr)    #provides plots and some statistics
library(ggtext)    #allows HTML-tags in text (element_markdown)
library(svglite)   #for saving images in svg format
library(broom)     #simplifies statistical data (tidy())

setwd("/media/sf_Genomik-data/user_data/Zemlianski Viacheslav/Cut project/Manuscript/Fig2")

theme_set(theme_classic() + theme(text = element_text(family = "Arial")))

################################################################################
# Fig2A, nucleus cross-section area at t0

Fig2A_data <- read.csv(file = "Fig2A_data.csv", header = TRUE) %>%              # loading raw data
  mutate(sample = paste(genotype, treatment, sep = "_"),
         sample = factor(sample, levels = c("WT_YES", "WT_YES_N",
                                            "Dcbf11_YES", "Dcbf11_YES_N")),
         genotype = factor(genotype, levels = c("WT","Dcbf11")),
         treatment = factor(treatment, levels = c("YES","YES_N")))

# Fig2A_t_1 <- t.test(c(Fig2A_data$area[Fig2A_data$sample == "WT_YES"]),          #statistical tests
#                     y = c(Fig2A_data$area[Fig2A_data$sample == "Dcbf11_YES"]),
#                     alternative = "greater",
#                     var.equal = T)
# Fig2A_t_2 <- t.test(c(Fig2A_data$area[Fig2A_data$sample == "Dcbf11_YES_N"]),
#                     y = c(Fig2A_data$area[Fig2A_data$sample == "Dcbf11_YES"]),
#                     alternative = "greater",
#                     var.equal = T)
# Fig2A_stats <- rbind(tidy(Fig2A_t_1), tidy(Fig2A_t_2)) %>%
#   mutate(p.adj = p.adjust(p.value, method = "holm"))
# rm(Fig2A_t_1, Fig2A_t_2)

# Trying ANOVA + post-hoc tests instead of Student's T

Fig2A_stats <- aov(area ~ sample,                                               # statistical test
                    data = Fig2A_data) %>%
  glht(mcp(sample = c("WT_YES - Dcbf11_YES <= 0", 
                      "Dcbf11_YES_N - Dcbf11_YES <= 0",
                      "WT_YES - WT_YES_N <= 0"
                      ))) %>%
  tidy()

ggplot(data = Fig2A_data, aes(x = genotype, fill = treatment, y = area)) +      # plot initialization
  ggtitle(expression("t"[0])) +                                                 # plot title
  stat_boxplot(geom ='errorbar',                                                # adding boxes
               position = position_dodge(0.8),
               width = 0.3) + 
  geom_boxplot(position = position_dodge(0.8),
               outlier.shape=NA,
               key_glyph = draw_key_polygon,
               colour="black", width = 0.7) +
  geom_point(position=position_jitterdodge(0.1),                                # adding datapoints
             alpha = 0.2,
             key_glyph = draw_key_blank) +
  scale_fill_manual(values = c("white", "gray50"),                              # boxes colors and legend labels
                    labels = c("YES", 
                               expression("YES+NH"[4]*"Cl"))) +
  scale_x_discrete(labels = c("WT", "<i>\u0394cbf11</i>")) +                    # x-labels          
  scale_y_continuous(expand=c(0, 0), limits=c(8, 22),                           # y-axis scale and labels
                     breaks=c(seq(8, 22, 2))) +
  ylab("Nuclear cross-section area (\u03BCm<sup>2</sup>)") +                    # y-axis title
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black"),
        axis.title.y = element_markdown(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = c(0.8, 0.95),
        legend.title = element_blank(),
        legend.text = element_text(size = 9, colour = "black", 
                                   vjust = 0.5, hjust = 0),
        legend.key.height = unit(0.35, 'cm'),
        legend.key.width = unit(0.35, 'cm'),
        legend.background = element_blank(),
        plot.margin=unit(c(0.7, 0.2, 0, 0.2), "cm")                             # plot margins
  ) +
  geom_signif(xmin = c(0.8, 1.8, 0.8),                                          # adding errorbars
              xmax = c(1.8, 2.2, 1.2),
              y_position = c(19.5, 16, 17.3),
              annotations = ""
              ) +
  geom_text(data = Fig2A_stats,
            aes(x = c(1.3, 2, 1),
                y = c(19.5, 16, 17.3),
                fill = NULL,
                label = paste0("p=", format(adj.p.value, scientific=T, digits=2))),
            position = position_nudge(0, 0.4),
            size = 9/.pt,
            colour = symnum(Fig2A_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray"))
            )

ggsave("Fig2A.png", width = 8, height = 10, units = "cm", dpi = 600)            # saving the plot
ggsave("Fig2A.svg", width=8, height=10, units ="cm", dpi=600, fix_text_size=F)



################################################################################
# Fig2B, nucleus cross-section area at t-max

Fig2B_data <- read.csv(file = "Fig2B_data.csv", header = TRUE) %>%              # loading raw data
  mutate(sample = paste(genotype, treatment, sep = "_"),
         sample = factor(sample, levels = c("WT_YES", "WT_YES_N",
                                            "Dcbf11_YES", "Dcbf11_YES_N")),
         genotype = factor(genotype, levels = c("WT","Dcbf11")),
         treatment = factor(treatment, levels = c("YES","YES_N")))

# Fig2B_t_1 <- t.test(c(Fig2B_data$area[Fig2B_data$sample == "WT_YES"]),          #statistical tests
#                     y = c(Fig2B_data$area[Fig2B_data$sample == "Dcbf11_YES"]),
#                     alternative = "greater",
#                     var.equal = T)
# Fig2B_t_2 <- t.test(c(Fig2B_data$area[Fig2B_data$sample == "Dcbf11_YES_N"]),
#                     y = c(Fig2B_data$area[Fig2B_data$sample == "Dcbf11_YES"]),
#                     alternative = "greater",
#                     var.equal = T)
# Fig2B_stats <- rbind(tidy(Fig2B_t_1), tidy(Fig2B_t_2)) %>%
#   mutate(p.adj = p.adjust(p.value, method = "holm"))
# rm(Fig2B_t_1, Fig2B_t_2)

# Trying ANOVA + post-hoc tests instead of Student's T

Fig2B_stats <- aov(area ~ sample,                                               # statistical tests
                   data = Fig2B_data) %>%
  glht(mcp(sample = c("WT_YES - Dcbf11_YES <= 0", 
                      "Dcbf11_YES_N - Dcbf11_YES <= 0",
                      "WT_YES - WT_YES_N <= 0"))) %>%
  tidy()

ggplot(data = Fig2B_data, aes(x = genotype, fill = treatment, y = area)) +      # plot initialization
  ggtitle(expression("t"[max])) +                                               # plot title
  stat_boxplot(geom ='errorbar',                                                # adding boxes
               position = position_dodge(0.8),
               width = 0.3) +
  geom_boxplot(position = position_dodge(0.8),
               outlier.shape=NA,
               key_glyph = draw_key_polygon,
               colour="black", width = 0.7) +
  geom_point(position=position_jitterdodge(0.1),                                # adding datapoints
             alpha = 0.2,
             key_glyph = draw_key_blank) +
  scale_fill_manual(values = c("white", "gray50"),                              # boxes colors and legend labels
                    labels = c("YES", 
                               expression("YES+NH"[4]*"Cl"))) +
  scale_x_discrete(labels = c("WT", "<i>\u0394cbf11</i>")) +                    # x-labels          
  scale_y_continuous(expand=c(0, 0), limits=c(8, 22),                           # y-axis scale and labels
                     breaks=c(seq(8, 22, 2))) +
  ylab("Nuclear cross-section area (\u03BCm<sup>2</sup>)") +                    # y-axis title
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black", linewidth = 0.75),           # axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black"),
        axis.title.y = element_markdown(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "none",                                               # removing legend
        plot.margin=unit(c(0.7, 0.2, 0, 0.2), "cm")                             # plot margins
  ) +
  geom_signif(xmin = c(0.8, 1.8, 0.8),                                          # adding errorbars
              xmax = c(1.8, 2.2, 1.2),
              y_position = c(21.4, 20, 20),
              annotations = ""
  ) +
  geom_text(data = Fig2B_stats,
            aes(x = c(1.3, 2, 1),
                y = c(21.4, 20, 20),
                fill = NULL,
                label = paste0("p=", format(adj.p.value, scientific=T, digits=2))),
            position = position_nudge(0, 0.4),
            size = 9/.pt,
            colour = symnum(Fig2B_stats$adj.p.value,
                            cutpoints = c(0, 0.05, Inf),
                            symbols = c("black", "gray"))
  )

ggsave("Fig2B.png", width = 8, height = 10, units = "cm", dpi = 600)            # saving the plot
ggsave("Fig2B.svg", width=8, height=10, units="cm", dpi=600, fix_text_size=F)
