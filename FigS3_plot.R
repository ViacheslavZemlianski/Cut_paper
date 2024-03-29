library(ggplot2)   #basic ggplot
library(ggpubr)    #provides plots and some statistics
library(ggtext)    #allows HTML-tags in text (element_markdown)

setwd("/media/sf_Genomik-data/user_data/Zemlianski Viacheslav/Cut project/Manuscript/FigS3")

#FigS3, cbf11KO vs cbf11KO tor2ts

FigS3_data <- read.csv(file = "FigS3_data.csv", header = TRUE) %>%
  mutate(genotype = factor(genotype, levels = c("WT","Dcbf11","Dcbf11_tor2ts")),
         cat_freq = cat_freq*100)

FigS3_stats <- t.test(c(FigS3_data$cat_freq[FigS3_data$genotype == "Dcbf11"]),   #statistical tests
                     y = c(FigS3_data$cat_freq[FigS3_data$genotype == "Dcbf11_tor2ts"]),
                     alternative = "two.sided",
                     var.equal = T) %>%
  tidy()

ggplot(FigS3_data, aes(x = genotype, y = cat_freq, fill = genotype)) +
  geom_bar(stat="summary", fun = "mean", 
           na.rm = T, colour="black", width = 0.3) +
  scale_fill_manual(values = c("white","grey90","grey60")) +
  geom_point(size = 2, show.legend = F) +
  geom_errorbar(data = aggregate(cat_freq~genotype,                             #adding errorbars
                                 data = FigS3_data,
                                 FUN = function(FigS3_data)
                                   c(AVG = mean(FigS3_data, na.rm=T),
                                     SD = sd(FigS3_data, na.rm=T))),
                aes(y = cat_freq[,"AVG"],
                    ymin = pmax(cat_freq[,"AVG"]-cat_freq[,"SD"], 0),
                    ymax = cat_freq[,"AVG"]+cat_freq[,"SD"]),
                width=0.1,
                colour="black", linewidth=0.75
  ) +
  scale_x_discrete(labels = c("WT", 
                              "<i>\u0394cbf11</i>", 
                              "<i>\u0394cbf11<br>tor2-S</i>")) +
  ylab("Frequency of catastrophic mitosis (%)") +
  scale_y_continuous(expand=c(0,0), limits=c(0,18), breaks=c(0,3,6,9,12,15,18)) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.75),           #axes design
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 9, colour = "black",
                                       hjust = 0.5, vjust = 0.9),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        legend.position = "none",
        plot.margin=unit(c(0.2, 0.2, 0, 0.2), "cm")
  ) +
  geom_signif(comparisons = list(c("Dcbf11", "Dcbf11_tor2ts")),
              annotation = paste0("p=", format(FigS3_stats$p.value, 
                                                digits = 2)),
              textsize = 9/.pt)

ggsave("FigS3.png", width = 8, height = 8, units = "cm", dpi = 600)             #saving the plot
ggsave("FigS3.svg", width=8, height=8, units="cm", dpi=600, fix_text_size=F)
