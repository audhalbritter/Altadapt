### FIGURE S3
library("cowplot")

# Defien colours
sig <- rgb(red=0, green=0, blue=0, alpha=170, max=255)
non.sig <- rgb(red=140, green=140, blue=140, alpha=100, max=255)


# Combine data
dat5 <- dat4 %>% 
  mutate(p.value = ifelse(is.na(p.value), 5, p.value)) %>% 
  mutate(pch.nr = ifelse(p.value == 5, 17, 16)) %>% # shape of points, NA = +
  mutate(Pvalue = ifelse(p.value < 0.05, "< 0.05", ifelse(p.value == 5, "NA", "> 0.05"))) %>% # point colour
  mutate(SampleSize = as.numeric(cut(sample.size, c(0,5,10,20,40,52)))) # point size

# Settings
axis.dim <- theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 12),
                  axis.ticks = element_line(size = 1), 
                  plot.title = element_text(size = 14),
                  legend.title = element_text(size = 12))


#### HEIGHT ####
d1 <- dat5 %>% filter(var == "height") %>% filter(!is.na(breed))
yylim <- range(d1$slopes)

SlopePlot <- ggplot(d1, aes(x = r.temp, y = abs.slopes)) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, color = "black", linetype = "dashed", size = 0.7) +
  geom_point(aes(size = sample.size, colour = Pvalue, shape = Pvalue, alpha = Pvalue)) + 
  ylim(-0.2, yylim[2]) +
  scale_colour_manual(values = c("grey50","grey80", "black")) +
  scale_alpha_manual(values = c(0.8, 0.8, 1)) +
  scale_shape_manual(values = c(16, 16, 17)) +
  scale_size(name = "Sample size", labels = c(" ", "0 - 5", "6 - 10", "11 - 20", "21 - 40", "> 50")) +
  geom_hline(yintercept = 0, col = "grey") +
  labs(x = "", y = "Absolute regression coefficient") +
  panel_border(colour = "black", remove = FALSE) +
  axis.dim

# r.temp
h1 <- SlopePlot + labs(title = "Height")


#### BIOMASS #### 
d2 <- dat5 %>% filter(var == "biomass") %>% filter(!is.na(breed))
yylim <- c(-0.2, d2$slopes[2])

#r.temp
b1 <- SlopePlot %+% d2 + labs(x = "Temperature range in °C", y = " ", title = "Biomass")

#### PHENOLOGY #### 
d3 <- dat5 %>% filter(var == "phenology") %>% filter(!is.na(breed))
yylim <- c(-0.2, d3$slopes[2])

#r.temp
p1 <- SlopePlot %+% d3 + labs(x = "", y = " ", title = "Phenology")



#### MAIN FIGURE
Fig3SlopePlots <- plot_grid(h1 + theme(legend.position="none"), 
                            b1 + theme(legend.position="none"), 
                            p1 + theme(legend.position="none"), 
                            align = "hv",
                            hjust = -1,
                            nrow = 1)

# extract legend with 
SlopePlot2 <- ggplot(d1, aes(x = r.temp, y = slopes)) +
  geom_point(aes(size = sample.size), color = "grey50") + 
  ylim(yylim) +
  scale_size(name = "Sample size", labels = c(" ", "0 - 5", "6 - 10", "11 - 20", "21 - 40", "> 50")) +
  geom_hline(yintercept = 0, col = "grey") +
  labs(x = "Temperature range in °C", y = "Regression coefficient") +
  panel_border(colour = "black", remove = FALSE) +
  axis.dim

p6 <- SlopePlot2 + geom_point(aes(size = SampleSize), color = "grey50")
legend_t <- get_legend(p6 + theme(legend.position="right"))
Fig3MainSlopePlot <- plot_grid(Fig3SlopePlots, legend_t, ncol = 2, rel_widths = c(1,0.2))
ggsave("FinalFigures/Fig3MainSlopePlot.pdf", Fig3MainSlopePlot, height = 4, dpi = 300)

