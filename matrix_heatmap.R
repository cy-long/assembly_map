#nolint start
library(ggplot2)
library(viridis)
library(reshape2)
library(MetBrewer)

mat_heat_1 <- (matL- matR)/(matL+matR)/2
mat_heat_2 <- pmax(matL, matR)/matO

mat_heat <- mat_heat_2
c_ratio <- 1

# mat_heat <- log(mat_heat_2,10)
# c_ratio <- 0

# building data frame by melting matrix
molten_J <- melt(mat_heat)
names(molten_J) <- c("row", "column", "value")
# plot
ggplot(data = molten_J, aes(x = column, y = row, fill = value)) + 
  geom_tile() +
  coord_equal() +
  scale_fill_gradient2(low = met.brewer("Hiroshige", n = 100)[1], mid = "#FFFFFF", 
                       high = met.brewer("Hiroshige", n = 100)[80], midpoint = c_ratio, 
                       limits = c(min(mat_heat),max(mat_heat)), na.value = "grey97") +
  #scale_fill_viridis(option = "plasma", limits = c(min_j, max_j)) +
  scale_y_discrete(limits = rev(levels(molten_J$row)), expand = c(0, 0)) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(1, 'cm'),
        axis.text.x.top = element_text(size = 28, angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 28),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.x = element_blank(),
        #axis.ticks.y = element_blank()
        #title = element_text("Test")
  )
#nolint end