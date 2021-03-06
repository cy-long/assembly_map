#nolint start
library(ggplot2)
library(viridis)
library(reshape2)
library(MetBrewer)

mat_heat <- t(pathpr)
c_ratio <- 0.5

# Specify row/col names for transition matrix
names <- c("0")
for (s in 1:num) {
  names <- c(names, c(apply(sub_coms[[s]], 2, convert2names)))
}
colnames(mat_heat) = rownames(mat_heat) = paste("s",names,sep="")

# mat_heat <- log(mat_heat_2,10)
# c_ratio <- 0

# building data frame by melting matrix
molten_J <- melt(mat_heat)
names(molten_J) <- c("row", "column", "Pr")
# plot
ggplot(data = molten_J, aes(x = column, y = row, fill = Pr)) + 
  geom_tile() + 
  scale_fill_gradient2(
    low = met.brewer("Hiroshige", n = 100)[1],
    # mid = "#FFFFFF", 
    high = met.brewer("Hiroshige", n = 100)[80],
    midpoint = c_ratio, 
    # limits = c(-0.04,0.04),
    limits = c(0,1),
    na.value = "grey95"
    ) +
  coord_fixed(ratio = 8) + 
  scale_y_discrete(limits = rev(levels(molten_J$row)), expand = c(0, 0)) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        legend.position = "bottom",
        legend.title = element_text(),
        legend.text = element_text(size = 10, angle= 60,vjust = 0.5),
        legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(1, 'cm'),
        axis.text.x.top = element_text(size = 10, angle = 60, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.x = element_blank(),
        #axis.ticks.y = element_blank()
        # title = element_text("test")
  )

#nolint end