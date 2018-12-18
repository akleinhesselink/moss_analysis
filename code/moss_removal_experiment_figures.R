source('code/moss_removal_analysis.R')

library(scales)
library(grid)
library(gridExtra)

source('code/moss_theme.R')

#### plot parameters
xTitle = "Position on environmental gradient"
yTitle1 = 'Percent survival'
yTitle2 = 'Mean plant size (log g)\n'
yTitle3 = 'Inflorescences per plant (no.)\n'


label_df <- data.frame(label = c('a)', 
                                 'b)', 
                                 'c)', 
                                 'd)', 
                                 'e)', 
                                 'f)'), 
                       species = c('Bromus', 
                                   'Vulpia'), 
                       x_pos = 0.6, 
                       y_pos = 1, 
                       Treatment = 'Moss patch')


my_colors <- c('blue', 'black', 'red')


prop_success$measure <- 'success'
size$measure <- 'size'
infls$measure <- 'infls'

fig1 <- 
  ggplot(data = prop_success, aes( x = Stress, 
                                   y = result, 
                                   ymin = result - stand_err, 
                                   ymax = result + stand_err, group = Treatment, 
                                   color = Treatment, 
                                   shape = Treatment)) + 
  geom_point(position = position_dodge(width = 0.2), size = 4) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) +
  geom_text(data = label_df[1:2,], 
            aes( x = x_pos, y = y_pos, label = label, ymin = 0, ymax = 0, color = NULL), size = 5, show.legend = FALSE) + 
  scale_y_continuous(yTitle1, labels = percent_format()) + 
  xlab(xTitle) + 
  scale_shape_manual(values = c(19, 2, 15))  + 
  scale_color_manual(values = c('black', 'black', 'darkgray')) + 
  guides(color = guide_legend(override.aes = list(linetype = 0, size = 5))) +
  facet_grid( . ~ species) +   
  moss_theme
  

fig2 <- 
  ggplot(data = size, aes( x = Stress, y = result, ymin = result - stand_err, 
                                         ymax = result + stand_err, group = Treatment, 
                                         color = Treatment, shape = Treatment)) + 
  geom_point(position = position_dodge(width = 0.2), size = 4) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) + 
  geom_text(data = label_df[3:4,], aes( x = x_pos, y = y_pos*max(size$result) + 0.3, label = label, ymin = 0, ymax = 0), size = 5, vjust = 0.1, show.legend = FALSE) + 
  scale_y_continuous(yTitle2) + 
  xlab(xTitle) + 
  scale_shape_manual(values = c(19, 2, 15))  + 
  scale_color_manual(values = c('black', 'black', 'darkgray')) + 
  guides(color = guide_legend(override.aes = list(linetype = 0, size = 5))) +
  facet_grid( . ~ species) +   
  moss_theme

fig2

fig3 <- ggplot(data = infls, aes( x = Stress, y = result, ymin = result - stand_err, 
                                ymax = result + stand_err, group = Treatment, 
                                color = Treatment, shape = Treatment)) + 
  geom_point(position = position_dodge(width = 0.2), size = 4) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) + 
  geom_text(data = label_df[5:6,], 
            aes( x = x_pos, y = y_pos*max(infls$result) + 0.2, 
                 label = label, ymin = 0, ymax = 0), 
            size = 5, vjust = 0.1, show.legend = FALSE) + 
  scale_y_continuous(yTitle3) + 
  xlab(xTitle) + 
  scale_shape_manual(values = c(19, 2, 15))  + 
  scale_color_manual(values = c('black', 'black', 'darkgray')) + 
  guides(color = guide_legend(override.aes = list(linetype = 0, size = 5))) +
  facet_grid( . ~ species) +   
  moss_theme


fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}

top_panel <- 
  fig1 + 
  theme( 
    legend.position = 'none', 
    strip.background = element_blank(),
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    plot.margin = margin( 0, 11.3, 0, 4, unit = 'line'))

mid_panel <- 
  fig2 + 
  scale_y_continuous(yTitle2, labels = fmt_dcimals(1)) + 
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    strip.text = element_blank(), 
    plot.margin = margin( 0, 1, 0, 4, unit = 'line'))


bottom_panel <- 
  fig3 + 
  theme(legend.position = 'none', 
        strip.text = element_blank(), 
        plot.margin = margin( 0, 11.3, 1, 4, unit = 'line'), 
        axis.title.x = element_text(margin = margin(1, 0, 0, 0, unit = 'line')))  


plot_all <- 
  grid.arrange(arrangeGrob(top_panel, mid_panel, bottom_panel, 
                           ncol = 1, 
                           nrow = 3, 
                           heights = c(1, 0.8, 1.12) ))


ggsave(path = 'figures', 
       filename = 'all.plots.png', 
       plot = plot_all, 
       width = 8, 
       height = 9, 
       units = 'in', 
       dpi = 300)

