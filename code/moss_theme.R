##### make ggplot theme for moss paper 

moss_theme <- 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.key = element_blank(),
    legend.key.width = unit(3, 'line'),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.title.x = element_text(vjust = -0.5, size = 14),
    axis.title.y = element_text(vjust = 0, size = 14),
    axis.text = element_text(size = 12, color = 'black'),
    plot.margin = unit(c(1, 1, 1, 3), 'line'),
    strip.text = element_text(face = 'italic', size = 15), 
    axis.ticks.length = unit( 5, unit = 'pt') 
  )

xlab_distance <- 'Position towards NW on environmental gradient (m)'

