##### make ggplot theme for moss paper 

moss_theme = theme_set(theme_bw() + 
                       theme( panel.grid.major = element_blank(), 
                              legend.key = element_blank(), 
                              legend.key.width = unit(3, 'line'),
                              legend.text = element_text(size = 12),
                              legend.title = element_text(size = 12),
                              axis.title.x = element_text(vjust = -0.5, size= 14), 
                              axis.title.y = element_text(vjust = 0, size = 14), 
                              axis.text = element_text(size = 12), 
                              plot.margin = unit(c(1,1,1,3), 'line')))

xlab_distance = 'Position away from SE end of gradient (m)'
