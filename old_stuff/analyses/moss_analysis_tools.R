#### Moss analysis tools

counter = function(data_row){
  #### Counts total number of emergents
  n = 0
  n_total = 0
  for (i in 1:length(data_row)){
    ngerm = data_row[i] - n
    if(ngerm > 0){
      n_total = as.numeric(ngerm + n_total)
    }
    n = data_row[i]
  }
  return(n_total)
} 

emergence_over_time = function(df_all, df_low, df_high, rows = c(3:5)){ 
  #### make side by side plots for emergence over days, high and low stress
  par(mfrow = c(1,2))
  matplot(x = df_all$time, df_all[, 3:5], type = 'n', xlab = "Days since planting", ylab = "No. plants/plot")
  matlines(x = df_low$time, df_low[, 3:5], col = 1, pch = 1:3, lty = 1:3, type = 'b')
  title (main = "Vulpia low stress")
  matplot(x = df_all$time, df_all[, 3:5], type = 'n', xlab = "Days since planting", ylab = "No. plants/plot")
  matlines(x = df_high$time, df_high[, 3:5], col = 1, pch = 1:3, lty = 1:3, type = 'b')
  title (main = "Vulpia high stress")
  legend('bottomright', c("bare control", "moss covered", "moss removed"), pch = 1:3, lty = 1:3, bty = 'n') 
}

calc_interaction_plot_means = function(df, block = "block", t1 = 'stress', t2 = "treatment", response = "mass_per_plant_mg", lg = TRUE){ 
  
  plot_means = with(df, aggregate(log(df[,  response]), list(t1 = df[, t1], t2 = df[, t2]),  mean))  
  plot_means$se = with(df, aggregate(log(df[, response]), list(t1 = df[, t1], t2 = df[, t2]), function(x) sd(x)/sqrt(44 - 4)))[, 3]
  
  return(plot_means)
}


treatment_by_stress_plots = function(df, t1 = 'stress', t2 = 'treatment', y_var = 'mass_per_plant_mg', plot_means, lg = TRUE){ 
  
  par(mfrow = c(1,1), mgp = c(1.5,0.2,0),  las = 1, xpd = FALSE, oma = c(1,1,1,1), mar = c(3,3,1,0), fin = c(5,5), tck = 0.02, cex = 1)  
  
  
  interaction.plot(df[, t1], df[, t2], df[, y_var], xlab = "", bty = 'l',                 
                   ylab = expression(paste(italic(species), " biomass (log g)")), type = 'b', legend = FALSE,
                   lty = c(2,1), cex = 0.75, cex.axis = 0.75, cex.lab = 0.75, xaxt = 'n')
  
  axis(side = 1, labels = FALSE, tick = TRUE, lwd = 0, lwd.ticks = 1, tcl = 0.2, at = c(1,2))                 
  title(xlab = "stress level", line = 1,  xpd = TRUE, cex.lab = 0.9)
  text(1, -1.4, "high", xpd = TRUE, cex = 0.75)
  text(2, -1.4, "low", xpd = TRUE, cex = 0.75)
  
  #arrows(1, plot_means$x[1] + plot_means$se[1], 1, plot_means$x[1] - plot_means$se[1], angle = 90, code = 3,col = 1, length = 0.1)
  #arrows(2, plot_means$x[3] + plot_means$se[3], 2, plot_means$x[3] - plot_means$se[3], angle = 90, code = 3,col = 2, length = 0.1)
  #arrows(1, plot_means$x[2] + plot_means$se[2], 1, plot_means$x[2] - plot_means$se[2], angle = 90, code = 3,col = 3, length = 0.1)
  #arrows(2, plot_means$x[4] + plot_means$se[4], 2, plot_means$x[4] - plot_means$se[4], angle = 90, code = 3,col = 4, length = 0.1)
  

}





