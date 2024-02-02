# make png fig of summary...
# grouping specifies the name of the column which should be used to subset the data, can be used for different clusters for example 
multispec <- function(df, grouping, theme = 'aurora', bins = 1000, normalize = FALSE, params = NULL, vlines = c(5.5,8.5,30.5,50.5,66.5,78.5), scaled = FALSE, ncol = 4, nrow = 2){
  if (is.null(params)){
    p1 <- grep("FSC|LightLoss|SSC", colnames(df), value = TRUE)
    p1 <- grep("-A", p1, value = TRUE)
    p2 <- grep("\\d.*-A", colnames(df),value = TRUE)
    params <- c(p1,p2)
  }
  plots <- list()
  for (i in sort(unique(df[,grouping]))){
    tmp <- df[df[,grouping]==i,]
    plots[[i]] <- spectralplot_df(df = tmp, 
                                  theme = theme, 
                                  save = FALSE, 
                                  bins = bins, 
                                  normalize = normalize,
                                  params = params, 
                                  vlines = vlines,
                                  scaled = scaled,
                                  plot_title = i)
  }
  final_plot <- ggpubr::ggarrange(plotlist = plots, ncol = ncol, nrow = nrow)
  return(final_plot)
}