
PlotTheme <- list(FontColor = theme(plot.title = element_text(colour = "black", face = "bold", hjust = 0.5),
                                    axis.title = element_text(colour = "black"),
                                    axis.text = element_text(colour = "black"),
                                    legend.text = element_text(colour = "black")),
                  Color = rep(c("#1F77B4", "#E39802", "#2CA02C", "#B60A1C", "#9467BD",
                                "#C46487", "#7F7F7F", "#17BECF", "#AEC7E8", "#FF7F0E",
                                "#638B66", "#D62728", "#C5B0D5", "#F7B6D2", "#C7C7C7",
                                "#5FA2CE", "#FFBB78", "#98DF8A", "#FF9896", "#8C564B",
                                "#DBD4C5", "#E377C2", "#57606C", "#BCBD22", "#FFDD71",
                                "#C3CE3D", "#C49C94", "#9EDAE5"), 5),
                  Box = theme(axis.line = element_blank(),
                              panel.border = element_rect(fill = NA, colour = "black", linewidth = 1)))