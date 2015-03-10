#' Adding customized black and white theme to ggplot objects and controls the font size.
#'
#' The function takes a ggplot object as an input and modifies the theme of it
#'
#' @param inputplot A ggplot object.  \code{s}
#' @return A ggplot object with a custom theme added
ggplot_add_theme <- function(inputplot){
  inputplot <- inputplot + theme(axis.text = element_text(size=12),
                                axis.title=element_text(size=14,face="bold"),
                                title = element_text(size=16,face="bold"),
                                legend.title = element_text(size= 10, face="bold"),
                                legend.text = element_text(size = 8),
                                legend.key = element_rect(fill = NA,color = "white"),
                                panel.grid.major = element_line(colour = "grey80", linetype="dashed"),
                                panel.background = element_rect(fill = NA,color = "grey80"))
  return(inputplot)
}
