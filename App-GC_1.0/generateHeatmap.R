generateHeatmap <- function(aligned_data){
  # Alignment Heatmap 
  p <- ggplot(aligned_data, aes(x=run, y=mean_RT, key=mean_RT)) + geom_tile(aes(fill=aligned), color="grey") + 
    ggtitle("Peak Alignment") + xlab("") + ylab("Average RT") + labs(fill = "Alignment") +
    scale_fill_manual(values=c("#999999","darkred")) +
    theme_minimal() + theme(axis.text=element_text(size=8))
  ggplotly(p, tooltip = c("x","y"), source="heatmap") %>% layout(xaxis=list(fixedrange=TRUE)) %>% layout(yaxis=list(fixedrange=TRUE))
}
