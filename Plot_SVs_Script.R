### This code plots interesting regions of a reference genome. 
### In my case, disease-candidate SVs in an specific specie. 

# Load the "ggplot2" library.
library("ggplot2")

# Define the function with all the plots.
plot_SVs_distribution <- function(RefGenome, SVsFile, ZoomYesOrNot='Yes', xZoom=c(1,18), yZoom=c(1, 10e9)) {
  
  # Read the Genome file and create a factor object with the Chromosome's names column. 
  Lens <- read.table(RefGenome, header = FALSE, col.names = c("Chromosome", "Size"))
  Lens$Chromosome <- factor(Lens$Chromosome, levels = Lens$Chromosome)
  
  # Repeat it with the SVs file. 
  SVs <- read.table(SVsFile, header = FALSE, col.names = c("Chromosome", "Start","End","Type"))
  SVs$Chromosome <- factor(SVs$Chromosome, levels = Lens$Chromosome)
  SVs$SV_Type<- factor(SVs$Type, levels = c('INV','DEL','BND')) # In my case these were the only SV types presented.
  
    # Let's do the plots.
  # Define the plot of the chromosomes.
  genome_plot <- ggplot() +
    geom_bar(data = Lens, 
             aes(x = Chromosome, y = Size),
             stat='identity',
             fill='grey80',
             width=0.3) +
    theme_minimal() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Define the plot of the SVs. The "yend" parameter is multiplied by 1.01 because otherwise
  # the SVs are not shown (as the size relation between SVs and full chromosomes is so low) 
  variants <- geom_segment(data = SVs, 
                           aes(x=Chromosome, xend=Chromosome, y=Start, yend=(Start+End)*1.01),
                           color='blue',
                           linewidth = 5)
  
  # The final general plot combines chromosomes and SVs.             
  final_plot <- genome_plot + variants + 
    labs(title = "Disease causant candidate SVs", y = "Base pairs", x = "Chromosomes")
  
  # If your SVs are so concentrated in a certain region (that's my case), you can zoom in there.
  if (ZoomYesOrNot=="Yes") {
    # ¡! If you want to focus on X, Y or other non-numerical named chromosomes... ex/ In human genome, X and Y would be 22 and 23.
    variants4zoom <- geom_segment(data = SVs, 
                                  aes(x=Chromosome, xend=Chromosome, y=Start, yend=(Start+End)*1.001), 
                                  color='blue',
                                  linewidth = 50) 
    
    labels_SV <- geom_point(data = SVs, aes(x = as.numeric(Chromosome), y = Start, 
                                            shape = SV_Type,color=SV_Type), size = 2) 
                            #position=position_jitterdodge(dodge.width = 0.5)) -> useful if there are a lot of SVs in the cluster
    
    final_plot_zoomed <- genome_plot + variants4zoom + labels_SV +
      coord_cartesian(ylim = yZoom, xlim= xZoom) +
      scale_shape_manual(values = c(15, 17, 19, 18)) +  # Define shapes: Square, Triangle, Circle, Diamond
      scale_color_manual(values = c('red','black','orange')) +
      labs(title = "Disease causant candidate SVs", y = "Base pairs", x = "Chromosome", subtitle='Zoomed Chromosome')
    
    # And the general plot with the clustered zone labeled. 
    final_plot_square <- final_plot + annotate("rect", xmin = xZoom[1]-0.5, xmax = xZoom[1]+0.5, 
                                               ymin = yZoom[1], ymax = yZoom[2], alpha = 0.2, fill = "red") +
      labs(subtitle='Clustered SVs')
    
    return(list('GeneralPlot' = final_plot, 'SquaredPlot' = final_plot_square, 'ZoomedPlot' = final_plot_zoomed))
  }
  else{
    return(list('GeneralPlot' = final_plot))
  }
}

# Save the results (plots) of the function when it's applied to your dataset.
all_plots <- plot_SVs_distribution('chr_tab_len.txt','SVs.txt',"Yes",xZoom=c(17,17), yZoom=c(1.9e7, 3e7)) # For ex., focus on chr 1 (e3-e5)
print(all_plots)

# Save plots automatically
n <- 0
for (i_plot in all_plots) {
  n <- n+1
  ggsave(filename = paste("plot",as.character(n),'.jpeg',sep=''), plot = i_plot, width = 6, height = 4)
}
