### INFO: This script takes 5 pairs of files. These files are pairs of: 
# SVs detected in a comparison between genetic pools, sorted by p-value; 
# and a file with the number of statistically significant SVs after some tests. 
### FUNCTION: It generates a Venn Diagram to represent how many of these SVs are statistically significant in both comparisons and how many are private of one or the other comparison.

library("VennDiagram")

# Define a function to process data so then the venn.diagram function can correctly interprete it.
pre_venn_function <- function(SVs_file,stats_file) {

  # Read the file with SVs and the file with the number of statistically significant SVs.
  SVs <- read.table(SVs_file,header=FALSE,sep="\t")
  sig <- as.numeric(unlist(strsplit(readLines(stats_file,warn=FALSE)[3], "\t")))
  Bonf <- sig[4] # 4th number corresponds to the number of SVs afeter Bonferroni correction. 
  
  # Define all columns we are interested in. 
  chr <- SVs$V1
  chr_start <- SVs$V2
  chr2 <- SVs$V4
  pos2 <- SVs$V5
  len <- SVs$V6
  sv_type <- SVs$V9
  
  # For each SV, generate a new SV with only the necessary information to compare SVs.
  new_SVs <- list()
  for (i in 1:Bonf) {
    new_SVs[i] <- paste(chr[i],chr_start[i],chr2[i],pos2[i],len[i],sv_type[i],sep="_")
  }
  
  return (list(unlist(new_SVs),Bonf)) # We want the new SVs and the number of statistically significant SVs after Bonferroni.
}

# Define "args" to let the user specify arguments in the terminal
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)

# Define empty lists to add the results of the pre_venn_function function.
x <- list()
Bonfs <- list()

# For each pair of arguments (SVs and statistics file):
n <- 0
for (i in seq(1,n_args,by=2)){
  n <- n+1
  y <- pre_venn_function(args[[i]],args[[i+1]]) # Apply the pre-processig function.
  x[[n]] <- y[[1]] # Save results in the "x" list.
  Bonfs[[n]] <- y[[2]] # And save the Bonferroni number.
}

# Generate the Venn diagram with this data.
venn.diagram(x,
             filename="Venn_Diagram.png",
             category = c("A","B","C","D","E"), # Names of the circles in the diagram.
             fill = c("olivedrab1","royalblue1","chocolate1","hotpink","goldenrod2"), # Colors of the circles.
             lwd = 0.5, # Circles' perimeter width.
             cex = 1.5, # Font size for numbers.
             cat.cex = 1.5, # Font size for names.
             cat.default.pos = "outer", cat.dist = c(0.2,0.25,0.2,0.2,0.27), # To relocate the names of the circles
             height = 5000 , width = 5000, # Image size.
             margin = 0.11) # Margins of the image
