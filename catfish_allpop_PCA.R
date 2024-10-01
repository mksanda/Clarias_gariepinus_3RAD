#!/usr/bin/env Rscript

#
# 
# Modified code from Angel Rivera-Colon and Julian Catchen's Stacks workshop
#
#
library(adegenet)
library(ggplot2)
library(RColorBrewer)

# Specify inputs
dir <- '.'

genepop_f <- paste(dir, '/populations.snps.gen', sep='')
# Set color palette
cols <- brewer.pal(9, 'Paired')

# Load a GenePop file into a GenInd object, add population info 
genepop_to_genind <- function(genepop_f, pop.idx=1){
  genInd <- import2genind(genepop_f)
  pop(genInd) <- sapply(strsplit(indNames(genInd), '_'),
                        function(genInd){paste(genInd[2], genInd[1], sep='')})
  return(genInd)
}

# Generate a population sample table from the GenInd object
make_samples_tbl <- function(genInd){
  samples <- indNames(genInd)
  pops <- pop(genInd)
  #lib <- sapply(strsplit(indNames(genInd), '_'),
  #function(genInd){genInd[3]})
  samples.df <- data.frame(samples, pops)
  colnames(samples.df) <- c('Indv', 'Pop')
  return(samples.df)
}


# Run a PCA on the genInd
run_pca <- function(genInd){
  # Scale genInd
  genInd.scaled <- scaleGen(genInd, NA.method='mean')
  # Do PCA
  pca <- dudi.pca(genInd.scaled, scannf=FALSE, nf=10)
  return(pca)
}

get_pca_eigenvalues <- function(pca, max_eigen=10){
  if (max_eigen <= 2){
    stop('Error: max_eigen must be > 2')
  }
  # Initialize output vector
  eig <- vector("numeric", max_eigen)
  # Calculate proportions
  eigTotal <- sum(pca$eig)
  for (i in 1:max_eigen){
    pc <- pca$eig[i]/eigTotal # % explained by a PC
    eig[i] <- pc
  }
  return(eig)
}

# Add PCs to the samples dataframe
add_pcs_to_df <- function(samples.df, pca){
  samples.df$PC1 <- pca$li$Axis1
  samples.df$PC2 <- pca$li$Axis2
  samples.df$PC3 <- pca$li$Axis3
  samples.df$PC4 <- pca$li$Axis4
  samples.df$PC5 <- pca$li$Axis5
  samples.df$PC6 <- pca$li$Axis6
  samples.df$PC7 <- pca$li$Axis7
  samples.df$PC8 <- pca$li$Axis8
  samples.df$PC9 <- pca$li$Axis9
  samples.df$PC10 <- pca$li$Axis10
  return(samples.df)
}

# Make PCA scatter plot in ggplot2
plot_pca <- function(samples.df, eigenvalues, colors_palette,
                     axis_x=1, axis_y=2, outdir='.',
                     min_x=NA, max_x=NA, min_y=NA, max_y=NA){
  pcx <- paste('PC', axis_x, sep='')
  pcy <- paste('PC', axis_y, sep='')
  eigx <- eigenvalues[axis_x]
  eigy <- eigenvalues[axis_y]
  message(paste(pcx, eigx, pcy, eigy))
  # Check the values
  if (!(pcx %in% colnames(samples.df))){
    stop(paste('Error:', pcx, 'not in samples.df'))
  }
  if (!(pcy %in% colnames(samples.df))){
    stop(paste('Error:', pcy, 'not in samples.df'))
  }
  # Plot the scatter plot
  fig <- ggplot(samples.df, aes(x=samples.df[,pcx], y=samples.df[,pcy],
                                color=samples.df[,"Pop"])) +
    # Add origin lines
    geom_hline(yintercept=0, linetype="dashed", color = "gray") +
    geom_vline(xintercept=0, linetype="dashed", color = "gray") +
    # Add points
    geom_point(
      size=4,
      alpha=0.65) +
    scale_color_manual(values = c("#FFDB6D", "#00AFBB","#999999", "#2E4053",
                                  "#922B21", "#008000", "#0000FF","#CCCCFF", "#A36262")) +
    # Add titles
    labs(title='Clarias gariepinus PCA',
         x=paste('PC',axis_x,' ',format(round((eigx*100),2)),'%'),
         y=paste('PC',axis_y,' ',format(round((eigy*100),2)),'%')) +
    # Set Themes
    theme_light(base_size=16) +
    theme(panel.grid = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    # Centralize title
    theme(plot.title=element_text(hjust=0.5)) +
    # Legend
    labs(col="Pop ID")
  # Save plot
  f = paste(outdir, '/rad_pca.n', nrow(samples.df), '.pc', axis_x, '_pc', axis_y, '.pdf', sep='')
  ggsave(f, plot=fig, width=8, height=6)
  return(fig)
}


# Run code
genInd <- genepop_to_genind(genepop_f)
samples.df <- make_samples_tbl(genInd)

#Modify Pop name
samples.df["Pop"] <- c("f_CMC","f_CMC","f_CMC","f_CMC","f_CMC","f_CMC","f_CMC","f_CMC","f_CMC","f_CMC","f_CMC","f_ODC","f_ODC","f_ODC","f_ODC","f_ODC","f_ODC","f_ODC","f_ODC","f_ODC","f_ODC","f_ODC","f_ODC","f_LAC","f_LAC","f_LAC","f_LAC","f_LAC","f_LAC","f_LAC","f_LAC","f_LAC","f_LAC","f_LAC","f_LAC","f_SAC","f_SAC","f_SAC","f_SAC","f_SAC","f_SAC","f_SAC","f_SAC","f_SAC","f_SAC","f_SAC","w_BYC","w_BYC","w_BYC","w_BYC","w_BYC","w_BYC","w_BYC","w_BYC","w_LGC","w_LGC","w_LGC","w_LGC","w_LGC","w_LGC","w_LGC","w_LGC","w_LGC","w_LGC","w_LGC","w_LGC","w_DKAL","w_DKAL","w_DKAL","w_DKAL","w_DKAL","w_DKC","w_DKC","w_DKC","w_DKC","w_DKC","w_DKC","w_DKC","w_DKC","w_DKC","w_DKC","w_DKC","w_DKC","w_DKC","w_KDC","w_KDC","w_KDC","w_KDC","w_KDC","w_KDC","w_KDC","w_KDC","w_KDC","w_KDC")
# PCA
pca <- run_pca(genInd)
pca.eig_vals <- get_pca_eigenvalues(pca)
samples.df <- add_pcs_to_df(samples.df, pca)
pca.fig <- plot_pca(samples.df, pca.eig_vals, col, 1, 2, dir)
#pca.fig