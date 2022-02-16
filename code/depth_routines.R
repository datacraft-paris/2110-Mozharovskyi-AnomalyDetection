# Routines for curves' parametrization and box-plot colouring
# Author: Pavlo Mozharovskyi

# Function for re-parametrizing data by length
# All the data are parametrized on [0,1] with 'nArgs' ticks (uniformly)
paramLength <- function(fdata, nArgs){
  rez <- list()
  for (i in 1:length(fdata)){
    lengths <- c(0, cumsum(
      sqrt(rowSums((fdata[[i]]$vals[2:nrow(fdata[[i]]$vals),] - 
                      fdata[[i]]$vals[1:(nrow(fdata[[i]]$vals) - 1),])^2))))
    lengths <- lengths / lengths[length(lengths)]
    args <- 0:(nArgs - 1) / (nArgs - 1)
    rez[[i]] <- list(args = args, vals = matrix(0, nrow = nArgs, 
                                                ncol = ncol(fdata[[i]]$vals)))
    for (j in 1:ncol(fdata[[i]]$vals)){
      rez[[i]]$vals[,j] <- approx(lengths, fdata[[i]]$vals[,j], n = nArgs)$y
    }
  }
  return (rez)
}

# Function for re-parametrizing data by time/argument
# All the data are parametrized on [0,1] with 'nArgs' ticks (uniformly)
paramTime <- function(fdata, nArgs){
  rez <- list()
  for (i in 1:length(fdata)){
    times <- fdata[[i]]$args
    times <- times / times[length(times)]
    args <- 0:(nArgs - 1) / (nArgs - 1)
    rez[[i]] <- list(args = args, vals = matrix(0, nrow = nArgs, 
                                                ncol = ncol(fdata[[i]]$vals)))
    for (j in 1:ncol(fdata[[i]]$vals)){
      rez[[i]]$vals[,j] <- approx(times, fdata[[i]]$vals[,j], n = nArgs)$y
    }
  }
  return (rez)
}

# Return a vector of curve box-plot colors w.r.t. depths
getColorsBoxplot <- function(depths, portOutl = 0.2){
  #col1 <- c("red", "royalblue", "royalblue", "royalblue4")
  #col1 <- c("red", "royalblue", "royalblue", "royalblue4")
  alp1 <- c(0.5, 0.25, 0.75, 0.1)
  col1 <- c("order4", "order1", "order2", "order3")
  lty1 <- c(4, 3, 2, 1)
  ncolors1 <- length(depths)
  my_palette1 <- colorRampPalette(c("yellow", "red"))(n = ncolors1)
  depth.colors <- my_palette1
  depth.alphas <- rep(0, ncolors1)
  depth.ltys <- rep(0, ncolors1)
  depth.colors[rank(depths) <= floor(ncolors1 / 2)] <- col1[2]
  depth.colors[rank(depths) > floor(ncolors1 / 2)] <- col1[3]
  depth.colors[rank(depths) <= floor(ncolors1 * portOutl)] <- col1[1]
  depth.colors[which.max(depths)] <- col1[4]
  depth.alphas[rank(depths) <= floor(ncolors1 / 2)] <- alp1[2]
  depth.alphas[rank(depths) > floor(ncolors1 / 2)] <- alp1[3]
  depth.alphas[rank(depths) <= floor(ncolors1 * portOutl)] <- alp1[1]
  depth.alphas[which.max(depths)] <- alp1[4]
  depth.ltys[rank(depths) <= floor(ncolors1 / 2)] <- lty1[2]
  depth.ltys[rank(depths) > floor(ncolors1 / 2)] <- lty1[3]
  depth.ltys[rank(depths) <= floor(ncolors1 * portOutl)] <- lty1[1]
  depth.ltys[which.max(depths)] <- lty1[4]
  return (list(colors = depth.colors, alphas = depth.alphas, ltys = depth.ltys))
}
