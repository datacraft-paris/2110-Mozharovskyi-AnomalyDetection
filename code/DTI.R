################################################################################
# File:             DTI.R
# Created by:       Pierre Lafaye De Micheaux
# First released:   06.11.2018
# Last revised:     06.11.2018
# 
# Contains R-functions to read DTI data.
################################################################################

paste.dir.file <- function(file, dir = "") {
  nchar.file <- nchar(file)
  if (substr(file, nchar.file, nchar.file) == "/") stop("argument 'file' should not contain a trailing slash!")
  dir <- path.expand(dir)
  if ((substr(dir, 1, 1)) == ".") dir <- paste(getwd(), substr(dir, 2, nchar(dir) - 1), sep = "")
  nchar.dir <- nchar(dir)
  if ((nchar.dir !=0) & (substr(dir, nchar.dir, nchar.dir) != "/")) dir <- paste(dir, "/", sep = "")
  dirname.file <- paste(dirname(file), "/", sep = "")
  if (dirname(file) != ".") { # 'file' seems to contain the full path
    if ((nchar.dir != 0) & (dirname.file != dir)) stop("'dir' argument is different from file.path(file)! ")
    return(file)
  }
  file <- basename(file)
  res <- paste(dir, file, sep = "")
  return(res)
}


# To obtain the coordinates of a voxel in the Standard Space from its coordinates in the Native Space ...
Native.to.Standard <- function(coords, fileTNa2St = "0099A_Na2St.mat", fileNa = "0099A_FA_Na.nii",
                               fileSt = "0099A_FA_St.nii", voxStart0 = FALSE, fileTNa2St.dir = "", fileNa.dir = "", fileSt.dir = "") {
# coords should be a 3 x N matrix, N being the number of points to consider.
  TNa2St <- as.matrix(read.table(paste.dir.file(fileTNa2St, fileTNa2St.dir)))
  Lfa <- f.read.header(paste.dir.file(fileNa, fileNa.dir)) # Native Space.
  LfaSt <- f.read.header(paste.dir.file(fileSt, fileSt.dir)) # Standard Space.
  Scale2 <- Scale1 <- diag(4)
  # To go from Voxel Space to Image Space (source):
  diag(Scale1)[1:3] <- Lfa$pixdim[2:4]
  # To go from Image Space to Voxel Space (destination):
  diag(Scale2)[1:3] <- 1 / LfaSt$pixdim[2:4]
  if (!voxStart0) coords <- coords - 1
  coords <- as.matrix(coords) ; coords <- rbind(coords, rep(1, ncol(coords)))
  # We may flip (e.g., FSL has a Radiological convention and others have Neurological):
  if (orientation(Lfa) == 1) coords[1, ] <- Lfa$dim[2] - 1 - coords[1, ]
  # The coordinates in the Standard Space are:
  res <- Scale2 %*% TNa2St %*% Scale1 %*% coords
  if (orientation(LfaSt) == 1) res[1, ] <- LfaSt$dim[2] - 1 - res[1, ]
  if (!voxStart0) res <- res + 1
return(res[1:3, ])
}

# ... and back.
Standard.to.Native <- function(coords, fileTNa2St = "0099A_Na2St.mat", fileNa = "0099A_FA_Na.nii",
                               fileSt = "0099A_FA_St.nii", voxStart0 = FALSE, fileTNa2St.dir = "", fileNa.dir = "", fileSt.dir = "") {
# coords should be a 3 x N matrix, N being the number of points to consider.
  TNa2St.inv <- solve(as.matrix(read.table(paste.dir.file(fileTNa2St, fileTNa2St.dir))))
  Lfa <- f.read.header(paste.dir.file(fileNa, fileNa.dir)) # Native Space.
  LfaSt <- f.read.header(paste.dir.file(fileSt, fileSt.dir)) # Standard Space.
  Scale2 <- Scale1 <- diag(4)
  # To go from Image Space to Voxel Space (Native):
  diag(Scale1)[1:3] <- 1 / Lfa$pixdim[2:4]
  # To go from Voxel Space to Image Space (Standard):
  diag(Scale2)[1:3] <- LfaSt$pixdim[2:4]
  if (!voxStart0) coords <- coords - 1
  coords <- as.matrix(coords) ; coords <- rbind(coords, rep(1, ncol(coords)))
  # The coordinates in the Standard Space are:
  if (orientation(LfaSt) == 1) coords[1, ] <- LfaSt$dim[2] - 1 - coords[1, ]
  res <- Scale1 %*% TNa2St.inv %*% Scale2%*% coords
  # We may flip (e.g., FSL has a Radiological convention and others have Neurological):
  if (orientation(Lfa) == 1) res[1, ] <- Lfa$dim[2] - 1 - res[1, ] 
  if (!voxStart0) res <- res + 1
return(res[1:3, ])
}

# To read files similar to: 0099Avox_left.txt
# WARNING! : voxel indices start at 0 in files output by MRtrix software
f.read.fibers <- function (fibers.file = "0099Avox_left.txt", fibers.dir = "") {
  fibers.file <- paste.dir.file(fibers.file, fibers.dir)
   nbpts.total <- as.numeric(scan(file = fibers.file, what = "", nlines = 1, skip = 5)[2])
   fibers.coordinates <- matrix(scan(file = fibers.file, nlines = nbpts.total, skip = 6),
                                ncol = 3, byrow = TRUE)
   nbfibers <- as.numeric(scan(file = fibers.file, what = "", nlines = 1, 
                               skip = 7 + nbpts.total)[2])
   tmp <- scan(file = fibers.file, skip = 8 + nbpts.total, blank.lines.skip = FALSE)
   tmp <- tmp[-length(tmp)]
   vect.of.nbpts.per.fiber <- c(tmp[1], tmp[which(is.na(tmp)) + 1])
   if (length(vect.of.nbpts.per.fiber) != nbfibers) 
               stop("There is a problem in function f.read.fibers()!")
   # Let's call M the matrix 'indices.of.fibers'.
   # M[n, 1] gives the row index (in matrix 'fibers.coordinates') of 
   #         the lower extremity of the nth fiber (in the 2mm region);
   # M[n, 2] gives the row index (in matrix 'fibers.coordinates') of 
   #         the other extremity of the nth fiber:
   indices.of.fibers <- cbind(cumsum(c(1, vect.of.nbpts.per.fiber[-nbfibers])),
                           cumsum(vect.of.nbpts.per.fiber))
 return(list(nbfibers = nbfibers, vect.of.nbpts.per.fiber = vect.of.nbpts.per.fiber,
             fibers.coordinates = fibers.coordinates + 1, indices.of.fibers = indices.of.fibers))# We add +1 because fibers.coordinates are given with voxel indices starting at 0.
}

# fibers.coordinates is a Nx3 matrix that contains the coordinates of all
#   the fibers of the subject considered. Overall, there are N coordinates.

# indices.of.fibers is a nx2 matrix that contains the indices of the n fibers considered 
#     for a given subject. (See above for a description of this matrix M).


# So for example, the coordinates of the nth fiber are given by:
f.coords.one.fiber <- function(n, fibers.coordinates, indices.of.fibers) {
 return(fibers.coordinates[indices.of.fibers[n, 1]:indices.of.fibers[n, 2], ])
}


f.fibers.node <- function(fibers.coordinates, indices.of.fibers, parcellated.volume, lower = FALSE) {
# Returns a vector of the index node associated to each one of the n fibers.
 extrem <- if (lower) 1 else 2
 coord.extr <- fibers.coordinates[indices.of.fibers[, extrem], , drop = FALSE] 
 res <- parcellated.volume[ceiling(cbind(coord.extr, rep(1, nrow(coord.extr))))]
return(res)
}

f.portions.lengths.one.fiber <- function(n, fibers.coordinates, indices.of.fibers) {
 # Computes the d_h values (h = 1, ..., n_f - 1) for the nth fiber.
 mat <- f.coords.one.fiber(n,fibers.coordinates, indices.of.fibers)
 nf <- nrow(mat)
 res <- apply((mat[1:(nf - 1), ] - mat[2:nf, ]) ^ 2, FUN = function(x) {sqrt(sum(x))}, MARGIN = 1)
return(res)
}

f.portions.lengths.fibers <- function(fibers.coordinates, indices.of.fibers) {
 # Computes the d_h values (h = 1, ..., n_f - 1) for the all fibers.
 # Returns a list; nth element of the list is: f.portions.lengths.one.fiber(n, fibers.coordinates, indices.of.fibers)
 nf <- nrow(fibers.coordinates)
 tmp <- apply((fibers.coordinates[1:(nf - 1), ] - fibers.coordinates[2:nf, ]) ^ 2, FUN = function(x) {sqrt(sum(x))}, MARGIN = 1)
 indices.of.fibers[, 2] <- indices.of.fibers[, 2] - 1
 res <- apply(indices.of.fibers, MARGIN = 1, FUN = function(x){tmp[":"(x[1], x[2])]})
 return(res)
}



f.length.one.fiber <- function(n, fibers.coordinates, indices.of.fibers) {
 # Computes the length L_{i,j} given in the mathematical formula.
 res <- sum(f.portions.lengths.one.fiber(n, fibers.coordinates, indices.of.fibers))
return(res)
}


f.length.fibers <- function(fibers.coordinates, indices.of.fibers) {
 # Computes the length L_{i,j} given in the mathematical formula.
 # res[n] is: f.length.one.fiber(n, fibers.coordinates, indices.of.fibers)
  res <- sapply(f.portions.lengths.fibers(fibers.coordinates, indices.of.fibers), FUN = sum)
return(res)
}


dicho <- function(n, vec.l, portions.lengths.fiber) {
# Returns the index i such that lL is between c(0,tmp)[i] and c(0,tmp)[i + 1].
  small.func <- function(a, b, lL) {
    if (lL < tmp[1]) return(1)
    while (a < b - 1) {
      center <- as.integer((a + b) / 2)
      if (tmp[center] < lL) a <- (a + b) / 2 else b <- (a + b) / 2
     }
  return(center)
  }
 tmp <- cumsum(portions.lengths.fiber)
 nf <- length(tmp) + 1
 L <- tmp[nf - 1]
 vec.lL <- vec.l * L # We search lL.
 i <- rep(NA, length(vec.l))
 for (l in 1:length(vec.l)) {
  a <- 0
  b <- nf
  center <- small.func(a, b, vec.lL[l])
  i[l] <- center
 }
return(list(i = i, L = L)) # lL is between these two points.
}

# A few intermediary functions:
# http://math.stackexchange.com/questions/293116/rotating-one-3-vector-to-another
cross.prod <- function(u, v) {
# Cross product between two vectors (return the vector u x v).
# Warning: do not confuse with crossprod() in base R.
 return(c(u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3], u[1] * v[2] - u[2] * v[1]))
}

norm <- function(u) {
# Norm of vector u.
 return(sqrt(sum(u ^ 2)))
}

dot.prod <- function(u, v) {
# Returns the dot product u.v:
 return(sum(u * v))
}

rotate <- function(a, b, eps = 10 ^ -6) {
# We rotate vector a onto vector b.
 theta <- acos(dot.prod(a, b) / (norm(a) * norm(b)))
 if (is.nan(theta)) theta <- 0
 if ((pi - theta) < eps) {
  ei <- c(0, 0, 0)
  ei[which.min(abs(a))] <- 1
  x <- cross.prod(a, ei)
 } else x <- cross.prod(a, b)
 x <- x / norm(x)
 A <- matrix(c(0, -x[3], x[2], x[3], 0, -x[1], -x[2], x[1], 0), nrow = 3, ncol = 3, byrow = TRUE)
 if (theta < eps) R <- diag(3) else R <- diag(3) + sin(theta) * A + (1 - cos(theta)) * A %*% A
# The vector Ra should have same direction as b.
 return(as.vector(R %*% as.matrix(a)))
}

coord.interp <- function(n, vec.l, fibers.coordinates, indices.of.fibers) {
# Coordinates (interpolated between the two nearby points) of the point 
# at length lL along the nth fiber.
# Note: l should be in [0, 1].
 portions.lengths.fiber <- f.portions.lengths.one.fiber(n, fibers.coordinates, indices.of.fibers)
 tmp <- cumsum(portions.lengths.fiber)
 res.dicho <- dicho(n, vec.l, portions.lengths.fiber)
 vec.lL <- vec.l * res.dicho$L
 vec.i <- res.dicho$i
# Point at length lL along the fiber is at a distance of beta (on the right) 
# from the point at length tmp[i]:
 coords <- f.coords.one.fiber(n, fibers.coordinates, indices.of.fibers)
 res <- matrix(NA, nrow = 3, ncol = length(vec.l))
 for (l in 1:length(vec.l)) {
  beta <- vec.lL[l] - c(0,tmp)[vec.i[l]]
  res[, l] <-  coords[vec.i[l], ] + rotate(c(beta, 0, 0), coords[vec.i[l] + 1, ])
 }
return(res)
}

f.fa.values.unit.segment <- function(nfa, FA.native.file, fibers.coordinates, indices.of.fibers, FA.native.dir = "") {
  FA.native.file <- paste.dir.file(FA.native.file, FA.native.dir)
  vec.l <- seq(from = 0, to = 1, length = nfa)
  vol <- f.read.volume(FA.native.file)
  nbfibers <- nrow(indices.of.fibers)
  res <- matrix(NA, nrow = nbfibers, ncol = nfa)
  for (n in 1:nbfibers) {
    tmp <- coord.interp(n, vec.l, fibers.coordinates, indices.of.fibers)
    tmp <- t(round(tmp))
    res[n, ] <- vol[cbind(tmp, rep(1, nrow(tmp)))]
  }
  return(res)
}

f.write.FA.unit.segments <- function(nfa = 100,
                                     file.write = FALSE,
                                     fibers.file,
                                     FA.native.file,
                                     native.parcellated.dilated.destination.region.mask.file,
                                     native.source.region.mask.file,
                                     file.out = "FA_unit_segment.txt",
                                     dir.out = "Results/",
                                     fibers.dir = "",
                                     FA.native.dir = "",
                                     native.parcellated.dilated.destination.region.mask.dir = "",
                                     native.source.region.mask.dir = ""
                                     ) {
  
  res <- f.read.fibers(fibers.file, fibers.dir)

# because some fibers coordinates given by MRtrix were not correct (were outside the array), I had to do this:
  hdr <- f.read.header(paste.dir.file(FA.native.file, FA.native.dir))$dim[2:4]
  res$fibers.coordinates[, 1][res$fibers.coordinates[, 1] > hdr[1]] <- hdr[1]
  res$fibers.coordinates[, 2][res$fibers.coordinates[, 2] > hdr[2]] <- hdr[2]
  res$fibers.coordinates[, 3][res$fibers.coordinates[, 3] > hdr[3]] <- hdr[3]

  res$fibers.coordinates[, 1][res$fibers.coordinates[, 1] < 1] <- 1
  res$fibers.coordinates[, 2][res$fibers.coordinates[, 2] < 1] <- 1
  res$fibers.coordinates[, 3][res$fibers.coordinates[, 3] < 1] <- 1

  val <- fa.values.unit.segment <- f.fa.values.unit.segment(nfa, FA.native.file,res$fibers.coordinates, res$indices.of.fibers, FA.native.dir)
  native.parcellated.dilated.destination.region.mask.file <- paste.dir.file(native.parcellated.dilated.destination.region.mask.file,
                                                                            native.parcellated.dilated.destination.region.mask.dir)
  native.source.region.mask.file <- paste.dir.file(native.source.region.mask.file,native.source.region.mask.dir)
  parcellated.volume <- f.read.nifti.volume(native.parcellated.dilated.destination.region.mask.file)
  fibers.node.destination <- f.fibers.node(res$fibers.coordinates, res$indices.of.fibers, parcellated.volume, lower = FALSE)
  parcellated.volume <- f.read.nifti.volume(native.source.region.mask.file)
  fibers.node.source <- f.fibers.node(res$fibers.coordinates, res$indices.of.fibers, parcellated.volume, lower = TRUE)
  lengths <- rep(NA, length(fibers.node.destination))
  for (n in 1:length(fibers.node.destination)) lengths[n] <- f.length.one.fiber(n, res$fibers.coordinates, res$indices.of.fibers)
  val <- cbind(fibers.node.source, fibers.node.destination, lengths, val)
  if (file.write) {
    file.out <- paste.dir.file(file.out, dir.out)
    write.table(val, file.out, col.names = FALSE)
  }
  else return(val)
}


# Before using the function below, issue a:
# setwd("~/HeritabilityMotor/Data/NativeSpace/FAs/")
# subjectIDs <- unlist(strsplit(dir(),"_FA_Na.nii",fixed=TRUE))
# setwd("/home/pierre/HeritabilityMotor") 
# source("/home/pierre/HeritabilityMotor/Tools/DTI.R")
# f.compute.all.subjects.FA(subjectIDs=subjectIDs)
f.compute.all.subjects.FA <- function(K = 3, # 2^K nodes in the parcellation
                                      sphere.source = 6,  # size in mm
                                      sphere.destination = 4,  # size in mm
                                      nfa = 100,  # Number of points in the [0,1] segment
                                      subjectIDs = c("0099A"),
                                      standard.space.dir = "~/HeritabilityMotor/Data/StandardSpace",
                                      standard.space.file = "JHU-ICBM-FA-2mm.nii",
                                      standard.source.region.mask.dir = "~/HeritabilityMotor/Data/StandardSpace/Masks/", # Nifti mask dir of source region in Standard Space
                                      source.region.mask = "cst_seed",
                                      standard.destination.region.mask.dir = "~/HeritabilityMotor/Data/StandardSpace/Masks/", # Nifti mask dir of destination region in Standard Space
                                      destination.region.mask = "cst_cortex",
                                      native.source.region.mask.dir = "~/HeritabilityMotor/Data/NativeSpace/Masks/", # Nifti mask dir of source region in Native Space
                                      native.destination.region.mask.dir = "~/HeritabilityMotor/Data/NativeSpace/Masks/", # Nifti mask dir of destination region in Native Space
                                      FA.native.dir = "~/HeritabilityMotor/Data/NativeSpace/FAs/",
                                      Transform.Matrices.dir = "~/HeritabilityMotor/Data/AffineTransf/",
                                      Results.dir = "~/HeritabilityMotor/Results/",
                                      fibers.dir = "~/HeritabilityMotor/Data/NativeSpace/Fibers/",
                                      native.parcellated.dilated.destination.region.mask.dir = "~/HeritabilityMotor/Data/NativeSpace/Masks/",
                                      uniform.parcellate.dir = "~/HeritabilityMotor/Tools/",
                                      uniform.parcellate.tool = "uniform_parcellate.m",
                                      hemispheres = c("left", "right"),
				      process.mask = TRUE	
                                      ) {

    
# To be done only once for all subjects:
########################################
  if (process.mask) {

    for (hemisphere in hemispheres) {
# Dilation of Source mask in Standard Space:
      system(paste("fslmaths ", paste.dir.file(paste(source.region.mask, "_", hemisphere, ".nii", sep = ""),
                                               standard.source.region.mask.dir), " -kernel sphere ", sphere.source, " -dilM ",
                   paste.dir.file(paste("dilated_", source.region.mask, "_", hemisphere, ".nii", sep = ""),
                                  standard.source.region.mask.dir), sep = ""))
      system(paste("gunzip -f ", paste.dir.file(paste("dilated_", source.region.mask, "_", hemisphere, ".nii", ".gz", sep = ""),
                                                standard.source.region.mask.dir)))
    
# Dilation of Destination mask in Standard Space:
      system(paste("fslmaths ", paste.dir.file(paste(destination.region.mask, "_", hemisphere, ".nii", sep = ""),
                                               standard.destination.region.mask.dir), " -kernel sphere ", sphere.destination,
                   " -dilM ", paste.dir.file(paste("dilated_", destination.region.mask, "_", hemisphere, ".nii", sep = ""),
                                             standard.destination.region.mask.dir), sep = ""))
      system(paste("gunzip -f ", paste.dir.file(paste("dilated_", destination.region.mask, "_", hemisphere, ".nii", ".gz", sep = ""),
                                                standard.destination.region.mask.dir)))
    
# Parcellation of Destination Dilated mask in Standard Space:
      current.dir <- getwd()
      setwd(uniform.parcellate.dir)
      system(paste("./matlab_batcher.sh ", uniform.parcellate.tool, " \"K=", K, "\" \"Msk='",
                   paste.dir.file(paste("dilated_", destination.region.mask, "_", hemisphere, ".nii", sep = ""),
                                  standard.destination.region.mask.dir), "'\" \"Out='",
                   paste.dir.file(paste("parcellated_dilated_", destination.region.mask, "_", hemisphere, ".nii", sep = ""),
                                  standard.destination.region.mask.dir), "'\"", sep = ""))
      setwd(current.dir)
    }
  }

# Conmputation of the results for all subjects:
###############################################
  require(AnalyzeFMRI)
  for (subjectID in subjectIDs) {
    cat("\n"); cat("Processing subject "); cat(subjectID);cat("\n")
    ### Computing subj_Na2St.mat and subj_St2Na.mat in AffineTransf
    system(paste("flirt -in ", paste.dir.file(paste(subjectID, "_FA_Na.nii", sep = ""), FA.native.dir), " -ref ",
                 paste.dir.file(standard.space.file,standard.space.dir), " -omat ",
                 paste.dir.file(paste(subjectID, "_Na2St.mat", sep = ""), Transform.Matrices.dir),
                 " -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear", sep = ""))
    system(paste("convert_xfm -omat ", paste.dir.file(paste(subjectID, "_St2Na.mat", sep = ""), Transform.Matrices.dir),
                 " -inverse ", paste.dir.file(paste(subjectID, "_Na2St.mat", sep = ""), Transform.Matrices.dir), sep = ""))

    for (hemisphere in hemispheres) {
      
# Transformation of Source (2mm) Dilated mask from Standard to Native Space:
      system(paste("flirt -in ", paste.dir.file(paste("dilated_", source.region.mask, "_", hemisphere, ".nii", sep = ""),
                                                standard.source.region.mask.dir), " -ref ",
                   paste.dir.file(paste(subjectID, "_FA_Na.nii", sep = ""), FA.native.dir),
                   " -applyxfm -init ", paste.dir.file(paste(subjectID, "_St2Na.mat", sep = ""),
                                                       Transform.Matrices.dir), " -out ",
                   paste.dir.file(paste(subjectID, "_dilated_", source.region.mask, "_", hemisphere, ".nii", sep = ""),
                                  native.source.region.mask.dir), " -paddingsize 0.0 -interp nearestneighbour", sep = ""))
      system(paste("gunzip -f ", paste.dir.file(paste(subjectID, "_dilated_", source.region.mask, "_", hemisphere, ".nii", ".gz", sep = ""),
                                                native.source.region.mask.dir), sep = ""))

# Transformation of Destination Dilated Parcellated mask from Standard to Native Space:    
      system(paste("flirt -in ", paste.dir.file(paste("parcellated_dilated_", destination.region.mask, "_", hemisphere, ".nii", sep = ""),
                                                standard.destination.region.mask.dir), " -ref ",
                   paste.dir.file(paste(subjectID, "_FA_Na.nii", sep = ""), FA.native.dir), " -applyxfm -init ",
                   paste.dir.file(paste(subjectID, "_St2Na.mat", sep = ""), Transform.Matrices.dir), " -out ",
                   paste.dir.file(paste(subjectID, "_parcellated_dilated_", destination.region.mask, "_", hemisphere, ".nii", sep = ""),
                                  native.destination.region.mask.dir), " -paddingsize 0.0 -interp nearestneighbour", sep = ""))
      system(paste("gunzip -f ", paste.dir.file(paste(subjectID, "_parcellated_dilated_", destination.region.mask, "_", hemisphere,
                                                      ".nii", ".gz", sep = ""), native.destination.region.mask.dir), sep = ""))

# Compute (and write) the ****_FA_unit_segment.txt files
      f.write.FA.unit.segments(nfa,
                               file.write = TRUE,
                               fibers.file = paste(subjectID, "vox_", hemisphere, ".txt", sep = ""),
                               FA.native.file = paste(subjectID, "_FA_Na.nii", sep = ""),
                               native.parcellated.dilated.destination.region.mask.file = paste(subjectID, "_parcellated_dilated_", destination.region.mask, "_", hemisphere, ".nii", sep = ""),
                               native.source.region.mask.file=paste(subjectID, "_dilated_", source.region.mask, "_", hemisphere, ".nii", sep = ""),
                               file.out = paste(subjectID, "_FA_unit_segment_", hemisphere, ".txt", sep = ""),
                               dir.out = Results.dir,
                               fibers.dir = fibers.dir,
                               FA.native.dir = FA.native.dir,
                               native.parcellated.dilated.destination.region.mask.dir = native.parcellated.dilated.destination.region.mask.dir,
                               native.source.region.mask.dir = native.source.region.mask.dir
                               )
    }
  }
}


