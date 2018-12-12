### Created to find mean species from 36 rodent species using template from 3d-vs-ct-scanning project (Marcy et al. 2018). The mean species found with this script will be used to create a new template for this project characterizing rodent morphology.

library(geomorph)

# Read in data file and get it into format that geomorph recognizes
my_data <- read.csv(file.choose(), header = T) #choose the 18_01_30_AusRodentPhase1.csv file

# This file already has problem landmarks deleted, leaving 58 landmarks, 145 semilandmarks, and 86 patches for a total of 289 points

data <- my_data[, 2:(dim(my_data)[2])] #gets rid of LM names for geomorph formatting reasons
data<-t(data)
A <- arrayspecs(data, 289, 3) #convert to 3D array with 289 points and 3 dimensions

## Create tables of metadata
library(stringr)

# create CatNum and scan type metadata
foo <- my_data[0,]  # gets a list with photo file name
sp_names <- names(foo[-1])  # creates a list of species names
categories <- strsplit(sp_names, "_")  #separates important data in list
my_classifiers <- matrix(unlist(categories), ncol=2, byrow=T)  #reads list into matrix
colnames(my_classifiers) <- c("Genus", "Species")
sp_info <- as.data.frame(my_classifiers)  #converts to data frame so can index using $
sp_info$Both <- with(sp_info, interaction(Genus, Species))

# Create vectors of patch points and semi-landmarks
patches <- str_detect(u_pt_names, "PAT")
pat_num <- which(patches == TRUE)
sli_matrix <- read.csv(file.choose(), header = T) #choose smatrix.csv from methodology paper folder


## Analysis time

# Perform procrustes alignment
Y <- gpagen(A, Proj = TRUE, ProcD = TRUE, curves = sli_matrix, surfaces = pat_num)


# Create PCA
gp <- as.factor(sp_info$Genus)
col.gp <- rainbow(20)
names(col.gp) <- levels(gp)
col.gp <- col.gp[match(gp, names(col.gp))] #colors by Genus

pca <- plotTangentSpace(Y$coords, groups = col.gp, axis1 = 1, axis2 = 2, label = sp_info$Both, verbose = T)

#find mean specimen for template
findMeanSpec(Y$coords) 
