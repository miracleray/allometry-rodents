# Utility Functions
# Contains functions used throughout the .Rmd scripts
# Loaded at the top of each script

##########################
# WriteMetadata
##########################
WriteMetadata <- function(threeD.array, cols) {
        # Makes metadata table from specimen filenames for shape coordinates.
        #
        # Args:
        #    threeD.array: 3D array (p x k x n), where p is the number of landmarks, k is the dimension (3), and n is the number of specimens. Assumes the 1st column of landmark names has been removed. 
        #    cols: a character vector of column names of length n-1, where n is the number of underscores separating metadata information in the specimen filenames. Assumes filenames contain information in the same order and the appropriate names are given in this order in cols. 
        #
        # Returns: 
        #    A dataframe containing the metadata for each specimen in the same order as specimens in the 3D array of shape data.
        
        # Remove 'ind' that bilat.symmetry() appends to specimen names
        names <- gsub("ind", "", dimnames(threeD.array)[[3]])
        
        # Convert name vectors into data frame
        categories <- strsplit(names, "_") 
        my.classifiers <- matrix(unlist(categories), ncol = length(cols), byrow = T) 
        colnames(my.classifiers) <- cols
        sp.info <- as.data.frame(my.classifiers)
        
        return(sp.info)
}

##########################
# StartSlider
##########################
StartSlider <- function(cur.num){
        # Creates a _geomorph_ slider matrix with correct After and Before landmark numbers for semi-landmarks in the middle of their curve. 
        #
        # Args:
        #    cur.num: an integer vector of landmarks categorized as curve semi-landmarks
        #
        # Returns:
        #   A matrix ready for geomorph's gpagen() curves = argument. It will be missing appropriate Before and After numbers for fixed landmarks bounding every curve.
        Before <- (cur.num - 1)
        Slider <- cur.num
        After <- (cur.num + 1)
        matrix <- cbind(Before, Slider, After)
        return(matrix)
}

##########################
# FindLMnumbers
##########################
FindLMnumbers <- function(i, index, pt.names, key, matrix){
        # Finds row numbers in shape data for 2 named LMs bounding a curve. 
        #
        # Args:
        #    i: an integer corresponding to the row number of a curve name in the key argument fed to `FinishSlider()`
        #    index: an integer vector of curve points row numbers
        #    pt.names: a list of point names in the landmarking protocol
        #    matrix: a slider matrix with 3 columns to be edited, e.g. one created by the function `StartSlider()`
        #
        # Returns:
        #    A slider matrix with correct landmark numbers for LMs bounding one curve, given by row number, i. 
        
        bef.LM <- as.character(key$Before[i])  # Before landmark name
        
        # Checks if this landmark is bilateral and if so, appends the appropriate "R" or "L" so the correct landmark number is retrieved later
        if ((sum(str_detect(pt.names, bef.LM)) > 1)) {
                side <- str_sub(pt.names[min(index)], -3, -3)
                bef.LM <- paste(bef.LM, side)
        } 
        # Finds & inserts landmark number in slider matrix Before column
        LM.num <- which(str_detect(pt.names, bef.LM))
        matrix[min(index), 1] <- LM.num
        
        aft.LM <- as.character(key$After[i])  # After landmark name
        
        # Checks if this landmark is bilateral
        if (sum(str_detect(pt.names, aft.LM)) > 1) {
                side <- str_sub(pt.names[min(index)], -3, -3)
                aft.LM <- paste(aft.LM, side)
        } 
        # Finds & inserts landmark number in After column
        LM.num <- which(str_detect(pt.names, aft.LM))
        matrix[max(index), 3] <- LM.num
        
        return(matrix)
}

##########################
# FinishSlider
##########################
FinishSlider <- function(pt.names, key, matrix){
        # Corrects numbers for bounding LMs for all curves in a sliding matrix intended for gpagen()'s argument curves. 
        #
        # Args:
        #    pt.names: a character list of point names
        #    key: the user-generated matrix with column "Curve" listing unique curve names, "Before" and "After" listing bounding landmark names
        #    matrix: a 3-column matrix initiated by `StartSlider()`.
        #
        # Returns
        #    A completed slider matrix ready for geomorph's gpagen()
        
        # Cycles through unique curve names in the template
        for (i in 1:length(key$Curve)){  # i passed to Give_LMs()
                curve <- as.character(key$Curve[i]) 
                
                # Checks to make sure more than 1 sliding pt in curve
                # If so, gets the index numbers included in curve
                # Also extracts letter which might indicate L or R
                if (sum(str_detect(pt.names, curve) != 0)){
                        index <- which(str_detect(pt.names, curve))
                        test.LR <- str_sub(pt.names[min(index)], -3, -3)
                        
                        # Tests for and handles center-line curves
                        if (test.LR != "L" & test.LR != "R"){
                                matrix <- FindLMnumbers(i, index, pt.names, key, matrix)
                        } else { # Handles bilateral curves, expecting the coordinates in reverse alphabetical order that lists the R side points immediately before the L side points
                                midpt <- length(index)/2
                                R.index <- index[1:midpt]
                                L.index <- index[(midpt + 1):(midpt*2)]
                                matrix <- FindLMnumbers(i, R.index, pt.names, key, matrix)
                                matrix <- FindLMnumbers(i, L.index, pt.names, key, matrix)
                        }
                }
        }
        return(matrix)
}

##########################
# FindPairs
##########################
FindPairs <- function(pt.names){
        # Creates table of paired bilateral landmarks for bilat.symmetry().
        #
        # Args:
        #   pt.names: a character vector of landmark names.
        #
        # Returns:
        #   2 column data table of paired landmarks ready for geomorph's bilat.symmetry()'s land.pair argument.
        
        pairs <- data.table("Right" = numeric(), "Left" = numeric())
        
        # Removes R and L designations so pairs can be detected
        no.side.names <- gsub(" R", "", pt.names)
        no.side.names <- gsub(" L", "", no.side.names)
        
        # Checks if point has a pair and if so, their index #s are paired
        for(i in unique(no.side.names)){
                index <- str_detect(no.side.names, i)
                if (sum(index) == 2) {
                        new.pair <- which(index == TRUE)
                        pairs <- rbind(pairs, as.list(new.pair))
                } 
        }
        return(pairs)
}

##########################
# RepAbility
##########################
RepAbility <- function(coords, ids, n.Rep, print = TRUE, export = FALSE, filename = NULL) {
        # Calculates repeatability (R) for GMM studies.
        #
        # Args:
        #    coords: a 3D array (p x k X n) of shape coordinates.
        #    ids: a list of identifiers used to find replicates, e.g. CatNum.
        #    n.Rep: number of repetitions taken for each individual
        #    print: if TRUE, prints ANOVA and R to the console.
        #    export: if TRUE, exports ANOVA and R to .csv.
        #    filename: the filename used to save the .csv file.
        #
        # Returns:
        #    A table with the ANOVA and the value of R, repeatability.
        
        # Calculations from formulas 1-3 in Fruciano 2016
        r.gdf <- geomorph.data.frame(coords = coords, ind = factor(ids))
        rep.er <- procD.lm(coords ~ ind, data = r.gdf, iter = 999)
        S.sq.A <- ((rep.er$aov.table$MS[1] - rep.er$aov.table$MS[2]) / n.Rep)  # among-individuals variance component: 
        S.sq.W <- rep.er$aov.table$MS[2]  # within-individual variance
        R <- S.sq.A / (S.sq.W + S.sq.A)  # analogue of the intraclass correlation coeffiecent
        
        table <- rep.er$aov.table
        table$Repeatability <- R
        
        if (print) {
                print(rep.er$aov.table)
                cat("\n","Repeatability =", R)
        }
        if (export) {
                write.csv(table, file = paste(filename, ".csv", sep = ""))
        } else {
                return (table)
        }
}

##########################
# PlotByGroup
##########################
PlotByGroup <- function(metadata, column, color.key){
        # Matches colors or other plotting attributes to each specimen according to a grouping factor or a column number.
        #
        # Args:
        #   metadata: metadata table, often created with WriteMetadata().
        #   column: 1 string matching the column name with target groups.
        #   color.key: a vector of attributes listed in the same order as the unique group descriptors given by levels(as.factor(metadata$column))
        # 
        # Returns:
        #    A vector of colors the length of specimen number, with colors according to the group descriptor of each individual, ready for plot().
        
        if (is.numeric(column)) {
                col.num <- column
        } else {
                col.names <- unlist(dimnames(metadata)[[2]])
                col.num <- which(col.names == column)
        }
        
        grp <- as.factor(metadata[, col.num])
        names(color.key) <- sort(unique(grp))
        grp.col <- color.key[match(grp, names(color.key))]
        return(grp.col)
}

##########################
# MatchSpecShape
##########################
MatchSpecShape <- function(spec, info, shape){
        # Matches a specimen name to its 3D shape information.
        #
        # Args:
        #    spec: the specimen dimname provided by plotOutliers(), usually the filename of the specimen provided while landmarking.
        #    info: the metadata table which contains a column named CatNum.
        #    shape: 3D shape array in (p x k x n), where n = specimen number.
        #
        # Returns:
        #    The shape data for 1 specimen of name "spec" in format landmarks x coordinates (p x k).
        
        categories <- unlist(strsplit(names(spec), "_"))
        catnum <- categories[str_which(categories, "[A-Z][0-9]")]  # detects CatNum
        spec.shape <- shape[, , which(info$CatNum == catnum)]
        return(spec.shape)
}

##########################
# PlotPCA
##########################
PlotPCA <- function(pca, PCx, PCy, col.grp, pch.grp = 16, flip.axis1 = F, flip.axis2 = F) {
        # Plots PCAs with specimens colored (and optionally given different points) according to groups, reports PC axis variation in %; optional axis flipping.
        # 
        # Args:
        #    PCA: an object of class plotTangentSpace
        #    PCx: the PC intended for the x-axis. Default is PC1.
        #    PCy: the PC intended for the y-axis, usually x > y. Default is PC2.
        #    col.grp: a vector of colors ordered in the same way as specimens, usually made with PlotByGroup().
        #    pch.grp: an optional vector for point shapes, also usually made with PlotByGroup() function. Default is a filled circle.
        #    flip.axis1: If TRUE, reverses sign for all coordinates of PCx
        #    flip.axis2: If TRUE, reverses sign for all coordinates of PCy
        #
        # Returns:
        #    If return.PCA is TRUE, returns the pca object from plotTangentSpace(). If FALSE, returns a plot coloring the PCA by groups specified by col.grp and optionally pch.grp. 
        
        # Handle flipped axes, if there are any
        if (flip.axis1 == TRUE) {
                pca$pc.scores[, PCx] <- -(pca$pc.scores[, PCx])
        }
        
        if (flip.axis2 == TRUE) {
                pca$pc.scores[, PCy] <- -(pca$pc.scores[, PCy])
        }
        
        # Write x and y labels with proportion of variance for PCx and PCy
        PCs <- pca$pc.summary$importance
        PCx.per <- round(PCs[2, PCx] * 100, digits = 1)  # % with 1 decimal
        PCx.lab <- paste("PC", PCx, " (", PCx.per, "%)", sep = "")
        PCy.per <- round(PCs[2, PCy] * 100, digits = 1)
        PCy.lab <- paste("PC", PCy, " (", PCy.per, "%)", sep = "")
        
        PCA.plot <- plot(x = pca$pc.scores[, PCx],
                         y = pca$pc.scores[, PCy], 
                         xlab = PCx.lab, 
                         ylab = PCy.lab,
                         asp = TRUE,
                         col = col.grp, 
                         pch = pch.grp, 
                         bg = col.grp,
                         cex = 1.5,
                         cex.axis = 1.3, 
                         cex.lab = 1.3)
}

##########################
# DoPlotPCA
##########################
DoPlotPCA <- function(shape, PCx, PCy, col.grp, pch.grp = 16, return.PCA = F, flip.axis1 = F, flip.axis2 = F) {
        # Runs and plots PCAs with specimens colored (and optionally given different points) according to groups, reports PC axis variation in %; optional axis flipping.
        # 
        # Args:
        #    shape: a 3D array of shape coordinates in (p x k x n) format
        #    PCx: the PC intended for the x-axis.
        #    PCy: the PC intended for the y-axis, usually x > y. 
        #    col.grp: a vector of colors ordered in the same way as specimens, usually made with PlotByGroup().
        #    pch.grp: an optional vector for point shapes, also usually made with PlotByGroup() function. Default is a filled circle.
        #    return.PCA: If TRUE, returns the PCA data (run with groups set to col.grp) without a fancy plot. Default is FALSE.
        #    flip.axis1: If TRUE, reverses sign for all coordinates of PCx
        #    flip.axis2: If TRUE, reverses sign for all coordinates of PCy
        #
        # Returns:
        #    If return.PCA is TRUE, returns the pca object from plotTangentSpace(). If FALSE, returns a plot coloring the PCA by groups specified by col.grp and optionally pch.grp. 
        
        pca <- plotTangentSpace(shape, groups = col.grp, axis1 = PCx, axis2 = PCy, verbose = T)
        
        if (return.PCA == TRUE) {
                return(pca)
        }
        
        # Handle flipped axes, if there are any
        if (flip.axis1 == TRUE) {
                pca$pc.scores[, PCx] <- -(pca$pc.scores[, PCx])
        }
        
        if (flip.axis2 == TRUE) {
                pca$pc.scores[, PCy] <- -(pca$pc.scores[, PCy])
        }
        
        # Write x and y labels with proportion of variance for PCx and PCy
        PCs <- pca$pc.summary$importance
        PCx.per <- round(PCs[2, PCx] * 100, digits = 1)  # % with 1 decimal
        PCx.lab <- paste("PC", PCx, " (", PCx.per, "%)", sep = "")
        PCy.per <- round(PCs[2, PCy] * 100, digits = 1)
        PCy.lab <- paste("PC", PCy, " (", PCy.per, "%)", sep = "")
        
        PCA.plot <- plot(x = pca$pc.scores[, PCx],
                         y = pca$pc.scores[, PCy], 
                         xlab = PCx.lab, 
                         ylab = PCy.lab,
                         asp = TRUE,
                         col = col.grp, 
                         pch = pch.grp, 
                         bg = col.grp,
                         cex = 1.5,
                         cex.axis = 1.3, 
                         cex.lab = 1.3)
}

##########################
# FoundInRegion
##########################
FoundInRegion <- function(spec.info, regions, inc.partial = FALSE) {
        # Returns a logical vector of which specimens occur in the region(s) of interest, in same order as spec.info table.
        #
        # Args:
        #    spec.info: table of specimen information which codes for presence (1), absence (0), and partial presence (0.5) in 7 columns for the regions listed below, in the order listed below. 
        #    regions: a list of numbers corresponding to the 7 regions defined by Breed & Ford 2006. Acceptable numbers are:
        #          1 = "Savannah"
        #          2 = "AridZone"
        #          3 = "NEWetForest"
        #          4 = "NEDryForest"
        #          5 = "Pilbara"
        #          6 = "SW"
        #          7 = "SE"
        #    inc.partial: if TRUE, includes species who occur in regions less frequently or only under certain conditions. Default is FALSE.
        #
        # Returns:
        #    A logical vector of which specimens occur in the region(s) of interest, in same order as spec.info table.
        
        n.spec <- dim(spec.info)[1]
        is.there <- vector(mode = "logical", length = n.spec)
        first.region <- which(colnames(spec.info) == "Savannah")
        
        test.num <- 1
        if (inc.partial == TRUE) {
                test.num <- 0.5  # now, test.num picks up partial presence
        }
        
        for (i in regions) {  # for each region of interest...
                region.col <- first.region + (i - 1)  # finds column index
                for (i in 1:n.spec)  # ...check if each specimen is there
                        if (spec.info[i, region.col] >= test.num) {
                                is.there[i] <- TRUE
                        }
        }
        return(is.there)
}

##########################
# PointOutDiffSpp
##########################
PointOutDiffSpp <- function(spec.info) {
        # Makes a vector of pch values where species of the same genus have different types of plot points.
        #
        # Args:
        #     spec.info: table of specimen information including columns named, "Genus" and "Species".
        #
        # Returns:
        #     A numeric vector of pch values which can be fed to PlotByGroup() such that species from the same genus have different points. 
        
        spec.info$Taxa <- paste(str_sub(spec.info$EGenus, 1, 3), str_sub(spec.info$Species, 1, 3), sep = "_")  # Mus musculus -> Mus_mus
        unique.taxa <- levels(as.factor(spec.info$Taxa))
        pch.taxa <- rep(21, length(unique.taxa))  # default pch = #21 circle
        
        # Give genera with multiple species different points for each spp
        n.spp <- 1  # variable to count species in a genus
        for (i in 1:(length(unique.taxa) - 1)) {
                
                if (str_sub(unique.taxa[i], 1, 3) == str_sub(unique.taxa[i + 1], 1, 3)) {
                        pch.for.taxa <- 21 + n.spp
                        if (pch.for.taxa > 25) {  # if genus has >5 spp
                                n.spp <- -14  # restart pch at 7
                        }
                        pch.taxa[i + 1] <- 21 + n.spp
                        n.spp <- n.spp + 1  # species counter
                        
                } else {
                        n.spp <- 1  # reset if no more species in that genus 
                }
        }
        return(pch.taxa)
}

##########################
# FindNatives
##########################
FindNatives <- function(spec.info, column, invasives) {
        # Creates a logical vector where TRUE indicates a native species and FALSE indicates an invasive species.
        #
        # Args:
        #    spec.info: table of specimen information that includes a column where the list invasives will distinguish the invasive species.
        #    column: the column name in spec.info to be used to ID invasives.
        #    invasives: a list of identifiers (such as species or genus names) for invasive species found in the column above.
        #
        # Returns:
        #    A logical vector in same order as spec.info table where TRUE indicates a native and FALSE indicates an invasive.
        
        n.spec <- dim(spec.info)[1]
        is.invasive <- vector(mode = "logical", length = n.spec)
        test.col <- which(colnames(spec.info) == column)
        
        for (i in 1:n.spec) {  # for each specimen...
                for (k in 1:length(invasives)) {  # ... test each invasive name
                        if (spec.info[i, test.col] == invasives[k]) {
                                is.invasive[i] <- TRUE
                        }
                }
        }
        is.native <- !is.invasive
        return(is.native)
}

##########################
# MatchTips
##########################

MatchTips <- function(tree, vector, verbose = TRUE) {
        # Creates a vector of indices to match other data to the tree with no duplicated/replicated entries of the tree. Also returns a list of specimens not found in the tree and how many times replicated specimens were replicated.  
        #
        # Args:
        #    tree: a phylogenetic tree, like that created by rtree()
        #    vector: a vector of names from another dataset, such as shape data.
        #    verbose: if TRUE, prints progress as function creates 3 returns.
        #
        # Returns:
        #    A list with 3 items: a vector with no replicated specimens of indices to match other data to the tree, specimens missing from the tree, and a list of replicated specimens.
        
        matched <- match(vector, tree$tip.label)
        
        # Test if NAs exist
        if (any(is.na(matched))){
                vector.nas <- which(is.na(matched))
                
                # Remove NA from matched list
                matched <- matched[-vector.nas]
                
                if (verbose) {
                        cat(paste0("These are not in the tree: \n\t"))
                        cat(paste0(vector[vector.nas], "\n\t"))
                        cat("\n")
                }
        } else {
                vector.nas <- NULL
        }
        
        # Check for cases of more than 1
        if (length(unique(matched)) != length(matched)) {
                count <- table(tree$tip.label[matched])
                reps <- which(count != 1)
                reps.num <- count[reps]
                reps.names <- names(reps.num)
                
                matched <- unique(matched)
                
                if (verbose) {
                        cat(paste0("Here are the replicated tips: \n"))
                        for (i in 1:length(reps)) {
                                cat(paste0(reps.names[i], "(", reps.num[i], ") " ))
                        }
                        cat("\n")
                }
        } else {
                reps.num <- NULL
        }
        
        return(list("matched" = matched, "NAs" = vector.nas, "reps" = reps.num))
}

##########################
# ReadableTable
##########################

ReadableTable <- function(table, morph.dis) {
        # Replaces numbers in a pvalue matrix returned by morphol.disparity with "more" or "less" depending on whether the column's species has significantly more or less morphological disparity than the row species. 
        #
        # Args:
        #    table: the pvalue table returned by morphol.disparity().
        #    morph.dis: the procrustes disparity table returned by morphol.disparity().
        #
        # Returns:
        #    A table of the same dimensions as the argument table with cells containing "" for non-significant values, and more or less for values where significance was found. 
        
        for (c in 1:dim(table)[1]) { # loop through cols then rows
                for (r in 1:dim(table)[2]) {
                        if (table[r, c] < 0.05) {  # if significant, test if column specie's disparity is more or less than the row species' disparity
                                if (results.morph[c] > results.morph[r]) {
                                        table[r, c] <- "more"  
                                } else (table[r, c] <- "less")
                        } else (table[r, c] <- NA)
                }
        }
        return(table)
}