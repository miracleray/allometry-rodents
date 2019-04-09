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
        
        pairs <- NULL
        
        # Removes R and L designations so pairs can be detected
        no.side.names <- gsub(" R", "", pt.names)
        no.side.names <- gsub(" L", "", no.side.names)
        
        # Checks if point has a pair and if so, their index #s are paired
        for(i in unique(no.side.names)){
                index <- which(no.side.names == i)
                if (length(index) == 2) {
                        pairs <- rbind(pairs, t(index))
                }
        }
        colnames(pairs) <- c("Right", "Left")
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
                         main = "Principal Component Analysis",
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

##########################
# FindRateDiff
##########################

FindRateDiff <- function(test) {
        # Determines if rate of group of interest is greater or less than other group.
        #
        # Args:
        #    test: object of class "evolrate" produced by geomorph's compare.evol.rates function with 2 groups and group 1 defined as the group of interest. 
        #
        # Returns: 
        #    A string of "Faster" or "Slower" based on test outcome.
        difference <- test$sigma.d.gp[1] - test$sigma.d.gp[2]
        if (difference > 0) {
                return("Faster")
        } else {
                "Slower"
        }
}

##########################
# plotGMphyloMS
##########################

plotGMPhyloMS<-function(phy,A,tip.labels=TRUE,node.labels=TRUE,ancStates=TRUE, xaxis=1, yaxis=2, zaxis=NULL, plot.param = list(), shadow=FALSE){
        if(any(is.na(A))==T){
                stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
        if (length(dim(A))==3){ 
                if(is.null(dimnames(A)[[3]])){
                        stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")  }
                x<-two.d.array(A)}
        if (length(dim(A))==2){ 
                if(is.null(dimnames(A)[[1]])){
                        stop("Data matrix does not include taxa names as dimnames for rows.")  }
                x<-A }
        if (!inherits(phy, "phylo"))
                stop("tree must be of class 'phylo.'")
        phy <- reorder.phylo(phy)
        N<-length(phy$tip.label)
        Nnode <- phy$Nnode
        if(!(N - phy$Nnode + (tabulate(phy$edge[, 1])[N + 1] >= 2)) == 2)
                stop("tree is not fully bifurcating (consider 'multi2di' in ape.")
        
        if(N!=dim(x)[1]){
                stop("Number of taxa in data matrix and tree are not not equal.")  }
        if(length(match(rownames(x), phy$tip.label))!=N) 
                stop("Data matrix missing some taxa present on the tree.")
        if(length(match(phy$tip.label,rownames(x)))!=N) 
                stop("Tree missing some taxa in the data matrix.")
        x<-x[phy$tip.label, ]  
        anc.states <- anc.Bayes(phy, x)
        all.data<-rbind(x,anc.states)  
        pcdata<-prcomp(all.data)$x 
        pcdata<-pcdata-matrix(rep(pcdata[(N+1),],nrow(pcdata)), nrow=nrow(pcdata),byrow=T)  #phylogenetic mean adjustment
        #plotting  
        p.p <- plot.param
        if(is.null(p.p$t.bg)) p.p$t.bg="black" ; if(is.null(p.p$t.pch)) p.p$t.pch=21
        if(is.null(p.p$t.cex)) p.p$t.cex=2 ; if(is.null(p.p$n.bg)) p.p$n.bg="white"
        if(is.null(p.p$n.pch)) p.p$n.pch=21 ; if(is.null(p.p$n.cex)) p.p$n.cex=1.25
        if(is.null(p.p$l.col)) p.p$l.col="black" ; if(is.null(p.p$lwd)) p.p$lwd=3
        if(is.null(p.p$txt.adj)) p.p$txt.adj=c(-.1,-.1) ; if(is.null(p.p$txt.col)) p.p$txt.col="black"
        if(is.null(p.p$txt.cex)) p.p$txt.cex=1 ; if(is.null(p.p$n.txt.adj)) p.p$n.txt.adj=c(-.1,-.1) 
        if(is.null(p.p$n.txt.col)) p.p$n.txt.col="black" ; if(is.null(p.p$n.txt.cex)) p.p$n.txt.cex=0.6
        limits = function(x,s){ 
                r = range(x)
                rc=scale(r,scale=F)
                l=mean(r)+s*rc}
        # regular 2D phylomorphospace
        if(is.null(zaxis)){
                if(tip.labels==TRUE){
                        plot(pcdata[,xaxis],pcdata[,yaxis],type="n",xlim=limits(pcdata[,xaxis],1.5),ylim=limits(pcdata[,yaxis],1.5),asp=1,
                             xlab = paste("PC", xaxis), ylab = paste("PC", yaxis)) }
                if(tip.labels==FALSE) {
                        plot(pcdata[,xaxis],pcdata[,yaxis],type="n",asp=1, xlab = paste("PC", xaxis), ylab = paste("PC", yaxis)) }
                for (i in 1:nrow(phy$edge)){
                        lines(pcdata[(phy$edge[i,]),xaxis],pcdata[(phy$edge[i,]),yaxis],type="l",col=p.p$l.col,lwd=p.p$lwd)
                }
                points(pcdata[1:N,xaxis], pcdata[1:N,yaxis],pch=p.p$t.pch, bg=p.p$t.bg, cex=p.p$t.cex)
                points(pcdata[(N+1):nrow(pcdata),xaxis], pcdata[(N+1):nrow(pcdata),yaxis],pch=p.p$n.pch, bg=p.p$n.bg, cex=p.p$n.cex)
                if(tip.labels==TRUE){
                        text(pcdata[1:N,xaxis],pcdata[1:N,yaxis],rownames(pcdata)[1:N],
                             col=p.p$txt.col,cex=p.p$txt.cex,adj=p.p$txt.adj)}
                if(node.labels==TRUE){
                        text(pcdata[(N + 1):nrow(pcdata),xaxis],pcdata[(N + 1):nrow(pcdata),yaxis],rownames(pcdata)[(N + 1):nrow(pcdata)],
                             col=p.p$n.txt.col,cex=p.p$n.txt.cex,adj=p.p$n.txt.adj)}
        }
        # 3d phylomorphospace in rgl
        if(is.numeric(zaxis)){
                plot3d(pcdata[1:N,xaxis], pcdata[1:N,yaxis], pcdata[1:N,zaxis],type="s",xlim=limits(pcdata[,xaxis],1.5),
                       ylim=limits(pcdata[,yaxis],1.5), zlim=limits(pcdata[,zaxis],1.5), asp=1,
                       xlab= paste("PC",xaxis), ylab= paste("PC",yaxis), zlab=paste("PC",zaxis),
                       col= p.p$t.bg, size=p.p$t.cex)
                if(p.p$n.bg == "white"){ p.p$n.bg <- "grey"}
                points3d(pcdata[(N + 1):nrow(pcdata),xaxis], pcdata[(N + 1):nrow(pcdata),yaxis], 
                         pcdata[(N + 1):nrow(pcdata),zaxis], 
                         col= p.p$n.bg, size=p.p$n.cex*4)
                for (i in 1:nrow(phy$edge)) {
                        lines3d(pcdata[(phy$edge[i, ]), xaxis], pcdata[(phy$edge[i, ]), yaxis],pcdata[(phy$edge[i, ]), zaxis], 
                                col=p.p$l.col, lwd=p.p$lwd)}
                if(tip.labels==TRUE){
                        text3d(pcdata[1:N,xaxis],pcdata[1:N,yaxis],pcdata[1:N,zaxis],rownames(pcdata)[1:N],
                               col=p.p$txt.col,cex=p.p$txt.cex,adj=p.p$txt.adj) }
                if(node.labels==TRUE){
                        text3d(pcdata[(N + 1):nrow(pcdata),xaxis],pcdata[(N + 1):nrow(pcdata),yaxis],
                               pcdata[(N + 1):nrow(pcdata),zaxis],rownames(pcdata)[(N + 1):nrow(pcdata)],
                               col=p.p$n.txt.col,cex=p.p$n.txt.cex,adj=p.p$n.txt.adj) }
        }
        # 3d phylomorphospace in rgl with time on Z-axis
        if(is.character(zaxis)){
                zaxis <- getNodeDepth(phy)
                zaxis <- abs(zaxis - max(zaxis))
                view3d(phi=90, fov=30)
                plot3d(pcdata[,xaxis],pcdata[,yaxis],zaxis,type="n",xlim=limits(pcdata[,xaxis],1.5),
                       ylim=limits(pcdata[,yaxis],1.5),
                       zlim=c(max(zaxis), min(zaxis)),
                       asp=c(1,1,1),
                       xlab= paste("PC",xaxis), ylab= paste("PC",yaxis), zlab="Time")
                points3d(pcdata[1:N,xaxis], pcdata[1:N,yaxis], zaxis[1:N],
                         col= p.p$t.bg, size=p.p$t.cex*4)
                if(p.p$n.bg == "white"){ p.p$n.bg <- "grey"}
                points3d(pcdata[(N + 1):nrow(pcdata),xaxis], pcdata[(N + 1):nrow(pcdata),yaxis], 
                         zaxis[(N + 1):nrow(pcdata)], 
                         col= p.p$n.bg, size=p.p$n.cex*4)
                for (i in 1:nrow(phy$edge)) {
                        lines3d(pcdata[(phy$edge[i, ]), xaxis], pcdata[(phy$edge[i, ]), yaxis],zaxis[(phy$edge[i, ])], 
                                col=p.p$l.col, lwd=p.p$lwd)}
                if(tip.labels==TRUE){
                        text3d(pcdata[1:N,xaxis],pcdata[1:N,yaxis],zaxis[1:N],rownames(pcdata)[1:N],
                               col=p.p$txt.col,cex=p.p$txt.cex,adj=p.p$txt.adj) }
                if(node.labels==TRUE){
                        text3d(pcdata[(N + 1):nrow(pcdata),xaxis],pcdata[(N + 1):nrow(pcdata),yaxis],
                               zaxis[(N + 1):nrow(pcdata)],rownames(pcdata)[(N + 1):nrow(pcdata)],
                               col=p.p$n.txt.col,cex=p.p$n.txt.cex,adj=p.p$n.txt.adj) }
                if(shadow==TRUE){
                        #plot shadow version at base
                        points3d(pcdata[1:N,xaxis], pcdata[1:N,yaxis], max(zaxis),
                                 col= p.p$t.bg, size=p.p$t.cex*4, alpha = 0.5)
                        points3d(pcdata[(N + 1):nrow(pcdata),xaxis], pcdata[(N + 1):nrow(pcdata),yaxis], 
                                 max(zaxis), col= p.p$n.bg, size=p.p$n.cex*4,alpha = 0.5)
                        for (i in 1:nrow(phy$edge)) {
                                lines3d(pcdata[(phy$edge[i, ]), xaxis], pcdata[(phy$edge[i, ]), yaxis],max(zaxis), 
                                        col=p.p$l.col, lwd=p.p$lwd, alpha = 0.5)}
                }
        }
        if(ancStates==TRUE){ return(anc.states)  }
}
