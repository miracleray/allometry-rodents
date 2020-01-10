# Function modified from phytools's plotTraitbyBranch assuming tip.traits is a vector = the number of tips.

tip.traits <- info.means$MeanCsize
tree <- aus.tree

phytools.branch.colors <- function(tree, tip.traits, col.fun) {
        
        # Estimate ancestral states by phytools assuming mode = tips
        traits <- c(tip.traits[tree$tip.label], fastAnc(tree, tip.traits))
        names(traits)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)
        trait.matrix <- matrix(traits[tree$edge], nrow(tree$edge), 2)
        traits <-rowMeans(trait.matrix)
        
        # Color breaks by phytools
        col <- gray(1000:1/1000) # colors
        tol <- 1e-6
        xlims <- range(traits) + c(-tol, tol)
        breaks <- 0:1000/1000 * (xlims[2] - xlims[1]) + xlims[1]
        
        # Select the right color for the one trait
        color.selector <- function(one.trait, col, breaks) {
                break.index <- 1  # initialise the break index
                while(one.trait >= breaks[break.index] && one.trait > breaks[break.index+1]) {
                        # Increment index
                        break.index <- break.index + 1
                }
                return(col[break.index])  # return the corresponding color
        }
        
        # Select & return the color for all traits
        branch.colors <- sapply(traits, color.selector, col, breaks)
        return(darken(branch.colors))
}