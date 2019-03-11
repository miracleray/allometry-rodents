# Australian rodent skull allometry (Chapter 3 of Thesis) - data and code
Code author: Ariel E. Marcy

To cite the paper and/or code:
> Coming soonish

As of March 2018, this is still a work in progress.

*All of the scripts are in RMarkdown format (.Rmd), which can be opened in RStudio. There, you can edit and run code chunks as normal, or you can click the Knit button to create HTML versions with both code and output. After cloning this repo, remember to either set your working directory to the allometry-rodents folder on your computer, or open an RStudio project from that folder.*

## Data
**Landmarking data:**
* 3D surface scanned meshes for all skulls in the study will be available via [MorphoSource](https://www.morphosource.org/)
* [Raw_Coordinates.csv](Data/Raw/3D_coords.csv) - the shape coordinates from landmarking 3D  skulls in Viewbox 

**Museum metadata provided by Curators:**
* [Australian Museum specimens](/Data/Raw/AM_muridae_skulls.csv)
* [Melbourne Museum specimens](/Data/Raw/MV_muridae_skulls.csv)
* [Queensland Museum specimens](/Data/Raw/QM_muridae_skulls.csv)
* [South Australian Museum specimens](/Data/Raw/SAM_muridae_skulls.csv)

**Ecological metadata:**
* [Trait data from Breed & Ford 2007](/Data/Processed/in_ex_traits.csv)

If you use these data, please cite the original authors:
> Breed B & Ford F. 2007. Native Mice and Rats. Australian Natural History Series, CSIRO Publishing: Colling-wood, Victoria, Australia, 185 pp. ISBN 978-0-6430-9166-5.

**Phylogenetic data:**
* [Fossil calibrated ultrametric tree modified from Smissen & Rowe 2018](/Data/Processed/Marcy-BEAST01.con.tre)

If you use these data, please cite this paper and the original authors:
> Smissen PJ & Rowe KC. 2018. Repeated biome transitions in the evolution of Australian Rodents. Molecular Phylogenetics and Evolution. 128:182–191. doi: 10.1016/j.ympev.2018.07.015.
    
## Analyses
The analysis workflow is broken down into smaller scripts explained below. Each script loads data created by the script before, so this workflow requires you to run the scripts in order. The intermediate data -- stored as .rda files in the /Data/Processed folder -- are too large to upload to Github. 

* **01-extract-data-for-analyses.Rmd** Extracts both 3D coordinate data as well as the metadata data from Viewbox and prepares it for analysis in `geomorph`. Separates coordinate data into big and small patch protocol datasets. Runs GPA with bilateral symmetry, merges the symmetric shape component with centroid size, and calculates asymmetric variation for both patch datasets.
* **02-calculate-user-error.Rmd** Allows users to view outliers and find major landmarking errors. Takes out replicated specimens from the shape data, finds their duplicates, and calculates user error based on repeatability for both patch datasets.
* **03-measure-static-allometry.Rmd** Tests correlation of PC1 and PC2 with centroid size. Uses `geomorph`'s `procD.allometry()` and `advanced.procD.lm()` to test for significant differences in allometric slopes and intercepts by genus and by clade. Plots allometric relationships using the Common Allometric Component. Plots "size-less" residuals of static allometry. **Generates Figure 1**.
* **04-plot-evolutionary-allometry.Rmd** Finds mean shape and centroid size in preparation for evolutionary allometry and phylogenetic analyses. Plots PCA of shape data. **Generates Figure 2**.
* **05-plot-heatmaps.Rmd** Plot heatmaps showing a statistically rigorous visualization of how landmarks vary from PC1 mininum to PC1 maximum. Ilustrates the allometric changes in Australian rodents. **Generates Figure 3**.
* **06-measure-evolutionary-allometry** Organizes the tree, mean shape data, and metadata in preparation for phylogenetic analyses of shape using the 2018 Smissen & Rowe tree. Runs phylogenetic analyses including evolutionary signal, and phylogenetically-informed `procD.pgls()` linear model for evolutionary allometry. **Generates Table 1**.
* **06-find-evolutionary-rates.Rmd** Runs `geomorph`'s `compare.evol.rates()` function on species with locomotion and dietary specializations to see if they evolved more slowly or more quickly than sister taxa. **Generates Table 2**.

* **Sup01-combine-metadata-sex-info.Rmd** This supplementary script reattaches sex information and museum of origin to the specimen metadata using the 4 museums' original metadata.  

**Custom functions in the utility file:** The analyses call custom functions that are defined in the ..Data/Functions/utilities.R file.
