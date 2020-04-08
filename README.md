# Australian rodent skull allometry (Chapter 3 of Thesis) - data and code
Code authors: Ariel E. Marcy, Thomas Guillerme, Emma Sherratt

To cite the paper and/or code:
> Coming soonish

As of April 2020, this is still a work in progress. Relies on `R` (v. 3.6.1), `geomorph` (v. 3.1.3), and `landvR` (v. 0.4).

*All scripts are in RMarkdown format (.Rmd) and can be opened in RStudio. There, you can edit and run code chunks as normal or use the Knit button to create HTML versions with both code and output. After cloning this repo, remember to either set your working directory to the allometry-rodents folder on your computer or open an RStudio project from that folder.*

## Data
**Landmarking data:**
* [MorphoSource Project 561](https://www.morphosource.org/MyProjects/Dashboard/dashboard/select_project_id/561) publically provides 3D meshes for all surface scanned crania landmarked in the study.
* [Raw_Coordinates.csv](Data/Raw/Raw_Coord_Data.csv) provides the shape coordinates from landmarked 3D skulls, exported from Viewbox.

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
**The first three scripts prepare the data for analysis and plotting**, the intermediate data generated by script 03 and 04 are stored in .rda files and used in later scripts.

* [01-extract-data-for-analyses.Rmd](/Analysis/01-extract-data-for-analyses.Rmd) Extracts both 3D coordinate data and metadata from Viewbox. Prepares them for analysis in `geomorph` by running GPA with bilateral symmetry and then merging the symmetric shape component with centroid size for each specimen.
* [02-calculate-user-error.Rmd](/Analysis/02-calculate-user-error.Rmd) Allows users to view outliers and find major landmarking errors. Calculates user error based on repeatability.
* [03-prepare-data.Rmd](/Analysis/03-prepare-data.Rmd) Creates the metadata and graphing vectors needed for all subsequent analyses. Prepares mean shape data and the ultrametric phylogenetic tree for analyses that require a phylogeny. 

**The next three scripts produce the main results:**

* [04-measure-allometry.Rmd](/Analysis/04-measure-allometry.Rmd) Tests static and evolutionary allometry with `geomorph`'s `procD.lm`'s ANCOVA function. `RRPP`'s `pairwise` function. This includes Procrustes ANOVAs and pairwise homogeneity of slopes (HOS) tests. Tests whether the HOS result is distinct from null (global slope) expectancy with a randomization procedure. **Generates Table 1 and Table A4**.
* [05-phylogenetic-rarefaction.Rmd](/Analysis/05-phylogenetic-rarefaction.Rmd) Uses phylogenetic rarefaction to tests for impact of species-level sample sizes on allometric findings with rarifaction across every clade in the phylogeny using `landvR`. **Generates Figure 3 and Figure A2**.
* [06-measure-evo-rates.Rmd](/Analysis/06-measure-evo-rates.Rmd) Runs `geomorph`'s `compare.evol.rates()` function on species with locomotion and dietary specializations to see if they evolved more slowly or more quickly than sister taxa. **Generates Table A6**.

**The next two scripts produce the main multi-panel figures:** 

* [07-plot-heatmaps.Rmd](/Analysis/07-plot-heatmaps.Rmd) Plots landmark heatmaps showing statistically rigorous visualizations of allometric shape variation from the mean (consensus) shape to either the PC1 mininum or to PC1 maximum (dorsal and lateral views). **Generates the four images for Figure 4C-F**.
* [08-plot-allometry.Rmd](/Analysis/08-plot-allometry.Rmd) Plots static allometry with both Regression Scores and Predicted Values. Plots evolutionary allometry of mean shape data, a PCA of the mean shape data, combined with four heatmap views of the crania shape changes across PC1. **Generates Figure 1, Figure 2, and Figure 4**.

### Custom functions 
* Some analyses call custom functions, which are defined in the [..Data/Functions/utilities.R](/Data/Functions/utilities.R) file.

* **Sup01-combine-metadata-sex-info.Rmd** This supplementary script reattaches sex information and museum of origin to the specimen metadata using the four museums' original metadata. **Generates Table A1**
