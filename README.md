# Patch colour pattern space generation and analysis pipeline

## Scripts

### 2_Patches/2_Scripts/

These scripts should be run in the order in which they are numbered (excluding any deprecated or optional scripts).

#### 00_Patch_number_specimens.R
Counts the number of species for which we have matched male/female colour data, and calculates percentage coverage of species/genera

*Inputs:* Processed pixel data (RDS file); taxonomic data (CSV file)

#### 01_Colourspace_mapping.R
Aggregates patch data within species/sex and maps to various colourspaces (e.g. CIELAB, USML, TCS, etc)

*Input:* Processed pixel data (RDS File)
*Output:* Patch data in different colourspaces (RDS file)

#### 02_Patch_Analyse_features.R
Performs PCA to get colour pattern spaces

*Input:* Patch data in different colourspaces (RDS file, from 01_Colourspace_mapping.R)
*Output:* Reorganised patch data in different colourspaces (RDS file); Colour pattern spaces (RDS file)

#### 02a_Patch_Phylogenetic_signal.R
Calculates Pagel's Lambda for each PC of colour pattern space

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); first 100 Hackett backbone full trees from birdtree.org (TRE file)

#### 02b_Patch_Umap_iterations.R
Performs various iterations of UMAP with different parameters to find best for data visualisation

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); taxonomic data (CSV file)
*Outputs:* Various UMAP objects (RDS files)

#### 02c_Patch_Density_plots.R
NOT NECESSARY TO RUN - plots various groups/traits on the UMAP (e.g. males vs females, passerines vs non-passerines, dichromatism), density plots of each PC axis, and interactive UMAP plot which is useful for finding the names of specific species in the UMAP

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); canonical UMAP (RDS file, from 02b_Patch_Umap_iterations.R)

#### 02d_Patch_3d_umap.R
NOT NECESSARY TO RUN - plots an interactive 3D UMAP

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); taxonomic data (CSV file)

#### 02e_Patch_colour_grids
DEPRECATED - do not run

#### 02ei_Patch_colour_grids_generate.R
Generates patch colour grids for all species for visualisation purposes. Run once only.

*Inputs:* Patch data in sRGB colourspace (RDS file, from 01_Colourspace_mapping.R)
*Outputs:* Patch colour grids for each species (PNG files)

#### 02eii_Patch_colour_grids_plot.R
Plot patch colour grids on PCA and UMAP colour pattern spaces.

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); canonical UMAP (RDS file, from 02b_Patch_Umap_iterations.R); Patch colour grids for each species (PNG files, from 02ei_Patch_colour_grids_generate.R)
*Outputs:* Colour pattern space grid plots (PNG files)

#### 02eii_Patch_colour_grids_plot_plotsforposter.R
DEPRECATED - do not run

#### 02f_Patch_space_dimensionality.R
Assess effective dimensionality of colour pattern space (linear methods only).

*Inputs:* Reorganised patch data in different colourspaces (RDS file, from 02_Patch_Analyse_features.R); Colour pattern space (RDS file, from 02_Patch_Analyse_features.R)

#### 02g_Patch_space_reconstruction.R
DEPRECATED - do not run

#### 02h_Patch_modularity.R
Assess modularity of patch data (i.e., do patches evolve independently or as integrated modules).

*Inputs:* Reorganised patch data in different colourspaces (RDS file, from 02_Patch_Analyse_features.R)

#### 03_Patch_Diversity_measures.R
DEPRECATED - do not run
Calculates various global-level diversity metrics for the entire dataset

#### 03a_Patch_Diversity_phylogeny.R
The important part of this script is the "Calculate within-group diversity" snippet (i.e., the first snippet). This calculates within-taxonomic-group diversity (mean distance to centroid) for each taxonomic subgroup and also the proportion of variance explained by taxonomic subgroup membership.
The rest of the script is deprecated.

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); taxonomic data (CSV file); first 100 Hackett backbone full trees from birdtree.org (TRE file)
*Outputs:* Within-group diversity boxplots and phylogenetic tree plots (PNG, SVG files)

#### 03b_Patch_sexual_dichromatism.R
NO NEED TO RUN - for interest only
Calculate dichromatism for each species (Euclidean distance between male and female in colour pattern space) and plot on phylogeny and colour pattern space.

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); canonical UMAP (RDS file, from 02b_Patch_Umap_iterations.R); first 100 Hackett backbone full trees from birdtree.org (TRE file)
*Outputs:* Plots (PNG files)

#### 03c_Patch_Diversity_extinction_risk.R
NO NEED TO RUN - for interest only
Boxplot of colour diversity by IUCN category, perform ANOVA to check IUCN category and diversity relation, IUCN categories highlighted

*Inputs:* Patch diversity metric data (RDS file, from 03_Patch_Diversity_measures.R)

#### 03d_Patch_diversity_SES.R
Calculate Standard Effect Size (SES) of change in diversity metric as species at different IUCN threat levels are sequentially lost

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); IUCN Red List data (CSV file, from patch-pipeline/3_SharedScripts/match_taxonomies_neognath.R)
*Outputs:* Separate results tables for males, females, and combined sexes (CSV files); Results plots (PNG, SVG files)

#### 03e_Patch_within_taxo_group_diversity_SES.R
NO NEED TO RUN - for interest only.
Calculate SES change in within-taxonomic-group diversity (mean distance to centroid) with the loss of threatened species for each taxonomic subgroup.

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); IUCN Red List data (CSV file, from patch-pipeline/3_SharedScripts/match_taxonomies_neognath.R); taxonomic data (CSV file); first 100 Hackett backbone full trees from birdtree.org (TRE file)
*Outputs:* Plots of SES change (PNG file)

#### 04_Calculate_diversity_metrics.R
DEPRECATED - do not run
Calculates various global-level diversity metrics for the entire dataset and outputs as CSV files with TXT metadata files.

#### 04a_Patch_spatial_mapping.R
DEPRECATED - do not run.
Maps grid cell assemblage diversity spatially, but only relative to *global* space centroid (i.e., more a measure of how unusual bird assemblages are relative to the global assemblage, rather than a measure of within-assemblage diversity).

#### 04b_Patch_diversity_bioregions.R
Calculate and plot on map bioregion assemblage colour pattern diversity. Option to calculate mean distance to global centroid of entire space, or local centroid of bioregion assemblage. Local centroid should be used for analyses to go in the manuscript.

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); Bioregion shape files (SHP files); Species range map Presence-Absence Matrices (RDS file)
*Outputs:* Bioregion diversity values (SHP files); Maps of bioregion-level diversity (PNG files)

#### 04b_Patch_diversity_bioregions_plotsforposter.R
DEPRECATED - do not run.
Outdated version of bioregions-level diversity mapping that allows for custom fonts when plotting.

#### 04c_Patch_diversity_regional_SES.R
DEPRECATED - do not run.
Only calculates SES diversity change relative to *global* centroid - not assemblage centroid.

#### 04c_Patch_diversity_regional_SES_local_dev.R
Calculates SES diversity change with the loss of threatened species within each bioregion assemblage.
Note that the name of this is misleading - it can calculate SES diversity change relative to *global* and *local* (i.e. assemblage) centroid. For the manuscript, *local* centroid should be used.

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); Bioregion shape files (SHP files); Species range map Presence-Absence Matrices (RDS file); IUCN Red List data (CSV file, from patch-pipeline/3_SharedScripts/match_taxonomies_neognath.R)
*Outputs:* Bioregion SES data for plotting (RDS file, SHP file, CSV file); Full bioregion SES data for manuscript tables (CSV file); Plots of bioregion SES diversity change (PNG file)

#### 04d_Patch_spatial_mapping_assemblage_level.R
Maps grid cell assemblage diversity spatially, with options for diversity relative to global or local centroid (local centroid should be used for manuscript figures).

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); Species range map Presence-Absence Matrices (RDS file)
*Outputs:* Raster of grid cell assemblage diversity (RDS file); Species richness mask (RDS file); Plots of grid cell assemblage diversity (PNG file)

#### 05a_Patch_aesthetic_attract_SES.R
Calculate SES change in overall avian aesthetic value as threatened species are lost. Includes code to match iratebirds (Clements/eBird) taxonomy to Jetz taxonomy.

*Inputs:* [Santangeli et al 2023](https://doi.org/10.1038/s44185-023-00026-2) iratebirds aesthetic value dataset (CSV file); McTavish-Jetz taxonomy crosswalk (CSV file); IUCN Red List data (CSV file, from patch-pipeline/3_SharedScripts/match_taxonomies_neognath.R); Colour pattern space (RDS file, from 02_Patch_Analyse_features.R)
*Outputs:* Separate results tables for males, females, and combined sexes (CSV files); Results plots (PNG, SVG files)