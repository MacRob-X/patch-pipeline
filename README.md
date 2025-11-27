# patch-pipeline

## Scripts for Chris

### 2_Patches/2_Scripts/

#### 00_Patch_number_specimens.R
Counts the number of species for which we have matched male/female colour data, and calculates percentage coverage of species/genera\

*Inputs:* Processed pixel data (RDS file); taxonomic data (CSV file)\

#### 01_Colourspace_mapping.R
Aggregates patch data within species/sex and maps to various colourspaces (e.g. CIELAB, USML, TCS, etc)\

*Input:* Processed pixel data (RDS File)\
*Output:* Patch data in different colourspaces (RDS file)\

#### 02_Patch_Analyse_features.R
Performs PCA to get colour pattern spaces\

*Input:* Patch data in different colourspaces (RDS file, from 01_Colourspace_mapping.R)\
*Output:* Reorganised patch data in different colourspaces (RDS file); Colour pattern spaces (RDS file)\

#### 02a_Patch_Phylogenetic_signal.R
Calculates Pagel's Lambda for each PC of colour pattern space\

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); first 100 Hackett backbone full trees from birdtree.org (TRE file)\

#### 02b_Patch_Umap_iterations.R
Performs various iterations of UMAP with different parameters to find best for data visualisation\

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); taxonomic data (CSV file)\
*Outputs:* Various UMAP objects (RDS files)\

#### 02c_Patch_Density_plots.R
NOT NECESSARY TO RUN - plots various groups/traits on the UMAP (e.g. males vs females, passerines vs non-passerines, dichromatism), density plots of each PC axis, and interactive UMAP plot which is useful for finding the names of specific species in the UMAP\

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); canonical UMAP (RDS file, from 02b_Patch_Umap_iterations.R)\

### 02d_Patch_3d_umap.R
NOT NECESSARY TO RUN - plots an interactive 3D UMAP\

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); taxonomic data (CSV file)\

### 02e_Patch_colour_grids
DEPRECATED - do not run\

### 02ei_Patch_colour_grids_generate.R
Generates patch colour grids for all species for visualisation purposes. Run once only.\

*Inputs:* Patch data in sRGB colourspace (RDS file, from 01_Colourspace_mapping.R)\
*Outputs:* Patch colour grids for each species (PNG files)\

### 02eii_Patch_colour_grids_plot.R
Plot patch colour grids on PCA and UMAP colour pattern spaces.\

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); canonical UMAP (RDS file, from 02b_Patch_Umap_iterations.R); Patch colour grids for each species (PNG files, from 02ei_Patch_colour_grids_generate.R)\
*Outputs:* Colour pattern space grid plots (PNG files)

### 02eii_Patch_colour_grids_plot_plotsforposter.R
DEPRECATED - do not run\

### 02f_Patch_space_dimensionality.R
Assess effective dimensionality of colour pattern space (linear methods only).\

*Inputs:* Reorganised patch data in different colourspaces (RDS file, from 02_Patch_Analyse_features.R); Colour pattern space (RDS file, from 02_Patch_Analyse_features.R)\

### 02g_Patch_space_reconstruction.R
DEPRECATED - do not run\

### 02h_Patch_modularity.R
Assess modularity of patch data (i.e., do patches evolve independently or as integrated modules).\

*Inputs:* Reorganised patch data in different colourspaces (RDS file, from 02_Patch_Analyse_features.R)

### 03_Patch_Diversity_measures.R
DEPRECATED - do not run\
Calculates various global-level diversity metrics for the entire dataset\

### 03a_Patch_Diversity_phylogeny.R
The important part of this script is the "Calculate within-group diversity" snippet (i.e., the first snippet). This calculates within-taxonomic-group diversity (mean distance to centroid) for each taxonomic subgroup and also the proportion of variance explained by taxonomic subgroup membership.\
The rest of the script is deprecated.\

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); taxonomic data (CSV file); first 100 Hackett backbone full trees from birdtree.org (TRE file)\
*Outputs:* Within-group diversity boxplots and phylogenetic tree plots (PNG, SVG files)

### 03b_Patch_sexual_dichromatism.R
NO NEED TO RUN - for interest only\
Calculate dichromatism for each species (Euclidean distance between male and female in colour pattern space) and plot on phylogeny and colour pattern space.\

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); canonical UMAP (RDS file, from 02b_Patch_Umap_iterations.R); first 100 Hackett backbone full trees from birdtree.org (TRE file)\
*Outputs:* Plots (PNG files)

### 03c_Patch_Diversity_extinction_risk.R
NO NEED TO RUN - for interest only\
Boxplot of colour diversity by IUCN category, perform ANOVA to check IUCN category and diversity relation, IUCN categories highlighted\

*Inputs:* Patch diversity metric data (RDS file, from 03_Patch_Diversity_measures.R)

### 03d_Patch_diversity_SES.R
Calculate Standard Effect Size (SES) of change in diversity metric as species at different IUCN threat levels are sequentially lost\

*Inputs:* Colour pattern space (RDS file, from 02_Patch_Analyse_features.R); IUCN Red List data (CSV file, from patch-pipeline/3_SharedScripts/match_taxonomies_neognath.R)
*Outputs:* Separate results tables for males, females, and combined sexes (CSV files); Results plots (PNG, SVG files)

### 03e_Patch_within_taxo_group_diversity_SES.R
Calculate within-taxonomic-group diversity (mean distance to centroid) for each taxonomic subgroup