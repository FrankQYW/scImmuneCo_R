# scImmuneCo
## About
We developed scImmuneCo, a comprehensive resource of cell-type-specific transcriptional modules derived from single-cell RNA sequencing (scRNA-seq) data. Our analytical framework, applied to 17 immunological conditions spanning autoimmunity, immunodeficiency, infection, and hematologic malignancies, identified 873 co-expression modules across seven major immune cell types. 

ScImmuneCo represents a significant methodological advance by providing stable, reusable modules that overcome single-dataset limitations; resolving cell-type-specific functional programs lost in bulk analyses; capturing transitional cellular states often missed by conventional clustering. 

<img src="./man/figures/Figure1.png" width="100%" style="display: block; margin: auto;" />

## Intallation
``` r
library(devtools)
devtools::install_github("FrankQYW/scImmuneCo_R")
```


## Usage
### Preprocess of the data
A preprocessed PBMC single cell RNA-seq in Seurat format is required to run the scImmuneCo package. Users are suggested to use Azimuth to map the PBMC data to the reference cell types. 


``` r
library(Azimuth)
seurat_object <- RunAzimuth(seurat_object, reference = "pbmcref")
```

As a result, the meta_data of seurat_object will have a 'predicted.celltype.l1' column contain the cell type:


``` r
unique(seurat_object@meta.data$predicted.celltype.l1)

[1] "CD8 T"   "CD4 T"   "Mono"    "other T" "DC"      "B"       "NK"      "other" 
``` 

Users could also manually annotate the single cell data as long as the name of the cell type aligned with the ones in Azimuth reference. 



### Identify Differentially Expressed Modules Between Conditions
We provide a convenient one line function to test modules in all cell types between two conditions.   


``` r
res <- do_dem(
  seurat_object,
  sample_column = "batch",
  gmt = module_info,
  condition_column = "condition",
  condition1 = "SjS",
  condition2 = "HC"
)



The users can also choose a specific cell type to test the modules differene between two conditions. 


``` r
gsva_matrix <- gsva_cell_type(seurat_object, cell = 'CD4 T', sample_column = 'batch', 
          cell_column = "predicted.celltype.l1",)

res <- compare_condition(seurat_object, gsva = a, 
                       sample_column = 'batch', 
                       condition_column = 'condition', 
                       condition1 = 'SjS', 
                       condition2 = 'HC')
``` 



``` 
The overall results can be visualized by 

``` r
draw_volcano(res, p_cutoff = 0.05, logFC_cutoff = 0.25)

``` 

<img src="./man/figures/volcano.png" width="100%" style="display: block; margin: auto;" />


Can also visualize the differential expressed modules on meta module level in each cell type

``` r
Mono_res <- filter(res, res$cell_type == 'Mono')


``` 
<img src="./man/figures/meta_module_volcano.png" width="100%" style="display: block; margin: auto;" />
