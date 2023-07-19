# CopyClust

A Copy-Number-Based Integrative Cluster Classification Scheme for Breast Cancer Samples

## Installation

The CopyClust package can be installed in R:

```r
 devtools::install_github("camyoung54/CopyClust")
```
    
## Documentation

This package serves as an alternative to the [iC10 package]() developed by Rueda et al. for the classification of breast cancer samples into integrative clusters when gene expression data is unavailable. CopyClust implements an XGBoost-based classifier trained on the copy number profiles of the [METABRIC cohort](https://www.nature.com/articles/nature10983) to predict integrative cluster label based on copy number data alone. The iC10 package can be downloaded from CRAN:

```R
install.packages("iC10")
library(iC10)
```


## Authors

- [Cameron Young](https://www.github.com/camyoung54)


## License

[MIT](https://choosealicense.com/licenses/mit/)


## References

 - [Curtis et al., The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups, (2012).](https://www.nature.com/articles/nature10983)
 - [Ali et al., Genome-driven integrated classification of breast cancer validated in over 7,500 samples, (2014).](https://link.springer.com/article/10.1186/s13059-014-0431-1)

