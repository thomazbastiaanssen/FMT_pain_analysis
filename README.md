# Pathological arm

    ## [1] "Using the following formula: x ~ Treatment + Timepoint + (1 | animal_ID) + Treatment:Timepoint"
    ## [1] "Adjusting for FDR using Benjamini & Hochberg's procedure."
    ## [1] "Using the following formula: x ~ Treatment + Timepoint + (1 | animal_ID) + Treatment:Timepoint"
    ## [1] "Adjusting for FDR using Benjamini & Hochberg's procedure."
    ## [1] "Using the following formula: x ~ Treatment + Timepoint + (1 | animal_ID) + Treatment:Timepoint"
    ## [1] "Adjusting for FDR using Benjamini & Hochberg's procedure."
    ## [1] "Operating in interaction mode"
    ## [1] "94 were matched between table 1 and the columns of the adjacency matrix"
    ## [1] "786 were matched between table 2 and the rows of the adjacency matrix"
    ## [1] "Running annotation-based correlations"
    ## [1] "Running correlations for the following groups: All, abx + FMTCTR, abx + FMTDNBS"
    ## [1] "Fitting models for differential correlation testing"
    ## [1] "Model type:lm"
    ## [1] "Adjusting p-values using Benjamini & Hochberg's procedure."
    ## [1] "Using theoretical distribution."

``` r
library(patchwork)
```

``` r
(p_alpha + p_beta + plot_layout(widths = c(3, 2)) ) / p_genus + 
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = "i")
```

![](README_files/figure-gfm/plot_path-1.png)<!-- -->

# Therapeutic arm

    ## [1] "Using the following formula: x ~ Treatment + Timepoint + (1 | animal_ID) + Treatment:Timepoint"
    ## [1] "Adjusting for FDR using Benjamini & Hochberg's procedure."
    ## [1] "Using the following formula: x ~ Treatment + Timepoint + (1 | animal_ID) + Treatment:Timepoint"
    ## [1] "Adjusting for FDR using Benjamini & Hochberg's procedure."
    ## [1] "Using the following formula: x ~ Treatment + Timepoint + (1 | animal_ID) + Treatment:Timepoint"
    ## [1] "Adjusting for FDR using Benjamini & Hochberg's procedure."
    ## [1] "Operating in interaction mode"
    ## [1] "94 were matched between table 1 and the columns of the adjacency matrix"
    ## [1] "786 were matched between table 2 and the rows of the adjacency matrix"
    ## [1] "Running annotation-based correlations"
    ## [1] "Running correlations for the following groups: All, DNBS + FMTCTR, DNBS + FMTDNBS"
    ## [1] "Fitting models for differential correlation testing"
    ## [1] "Model type:lm"
    ## [1] "Adjusting p-values using Benjamini & Hochberg's procedure."
    ## [1] "Using theoretical distribution."

``` r
(t_alpha + t_beta + plot_layout(widths = c(3, 2)) ) / t_GBM + 
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = "i")
```

![](README_files/figure-gfm/plot_ther-1.png)<!-- -->

``` r
# 
# t_alpha
# t_beta
# t_GBM
# t_GMM
```
