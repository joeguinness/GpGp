
# GpGp

GpGp is an R package for fast approximate Gaussian process computation. 
The package includes implementations of the Vecchia's (1988) original 
approximation, as well as several updates to it, including the reordered 
and grouped versions of the approximation outlined in Guinness (2018).

## Installing

The package can be installed from CRAN with the usual R command

```{r}
install.packages("GpGp")
```

or directly from Github for the latest version

```{r}
devtools::install_github("joeguinness/GpGp")
```

## Basic Use

The main function for fitting models is called 'fit_model', and the
main function for doing predictions is called 'predictions'.

See this youtube video for a tutorial:
https://www.youtube.com/watch?v=phyB4n0CDWg&t=4s
