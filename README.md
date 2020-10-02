
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

We always recommend using multithreaded linear algebra libraries
in R, but for this package in particular, using multithreaded libraries
can have a big impact on performance. On a Mac, there is a very simple
way to [link to the Apple Accelerate Framework](https://gist.github.com/nicebread/6920c8287d7bffb03007).
On PC and Linux, it's more complicated, but you can use 
[Microsoft R Open](https://mran.microsoft.com/open) instead, which comes automatically with multithreaded libraries.

## Basic Use

The main function for fitting models is called 'fit_model', and the
main function for doing predictions is called 'predictions'.

See this youtube video for a tutorial:
https://www.youtube.com/watch?v=phyB4n0CDWg&t=4s
