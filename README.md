# SSP
Repository of the project Simulated Sampling Procedure for Community Ecology.

This repo include the five R scripts, as well as the data sets of the four study case used in the manuscript: 

A simulation-based protocol to estimate sampling effort in studies of ecological communities, by Edlin J. Guerra-Castro, Nuno 
Simoes, Juan Cruz-Motta and Maite Mascaro.


## Install package

```r
# install.pacakges("devtools") # install if needed
devtools::install_github("edlinguerra/SSP")
```

## Use package

See article [intro](articles/intro.html).

## Develop into R package

This [Quick Intro to Package Development with devtools](http://ucsb-bren.github.io/env-info/wk07_package.html) is now a bit dated. Best to check out:

- cheat sheet: [Package Development](https://github.com/rstudio/cheatsheets/raw/master/package-development.pdf). Note that many of the functions that were in `devtools` have been ported to the `usethis` R package.

- book: [R packages](http://r-pkgs.had.co.nz/) by Hadley Wickham

### `create_package()`

```r
usethis::create_package("../SSP")
```

### functions in R/*.R

Move the `Functions/*.R` into `R/*.R` where R expects to read in functions.

### `document()`

Self-document functions by placing cursor inside function and using RStudio menu Code > "Insert Roxygen Skeleton". It will insert a template for title, parameters and return object description. For more details see:

- [Object documentation · R packages](http://r-pkgs.had.co.nz/man.html)

Use Code > "Reflow Comment" to wrap comment text for readability. After updating any documentation, run this to update R documentation in `man/` folder.

```r
devtools::document()
```

### use libaries

For installing dependencies:

```r
usethis::use_package("vegan")
usethis::use_package("ape")
usethis::use_package("ggplot2")
```

### create vignette

Show how functions work in a vignette, eg:

```r
usethis::use_vignette("intro")
```

which places the vignette into folder `vignettes/`

### Update DESCRIPTION

Customize the DESCRIPTION, like:

```
Package: SSP
Title: Simulated Sampling Procedure for Community Ecology
Authors@R: 
    person(given = "Edlin",
           family = "José Guerra Castro",
           role = c("aut", "cre"),
           email = "edlinguerra@gmail.com")
License: MIT
URL: https://github.com/edlinguerra/SSP
BugReports: https://github.com/edlinguerra/SSP/issues
```

TODO:

```
Version: 0.0.0.9000
Description: What the package does (one paragraph).
```
### Create/update website

[Introduction to pkgdown • pkgdown](https://pkgdown.r-lib.org/articles/pkgdown.html)

```r
pkgdown::build_site()
```

In your Github repository Settings, turn on `docs/` for Github Pages, per [Configuring a publishing source for GitHub Pages - GitHub Help](https://help.github.com/en/articles/configuring-a-publishing-source-for-github-pages).

### TODO: create `data/`

Using data and scripts in `raw-data/`, create objects, eg `x`, and at end of script load into `data/`, eg with `devtools::use_data(x)`.

[Data · R packages](http://r-pkgs.had.co.nz/data.html)

### load locally

```r
devtools::load_all()
```


