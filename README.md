# RDHonest-vStata

This Stata package calculates honest and nearly-optimal one- and two-sided confidence intervals in fuzzy and sharp regression discontinuity designs based on local linear regression, see [RDHonest](https://github.com/kolesarm/RDHonest) package for an R version.

See the [vignette](https://tbarmstr.github.io/RDHonest-vStata/vignette/rdhonest_stata_vignette.html) for a description of the package, and the package manual available through `help rdhonest` once the package is installed.

## Installation

### Install via SSC

The released version of `rdhonest` is coming soon.

### Install development version via github

To install the package

```stata
* Remove program if it existed previously
capture ado uninstall rdhonest
* Install most up-to-date version
net install rdhonest, from("https://raw.githubusercontent.com/tbarmstr/RDHonest-vStata/master/")
```

### Manual Installation

To download the development version of these packages from GitHub, download the files
[`rdhonest.ado`](rdhonest.ado), [`rdhonest.sthlp`](rdhonest.sthlp), and put them into Stata's personal `ado` directory,
typically

- `c:\ado\personal` on Windows
- `~/Documents/Stata/ado` on a Mac
- `~/ado` on Linux
  
  Or run `sysdir` in Stata to find out local directories.

For more information on how to use personal ado files, please refer to [Stata Official FAQ](https://www.stata.com/support/faqs/programming/personal-ado-directory/).

## Structure

- the master folder contains the up-to-date version of the Stata package, including
  - codes: [`rdhonest.ado`](rdhonest.ado)
  - help file: [`rdhonest.sthlp`](rdhonest.sthlp)
  - Stata installation files: [`rdhonest.pkg`](rdhonest.pkg) and [`stata.toc`](stata.toc)

- subfolder [`tests`](https://github.com/tbarmstr/RDHonest-vStata/tree/main/tests) contains testing documentations:
  - testing script: [`rdhonest_test.do`](tests/rdhonest_test.do)
  - log file [`rdhonest_test.log`](tests/rdhonest_test.log)

- subfolder [`data`](https://github.com/tbarmstr/RDHonest-vStata/tree/main/data) contains 5 exemplary data sets, including
  - [`cghs.dta`](data/cghs.dta): from Oreopoulos (2006)
  - [`headst.dta`](data/headst.dta): from Ludwig and Miller (2007)
  - [`lee08.dta`](data/lee08.dta): from Lee (2008)
  - [`rcp.dta`](data/rcp.dta): from Battistin, Brugiavini, Rettore and Weber (2009)
  - [`rebp.dta`](data/rebp.dta): from Lalive (2008)
  
- subfolder [`vignette`](https://github.com/tbarmstr/RDHonest-vStata/tree/main/vignette) contains the vignette

## Bug reporting and Questions

Please open issues and leave feedback, use click on the [`Issues`](https://github.com/tbarmstr/RDHonest-vStata/issues) tab.

## References

- Main references:

  Armstrong, T. B., & Kolesár, M. (2018). Optimal inference in a class of regression models. *Econometrica*, 86, 655-683. [doi:10.3982/ECTA14434](https://doi.org/10.3982/ECTA14434).

  Armstrong, T. B., & Kolesár, M. (2020). Simple and honest confidence intervals in nonparametric regression. *Quantitative Economics*, 11(1), 1-39. [doi:10.3982/QE1199](https://doi.org/10.3982/QE1199)

  Kolesár, M., & Rothe, C. (2018). Inference in regression discontinuity designs with a discrete running variable. *American Economic Review*, 108, 2277–2304. [doi:10.1257/aer.20160945](https://www.aeaweb.org/articles?id=10.1257/aer.20160945)

- Data references:
  - [`cghs.dta`](data/cghs.dta): Oreopoulos, P. (2006). Estimating average and local average treatment effects of education when compulsory schooling laws really matter. *American Economic Review*, 96(1), 152-175. [doi: 10.1257/000282806776157641](https://www.aeaweb.org/articles?id=10.1257/000282806776157641)
  - [`headst.dta`](data/headst.dta): Ludwig, J., & Miller, D. L. (2007). Does Head Start improve children's life chances? Evidence from a regression discontinuity design. *The Quarterly journal of economics*, 122(1), 159-208. [doi:10.1162/qjec.122.1.159](https://doi.org/10.1162/qjec.122.1.159)
  - [`lee08.dta`](data/lee08.dta): Lee, D. S. (2008). Randomized experiments from non-random selection in US House elections. *Journal of Econometrics*, 142(2), 675-697. [doi:10.1016/j.jeconom.2007.05.004](https://doi.org/10.1016/j.jeconom.2007.05.004)
  - [`rcp.dta`](data/rcp.dta): Battistin, E., Brugiavini, A., Rettore, E., & Weber, G. (2009). The retirement consumption puzzle: evidence from a regression discontinuity approach. *American Economic Review*, 99(5), 2209-26. [doi:10.1257/aer.99.5.2209](https://www.aeaweb.org/articles?id=10.1257/aer.99.5.2209)
  - [`rebp.dta`](data/rebp.dta): Lalive, R. (2008). How do extended benefits affect unemployment duration? A regression discontinuity approach. *Journal of econometrics*, 142(2), 785-806. [doi:10.1016/j.jeconom.2007.05.013](https://doi.org/10.1016/j.jeconom.2007.05.013)
  