InCell
===========

InCell is a Shiny application for the analysis of InCell data.  


Run
===========
Dependencies: `shiny`, `shinyIncubator`, `shinythings`, `plotrix`.

```s
if(!require('devtools')) install.packages('devtools')
install.packages(c('shiny', 'plotrix'), repos="http://cran.rstudio.com")
devtools::install_github('dcurrier/shinythings')
devtools::install_gitbuh('rstudio/shiny-incubator')
```

Run the app:
```s
shiny::runGitHub("InCell", "dcurrier", launch.browser=T)
```