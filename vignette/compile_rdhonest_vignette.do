

* This .do file compiles the Stata Markdown file rdhonest_stata_vignette.stmd to html by first calling markstat to compile to markdown, and then calling pandoc to convert to html.  The extra step is necessary, since the markstat html option doesn't allow arbitrary command line options to be passed to pandoc.
*  - requires the Stata Markdown package (https://data.princeton.edu/stata/markdown, https://ideas.repec.org/c/boc/bocode/s458401.html), pandoc and the pandoc-secnos extension for pandoc (https://github.com/tomduck/pandoc-secnos)
*  - markstat_css_style_chunk.html contains the markstat.css file between <style> and </style>.  Calling with -H markstat_css_style_chunk.html instead of -c markstat.css mimics the markstat command behavior of calling pandoc without css and then including the css code in the file directly (pandoc adds some style information only when called without the css option).


markstat using rdhonest_stata_vignette.stmd, markdown

shell pandoc rdhonest_stata_vignette.md -s -o rdhonest_stata_vignette.html --filter pandoc-secnos --citeproc --mathjax --no-highlight -H markstat_css_style_chunk.html --number-sections

