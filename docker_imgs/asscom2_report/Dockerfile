FROM rocker/tidyverse
RUN install2.r prettydoc DT ape
#RUN ln -s /usr/local/lib/R/site-library/littler/examples/render.r /usr/local/bin

# Encapsulate the latest report R-markdown template:
#RUN cd ~ && pwd && wget --no-check-certificate https://raw.githubusercontent.com/cmkobel/assemblycomparator2/master/scripts/genomes_to_report_v2.Rmd
RUN cd /usr/local/bin && pwd && wget --no-check-certificate https://raw.githubusercontent.com/cmkobel/assemblycomparator2/master/scripts/genomes_to_report_v2.Rmd



