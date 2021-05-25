FROM rocker/tidyverse
RUN install2.r prettydoc DT
#RUN ln -s /usr/local/lib/R/site-library/littler/examples/render.r /usr/local/bin

# Encapsulate the latest report R-markdown template:
RUN cd ~ && wget --no-check-certificate https://raw.githubusercontent.com/cmkobel/assemblycomparator2/master/scripts/genomes_to_report_v2.Rmd
#RUN mv genomes_to_report_v2.Rmd /usr/local/bin
