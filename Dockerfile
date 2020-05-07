#    This file is part of iss-microglia.
#    Copyright (C) 2020  Emir Turkes
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Emir Turkes can be contacted at emir.turkes@eturkes.com

FROM rocker/rstudio:3.6.3

LABEL maintainer="Emir Turkes emir.turkes@eturkes.com"

COPY user-settings /home/rstudio/.rstudio/monitored/user-settings/

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
    zlib1g-dev \
 && Rscript \
    -e "install.packages('conflicted')" \
    -e "install.packages('peakRAM')" \
    -e "install.packages('rmarkdown')" \
    -e "install.packages('rprojroot')" \
    -e "install.packages('knitr')" \
    -e "install.packages('data.table')" \
    -e "install.packages('DT')" \
    -e "install.packages('plotly')" \
    -e "install.packages('BiocManager')" \
    -e "BiocManager::install('scater')" \
    -e "BiocManager::install('scran')" \
 && apt-get clean \
 && rm -Rf \
    /var/lib/apt/lists/ \
    /tmp/downloaded_packages/ \
    /tmp/*.rds
