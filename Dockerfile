################################
# Docker file for App-gc version 1.0
# Base on Ubuntu 20.04 and yufree/xcmsrocker
################################

FROM yufree/xcmsrocker:latest

#Developer/maintainer
MAINTAINER Apiwat Sangphukieo

ADD / /working_directory

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget
	
RUN R -e 'BiocManager::install(c( "GCalignR","filesstrings","msm","expm","pracma","doSNOW","erah","shinythemes","shinyFiles","shiny","plotly","htmlwidgets","fontawesome"))'



# Inspired by suchja/wine
# Install some tools required for creating the image
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		curl \
		unzip \
		ca-certificates

# Install wine and related packages
RUN dpkg --add-architecture i386 \
		&& apt-get update \
		&& apt-get install -y --no-install-recommends \
				wine \
				wine32 \
		&& rm -rf /var/lib/apt/lists/*

# Use the latest version of winetricks
RUN curl -SL 'https://raw.githubusercontent.com/Winetricks/winetricks/master/src/winetricks' -o /usr/local/bin/winetricks \
		&& chmod +x /usr/local/bin/winetricks

# Get latest version of mono for wine
RUN mkdir -p /usr/share/wine/mono \
	&& curl -SL 'http://sourceforge.net/projects/wine/files/Wine%20Mono/$WINE_MONO_VERSION/wine-mono-$WINE_MONO_VERSION.msi/download' -o /usr/share/wine/mono/wine-mono-$WINE_MONO_VERSION.msi \
	&& chmod +x /usr/share/wine/mono/wine-mono-$WINE_MONO_VERSION.msi

# Get MSPepSearch 32x version and extract
RUN curl -LO https://chemdata.nist.gov/download/peptide_library/libraries/software/current_releases/MSPepSearch/2019_02_22_MSPepSearch_x32.zip

RUN unzip 2019_02_22_MSPepSearch_x32.zip


# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/ && \
    chown shiny:shiny /var/lib/shiny-server

WORKDIR /working_directory

EXPOSE 3838

#COPY shiny-server.sh /usr/bin/shiny-server.sh

RUN mkdir -p /var/log/shiny-server && \
	chown shiny.shiny /var/log/shiny-server 

CMD ["if [ "$APPLICATION_LOGS_TO_STDOUT" != "false" ]; then ; exec xtail /var/log/shiny-server/ & ;fi"]