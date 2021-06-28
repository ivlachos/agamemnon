# image: AGAMEMNON
# This dockerfile is based on the one created by
# Titus Brown (available at https://github.com/ctb/2015-docker-building/blob/master/salmon/Dockerfile)
#FROM insilicodb/bioconductor 


# Example shiny app docker file
# https://blog.sellorm.com/2021/04/25/shiny-app-in-docker/

# get shiny server and R from the rocker project
#FROM rocker/shiny:4.0.5

# system libraries
# Try to only install system libraries you actually need
# Package Manager is a good resource to help discover system deps
#RUN apt-get update && apt-get install -y \
#    libcurl4-gnutls-dev \
#    libssl-dev
  

# install R packages required 
# Change the packages list to suit your needs
#RUN R -e 'install.packages(c(\
#              "shiny", \
#              "shinydashboard", \
#              "ggplot2" \
#            ), \
#            repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23"\
#          )'


# copy the app directory into the image
#COPY ./shiny-app/* /srv/shiny-server/

# run app
#CMD ["/usr/bin/shiny-server"]


FROM ubuntu:16.04

#COPY --from=0 / /

ENV PACKAGES git gcc make g++  liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf vim wget time bzip2 \
	libcurl4-gnutls-dev libssl-dev

ENV AGAMEMNON_VERSION 0.1.0

# AGAMEMNON will be available in /home/agamemnon
# pufferfish executable will be available in /home/pufferfish/build/src


WORKDIR /home


ADD https://cmake.org/files/v3.9/cmake-3.9.5-Linux-x86_64.sh /cmake-3.9.5-Linux-x86_64.sh
RUN mkdir /opt/cmake && \
	sh /cmake-3.9.5-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
	ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean


#RUN apt install apt-transport-https software-properties-common && \
#	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
#	add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' && \
#	apt update && \
#	apt install r-base && \
#	R --version

# install R packages required 
# Change the packages list to suit your needs
#RUN R -e 'install.packages(c(\
#              "shiny", \
#              "shinydashboard", \
#              "ggplot2" \
#            ), \
#            repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23"\
#          )'




# python and required packages
RUN yes | apt-get install software-properties-common && \
	add-apt-repository ppa:deadsnakes/ppa && \
	apt-get update

RUN yes | apt-get install python3.6
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	printf "\nyes\n\nyes\nyes\n" | bash Miniconda3-latest-Linux-x86_64.sh

RUN echo ${PATH} && \
	yes | /root/miniconda3/bin/conda install scipy && \
	yes | /root/miniconda3/bin/conda install pandas && \
	yes | /root/miniconda3/bin/conda install scikit-learn


# system libraries
# Try to only install system libraries you actually need
# Package Manager is a good resource to help discover system deps
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev

  
# from https://blog.sellorm.com/2021/04/25/shiny-app-in-docker/
# install R packages required 
# Change the packages list to suit your needs
#RUN R -e 'install.packages(c(\
#              "shiny", \
#              "shinydashboard", \
#              "ggplot2" \
#            ), \
#            repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23"\
#          )'


# copy the app directory into the image
#COPY ./shiny-app/* /srv/shiny-server/

# run app
#CMD ["/usr/bin/shiny-server"]


# snakemake
#ENV BASH_ENV ~/.bashrc
SHELL ["/bin/bash", "-c"]

ENV PATH=/root/miniconda3/bin:$PATH

RUN conda init bash
RUN conda install -c conda-forge r-base

RUN conda install -n base -c conda-forge mamba && \
	conda install -n base libiconv && \
	echo "conda activate base" >> ~/.bashrc && \
	mamba create -c conda-forge -c bioconda -n snakemake snakemake

RUN echo "conda activate snakemake" >> ~/.bashrc
#RUN conda activate snakemake && \
#	snakemake --help


RUN git clone https://github.com/ivlachos/agamemnon.git && \
	chmod +x agamemnon/scripts/*

#RUN wget https://github.com/ivlachos/agamemnon/releases/download/AGAMEMNON-v0.1.0/AGAMEMNON-v0.1.0.tar.gz && \
#	tar -xvzf AGAMEMNON-v0.1.0.tar.gz && \
#	mv AGAMEMNON-v0.1.0 agamemnon && \
#	rm -r agamemnon/binaries

RUN git clone https://github.com/COMBINE-lab/pufferfish.git

# Pufferfish && Cedar
RUN cd pufferfish && \	
	git checkout tags/v1.1.0 && \
	mkdir build && \
	cd build && \
	cmake ../ && \
	make 

# HISAT2
RUN wget https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download && \
	mv download hisat2-2.1.0-Linux_x86_64.zip && \
	unzip hisat2-2.1.0-Linux_x86_64.zip 



# move pufferfish to agamemnon required hierarchy
RUN mkdir -p /home/agamemnon/binaries/pufferfish && \
	cp /home/pufferfish/build/src/pufferfish /home/agamemnon/binaries/pufferfish/ && \
	chmod +x /home/agamemnon/binaries/pufferfish/pufferfish 

# move cedar to agamemnon required hierarchy
RUN mkdir -p /home/agamemnon/binaries/cedar && \
	cp /home/pufferfish/build/src/cedar /home/agamemnon/binaries/cedar/ && \
	chmod +x /home/agamemnon/binaries/cedar/cedar

RUN mv /home/hisat2-2.1.0 /home/agamemnon/binaries/hisat2
#	chmod +x /home/agamemnon/binaries/hisat2/


RUN ldconfig
