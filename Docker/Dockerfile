### image: AGAMEMNON v-0.1.0 ###
### This dockerfile is based on the one created by
### Titus Brown (available at https://github.com/ctb/2015-docker-building/blob/master/salmon/Dockerfile) ###
FROM ubuntu:16.04

ENV PACKAGES git gcc make g++ libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf vim wget time bzip2 \
    nano unzip libpthread-stubs0-dev

ENV PUFFERFISH_VERSION 0.9.0

### don't modify things below here for version updates etc. ###

WORKDIR /home


RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean


RUN wget https://cmake.org/files/v3.9/cmake-3.9.5-Linux-x86_64.sh && \
	mkdir /opt/cmake && \
	sh cmake-3.9.5-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
	ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake && \
	rm cmake-3.9.5-Linux-x86_64.sh

   
RUN curl -o /usr/local/bin/jq http://stedolan.github.io/jq/download/linux64/jq && \
	chmod +x /usr/local/bin/jq


# python and required packages
RUN yes | apt-get install software-properties-common && \
	add-apt-repository ppa:deadsnakes/ppa && \
	apt-get update
RUN yes | apt-get install python3.6

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	printf "\nyes\n\nyes\n" | bash Miniconda3-latest-Linux-x86_64.sh

RUN yes | /root/miniconda3/bin/conda install -c bioconda -c conda-forge snakemake && \
	rm Miniconda3-latest-Linux-x86_64.sh


### TwoPaCo installation ###
RUN wget https://github.com/01org/tbb/releases/download/2018_U4/tbb2018_20180411oss_lin.tgz && \
	tar zxf tbb2018_20180411oss_lin.tgz && \
	rm tbb2018_20180411oss_lin.tgz && \
	cd tbb2018_20180411oss && \
	cp -r include /usr/local && \
	cp -r lib/intel64/gcc4.7/* /usr/local/lib && \
	cd .. && \
	rm -rf tbb2018_20180411oss

RUN git clone https://github.com/fataltes/TwoPaCo.git

RUN cd TwoPaCo && \
	git checkout pufferize && \
	mkdir build && \
	cd build && \
	cmake ../src && \
	make

RUN git clone https://github.com/COMBINE-lab/pufferfish.git


### sdsl-lite installation ###
RUN	cd pufferfish && \
	git clone https://github.com/simongog/sdsl-lite.git && \
	cd sdsl-lite && \ 
	./install.sh ../


### Pufferfish installation ###
RUN cd pufferfish && \
	git checkout AGAMEMNON && \
	mkdir build && \
	cd build && \
	cmake ../ && \
	make


### HISAT2 binaries ###
RUN wget http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-Linux_x86_64.zip && \
	unzip hisat2-2.1.0-Linux_x86_64.zip && \
	rm -rf hisat2-2.1.0-Linux_x86_64.zip


### AGAMEMNON github repo ###
RUN mkdir AGAMEMNON && \
	cd AGAMEMNON && \
	git clone  && \
	cd agamemnon && \
	git checkout develop && \
	cp -a /home/AGAMEMNON/agamemnon/. /home/AGAMEMNON/ && \
	rm -rf /home/AGAMEMNON/agamemnon


### Create directories - Move files ###
RUN mkdir /home/AGAMEMNON/binaries && \
	mkdir /home/AGAMEMNON/binaries/pufferfish && \
	mkdir /home/AGAMEMNON/binaries/TwoPaco && \
	mkdir /home/AGAMEMNON/binaries/hisat2-2.1.0 && \
	cp -a /home/pufferfish/. /home/AGAMEMNON/binaries/pufferfish && \
	cp -a /home/TwoPaCo/. /home/AGAMEMNON/binaries/TwoPaCo && \
	cp -a /home/hisat2-2.1.0/. /home/AGAMEMNON/binaries/hisat2-2.1.0 && \
	rm -rf pufferfish  && \
	rm -rf TwoPaCo  && \
	rm -rf hisat2-2.1.0


### Give permission to scripts ###
RUN chmod +x /home/AGAMEMNON/scripts/downloadTaxonomy.sh && \
	chmod +x /home/AGAMEMNON/scripts/editYaml.py && \
	chmod +x /home/AGAMEMNON/scripts/getTaxonomy.py && \
	chmod +x /home/AGAMEMNON/scripts/indexSize.py && \
	chmod +x /home/AGAMEMNON/scripts/single_cell-demultiplex.py && \
	chmod +x /home/AGAMEMNON/scripts/accessions2taxids.py && \
	chmod +x /home/AGAMEMNON/scripts/taxonomy.sh

RUN sed -i 's#"pufferfish":"/home/pufferfish",#"pufferfish":"/home/AGAMEMNON/binaries/pufferfish",#g' /home/AGAMEMNON/binaries/pufferfish/config.json && \
	sed -i 's#"twopaco":"/home/TwoPaCo",#"twopaco":"/home/AGAMEMNON/binaries/TwoPaCo",#g' /home/AGAMEMNON/binaries/pufferfish/config.json

RUN ldconfig

### END ###
