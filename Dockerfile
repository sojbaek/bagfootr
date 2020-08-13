FROM rocker/verse:latest

MAINTAINER   Songjoon Baek <sojbaek@hotmail.com>

WORKDIR /usr/src
RUN apt-get update && apt-get -y install tcl8.6-dev tk8.6-dev
RUN apt-get install -qqy x11-apps

#Firefox is to be installed to make sure X11 graphics libraries will be installed
#Samtools
RUN apt-get update && apt-get -y upgrade && \
	apt-get install -y build-essential wget \
		libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
	

RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar jxf samtools-1.9.tar.bz2 && \
	rm samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure --prefix $(pwd) && \
	make
ENV PATH=${PATH}:/usr/src/samtools-1.9 

# required R packages
RUN R -e "install.packages('hash')"
RUN R -e "install.packages('parallel')"
RUN R -e "install.packages('digest')"
RUN R -e "install.packages('data.table')"
RUN R -e "install.packages('Cairo')"
RUN R -e "install.packages('aplpack')"
RUN R -e "install.packages('devtools')"

RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm9')"
RUN R -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')"
RUN R -e "BiocManager::install('BSgenome.Rnorvegicus.UCSC.rn5')"
RUN R -e "BiocManager::install('BSgenome.Rnorvegicus.UCSC.rn6')"
RUN R -e "BiocManager::install('Rsamtools')"
RUN R -e "BiocManager::install('GenomeInfoDb')"
RUN R -e "BiocManager::install('GenomicRanges')"
RUN R -e "BiocManager::install('IRanges')"


# copy BagFoot related files

WORKDIR /home/bagfoot
RUN wget https://sourceforge.net/projects/bagfootr/files/example_data.tar.gz
RUN tar -xvf /home/bagfoot/example_data.tar.gz

RUN wget https://sourceforge.net/projects/bagfootr/files/bagfoot_prep_example.R
RUN wget https://sourceforge.net/projects/bagfootr/files/bagfoot_run_example.R


# bagfoot installation
RUN R -e "devtools::install_github('sojbaek/bagfootr')"


RUN rm /home/bagfoot/example_data.tar.gz
COPY Dockerfile /opt

#  The BaGFoot requires at least 16GB of memory.  To set the memory requirement of the Docker contatiner
#  See  https://stackoverflow.com/questions/44533319/how-to-assign-more-memory-to-docker-container/44533437
#  For a Mac system, run "XQuartz" for a X-windows support which is requires to run a plot function  
#	
#  To build:
#   sudo docker build -t bagfoot .
#
#  To run BaGFoot in a Linux Docker envrionment (https://skandhurkat.com/post/x-forwarding-on-docker/)
#    #The below lines stes the path to the .Xauthority file
#    XAUTH=$HOME/.Xauthority
#    touch $XAUTH
#    # Start the docker image
#    sudo docker run --tty --interactive --network=host --env DISPLAY=$DISPLAY --volume $XAUTH:/root/.Xauthority bagfoot bash
#
#  To run in a Mac envrionment (https://forums.docker.com/t/x11-forwarding-issues-with-release-10/11252/3 )
#  Run XQuartz
# 	a. Update preferences ‘Security’ tab - turn on 'Allow connection from network clients’
# 	b. Restart XQuartz and then check to see that it is listening on port 6000:
#    lsof -i :6000
#
#  	Get your local machine’s IP:
#   	 IP=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}') && echo "My IP is: $IP"
#  	Allow the local machine to talk to XQuartz
# 	 xhost + ${IP} 
# 	 Run your docker host:
#
#    sudo docker run -ti --rm  -e DISPLAY=${IP}:0 -e XAUTHORITY=/.Xauthority -v /tmp/.X11-unix:/tmp/.X11-unix -v ~/.Xauthority:/.Xauthority   bagfoot bash
#    In R, type 'source("bagfoot_run_example.R")' to run the example program.
