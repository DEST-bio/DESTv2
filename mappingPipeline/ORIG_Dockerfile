FROM ubuntu:18.04

# packaging dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        apt-utils \
        build-essential \
        autotools-dev \
        autoconf \
        automake \
        software-properties-common \
        libncurses5-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        openjdk-8-jdk \
        maven \
        wget \
        sed \
        gzip \
        perl \
        gnupg2 \
        libopenblas-dev \
        libbz2-dev \
        liblzma-dev \
        zlib1g-dev \
        bash-completion \
        vim \
        less \
        curl \
        git \
        pigz \
        ocaml-nox

# install Python
RUN apt-get install -y --no-install-recommends \
        python3-dev \
        python3-setuptools \
        python3-pip \
        python2.7-dev \
        python-setuptools \
        python-pip \
        python-pysam && \
        pip install setuptools --upgrade && \
        pip install pip --upgrade && \
        pip3 install setuptools --upgrade && \
        pip3 install pip --upgrade


# install R
ARG DEBIAN_FRONTEND=noninteractive
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
#    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' && \
#    apt install -y --no-install-recommends \
#        apt-transport-https software-properties-common \
RUN apt install -y --no-install-recommends \
    r-base \
    r-recommended && \
    rm -rf /var/lib/apt/lists/*

#Have not used
# install bam2fastq
#ARG BAM2FASTQ_VERSION=1.1.0
#WORKDIR /opt
#RUN wget https://gsl.hudsonalpha.org/static/software/bam2fastq-${BAM2FASTQ_VERSION}.tgz && \
#    tar xzf bam2fastq-${BAM2FASTQ_VERSION}.tgz && \
#    mkdir bam2fastq && \
#    mv bam2fastq-${BAM2FASTQ_VERSION} bam2fastq/${BAM2FASTQ_VERSION} && \
#    cd bam2fastq/${BAM2FASTQ_VERSION} && \
#    make && \
#    cd ../.. && \
#    rm bam2fastq-${BAM2FASTQ_VERSION}.tgz
#ENV PATH="/opt/bam2fastq/${BAM2FASTQ_VERSION}:${PATH}"


# install bwa
ARG BWA_VERSION=0.7.15
WORKDIR /opt
RUN echo $BWA_VERSION && \
    mkdir bwa && \
    cd bwa && \
    wget https://github.com/lh3/bwa/archive/v${BWA_VERSION}.tar.gz && \
    tar xzf v${BWA_VERSION}.tar.gz && \
    mv bwa-${BWA_VERSION} ${BWA_VERSION} && \
    cd ${BWA_VERSION} && \
    make && \
    cd .. && \
    rm v${BWA_VERSION}.tar.gz
ENV PATH="/opt/bwa/${BWA_VERSION}:${PATH}"


# install HTSlib (may not need all this)
ARG HTSLIB_VERSION=1.9
WORKDIR /opt
RUN wget https://github.com/samtools/htslib/archive/${HTSLIB_VERSION}.tar.gz && \
    tar xzf ${HTSLIB_VERSION}.tar.gz && \
    cd htslib-${HTSLIB_VERSION} && \
    autoheader && \
    autoreconf && \
    ./configure --prefix /opt/htslib/${HTSLIB_VERSION} && \
    make && \
    make install && \
    cd .. && \
    rm ${HTSLIB_VERSION}.tar.gz && \
    rm -rf htslib-${HTSLIB_VERSION}
ENV PATH="/opt/htslib/${HTSLIB_VERSION}/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/htslib/${HTSLIB_VERSION}/lib:${LD_LIBRARY_PATH}"


# install SAMTools
# LDFLAGS="-Wl,-R$/opt/htslib/${HTSLIB_VERSION}/lib && \

ARG SAMTOOLS_VERSION=1.9
WORKDIR /opt
RUN wget https://github.com/samtools/samtools/archive/${SAMTOOLS_VERSION}.tar.gz && \
    tar xzf ${SAMTOOLS_VERSION}.tar.gz && \
    cd samtools-${SAMTOOLS_VERSION} && \
    autoheader && \
    autoreconf -Wno-syntax && \
    ./configure \
       --prefix /opt/samtools/${SAMTOOLS_VERSION} \
       --with-htslib=/opt/htslib/${HTSLIB_VERSION} \
       LDFLAGS="-Wl,-R/opt/htslib/${HTSLIB_VERSION}/lib" && \
    make && \
    make install && \
    cd /opt && \
    rm -rf ${SAMTOOLS_VERSION}.tar.gz && \
    rm -rf samtools-${SAMTOOLS_VERSION}
ENV PATH="/opt/samtools/${SAMTOOLS_VERSION}/bin:${PATH}"


# install Cutadapt
ARG CUTADAPT_VERSION=2.3
RUN python3 -m pip install cutadapt==${CUTADAPT_VERSION}

ENV JRE_HOME=${JAVA_HOME}/jre

# install Picard Tool
ARG PICARD_VERSION=2.20.2
WORKDIR /opt
RUN mkdir -p /opt/picard/${PICARD_VERSION} && \
    cd /opt/picard/${PICARD_VERSION} && \
    wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar
ENV CLASSPATH="/opt/picard/${PICARD_VERSION}/picard.jar:$CLASSPATH"
ENV PICARD="/opt/picard/${PICARD_VERSION}/picard.jar"

#install FastQC
ARG FASTQC_VERSION=0.11.9
WORKDIR /opt
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip && \
    unzip fastqc_v${FASTQC_VERSION}.zip && \
    rm fastqc_v${FASTQC_VERSION}.zip && \
    cd FastQC && \
    chmod 755 fastqc
WORKDIR /opt
ENV PATH="/opt/FastQC:${PATH}"

#install Trimmomatic
#ARG TRIM_VERSION=0.36
#WORKDIR /opt
#RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIM_VERSION}.zip && \
#    unzip Trimmomatic-${TRIM_VERSION}.zip && \
#    rm Trimmomatic-${TRIM_VERSION}.zip && \
#    cd Trimmomatic-${TRIM_VERSION}
#WORKDIR /opt
#ENV TRIMMOMATIC="/opt/Trimmomatic-${TRIM_VERSION}/trimmomatic-${TRIM_VERSION}.jar"

#install BBmap

ARG BBMAP_VERSION=38.80

WORKDIR /opt
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz/download -O BBMap_${BBMAP_VERSION}.tar.gz && \
    tar -xvzf BBMap_${BBMAP_VERSION}.tar.gz && \
    rm BBMap_${BBMAP_VERSION}.tar.gz
ENV PATH="/opt/bbmap:${PATH}"

# install GenomeAnalysisTK

#ARG GATK_VERSION=3.8.1 Changing this will only work for version 4 of GATK, which is the RUN command associated with the github releases link commented out below

WORKDIR /opt
#RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
#    unzip gatk-${GATK_VERSION}.zip && \
#    rm gatk-${GATK_VERSION}.zip && \
#    mv gatk-${GATK_VERSION} gatk
RUN wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 && \
    tar -xvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 && \
    mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef gatk && \
    rm GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
ENV CLASSPATH="/opt/gatk/GenomeAnalysisTK.jar:$CLASSPATH"
ENV GATK="/opt/gatk/GenomeAnalysisTK.jar"

ENV PATH="/usr/lib/jvm/java-1.8.0-openjdk-amd64/bin:${PATH}"

#Install snape-pooled

WORKDIR /opt
RUN git clone https://github.com/EmanueleRaineri/snape-pooled.git && \
    cd snape-pooled && \
    make snape-pooled

ENV PATH="/opt/snape-pooled:${PATH}"

#Clone the repo

RUN git clone https://github.com/DEST-bio/DEST_freeze1.git

#Setup the hologenome
RUN /opt/DEST_freeze1/mappingPipeline/scripts/hologene_gen.sh
RUN bwa index /opt/hologenome/holo_dmel_6.12.fa && \
    samtools faidx /opt/hologenome/raw/D_melanogaster_r6.12.fasta && \
    samtools faidx /opt/hologenome/holo_dmel_6.12.fa && \
    java -jar $PICARD CreateSequenceDictionary \
              R=/opt/hologenome/holo_dmel_6.12.fa \
              O=/opt/hologenome/holo_dmel_6.12.dict

ARG CACHE_BUST=06:03:2020:1800

RUN python3 /opt/DEST_freeze1/mappingPipeline/scripts/PickleRef.py \
    --ref /opt/hologenome/raw/D_melanogaster_r6.12.fasta  \
    --output /opt/hologenome/raw/D_melanogaster_r6.12.fasta.pickled > pickled_ref.out

WORKDIR /opt
RUN rm -rf /opt/DEST_freeze1/ && \
    git clone https://github.com/DEST-bio/DEST_freeze1.git

ENV PATH="/opt/DEST_freeze1/mappingPipeline/scripts:${PATH}"

#Command to execute on "docker run" of this image
ENTRYPOINT ["/bin/bash", "/opt/DEST_freeze1/mappingPipeline/scripts/fq_to_sync_pipeline.sh"]

#Default parameters passed to the above function call. Can be overridden by supplying parameters after image name in "docker run"
CMD ["--help"]
