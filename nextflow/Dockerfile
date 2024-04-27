FROM debian:bullseye-slim

LABEL image.author.name "Yuntao Yang"
LABEL image.author.email "uestcyangyuntao@gmail.com"

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    pigz \
    python3 \
    python3-pip \
    make \
    vim-common \
    default-jre \
    && rm -rf /var/lib/apt/lists/*

# Install Cutadapt
RUN pip3 install cutadapt

# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip \
    && unzip fastqc_v0.12.1.zip \
    && chmod +x FastQC/fastqc \
    && mv FastQC /opt/ \
    && ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc \
    && rm fastqc_v0.12.1.zip

# Install Trim Galore!
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz \
	&& tar -xzvf 0.6.10.tar.gz \
	&& mv TrimGalore-0.6.10/trim_galore /usr/local/bin/ \
	&& rm -rf TrimGalore-0.6.10 0.6.10.tar.gz

# Install MultiQC
RUN pip3 install multiqc

# Install STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz \
    && tar -xzvf 2.7.11b.tar.gz \
	&& cd STAR-2.7.11b/source \
    	&& make STAR \
	&& cp STAR /usr/local/bin/ \
	&& cd ../../ \
	&& rm -rf 2.7.11b.tar.gz STAR-2.7.11b/
