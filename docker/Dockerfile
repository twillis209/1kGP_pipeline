FROM rocker/r-ver:latest

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    git \
    build-essential \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error \
    -n 4 \
    data.table \
    R.utils \
    magrittr \
    devtools \
    stringr \
    remotes

RUN wget -O /usr/local/bin/plink.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20240818.zip
RUN unzip /usr/local/bin/plink.zip -d /usr/local/bin

RUN wget -O /usr/local/bin/plink2.zip https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20241007.zip
RUN unzip /usr/local/bin/plink2.zip -d /usr/local/bin

ENTRYPOINT ["/bin/bash"]