FROM ubuntu:latest

ENV DEBIAN_FRONTEND noninteractive

ADD https://bootstrap.pypa.io/get-pip.py /tmp/get-pip.py
ADD https://raw.githubusercontent.com/dceoy/print-github-tags/master/print-github-tags /usr/local/bin/print-github-tags
ADD . /tmp/tmber

RUN set -e \
      && ln -sf bash /bin/sh \
      && ln -s python3 /usr/bin/python

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        autoconf ca-certificates curl g++ gcc gnuplot libbz2-dev \
        libcurl4-gnutls-dev liblzma-dev libncurses5-dev libssl-dev libz-dev \
        make pbzip2 pigz pkg-config python3 python3-distutils \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -eo pipefail \
      && chmod +x /usr/local/bin/print-github-tags \
      && print-github-tags --release --latest samtools/htslib \
        | xargs -I{} curl -SL \
          https://github.com/samtools/htslib/releases/download/{}/htslib-{}.tar.bz2 \
          -o /tmp/htslib.tar.bz2 \
      && tar xvf /tmp/htslib.tar.bz2 -C /usr/local/src --remove-files \
      && mv /usr/local/src/htslib-* /usr/local/src/htslib \
      && cd /usr/local/src/htslib \
      && autoheader \
      && autoconf \
      && ./configure \
      && make \
      && make install

RUN set -eo pipefail \
      && print-github-tags --release --latest --tar arq5x/bedtools2 \
        | xargs -i curl -SL {} -o /tmp/bedtools2.tar.gz \
      && tar xvf /tmp/bedtools2.tar.gz -C /usr/local/src --remove-files \
      && mv /usr/local/src/bedtools2-* /usr/local/src/bedtools2 \
      && cd /usr/local/src/bedtools2 \
      && make \
      && make install

RUN set -e \
      && /usr/bin/python3 /tmp/get-pip.py \
      && pip install -U --no-cache-dir pip /tmp/tmber \
      && rm -rf /tmp/get-pip.py /tmp/tmber

ENTRYPOINT ["/usr/local/bin/tmber"]
