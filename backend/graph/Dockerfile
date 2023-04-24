FROM ubuntu:20.04

# Put your own Dockfile commands here

WORKDIR /usr/src/app

ENV TZ="America/NewYork" \
    DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get -qq upgrade \
  && apt-get install -y build-essential \
  && apt-get install -y wget \
  && apt-get -q install -y libpq-dev \
  && apt-get install -y libc6-dev \
  && apt-get install libpthread-stubs0-dev \
  && apt-get -qq install libtcmalloc-minimal4 \
  && ln -s /usr/lib/libtcmalloc_minimal.so.4 /usr/lib/libtcmalloc_minimal.so \
  && mkdir dependencies \
  && wget https://github.com/Kitware/CMake/releases/download/v3.26.3/cmake-3.26.3-linux-x86_64.sh \
      -q -O /tmp/cmake-install.sh \
      && chmod u+x /tmp/cmake-install.sh \
      && mkdir /opt/cmake-3.26.3 \
      && /tmp/cmake-install.sh --skip-license --prefix=/opt/cmake-3.26.3 \
      && rm /tmp/cmake-install.sh \
      && ln -s /opt/cmake-3.26.3/bin/* /usr/local/bin \
  && wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz \
      && tar xvf boost_1_81_0.tar.gz --directory dependencies \
      && rm boost_1_81_0.tar.gz \
  && wget https://github.com/CrowCpp/Crow/releases/download/v1.0%2B5/crow-v1.0+5.tar.gz \
      && mkdir dependencies/crow \
      && tar xvf "crow-v1.0+5.tar.gz" --directory dependencies/crow \
      && rm "crow-v1.0+5.tar.gz" \
#  && apt-get update \
  && apt-get install -y ca-certificates lsb-core \
  && sh -c 'echo "deb http://apt.postgresql.org/pub/repos/apt $(lsb_release -cs)-pgdg main" > /etc/apt/sources.list.d/pgdg.list' \
  && wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | apt-key add - \
  && apt-get update \
  && apt-get install -y postgresql postgresql-contrib \
  && wget https://github.com/jtv/libpqxx/archive/refs/tags/6.4.5.tar.gz \
      && tar xvf 6.4.5.tar.gz --directory dependencies \
      && rm 6.4.5.tar.gz \
  && rm -rf /var/lib/apt/lists/*

ADD include ./include
COPY ./main.cpp ./
COPY ./CMakeLists.txt ./

RUN echo $(ls -1 .)
RUN cmake .
RUN make

CMD ["./backend", "8080", "osm-northeast"]