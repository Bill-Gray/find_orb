FROM ubuntu:latest
MAINTAINER Bill Gray 

SHELL ["/bin/bash", "-c"]
ENV PATH="${PATH}:~/bin"

# Update applications and install OS-level dependencies
RUN apt-get update \
	&& apt-get upgrade -y \
	&& apt-get install g++ make wget libncurses5-dev libcurl4-openssl-dev git -y

# Download and install find_orb and dependencies
RUN mkdir software && cd software \
	&& git clone https://github.com/Bill-Gray/find_orb.git \
	&& cd find_orb \
	&& /bin/bash DOWNLOAD.sh -d .. \
	&& /bin/bash INSTALL.sh -d .. -u
