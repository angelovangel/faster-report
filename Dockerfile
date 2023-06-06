FROM rstudio/r-base:4.2.2-focal
# fixed to match renv

LABEL author="Angel Angelov <aangeloo@gmail.com>"
LABEL description="Docker image containing the requirements for faster-report"

RUN apt-get update && apt-get install -y \
ksh procps libxt-dev libssl-dev libxml2-dev libfontconfig1-dev
# libxt-dev is required to solve the segfault error caused by cairoVersion() in R

# setup faster and fastkmers for linux
RUN wget -P bin https://github.com/angelovangel/faster/releases/download/v0.1.4/x86_64_linux_faster && \
mv bin/x86_64_linux_faster bin/faster && \
chmod 755 bin/faster

RUN wget -P bin https://github.com/angelovangel/fastkmers/releases/download/v0.1.3/fastkmers && \
chmod 755 bin/fastkmers

# setup renv
ENV RENV_VERSION 0.15.2
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

RUN mkdir -p renv
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.dcf renv/settings.dcf
RUN R --vanilla -s -e 'renv::restore()'

# the only purpose of this container is to run the report generation
COPY faster-report.R /bin/faster-report.R
RUN chmod 755 /bin/faster-report.R
ENTRYPOINT [ "faster-report.R" ]
