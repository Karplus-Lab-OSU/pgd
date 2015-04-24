FROM centos:7

MAINTAINER OSU OSL support@osuosl.org

EXPOSE 8000

# Dockerfiles do not support "here documents"
RUN echo "[osuosl]" > /etc/yum.repos.d/osuosl.repo
RUN echo "name=OSUOSL Repository 6 - x86_64" >> /etc/yum.repos.d/osuosl.repo
RUN echo "baseurl=http://ftp.osuosl.org/pub/osl/repos/yum/6/x86_64" >> /etc/yum.repos.d/osuosl.repo
RUN echo "enabled=1" >> /etc/yum.repos.d/osuosl.repo
RUN echo "gpgcheck=0" >> /etc/yum.repos.d/osuosl.repo

RUN yum -y update && yum -y install \
    epel-release

RUN yum -y update && yum -y install \
  bzip2 \
  cairo-devel \
  gcc \
  gcc-c++ \
  libffi \
  libffi-devel \
  mysql \
  mysql-devel \
  nodejs \
  npm \
  osuosl-dssp \
  python-devel \
  python-setuptools \
  tar

RUN npm -g install phantomjs

RUN easy_install pip

# Copy and configure pgd
WORKDIR /opt/pgd
RUN mkdir /opt/pgd/media /opt/pgd/static
# Copy requirements.txt separately for better caching
COPY ./requirements.txt /opt/pgd/requirements.txt
RUN pip install -r requirements.txt
# NB: copying the settings file is not a good idea when using volumes!
COPY . /opt/pgd/
RUN cp /opt/pgd/pgd/settings.py.dist /opt/pgd/pgd/settings.py

CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]
