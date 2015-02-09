FROM centos:7

MAINTAINER OSU OSL support@osuosl.org

EXPOSE 8000

RUN yum install -y \
  cairo \
  gcc \
  gcc-c++ \
  git \
  libcairo-devel \
  libffi \
  libffi-devel \
  mysql \
  mysql-devel \
  pycairo \
  python-setuptools \
  python-devel

RUN easy_install pip

# Copy and configure pgd
WORKDIR /opt/pgd
# Copy requirements.txt separately for better caching
COPY ./requirements.txt /opt/pgd/requirements.txt
RUN pip install -r requirements.txt
# NB: copying the settings file is not a good idea when using volumes!
COPY . /opt/pgd/
RUN cp /opt/pgd/settings.py.dist /opt/pgd/settings.py

CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]
