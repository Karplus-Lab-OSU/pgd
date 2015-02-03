FROM centos:7

MAINTAINER OSU OSL support@osuosl.org

EXPOSE 8000

#Sane configuration defaults
ENV DATABASE_ENGINE mysql
ENV DATABASE_NAME pgd
ENV DATABASE_USER pgd
ENV DATABASE_PASSWORD pgd
ENV DATABASE_HOST mysql

RUN yum install -y mysql mysql-devel gcc gcc-c++ git libcairo-devel cairo pycairo python-setuptools python-devel libffi libffi-devel

RUN easy_install pip

# Copy and configure pgd
WORKDIR /opt/pgd
# Copy requirements.txt separately for better caching
COPY ./requirements.txt /opt/pgd/requirements.txt
RUN pip install -r requirements.txt
COPY . /opt/pgd
RUN cp ./settings.py.dist ./settings.py

CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]
