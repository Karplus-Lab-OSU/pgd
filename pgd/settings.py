# Django settings for pgd project.

from decouple import config
import os
PROJECT_ROOT = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))

DEBUG = config('DEBUG', default=False, cast=bool)
TEMPLATE_DEBUG = DEBUG

ADMINS = (
    # ('Your Name', 'your_email@domain.com'),
)

MANAGERS = ADMINS

# At this time only MySQL is supported.  When other databases are
# supported, these settings will be refactored.
import dj_database_url
DATABASES = {
    'default': dj_database_url.config(default=config('DATABASE_URL', ''))
}
import sys
if 'test' in sys.argv:
    DATABASES['default']['OPTIONS'] = {'init_command': 'SET foreign_key_checks=0'}

# Warning: this value is insecure!
ALLOWED_HOSTS = config('ALLOWED_HOSTS', default='*')

# prefix used for the site.  ie. http://myhost.com/<SITE_ROOT>/
# for the django standalone server this should be ""
# for apache this is the url the site is mapped to, probably /pgd
SITE_ROOT = config('SITE_ROOT', default='')

# absolute path to the docroot of this site
# DOC_ROOT = config('DOC_ROOT', default='')
DOC_ROOT = PROJECT_ROOT

# Google analytics ID. Enter the full id, as in: UA-xxxxxx-x
GOOGLE_ID = config('GOOGLE_ID', default='UA-xxxxxx-x')

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# If running in a Windows environment this must be set to the same as your
# system time zone.
TIME_ZONE = config('TIME_ZONE', default='America/Chicago')

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = config('LANGUAGE_CODE', default='en-us')

SITE_ID = config('SITE_ID', default=1, cast=int)

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = config('USE_I18N', default=True, cast=bool)

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale.
USE_L10N = True

# If you set this to False, Django will not use timezone-aware datetimes.
USE_TZ = True

# Absolute path to the directory that holds media.
# Example: "/home/media/media.lawrence.com/"
MEDIA_ROOT = config('MEDIA_ROOT', default='%s/media' % DOC_ROOT)

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash if there is a path component (optional in other cases).
# Examples: "http://media.lawrence.com", "http://example.com/media/"
MEDIA_URL = config('MEDIA_URL', default='%s/%s/' % (SITE_ROOT, 'media'))

# URL prefix for admin media -- CSS, JavaScript and images. Make sure to use a
# trailing slash.
# Examples: "http://foo.com/media/", "/media/".
ADMIN_MEDIA_PREFIX = config('ADMIN_MEDIA_PREFIX', default='/media/')

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/home/media/media.lawrence.com/static/"
STATIC_ROOT = config('STATIC_ROOT', default='%s/static' % DOC_ROOT)

# URL prefix for static files.
# Example: "http://media.lawrence.com/static/"
STATIC_URL = config('STATIC_URL', default='%s/%s/' % (SITE_ROOT, 'static'))

# Additional locations of static files
STATICFILES_DIRS = (
    # Put strings here, like "/home/html/static" or "C:/www/django/static".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    'django.contrib.staticfiles.finders.DefaultStorageFinder',
)

# Make this unique, and don't share it with anybody.
SECRET_KEY = config('SECRET_KEY')

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
)

TEMPLATE_CONTEXT_PROCESSORS = (
    'django.contrib.auth.context_processors.auth',
    'django.core.context_processors.debug',
    'django.core.context_processors.media',
    'context_processors.PGDContextProcessor',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
)

ROOT_URLCONF = config('ROOT_URLCONF', default='pgd.urls')

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'pgd.wsgi.application'

TEMPLATE_DIRS = (
    '%s/templates' % DOC_ROOT
)

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.staticfiles',
    'django.contrib.admin',
    'registration',
    'pgd_core',
    'pgd_search',
    'pgd_splicer',
)

# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error when DEBUG=False.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        }
    },
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'filters': ['require_debug_false'],
            'class': 'django.utils.log.AdminEmailHandler'
        },
    },
    'loggers': {
        'django.request': {
            'handlers': ['mail_admins'],
            'level': 'ERROR',
            'propagate': True,
        },
    }
}

# PGD Specific settings
QUERY_LIMIT = config('QUERY_LIMIT', default=50000000, cast=int)
SEGMENT_SIZE = config('SEGMENT_SIZE', default=10, cast=int)
DATA_VERSION = config('DATA_VERSION', default='testing')
PGD_VERSION = config('PGD_VERSION', default='1.0.2')

# Django registration
ACCOUNT_ACTIVATION_DAYS = config('ACCOUNT_ACTIVATION_DAYS', default=5, cast=int)

# Email configuration
EMAIL_BACKEND = config('EMAIL_BACKEND', default='django.core.mail.backends.smtp.EmailBackend')
EMAIL_HOST = config('EMAIL_HOST', default='smtp.osuosl.org')
EMAIL_PORT = config('EMAIL_PORT', default='25')
DEFAULT_FROM_EMAIL = config('DEFAULT_FROM_EMAIL', default='registration@pgd.science.oregonstate.edu')
SERVER_EMAIL = config('SERVER_EMAIL', default='pgd@pgd.science.oregonstate.edu')
LOGIN_REDIRECT_URL = '%s/search/' % SITE_ROOT

# FTP settings
PDB_FTP_HOST = config('PDB_FTP_HOST', default='ftp.wwpdb.org')
PDB_REMOTE_DIR = config('PDB_REMOTE_DIR', default='/pub/pdb/data/structures/divided/pdb')
PDB_LOCAL_DIR = config('PDB_LOCAL_DIR', default='./pdb')
PDB_TMP_DIR = config('PDB_TMP_DIR', default='./tmp')
