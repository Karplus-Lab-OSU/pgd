
def PGDContextProcessor(request):
    from django.conf import settings
    return {
        'MEDIA_URL': settings.MEDIA_URL,
        'SITE_ROOT': settings.SITE_ROOT,
        'session':request.session
    }
