
def PGDContextProcessor(request):
    from django.conf import settings

    template_dict = {
        'MEDIA_URL': settings.MEDIA_URL,
        'SITE_ROOT': settings.SITE_ROOT,
        'session':request.session,
    }

    try:
        search = request.session['search']
        template_dict['count'] = search.querySet().count()
    except KeyError:
        # if search hasn't been completed yet discard error
        pass

    return template_dict
