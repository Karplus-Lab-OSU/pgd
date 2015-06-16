from django.views.generic import TemplateView
from django.conf import settings


class ExtraContextTemplateView(TemplateView):
    extra_context = {
        'SITE_ROOT': settings.SITE_ROOT,
        'PGD_VERSION': settings.PGD_VERSION,
        'ROOT': settings.SITE_ROOT,
        'DATA_VERSION': settings.DATA_VERSION,
        'GOOGLE_ID': settings.GOOGLE_ID
    }

    def get_context_data(self, *args, **kwargs):
        context = super(ExtraContextTemplateView, self).get_context_data(*args, **kwargs)
        context.update(self.extra_context)
        return context


class ReferencesView(ExtraContextTemplateView):
    template_name = "references.html"


class ContactUsView(ExtraContextTemplateView):
    template_name = "contactus.html"


class NewsView(ExtraContextTemplateView):
    template_name = "news.html"


class WelcomeView(ExtraContextTemplateView):
    template_name = "welcome.html"
