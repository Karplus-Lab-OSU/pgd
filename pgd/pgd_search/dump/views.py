from django.http import HttpResponse
from DataDump import dumpSearch

"""
render the results of the search as a TSV (tab separated file)
and return it to the user as a download
"""
def dataDump(request):
    #from DataDump import dumpSearch
    response = HttpResponse(mimetype="text/tab-separated-values")
    response['Content-Disposition'] = 'attachment; filename="data.tsv"'

    dumpSearch(request.session['search'], response)

    return response
