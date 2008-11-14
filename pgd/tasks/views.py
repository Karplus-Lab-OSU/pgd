from django.http import HttpResponse
from django.template import RequestContext, Context, loader
from django.conf import settings
from django.utils import simplejson

from taskserver import TaskClient

taskClient = TaskClient()

def task_processor(request):

    global taskClient

    if taskClient == None:
        taskClient = TaskClient()

    if not taskClient.connected:
        taskClient.connect(5800)

    return {}


def showtasks(request):

    task_processor(request)
    json_list = taskClient.listTasks()
    tasks = simplejson.loads(json_list)

    t = loader.get_template('tasks.html')
    c = RequestContext(request, {
        'MEDIA_URL': settings.MEDIA_URL,
        'tasks': tasks
    }, [task_processor])

    return HttpResponse(t.render(c))


def taskprogress(request):

    c = RequestContext(request, {
        'MEDIA_URL': settings.MEDIA_URL
    }, [task_processor])

    return HttpResponse(taskClient.progress(), mimetype='application/javascript')


def starttask(request):

    key = request.POST['key']

    c = RequestContext(request, {
        'MEDIA_URL': settings.MEDIA_URL
    }, [task_processor])

    return HttpResponse(taskClient.start(key), mimetype='application/javascript')


