from pyamf.remoting.client import RemotingService

gw = RemotingService('http://127.0.0.1:8000/search/gateway/')
service = gw.getService('myservice')

print service.echo('Hello World!')
