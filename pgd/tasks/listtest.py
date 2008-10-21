def test(list):
    list.append('bar')

list = []
list.append('foo')

test(list)

for item in list:
    print item

