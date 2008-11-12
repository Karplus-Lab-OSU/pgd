# subscripter is a way to make a property subscriptable
# this allows you to reference properties by index
class Subscripter():
    def __init__(self, key, parent):
        self.key = key
        self.parent = parent
        #add this instance to the parent. doing this here
        #makes defining subscriptor instance simpler because
        #you only need to specify the key once
        parent.__dict__[key] = self

    def __getitem__(self, i):
        return self.parent.__dict__['%s_%d' % (self.key, i)]

    def __setitem__(self, i, val):
        self.parent.__dict__['%s_%d' % (self.key, i)] = val

# Base foo contains functions we'd like to add in a normal way
# It can also contain normal Django fields
class BaseFoo(models.Model):

    normal_field = models.IntegerField()

    # this won't get called when the object is constructed
    # something about multiple inheritence in python
    # it must be called manually
    def __init__(self):
        models.Model.__init__(self)
        Subscripter('bar', self)

    # make class abstract
    class Meta:
        abstract = True


#add __module__ to the dictionary, used by django to do its magic
#this is the path to the file containing your module
mydict = {'__module__':'pgd.tasks.models'}

#dynamically add some django fields to the dictionary
for i in range(5):
    mydict['bar_%d' % i] = models.IntegerField()


#Call Types to dynamically create the Foo type inheriting from the base class. 
Foo = type('Foo', (BaseFoo,), mydict)
