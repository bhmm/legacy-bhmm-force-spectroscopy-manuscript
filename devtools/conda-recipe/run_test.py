import os
import nose
import bhmm


path = bhmm.__path__[0]
path = os.path.join(path, '..')
print ("executing nose from new path:", path)

os.chdir(path)
args = ['bhmm', '--nocapture', '--verbosity=2', '--with-doctest']
nose.main(args)
