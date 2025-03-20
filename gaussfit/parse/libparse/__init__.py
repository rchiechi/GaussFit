import platform
from .dorect import dorect
from .findtraces import findtraces
from .dodjdv import dodjdv
from .findmin import findmin, old_findmin
from .dummies import dummyListener, dummyPopen
from .doconductance import doconductance
from .doslm import doslm
from .readfiles import readfiles

if platform.system() in ('Linux', 'Darwin', 'Windows'):
    from .dolag import doLagMultiprocess as doLag
    from .findsegments import findSegmentsMultiprocess as findSegments
else:
    from .dolag import doLag
    from .findsegments import findSegments
