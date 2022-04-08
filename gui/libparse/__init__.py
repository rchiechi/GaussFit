import platform
from .ParseThread import ParseThread
from .GUIParse import GUIParse
from .GUIPlot import GUIPlot
# if platform.system() == 'Linux':
#     from .GUIPlot import GUIPlotMultiprocessing as GUIPlot
# else:
#     from .GUIPlot import GUIPlot
