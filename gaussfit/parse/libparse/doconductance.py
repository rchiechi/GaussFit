import asyncio
from gaussfit.parse.libparse.util import throwimportwarning

try:
    from scipy.stats import linregress, gmean
    import numpy as np
except ImportError as msg:
    throwimportwarning(msg)


async def dolag(xy, que):
    '''
    Make a lag plot of Y
    '''
    lag = {}
    logger = logging.getLogger(__package__+".dolag")
    logger.addHandler(QueueHandler(que))
    results = []
    
    try:
        async with asyncio.TaskGroup() as tg:
            tasks = [tg.create_task(_lagx(x, group, logger)) for x, group in xy]
        for task in tasks:
            results.append(task.result())
        for result_dict in results:
            lag.update(result_dict)  # Merge dictionaries
        
    except ExceptionGroup as eg:
        errors = set()
        for exc in eg.exceptions:
            set.add(str(exc))
        logger.warning(f"Lag calculation encountered exceptions: {','.join(errors)}")
    logger.info("Lag done.")
    return lag

async def doconductance(self):
    '''
    Find the conductance using a linear regression on the first four data points.
    '''
    self.logger.info("* * * * * * Computing G * * * * * * * *")
    self.loghandler.flush()
    tossed = 0
    voltages = []
    # Make a list of V = 0 the three values of V above and below
    for _i in range(-3, 4):
        _idx = self.avg.loc[0].index.tolist().index(0) + _i
        voltages.append(self.avg.loc[0].index.tolist()[_idx])
    

    async def _computeG(trace):
        self.SLM['G'][trace] = np.nan
        if self.opts.skipohmic and trace in self.ohmic:
            tossed += 1
            return
        _Y = []
        for _v in voltages:
            try:
                _Y.append(self.avg.loc[trace]['J'][_v])
            except KeyError:
                return
        try:
            _fit = linregress(voltages, _Y)
        except ValueError:
            self.logger.warning("Cannot compute conductance (probably because of unequal voltage steps.)")
            return
        self.logger.debug(f"G:{_fit.slope:.2E} (R={_fit.rvalue:.2f})")
        if _fit.rvalue ** 2 < self.opts.minr:
            self.logger.warn("Tossing G-value with R < %s", self.opts.minr)
            return
        if _fit.slope > self.opts.maxG:
            self.logger.warn(f"Tossing ridiculous G-value: {_fit.slope:.2E} > {self.opts.maxG}")
        self.G[trace] = _fit.slope
        self.SLM['G'][trace] = _fit.slope
    


    results = []
    
    try:
        async with asyncio.TaskGroup() as tg:
            tasks = [tg.create_task(_computeG(trace)) for trace in self.avg.index.levels[0]]
        for task in tasks:
            results.append(task.result())
        # for result_dict in results:
        #     lag.update(result_dict)  # Merge dictionaries
        
    except ExceptionGroup as eg:
        errors = set()
        for exc in eg.exceptions:
            set.add(str(exc))
        logger.warning(f"G calculation encountered exceptions: {','.join(errors)}")
    
    Gavg = gmean(list(self.G.values()), nan_policy='omit')
    self.logger.info("Average conductance: %.2E", Gavg)
    return Gavg
