import numpy as np
from collections import OrderedDict
import scipy.interpolate  # type: ignore
import logging
from logging.handlers import QueueHandler
from gaussfit.parse.libparse.dohistogram import dohistogram  # type: ignore
# import multiprocessing
# import concurrent.futures
# import pandas as pd

def dodjdv(conn, opts, que, df, avg):
    try:
        conn.put(_dodjdv(opts, que, df, avg))
    except Exception as e:
        conn.put(e)


# def process_trace(args):
#     """
#     Processes a single trace (designed for ThreadPoolExecutor).
#     """
#     trace, opts, linx, vfilterneg, vfilterpos, avg = args
#     logger = logging.getLogger(__package__ + ".dodjdv")
# 
#     try:
#         avg_trace = avg.loc[trace]
#         spl = scipy.interpolate.UnivariateSpline(avg_trace.index, avg_trace['J'], k=5, s=opts.smooth)
#         dd = scipy.interpolate.UnivariateSpline(avg_trace.index, avg_trace['J'], k=5, s=None).derivative(2)
#         spldd = dd(vfilterpos) + (-1 * dd(vfilterneg))
# 
#         ohmic_flag = len(spldd[spldd < 0]) > 0
#         if ohmic_flag and opts.skipohmic:
#             return trace, None, None, None, None, None, ohmic_flag, 0, 0
# 
#         filtered_trace = [(row[0], spl(row[0]), row[1].J) for row in avg_trace.iterrows()]
# 
#         spls_trace, splhists_trace, spls_norm_trace, spl_normhists_trace = {}, {'spl': [], 'hist': {}}, {}, {'spl': [], 'hist': {}}
#         ndc_cut_trace, ndc_tot_trace = 0, 0
# 
#         for x in linx:
#             try:
#                 d = spl.derivatives(x)
#                 if np.isnan(d[opts.heatmapd]):
#                     logger.warning("Got NaN computing dJ/dV")
#                     continue
#                 spls_trace[x] = d[opts.heatmapd]
#                 splhists_trace['spl'].append(np.log10(abs(d[opts.heatmapd])))
#                 ndc = d[1] * (x / spl(x))
#                 spls_norm_trace[x] = ndc
#                 ndc_tot_trace += 1
#                 if 0.0 < ndc < 10.0:
#                     spl_normhists_trace['spl'].append(ndc)
#                 else:
#                     ndc_cut_trace += 1
#             except (TypeError, ValueError) as msg:
#                 logger.warning('Error computing derivative: %s', str(msg))
#                 continue
# 
#         return trace, spls_trace, splhists_trace, spls_norm_trace, spl_normhists_trace, filtered_trace, ohmic_flag, ndc_cut_trace, ndc_tot_trace
# 
#     except Exception as e:
#         logger.exception(f"Error processing trace {trace}: {e}") #trace available!
#         return trace, None, None, None, None, None, False, 0, 0
# 
# 
# 
# def new_dodjdv(opts, que, df, avg):
#     logger = logging.getLogger(__package__ + ".dodjdv")
#     logger.addHandler(QueueHandler(que)) #Correct.
#     logger.info("* * * * * * Computing dY/dX  * * * * * * * *")
# 
#     linx = np.linspace(df.V.min(), df.V.max(), 200)
#     vfilterneg, vfilterpos = (linx[-1 * opts.vcutoff < linx < 0], linx[0 < linx < opts.vcutoff]) if opts.vcutoff > 0 else (linx[linx < 0], linx[linx > 0])
# 
#     with concurrent.futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
#         futures = [executor.submit(process_trace, (trace, opts, linx, vfilterneg, vfilterpos, avg)) for trace in avg.index.levels[0]]
# 
#         spls, spls_norm, splhists, spl_normhists = OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict()
#         filtered, ohmic = [('Potential', 'Fit', 'Y')], []
#         ndc_cut_total, ndc_tot_total = 0, 0
# 
#         for future in concurrent.futures.as_completed(futures):
#             trace, spls_trace, splhists_trace, spls_norm_trace, spl_normhists_trace, filtered_trace, ohmic_flag, ndc_cut, ndc_tot = future.result()
# 
#             if ohmic_flag:
#                 ohmic.append(trace)
#             if spls_trace:  # Check for None
#                for x in linx:
#                    if x in spls_trace:
#                        if x not in spls:
#                            spls[x] = []
#                        spls[x].append(spls_trace[x])
# 
#                    if x in spls_norm_trace:
#                        if x not in spls_norm:
#                            spls_norm[x] = []
#                        spls_norm[x].append(spls_norm_trace[x])
# 
#                    if x in splhists_trace:
#                         if x not in splhists:
#                             splhists[x] = {'spl': [], 'hist':{}}
#                         splhists[x]['spl'].extend(splhists_trace['spl'])
# 
# 
#                    if x in spl_normhists_trace:
#                         if x not in spl_normhists:
#                             spl_normhists[x] = {'spl': [], 'hist':{}}
#                         spl_normhists[x]['spl'].extend(spl_normhists_trace['spl'])
# 
#             if filtered_trace:  # Check for None
#                 filtered.extend(filtered_trace)
# 
#             ndc_cut_total += ndc_cut
#             ndc_tot_total += ndc_tot
# 
# 
#     logger.info("Non-tunneling traces: %s (out of %0d)", len(ohmic), len(avg.index.levels[0]))
#     if len(ohmic) == len(avg.index.levels[0]) and opts.skipohmic:
#         logger.error("You have elected to skip all traces: disable skip non-ohmic and re-parse!")
#         return ohmic, {}, {}, {}, {}, {}
#     if ndc_tot_total:
#         logger.info("NDC values not between 0 and 10: %s (%0.2f%%)", ndc_cut_total, (ndc_cut_total / ndc_tot_total) * 100)
#     for x in splhists:
#         splhists[x]['hist'] = dohistogram(np.array(splhists[x]['spl']), label='DJDV', que=que)
#     for x in spl_normhists:
#         spl_normhists[x]['hist'] = dohistogram(np.array(spl_normhists[x]['spl']), label='NDC', que=que)
# 
#     logger.info("dJdV complete.")
#     return ohmic, spls, splhists, spls_norm, spl_normhists, filtered


def _dodjdv(opts, que, df, avg):
    '''
    Fit a spline function to X/Y data,
    compute dY/dX and normalize.
    '''
    logger = logging.getLogger(__package__+".dodjdv")
    logger.addHandler(QueueHandler(que))
    logger.info("* * * * * * Computing dY/dX  * * * * * * * *")
    ohmic = []
    linx = np.linspace(df.V.min(), df.V.max(), 200)
    if opts.vcutoff > 0:
        logger.debug('Using %s cutoff for dj/dv', opts.vcutoff)
        vfilterneg, vfilterpos = np.linspace(-1*opts.vcutoff, 0, 200), np.linspace(0, opts.vcutoff.max(), 200)
    else:
        vfilterneg, vfilterpos = np.linspace(df.V.min(), 0, 200), np.linspace(0, df.V.max(), 200)
    if opts.vcutoff > 0:
        vfilterneg, vfilterpos = linx[-1*opts.vcutoff < linx < 0], linx[0 < linx < opts.vcutoff]
    else:
        vfilterneg, vfilterpos = linx[linx < 0], linx[linx > 0]

    spls = OrderedDict()
    spls_norm = OrderedDict()
    splhists = OrderedDict()
    spl_normhists = OrderedDict()
    ndc_cut, ndc_tot = 0, 0
    filtered = [('Potential', 'Fit', 'Y')]
    for x in linx:
        spls[x] = []
        splhists[x] = {'spl': [], 'hist': {}}
        spls_norm[x] = []
        spl_normhists[x] = {'spl': [], 'hist': {}}
    for trace in avg.index.levels[0]:
        try:
            spl = scipy.interpolate.UnivariateSpline(
                avg.loc[trace].index, avg.loc[trace]['J'], k=5, s=opts.smooth)
            dd = scipy.interpolate.UnivariateSpline(
                avg.loc[trace].index, avg.loc[trace]['J'], k=5, s=None).derivative(2)
        except Exception as msg:
            # TODO: Figure out how to catch all the various scipy errors here
            logger.error('Error in derivative calulation: %s', str(msg))
            continue
        try:
            spldd = dd(vfilterpos)  # Compute d2J/dV2
            spldd += -1*dd(vfilterneg)  # Compute d2J/dV2
        except ValueError as msg:
            logger.error('Error computing second derivative: %s', str(msg))
            continue
        if len(spldd[spldd < 0]):
            # record in the index where dY/dX is < 0 within vcutoff range
            ohmic.append(trace)
            if opts.skipohmic:
                continue
        else:
            for row in avg.loc[trace].iterrows():
                # filtered is a list containing only "clean" traces
                filtered.append((row[0], spl(row[0]), row[1].J))
        err = None
        for x in spls:
            try:
                d = spl.derivatives(x)
            except (TypeError, ValueError) as msg:
                err = str(msg)
                logger.warning('Error computing derivative: %s', str(msg))
                continue
            if np.isnan(d[opts.heatmapd]):
                logger.warning("Got NaN computing dJ/dV")
                continue
            spls[x].append(d[opts.heatmapd])
            splhists[x]['spl'].append(np.log10(abs(d[opts.heatmapd])))
            ndc = d[1] * (x/spl(x))
            spls_norm[x].append(ndc)
            ndc_tot += 1
            if 0.0 < ndc < 10.0:
                # ndc values beyond this range can safely be considered artifacts of the numerical derivative
                spl_normhists[x]['spl'].append(ndc)
            else:
                # but we keep track of how many data points we toss as a sanity check
                ndc_cut += 1
        if err:
            logger.error("Error while computing derivative: %s", str(err))

    logger.info("Non-tunneling traces: %s (out of %0d)",
                     len(ohmic), len(avg.index.levels[0]))
    if (len(ohmic) == len(avg.index.levels[0])) and opts.skipohmic:
        logger.error("You have elected to skip all traces: disable skip non-ohmic and re-parse!")
        return ohmic, {}, {}, {}, {}, {}
    if ndc_tot:
        logger.info("NDC values not between 0 and 10: %s (%0.2f%%)", ndc_cut, (ndc_cut/ndc_tot)*100)
    for x in splhists:
        splhists[x]['hist'] = dohistogram(np.array(splhists[x]['spl']), label='DJDV', que=que)
        spl_normhists[x]['hist'] = dohistogram(np.array(spl_normhists[x]['spl']), label='NDC', que=que)
    logger.info("dJdV complete.")
    return ohmic, spls, splhists, spls_norm, spl_normhists, filtered
