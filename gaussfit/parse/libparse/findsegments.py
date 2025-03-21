import numpy as np
import logging
from logging.handlers import QueueHandler
from gaussfit.parse.libparse.dohistogram import dohistogram

def findSegments(conn, opts, que, df):
    '''
    Break out each trace into four segments of
    0V -> Vmax, Vmax -> 0, 0V -> Vmin, Vmin -> 0V.
    '''
    logger = logging.getLogger(__package__+".findsegments")
    logger.addHandler(QueueHandler(que))
    logger.info("* * * * * * Finding segments   * * * * * * * *")
    # TODO set num_segments in opts
    # NOTE this is a crude hack because I forgot how Pandas works
    if opts.ycol < 0:
        logger.warning("Parsing segments when all columns are parsed may produce weird results!")
    try:
        if df.V.value_counts()[0] % int(opts.segments-1) != 0:
            logger.warning("Dataset does not seem to have %s segments.", int(opts.segments))
    except KeyError:
        logger.warning("Could not segment data by 0's.")
        return {}
    error = False
    logger.info("Breaking out traces by segments of 0V -> Vmin/max.")
    segments = {}
    segments_combined = {}
    segments_combined_nofirst = {}
    nofirsttrace = {}
    max_column_width = 0
    for _fn in df.index.levels[0]:
        _seg = None
        _trace = None
        if opts.ycol > 0:
            _n_traces = int(df.V.value_counts()[0] / 3)
            try:
                if df.loc[_fn]['V'][0] != 0.0:
                    logger.warning("J/V didn't start at 0V for %s", _fn)
            except KeyError:
                logger.error("\n\n\t>>>Cannot parse %s<<<\n\tCheck that J/V curves are not trucated or remove file and parse again.\n", _fn)
                continue
        else:
            _n_traces = len(df.index.levels[0])

        _idx = [_idx for _idx in df.loc[_fn].index]
        for _n, _i in enumerate(_idx):
            J = df.loc[_fn]['J'][_i]
            V = df.loc[_fn]['V'][_i]
            _Vmax = df.loc[_fn]['V'].max()
            _Vmin = df.loc[_fn]['V'].min()
            _last_V, _next_V = None, None

            if _n > 0:
                _last_V = df.loc[_fn]['V'][_idx[_n-1]]
            if _n+1 < len(_idx):
                _next_V = df.loc[_fn]['V'][_idx[_n+1]]

            if _last_V is None:  # First trace
                _seg = 0
                _trace = 0
            elif _next_V is None:  # Last trace
                pass
            elif V == 0 and _next_V != 0:  # Crossing zero
                _seg += 1
            elif V in (_Vmax, _Vmin) and _next_V not in (_Vmax, _Vmin) and V != 0:  # Turnaround at ends
                _seg += 1

            if _seg in (-1, opts.segments):  # New trace starts with new segment
                # print("_seg is 4\nLast V: %s\nThis V: %s\nNext V: %s\n" % (_last_V, V, _next_V))
                _seg = 0
                _trace += 1

            if _trace > _n_traces:
                logger.warning("Parsing trace %s, when there should only be %s",
                               _trace, _n_traces)

            if _seg not in segments:
                segments[_seg] = {}
                segments_combined[_seg] = {}
                segments_combined_nofirst[_seg] = {}
            if _trace not in segments[_seg]:
                segments[_seg][_trace] = {}

            if V not in segments[_seg][_trace]:
                segments[_seg][_trace][V] = []
                segments_combined[_seg][V] = []
                segments_combined_nofirst[_seg][V] = []
            _last_V = V
            segments[_seg][_trace][V].append(J)
            segments_combined[_seg][V].append(J)

            if V not in nofirsttrace:
                nofirsttrace[V] = []
            if _trace > 0:
                nofirsttrace[V].append(J)
                segments_combined_nofirst[_seg][V].append(J)
            if V != 0:
                if len(nofirsttrace[V]) > max_column_width:
                    max_column_width = len(nofirsttrace[V])
    if len(segments.keys()) != opts.segments:
        error = True
        logger.error('Expected %i segments, but found %i!', opts.segments, len(segments.keys()))
    else:
        logger.info('Found %s segments.', len(segments.keys()))

    segmenthists = {}
    segmenthists_nofirst = {}
    for _seg in segments:
        if _seg not in segmenthists:
            segmenthists[_seg] = {'combined': {}}
            segmenthists_nofirst[_seg] = {'combined': {}}
        for _trace in segments[_seg]:
            logger.debug('Segment: %s, Trace: %s', _seg, _trace)
            if _trace not in segmenthists[_seg]:
                segmenthists[_seg][_trace] = {}
                if _trace > 0:
                    segmenthists_nofirst[_seg][_trace] = {}
            for _V in segments[_seg][_trace]:
                if _V not in segmenthists[_seg][_trace]:
                    segmenthists[_seg][_trace][_V] = {}
                segmenthists[_seg][_trace][_V] = dohistogram(
                                                             np.array([np.log10(abs(_j)) for _j in segments[_seg][_trace][_V]]),
                                                             label='segment', que=que, opts=opts)
                if _trace > 0:
                    segmenthists_nofirst[_seg][_trace][_V] = segmenthists[_seg][_trace][_V]
        for _V in segments_combined[_seg]:
            if _V not in segmenthists[_seg]['combined']:
                segmenthists[_seg]['combined'][_V] = {}
            segmenthists[_seg]['combined'][_V] = dohistogram(
                                                             np.array([np.log10(abs(_j)) for _j in segments_combined[_seg][_V]]),
                                                             label='segment', que=que, opts=opts)
        for _V in segments_combined_nofirst[_seg]:
            if _V not in segmenthists_nofirst[_seg]['combined']:
                segmenthists_nofirst[_seg]['combined'][_V] = {}
            segmenthists_nofirst[_seg]['combined'][_V] = dohistogram(
                                                                     np.array(
                                                                              [np.log10(abs(_j))
                                                                               for _j in segments_combined_nofirst[_seg][_V]]),
                                                                     label='segment', que=que, opts=opts)
    for _V in nofirsttrace:
        if _V == 0:
            if len(nofirsttrace[_V]) > max_column_width:
                nofirsttrace[_V] = np.zeros(max_column_width)
                logger.warning("Setting J = 0 for all V = 0 in nofirsttrace.")
        nofirsttrace[_V] = np.array(nofirsttrace[_V])
    logger.info("Findsegments done.")
    conn.put((error, segmenthists, segmenthists_nofirst, nofirsttrace))
