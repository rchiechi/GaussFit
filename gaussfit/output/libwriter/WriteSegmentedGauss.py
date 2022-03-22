import os
import csv
import logging
import numpy as np
from scipy.special import stdtrit

logger = logging.getLogger('output')


def WriteSegmentedGauss(self, key=None):
    '''Write histograms of values of J broken out by segment to catch
    hysteretic behavior without smearing it out.'''

    if not self.segments:
        logger.warning("No segments found.")
        return

    if key == 'nofirst':
        _segments = self.segments_nofirst
        _label = 'Segment_NoFirst'
    else:
        _segments = self.segments
        _label = 'Segment'

    for segment in _segments:
        rows = {}
        _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Gauss_%s_%s.txt" % (_label, str(segment+1)))
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = ["Potential (V)"]
            _maxtrace = 0
            for trace in _segments[segment]:
                # TODO: Fix this hack
                if not isinstance(trace, int):
                    continue
                _maxtrace += 1
                headers += ["Log|J|",
                            "Standard Deviation",
                            "Standard Error of the Mean",
                            "%s%% confidence interval" % (100*(1-self.opts.alpha))]
                for x in _segments[segment][trace]:
                    _hist = _segments[segment][trace][x]
                    if x not in rows:
                        rows[x] = []
                    rows[x].append("%0.4f" % _hist['mean'])
                    rows[x].append("%0.4f" % _hist['std'])
                    _sem = float(_hist['std'])/np.sqrt(self.opts.degfree - 1 or 1)
                    rows[x].append("%0.4f" % _sem)
                    _t_val = _sem * stdtrit(self.opts.degfree - 1 or 1, 1 - self.opts.alpha)
                    rows[x].append("%0.4f" % _t_val)

            writer.writerow(headers)
            _V = list(rows.keys())
            _V.sort()
            for x in _V:
                while len(rows[x]) < _maxtrace * 3:
                    rows[x] += ['-', '-', '-']
                    logger.warning('Filling columns for segment %i, V=%s to match %s traces.', segment, x, _maxtrace)
                writer.writerow(["%0.4f" % x]+rows[x])

    # TODO: Don't just repeat the whole code block
    for segment in _segments:
        rows = {}
        _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Gauss_%s_Combined_%s.txt" % (_label, str(segment+1)))
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = ["Potential (V)",
                       "Log|J|",
                       "Standard Deviation",
                       "Standard Error of the Mean",
                       "%s%% confidence interval" % (100*(1-self.opts.alpha))]
            for x in _segments[segment]['combined']:
                _hist = _segments[segment]['combined'][x]
                if x not in rows:
                    rows[x] = []
                rows[x].append("%0.4f" % _hist['mean'])
                rows[x].append("%0.4f" % _hist['std'])
                _sem = float(_hist['std'])/np.sqrt(self.opts.degfree - 1 or 1)
                rows[x].append("%0.4f" % _sem)
                _t_val = _sem * stdtrit(self.opts.degfree - 1 or 1, 1 - self.opts.alpha)
                rows[x].append("%0.4f" % _t_val)

            writer.writerow(headers)
            _V = list(rows.keys())
            _V.sort()
            for x in _V:
                writer.writerow(["%0.4f" % x]+rows[x])
