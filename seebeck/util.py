import os
import logging
import csv
import numpy as np

logger = logging.getLogger(__package__)


def get_header_offset(fn):
    i = 0
    end_of_header = 0
    with open(fn) as fh:
        for _l in fh:
            i += 1
            if _l.strip() == '***End_of_Header***':
                end_of_header = i
    # print(f'end_of_head: {end_of_header}')
    return end_of_header


def get_delimiter(fn):
    delim_map = {'Space': ' ', 'Tab': '\t'}
    delimiter = delim_map['Tab']
    with open(fn) as fh:
        for _l in fh:
            _s = _l.strip().split()
            if not _s:
                continue
            if _s[0] == 'Separator':
                _key = _s[-1]
                if _key in delim_map:
                    delimiter = delim_map[_key]
                else:
                    logging.error(f'Unknown delimiter {_key}, using default.')
    return delimiter


def get_raw_data(opts, logq, raw_data):
    if opts.truetemp:
        return __get_raw_data_true_temp(opts, logq, raw_data)
    else:
        return __get_raw_data_dTn(opts, logq, raw_data)


def __get_raw_data_true_temp(opts, logq, raw_data):
    logging.info('True temp not yet implemented.')
    return __get_raw_data_dTn(opts, logq, raw_data)


def __get_raw_data_dTn(opts, logq, raw_data):
    for dt in opts.dTn:
        raw_data[f'DT{dt}'] = {'data': [], 'labels': [], 'dt': 0.0}
    for fn in opts.in_files:
        logq.put('Parsing %s' % os.path.basename(fn))
        delimiter = get_delimiter(fn)
        header_offset = get_header_offset(fn)
        with open(fn, newline='') as fh:
            _key = 'DT0'
            for dt in opts.dTn:
                if f'dt{dt}' in fn.lower():
                    _key = f'DT{dt}'
            content = csv.reader(fh, delimiter=delimiter)
            for i in range(header_offset):
                next(content)
            if not raw_data[_key]['labels']:
                raw_data[_key]['labels'] = next(content)
            elif raw_data[_key]['labels'] != next(content):
                logq.put("Warning: column headers do not match between files!")
            for row in content:
                i = -1
                for __column in row:
                    i += 1  # Otherwise ValueError will skip incrementing i
                    try:
                        while i >= len(raw_data[_key]['data']):
                            raw_data[_key]['data'].append([])
                        if not __column:
                            raw_data[_key]['data'][i].append(np.NaN)
                        else:
                            if abs(float(__column)) > opts.cuttoff_to_toss and i < 6:
                                logq.put(f'[get_raw_data] Tossing ridiculous value {__column}')
                                __column = 0
                            raw_data[_key]['data'][i].append(float(__column))
                    except ValueError:
                        logq.put(f'Error converting data from {_key} to float: {__column}')

    for dt in opts.dTn:
        while len(raw_data[f'DT{dt}']['data']) < opts.col_to_parse + 1:
            raw_data[f'DT{dt}']['data'].append([])
            raw_data[f'DT{dt}']['labels'].append('')

    for _key in raw_data:
        if len(raw_data[_key]['data']) > 5:
            raw_data[_key]['dt'] = np.mean(raw_data[_key]['data'][5])  # TODO: make sure this is the dT column
        else:
            logq.put(f"Could not convert real ΔT for {_key}.")
        # for i in range(len(raw_data[_key]['data'])):
        #     raw_data[_key]['data'][i] = np.array(raw_data[_key]['data'][i])
        #     if i == 5:  # TODO: make sure this is the dT column
        #         raw_data[_key]['dt'] = np.mean(raw_data[_key]['data'][i])
