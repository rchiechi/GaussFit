import os
import logging
import csv

logger = logging.getLogger(__package__)


def get_header_offset(fn):
    end_of_header = 0
    with open(fn) as fh:
        _l = fh.readline().strip()
        while _l:
            if _l == '***End_of_Header***':
                end_of_header = fh.tell()
                break
            _l = fh.readline().strip()
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


def __get_raw_data_true_temp(opts):
    logging.info('True temp not yet implemented.')
    return __get_raw_data_dTn(opts)


def __get_raw_data_dTn(opts, logq, raw_data):
    for dt in opts.dTn:
        raw_data[f'DT{dt}'] = {'data': [], 'labels': []}
    # for i in range(2*len(opts.dTn)):
    #     raw_data.append([])
    #     # For each dT creat two lists in raw_data: one dV and one dT
    # dts_present = []
    # for fn in opts.in_files:
    #     for dt in opts.dTn:
    #         if f'dt{dt}' in fn.lower():
    #             if dt not in dts_present:
    #                 dts_present.append(dt)
    # if len(dts_present) != opts.dTn:
    #     logging.error('All ΔT values are not present in input file names.')
    #     return raw_data
    # headers = X_Value	Raw Voltage (uV)	Corrected Voltage (uV)	Top T	BottomT	Delta-T (Deg.C)	Seebeck (uV/K)	Elapsed Time (s)	Comment
    for fn in opts.in_files:
        delimiter = get_delimiter(fn)
        header_offset = get_header_offset(fn)
        with open(fn) as fh:
            _key = 'DT0'
            for dt in opts.dTn:
                if f'dt{dt}' in fn.lower():
                    _key = f'DT{dt}'
            fh.seek(header_offset)  # seeks to the end of the LabView header
            if not raw_data[_key]['labels']:
                raw_data[_key]['labels'] = fh.readline().split(delimiter)
            content = csv.reader(fh, delimiter=delimiter)
            for row in content:
                i = -1
                for __column in row:
                    i += 1  # Otherwise ValueError will skip incrementing i
                    try:
                        while i >= len(raw_data[_key]['data']):
                            raw_data[_key]['data'].append([])
                        raw_data[_key]['data'][i].append(float(__column))
                        # if float(row[1]) < 200 and float(row[1]) > -200:  # skips bad data
                        #     raw_data[i].append(float(row[1]))
                        #     if dictionary["TrueTemp"]:
                        #         raw_data[i + dictionary["dT_middle"]].append([float(row[5])])
                    except ValueError:
                        logq.put(f'Error converting data from {_key} to float.')
                        continue
    # raw data now looks like [[v1],[v2],[v3],..., [t1],[t2],[t3],...]
    # if not opts.truetemp:     # if there is no dT data collected, remove those empty lists
    #     del raw_data[dictionary["dT_middle"]:]
 