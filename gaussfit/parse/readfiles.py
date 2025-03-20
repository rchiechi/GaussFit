import asyncio
import logging
from pathlib import Path
import pandas as pd
from .libparse.util import async_getfilechecksum

async def dedupe_files(input_files, logger):
    results = {}
    async with asyncio.TaskGroup() as tg:
        tasks = {tg.create_task(async_getfilechecksum(fn)): fn for fn in input_files}
    
    checksums = {}
    for task, fn in tasks.items():
        try:
            checksums[task.result()] = Path(fn)
        except Exception as e:
            # Log or handle specific exceptions
            logger.error(f"Error parsing {fn}")
    
    return checksums

def parse_two_columns(file_paths, opts, logger):
    df = pd.DataFrame()
    frames = {}
    for checksum, f in file_paths.items():
        with f.open('rb') as fh:
            try:
                headers = fh.readline().split(bytes(opts.delim, encoding=opts.encoding))
                headers = list(map(lambda x: str(x, encoding=opts.encoding), headers))
            except UnicodeDecodeError:
                logger.warning("Encountered an illegal unicode character in headers.")
        try:
            if headers:
                x, y = headers[opts.xcol].strip(), headers[opts.ycol].strip()
                frames[f] = pd.read_csv(f, sep=opts.delim, encoding=opts.encoding,
                                        usecols=(x, y))[[x, y]]
                frames[f].rename(columns={x: 'V', y: 'J'}, inplace=True)
            elif opts.X > opts.Y:
                raise pd.errors.ParserError("xcol cannot be greater than ycol without column headers.")
            else:
                frames[f] = pd.read_csv(f, sep=opts.delim,
                                        usecols=(opts.xcol,
                                                 opts.ycol),
                                        names=('V', 'J'), header=0)
        except OSError as msg:
            logger.warning("Skipping %s because %s", f, str(msg))
        except pd.errors.ParserError as msg:
            logger.warning("Skipping malformatted %s because %s", f, str(msg))
    return pd.concat(frames)

def parse_all_columns(file_paths, logger):
    df = pd.DataFrame()
    frames = {}
    for checksum, f in file_paths.items():
        try:
            df = pd.read_csv(f, sep=opts.delim,
                              index_col=opts.xcol,
                              header=0,
                              error_bad_lines=False,
                              warn_bad_lines=False)
            i = 0
            for col in df:
                frames['%s_%.2d' % (f, i)] = pd.DataFrame({'V': df.index, 'J': df[col]})
                i += 1
        except OSError as msg:
            logger.warning("Skipping %s because %s", f, str(msg))
    return pd.concat(frames)

async def readfiles(opts, **kwargs):
    '''Walk through input files and parse
    them into attributes '''
    logger = kwargs.get('logger', logging.getLogger(__package__+".readfiles"))
    input_files = opts.in_files
    if not isinstance(input_files, list):
        input_files = [input_files]
    
    file_paths = await dedupe_files(input_files, logger)

    logger.debug('Parsing %s', ', '.join( fn.stem for fn in file_paths.values() ))
    if opts.ycol > -1:
        logger.info("Parsing two columns of data (X=%s, Y=%s).", opts.xcol + 1, opts.ycol + 1)
        df = parse_two_columns(file_paths, opts, logger)
    else:
        logger.info("Parsing all columns of data.")
        df = parse_all_columns(file_paths, logger)

    if opts.xrange > 0 and not df.empty:
        logger.info(f"Pruning x-axis to +/- {opts.xrange}.")
        df = df[df.V > -1 * abs(opts.xrange)]
        df = df[df.V < abs(opts.xrange)]
    return df