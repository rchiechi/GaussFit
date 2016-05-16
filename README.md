# GaussFit

## Synopsis
This is a set of python scripts to parse I/V data, perform statistical analyses, generate plots and output parsed data as CSV files.
Most functions and options are available via a `tkinter` gui that you can call with `GaussFit.py -G`.

The `Stats.py` part is very much a work in progress, but it should, at the very least compute p-values using eac file as a degree of freedom.
The `autonobs` command will try to guess time intervals and use those as degrees of freedom. The commands for `GaussFit.py` and `Stats.py`
are identical (with a few special options for `Stats.py`); launching with `-G` will get you a GUI in either case.

### Input file format
GaussFit will try to guess the range of your *I/V* sweeps, but it works best with `0 -> + -> 0 -> – ->0` (i.e., three zeros per trace).


### Caveats
Input files containing *I/V* sweeps with different voltage step-sizes will cause problems with Gaussian fits unless you have a ton of data.
Say you have 1,000 traces with 0.5 V steps and 100 with 0.25 V steps. You will end up with 1,100 points in the histogram at each 0.5 V
interval and only 100 at each 0.25 V interval, thus the reliatbility of the fits (and size of the error bars) will vary considerably.

## Dependencies
```
python 3.4+
pandas
colorama
numpy
scipy
matplotlib
```

### Usage
You can just call the GUI with `GaussFit.py -G` and work from there.

To quickly plot your data without writing files:
```
GaussFit.py -n (-G) <files>
```
To write your data quickly and without plotting:
```
GaussFit.py -o <output name> (-G) <files>
```
To do both:
```
GaussFit.py -o <output name> -p (-G) <files>
```
