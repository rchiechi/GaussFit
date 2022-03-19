# GaussFit

## Synopsis

This is a set of python scripts to parse raw I/V data, perform statistical analyses, generate plots and output parsed data as CSV files.
Most functions and options are available via a `tkinter` GUI that you can call with `GaussFit.py -G`. The outputs of `GaussFit.py` have been tested extensively and, while the code is a mess, the results can be considered reliable provided all warning messages are heeded.

The `Stats.py` part is very much a work in progress, but it should, at the very least compute p-values using each file as a degree of freedom.
The `autonobs` command will try to guess time intervals and use those as degrees of freedom, but it relies on metadata in the input files and has only been tested with output containing this metadata in a very specific format.

The commands for `GaussFit.py` and `Stats.py` are identical (with a few special options for `Stats.py`); launching with `-G` will get you a GUI in either case.

Although we have done some testing of the statistically functions, do not rely on this software without using some robust, external test of the veracity of the output!

### Input file format

GaussFit will try to guess the range of your _I/V_ sweeps, but it works best with `0 -> + -> 0 -> – ->0` (i.e., three zeros per trace). The default assumes an input file format that is specific to our setups, but you can specify which columns to parse. (Though it will always refer to the x-axis as voltage and the y-axis as current-density when it generates plots.) It will try to determine if the input data contains `0 -> + -> 0 - -> 0` traces or just linear sweeps (which it will call "AFM data"). It assembles full forward-reverse sweeps for computing rectification ratios, lag plots, computing `DI/DV` plots etc. and will warn if there are leftover data. It will also warn if the data are excessively noisy or flat by complaining about "non-tunneling traces". The defaults for the spline fits used to compute the first and second derivatives work well for our input data, but if your data have small (>0.1 V) step-sizes, you should tweak the settings.

### Outputs

GaussFit writes everything to text files in the "parsed" sub-directory for auditing. It displays plots of data using `matplotlib` as well as `gnuplot` input files. All transformations and calculations are performed on individual `I/V` sweeps before fitting to a Gaussian. Given a sufficiently large dataset, it can interpolate transition voltage values because they are derived from histograms of spline fits, but those should be compared against the standard method, which is only capable of finding transition voltages corresponding to real input bias values.

### Caveats

Input files containing _I/V_ sweeps with different voltage step-sizes will cause problems with Gaussian fits unless you have a ton of data. Say you have 1,000 traces with 0.5 V steps and 100 with 0.25 V steps. You will end up with 1,100 points in the histogram at each 0.5 V
interval and only 100 at each 0.25 V interval, thus the reliability of the fits (and size of the error bars) will vary considerably.

Currently GaussFit only fits single Gaussian (or Lorentzian) functions. Thus, bimodal data need to be separated beforehand to avoid variances.

GaussFit assumes that each input file contains data from one junction, thus, it takes the number of input files as the default number of degrees of freedom for all statistical calculations.

Many of the functions, such as computing differential conductance and rectification ratios, depend on the ability to separate the data into individual _I/V_ sweeps. This dependence can cause problems with hysteretic data. We start our _I/V_ sweeps at 0 V and then sweep to the maxima and back, meaning that each full sweep is two _I/V_ curves. The software should be able to pick up individual traces for any sweeps that begin and end at the same bias, but will almost certainly fail to parse anything else.

`GaussFit.py` is fairly verbose by default and will summarize data that it cannot parse, that can safely be considered numerical artifacts (e.g., when computing derivatives) or that fall outside of a cut-off range specified by the user. It is vital to pay attention to this output to ensure that the majority of the input data are actually being parsed into the final outputs.

## Dependencies

```
python 3.4+
pandas
colorama
numpy
scipy
matplotlib
appdirs
```

### Usage

You can just call the GUI with `GaussFit.py -G` and work from there.

To plot your data quickly without writing files:

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

A default configuration file will be created on first-run. To see where it is:

```
GaussFit.py -h
```

### Issues

We are not software developers. We offer this software as-is for anyone interested in parsing large amounts of _I/V_ data from large-area tunneling junctions, but it is tested quite specifically on EGaIn and CP-AFM data. If you have problems, you are on your own.
