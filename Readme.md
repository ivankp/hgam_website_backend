# HGam website binning
This repository contains backend programs for the
[HGam binnig page](https://github.com/ivankp/hgam_website/tree/main/binning).

## `binner`
This is the backend program that the server runs when requests are sent on the
binning page.

The program directly takes the request query string as a single argument.
The program reads the necessary binary input files and produces JSON output.

## `convert_mxaods`
This program converts HGam MxAOD ROOT files to binary float dumps of the
observables' values.
This file format simplifies the binning page operation, avoids dependency on
ROOT, and improves the speed of data reading.

To generate binary files for the observables from MxAODs,
run the program once for the data files,
and once for the MC files.
For example:
1. `./bin/convert_mxaods data mxaods/data22/data22_37.8ifb.root`
2. `./bin/convert_mxaods data mxaods/NoSys/PhPy8EG_PDF4LHC21_*.root`

The first argument is the output directory where the data files will be
written.
The subsequent agruments are names of input MxAOD files.

The data (as oposed to MC) files require the respective luminosity to be
specified, so that MC weights can be scaled to get the expected number of
events.
The program expects luminosity to be specified in the MxAOD file names,
e.g. `_37.8ifb`.
If the files on disk don't have the respective luminosities in their names,
one can make links to them, using a command such as
`ln -s data.root data_37.8ifb.root`.

In order to modify what observables the program outputs,
edit the [`convert_mxaods.cc`](./src/convert_mxaods.cc) source code.

Every variable requires `branch_reader`s for reco and truth values,
as well as `write()` function calls for the respective output files inside
the event loop.

For variables depending on a certain number of jets, take care to place the
respective `write()` calls or `VAR` macros in a place where `enough_jets` has
the correct value.

Note, that `VAR` is a preprocessor macro.
It avoids a lot of code duplication and easy to make mistakes.
The meaning of `VAR` is **redefined** several times in the source code.

## `head`
A small program to inspect headers of the `.dat` files.

# Compilation and setup instructions
1. Run `make`.
2. Copy the compiled `hgam_website_binning/bin/binner` program to
   `hgam_website/binning/binner`.

The `Makefile` is set up to compile the `binner` program statically.
This should enable the program to work with any Linux distribution,
avoiding library dependencies, and version compatibility problems.
