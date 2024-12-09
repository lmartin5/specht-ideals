# Specht Ideals

This code uses the `SpechtModule` package in `Macaulay2` to calculate examples 
of Hilbert functions, series, and polynomials of Specht ideals. The ultimate goal 
is to compute examples of equivariant Hilbert series, of which they are very few currently
known.

## Table of Contents

- [Pre-requisites](#pre-requisites)
- [Usage](#usage)
- [Contributions](#contributions)
- [Acknowledgements](#acknowledgements)
- [License](#license)

# Pre-requisites
<sup>[(Back to top)](#table-of-contents)</sup>

This project assumes that you have successfully installed `Macaulay2` on your machine.
The functions included in this repo require the `SpechtModule` and `Posets` packages.
To install the packages, run the commands

```sh
i1: installPackage "SpechtModule"
i2: installPackage "Posets"
```
while in a `M2` shell.

# Usage
<sup>[(Back to top)](#table-of-contents)</sup>

To use the functions included in the `specht-functions.m2` file, run the command
```sh
i1: load "specht-functions.m2"
```
while in a `M2` shell.

# Acknowledgements
<sup>[(Back to top)](#table-of-contents)</sup>

HERE

# License
<sup>[(Back to top)](#table-of-contents)</sup>

This project uses the GPL-3.0 license. See the file LICENSE for details.