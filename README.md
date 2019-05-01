# pySOS
`pySOS` is a python interface for the Successive Orders of Scattering (SOS) radiative transfer code. The SOS is a radiative transfer code developed by CNES

The `pySOS` interface aims to incorporate the creation of run scripts and parsing of output results for the SOS model. It also incorporates helpers to perform common tasks like calculating the radiance for a certain band instead of a wavelength or running the model for multiple wavelengths. 

This code was inspired by [py6S](https://github.com/robintw/Py6S) by Robin Wilson.

## Installation

The installation of the `pySOS` has two parts.

First, you need to install the SOS software package https://github.com/CNES/RadiativeTransferCode-SOS.

Second, you need to download the last version of the `pySOS` from [github](https://github.com/fnemina/pySOS/releases/latest).

Once downloaded decompress it, go to the folder containing the code and run

```bash
python setup.py install
```

To then check that software installed correctly

```python
# Load pyOSOAA module
import pySOS
# Run the test suite
pySOS.test()
```
the following output should appear at the end of the screen
```
SOS wrapper script by Francisco Nemiña
Inspired by Py6S wrapper by Robin Wilson
Using SOS located at /home/.../RT/SOS6.2
Running SOS using a set of test parameters
The results are:
Expected result: 0.346131
Actual result: 0.346131
#### Results agree PyOSOAA is working correctly

```
