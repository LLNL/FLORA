# Flora Equilibrium and Stability Code

### Documentation and Papers
    [User manual] (FLORA User's Manual UCID-20400 (April 1985).pdf)

    [Caponi, M.Z., Cohen, B.I., and Freis, R.P., Stabilization of Flute Modes by Finite Larmor Radius and Surface Effects, Phys. Fluids 30, 1410 (1987).](https://doi.org/10.1063/1.866254)

    [Dobrott, D., Kendrick, S., Freis, R.P., and Cohen, B.I., Magnetohydrodynamic Stability of Hot-Electron Stabilized Tandem Mirrors, Phys. Fluids 30, 2149 (1987).](https://doi.org/10.1063/1.866149)

    [Post, R. F., T. K. Fowler, R. Bulmer, J. Byers, D. Hua, and L. Tung. "Axisymmetric tandem mirrors: stabilization and confinement studies." Fusion science and technology 47, no. 1T (2005): 49-58.](https://www.osti.gov/servlets/purl/15014491)

### Compatibility:
Compatible with Python>3.5 and Gfortran>=10. Some compiler switches used are not available in Gfortran < 10. It's been tested with Gfortran 10.0.1, 11.2, 11.3, and 12.2 and Python 3.6-3.11 on Mac M1, Linux RHEL7 x86_64, and Pop 22.04 x86_64. 
#### Anaconda warning:
If using Python from an Anaconda installation pay close attention to the GCC version printed out by Anaconda when you run Python. If it reports a GCC(7.3.0) then you will need to install a non-Anaconda python as your Anaconda Python is incompatible with the Gfortran required for Flora. The test described below will generate a segmentation violation if Anaconda and Gfortran are incompatible.
### Building:
pip install forthon<br>
python setup.py build install # if you have permissions to install python packages <br>
python setup.py build install --user # to install in your user python area if you don't have needed permissions

### To run the code:
$ python<br>
\>>>import flora<br>
\>>>flora.glr.glrgen()

To setup case set variables with:<br>
from flora import glr<br>
glr.varname1 = ....<br>
glr.varname2 = ....<br>

For a list of variables use <br>
glr.varlist()

Because of the way the external objects are created a dir(glr) will not reveal the variable names.

### Examples:

There are a couple of cases in the examples directory. Run with:

python test1.py <br>
python test2.py <br>

Uncomment the call to plots1() for some graphics.

Compare with test1_ref.log and test2_ref.log. These logs are very old, compare for general agreement only. 

### Release 

Flora is released under an LGPL license.  For more details see the
NOTICE and LICENSE files.

``LLNL-CODE-431811``
------
--------
