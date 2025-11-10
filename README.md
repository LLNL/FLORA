# Flora Equilibrium and Stability Code

### Abstact:

FLORA solves, in a 2-D domain (radial and axial dimensions) with a specified azimuthal Fourier mode, for the linearized stability of a long, thin, axisymmetric plasma equilibrium in an applied magnetic field. Before the stability equation is solved, FLORA solves a set of simple equations for pressure balance that specify the equilibrium magnetic field and plasma pressure in the long-thin limit given a simplified description of the magnetic coils. It uses an initial-value method for the linear stability problem in which an equilibrium is given an initial perturbation to its magnetic field, and the temporal behavior of the perturbation is followed. The perturbation has been Fourier expanded in the azimuthal coordinate; each azimuthal mode must be examined separately. The complex partial differential equation of motion for the perturbed radial displacement of the field lines is solved as a coupled system of two real p.d.e.'s and the solution consists of two parts, the real part and the imaginary part. The system is solved by bringing the coupling terms in each equation to the right side and using an iterative technique. [Here](docs/Flora_archive_abstract.pdf) is the full abstract with references.

### FLORA cases from archival microfiche. Scanned and saved as pdf.  
* [Microfiche Readme] (docs/README_microfiche.md)
* [Microfiche set 1](docs/BOXU21_FICHE_SET1)  
* [Microfiche set 2](docs/BOXU21_FICHE_SET2)  
* [Microfiche set 3](docs/BOXU21_FICHE_SET3)

### Documentation and Papers:
[FLORA User Manual](docs/flora_manual.pdf)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Addendum: [Limitation of FLORA equilibrium package](docs/Comment_on_FLORA_Equilibrium_Solution.pdf)  

Post, R. F., T. K. Fowler, R. Bulmer, J. Byers, D. Hua, and L. Tung. "Axisymmetric tandem mirrors: stabilization and confinement studies." [Fusion science and technology 47, no. 1T (2005): 49-58.](docs/Axisymmetric_Tandem_Mirrors_Stabilization_and_Confinement_Studies.pdf)

Caponi, M.Z., Cohen, B.I., and Freis, R.P., Stabilization of Flute Modes by Finite Larmor Radius and Surface Effects, [Phys. Fluids 30, 1410 (1987).](docs/Caponi_et_al.pdf)
    
Dobrott, D., Kendrick, S., Freis, R.P., and Cohen, B.I., Magnetohydrodynamic Stability of Hot-Electron Stabilized Tandem Mirrors, [Phys. Fluids 30, 2149 (1987).](docs/Dobrott_et_al.pdf)

Bruce I. Cohen, Robert P. Freis, William A. Newcomb; Interchange, rotational, and ballooning stability of long‐thin axisymmetric systems with finite‐orbit effects. [Phys. Fluids 1 May 1986; 29 (5): 1558–1577.](docs/Cohen-Freis-Newcomb.pdf)

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
