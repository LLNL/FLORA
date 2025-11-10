# Microfiche Archive Files

A cache of microfiche film was found containing FLORA simulations done in the early 80s. These were scanned and the images saved in pdf files. The pdf files were then processed on an Ubuntu Linux system to extract text from the images. Each image was processed with the Tesseract (pytesseract) Open Source OCR software and output into a .txt file. This was done as a batch process and the files are crudely formatted. 

Note that the images are not the output of the FLORA Python module contained in this repository but were output either directly by Fortran77 code or by the Basis code version of FLORA. Basis is a deprecated system for development of scienctific codes developed at LLNL.

### Archival microfiche sets.   

* [Microfiche set 1](docs/BOXU21_FICHE_SET1)  
* [Microfiche set 2](docs/BOXU21_FICHE_SET2)  
* [Microfiche set 3](docs/BOXU21_FICHE_SET3)

### Using ChatGPT to improve the OCR output. 

<strong>Note that this option requires a ChatGPT account.</strong> 

The file 251012_EDISK70T_FICHE1.pdf, from set 1, was displayed. [This screengrab](screengrab.png) was captured on a mac of some text. The png file from the mac was dragged/dropped onto the [ChatGPT page](https://chatgpt.com). The prompt "OCR this image" produced the following output. ChatGPT also offers other output formats.
<pre>
aname      .  r0 1.295432948026-2434 r7 0.0.
m          .  2
bias       .  1.e-02
lmax       .  4
nmax       .  400
nddgo      .  20
f1         .  1.
fzx        .  1.
fj1        .  -1.
fjx        .  -1.
ex1        .  1.
fx1        .  0.
fps1       .  0.5
fu         .  0.5
fv         .  0.5
fz         .  0.5
kplotb     .  0
kzs        .  0
jfour      .  1
nfourp     .  1500
nfourmax   .  5
sf6        .  0.
sf8        .  0.
swg1       .  1.
swg2       .  1.
swg3       .  1.
swg4       .  1.
0

echarg     .  4.8e-10
lb         .  0.
rv         .  0.
xm         .  0.
rvl        .  3.95
zmax       .  100.
dt         .  1.e-07
theta0     .  0.
pfudge     .  0.5
ps#rel     .  0.6
z1rel      .  0.5638
z2rel      .  0.9638
z3rel      .  0.7338
nslosum    .  2.e+13
cold       .  1.e-02
z3rel      .  0.5038
z4rel      .  0.9638
ncenter    .  1.e+13
betslis1   .  5.e-02
betcent    .  5.e-02
rp         .  3.
ztrans     .  0.4
ltrans     .  0.1
bmg        .  100000.
rpl        .  1.7
bceng      .  1114.
eps0       .  1.e-06
ph0e1      .  0.1
ph0lng     .  0.1
kn         .  23
xpot       .  1.
ypot       .  2.
vpot       .  2.
betcen     .  5.e-02
betslis2   .  5.e-02
betrng     .  0.
ps1rel     .  1.e-02
ps2rel     .  1.05
p2evde     .  0.25
ps3rel     .  1.05
pxp1       .  2.
pxp2       .  2.
safe       .  0.9
pelf       .  8.9
p2ev.de    .  0.25
ps#rel     .  0.58
p2floor    .  0.
</pre>
