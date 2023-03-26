# Disorder Conductivity

Compute a spatial map of optical conductivity from eigenvalues and eigenvectors of a Hartree-Fock-Bogoliubov solution of a disordered Hubbard model.

Read the note on [Conductivity](Conductivity.md).

Results should follow the conventions of https://code.osu.edu/oneal.144/TechCode and https://code.osu.edu/oneal.144/TrivediGroup.


## Instructions

### How to build

```
python build.py
```


### How to run

To add list of bonds to the FITS file,
```
python addbonds.py output.fits
```

To compute paramagnetic susceptibility,
```
python diamagnetic.py output.fits
```

To compute paramagnetic susceptibility,
```
python paramagnetic.py output.fits --frequency 0.1 --broadening 0.01
```

To compute all,
```
python cgs.py output.fits --frequency 0.1 --broadening 0.01
```
