# ppInteractions

Plugin to be used with CRPropa 3 to handle proton-proton scattering. 
It models the interaction of high-energy cosmic rays with hydrogen-rich gas in astrophysical environments.


The interaction and spectra of secondary particles are parametrised following:

- Kelner, Aharonian, Bugayov. Phys. Rev. D 74 (2006) 034018. (see also erratum: PRD 79 (2009) 039901).

- Kafexhiu, Aharonian, Taylor, Vila. Phys. Rev. D 90 (2014) 123014.


Information about CRPropa 3 can be found [here](https://github.com/CRPropa/CRPropa3/) and in the paper:

- Alves Batista et al. J. Cosmol. Astropart. Phys. 05 (2016) 038. [arXiv:1603.07142](https://arxiv.org/abs/1603.07142)

- Alves Batista et al. J. Cosmol. Astropart. Phys. 09 (2022) 035. [arXiv:2208.00107](https://arxiv.org/abs/2208.00107)

## Notes
- The production of secondaries does not conserve energy at each interaction, only statistically.
- The `NucleusNucleusInteraction` module produces mesons (pions, eta). The `ParticleDecay` model performs the actual decays of these particles.

## To-do
- For now this module only performs proton-proton interactions, but it will soon be extended for nucleus-nucleus interactions.


## Installation procedure

To install this module, you will need to have CRPropa 3 installed. 
Go to https://github.com/CRPropa/CRPropa3/ and follow the instructions.
Note that you need at least CRPropa 3.2.1 installed.

Now proceed with the installation of this module.

1. Download the latest version of the code.
```
git clone https://github.com/rafaelab/ppInteractions.git
```

2. Navigate into the downloaded folder and create a folder called "build/".

3. Install the code with CMake and then make:
```
cmake ..
make && make install
```
To ensure that CRPropa and Python are correctly found, you might want to consider explicitly setting some of the flags manually, as below.
```
cmake .. \
    -DCMAKE_INSTALL_PREFIX=$PWD \
    -DCMAKE_CXX_COMPILER=clang++ \
    -DCRPropa_INCLUDE_DIR=$CRPropa_DIR/include \
    -DCRPropa_LIBRARY=$CRPropa_DIR/libcrpropa.dylib \
    -DPython_EXECUTABLE=python3.12 \
    -DPython_LIBRARY=$PYTHON_DIR/lib/libpython$PYTHON_VERSION_MAJOR.$PYTHON_VERSION_MINOR.dylib \
    -DPython_INCLUDE_DIRS=$PYTHON_DIR/include/python$PYTHON_VERSION_MAJOR.$PYTHON_VERSION_MINOR \
    -DPython_INSTALL_PACKAGE_DIR=$PWD \
    -DCRPropa_SWIG_PATH=$CRPropa_DIR/share/crpropa/swig_interface
```


4. If the code was properly compiled, you are (probably) ready to go!
Try the following:
```
python -c "from ppInteractions import *"
```
Make sure to add the path where ppInteractions.py is created to your PYTHONPATH, if you want to expose it globally.
Alternatively, this can be hard-coded in your python script.



## Disclaimer
This program is provided 'as is', without warranties of any kind. 
Please use your discernement to interpret the results obtained with it.
