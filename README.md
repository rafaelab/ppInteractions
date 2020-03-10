# ppInteractions

Plugin to be used with CRPropa 3 to handle proton-proton scattering. 
It models the interaction of high-energy cosmic rays with hydrogen-rich gas in astrophysical environments.


The interaction and spectra of secondary particles are parametrised following:

- Kelner, Aharonian, Bugayov. Phys. Rev. D 74 (2006) 034018. (see also erratum: PRD 79 (2009) 039901).

- Kafexhiu, Aharonian, Taylor, Vila. Phys. Rev. D 90 (2014) 123014.


Information about CRPropa 3 can be found [here](https://github.com/CRPropa/CRPropa3/) and in the paper:
Alves Batista et al. J. Cosmol. Astropart. Phys. 05 (2016) 038.

## Notes
- The production of secondaries does not conserve energy at each interaction, only statistically.
- The NucleusNucleusInteraction module produces mesons (pions, eta). The ParticleDecay model performs the actual decays of these particles.

## To-do
- For now this module only performs proton-proton interactions, but it will soon be extended for nucleus-nucleus interactions.

## Disclaimer
This plugin is still being tested and may contain bugs. 
Use it at your own risk.


## Installation procedure

To install this module, you will need to have CRPropa 3 installed. 
Go to https://github.com/CRPropa/CRPropa3/ and follow the instructions.

Now proceed with the installation of this module.

1. Download the latest version of the code.
```
git clone https://github.com/rafaelab/ppInteractions.git
```

2. Navigate into the downloaded folder and create a folder called "build/".

3. Install the code with CMake and then make:

```
cmake ..
make
```

4. If the code was properly compiled, you are ready to go!
Make sure to add the path where ppInteractions.py is created to your PYTHONPATH.
Alternatively, this can be hard-coded in your python script.

For consistency with CRPropa, the flag `-DCMAKE_CXX_FLAGS="-std=c++11"` may be required to force compilers to adopt C++11.