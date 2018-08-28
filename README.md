# parallelSEE: A parallel Monte Carlo ray-tracing code for efficient simulations of secondary electron emission (SEE) in arbitrarily complex geometries under electron irradiation

Codes based on "Calculation of secondary electron emission yields from low-energy electron deposition in tungsten surfaces" https://www.sciencedirect.com/science/article/pii/S0169433218311176  & "Monte Carlo Raytracing Method for Calculating Secondary Electron Emission from Micro-Architected Surfaces" https://arxiv.org/pdf/1806.00205.pdf by Hsing-Yin (Iris) Chang and Andrew Alvarado.

Also see "Monte Carlo modeling of low-energy electron-induced secondary electron emission yields in micro-architected boron nitride surfaces" by Hsing-Yin (Iris) Chang and Andrew Alvarado for a similar approach.

# Code overview
```insulator_hBN.py``` is the main program based on the discrete Monte Carlo model, which includes multiple elastic and inelastic scattering processes: 
(i) Mott’s theory for electron-atom interactions, (ii) Ritchie’s theory for electron-electron interactions, (iii) Fro ̈hlich’s perturbation theory for electron-phonon interactions, and (iv) Ganachaud and Mokrani’s semi-empirical model for electron-polaron interactions.

# Running the program
To run the program, simply execute ```makedirect.sh``` script.

Once it finishes, use ```SEYcalc.sh``` to create a file that has the SEE yield and standard error for each energy. ```SEYcalc``` uses another program ```avecalc``` and also ```applog.sh```. ```applog.sh``` appends output files from each directory into the comm directory. 

# Input files
Sample input files for hexagonal boron nitride (```DESCS_hBN.txt```, ```ELoss_hBN.txt```, ```ThetaEl_hBN.txt```) are provided in the sample_inputs/ folder.

See below for a list of the possible input parameters to these programs and what they do.

# Input parameters
* Primary electron energy [keV] 
* Incident angle [rad]
* Sampling times




