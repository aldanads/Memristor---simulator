

% Charge transport based on Schotty barrier 

A0 = 1.202E6; % Richardson constant (A*(m^-2)*(K^-2))
T = 300; % Temperature (K)
%Constante de Boltzmann (eV/K)
kb=8.6173324E-5;
density_Vs = 0.22; % (Vs/nm2)
density_Vs_cm2 = 0.22*1E14;
er = 3.1; % Dielectric constant MoS2
e0 = 8.854187817E-12; % Permittivity of free space (C^2 / (N*m^2))
e_charge = 1.60217663E-19; % Electron charge (C)

% Sangwan, Vinod K., Hong-Sub Lee, Hadallia Bergeron, Itamar Balla, Megan E. Beck, Kan-Sheng Chen, and Mark C. Hersam. 
% "Multi-terminal memtransistors from polycrystalline monolayer molybdenum disulfide." Nature 554, no. 7693 (2018): 500-504.
% Sangwan says that phi_bn (effective Schottky barrier) range between 80–125 meV in experiments
% They use 20 meV to 280 meV in simulations
phi_b0 = 0.385; % Schottky barrier (eV) not modulated  (fitting parameter)

% w is in the range of a few nanometers --> region with excess of dopants
% we are simulating this region
w = 1e-9;

% Search doping values in: Kim, In Soo, et al. 
% "Influence of stoichiometry on the optical and electrical properties of chemical vapor deposition derived MoS2." 
% ACS nano 8.10 (2014): 10551-10558.
% delta_n (cm^-3) --> excess of doping in region w
delta_n = density_Vs * 1E27;
% They have values of w*delta_n between 10E10 and 3.5E9 cm^-2 

% A constant
A = 10-5;
D = 5.5E-9;

doping_term = e_charge/(er*e0) * sqrt(w*delta_n/(4*pi));
