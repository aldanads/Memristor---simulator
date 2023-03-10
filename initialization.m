%-----------------------------------------------------------------%
% Simulator: Lateral memtransistors from monolayer single-crystal %
% molybdenum disulfide (MoS2)                                     %
% Starting Date: 2021/11/16                                       %
% Samuel Aldana Delgado                                           %
%-----------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Purpose: This function initializes the main parameters of the simulator
%%
%%------------------------------------------------------------------------
%% -------------------- Parameters ---------------------------------------
%% d(1) = dx --> dimension of the grid in x axis
%% d(2) = dy --> dimension of the grid in y axis
%% h(1) = hx --> Step in x axis
%% h(2) = hy --> Step in y axis
%% Grid_S --> 2D hexagonal grid with d(1) and d(2) dimension and h(1) and h(2) separation between grid points
%%       within the Grid_S matrix:
%%                                - 0 forbidden position, not part of the crystal
%%                                - 1 Sulfur atom, S
%%                                - 2 Sulfur vacancy, Vs
%% dose
%% dose_width
%% -----------------------------------------------------------------------



function [Grid_S,phy_const,ActE,Nparticles,Vs_ij,parameters,direction_vac,direction_phy_pl,direction_data,res_fit_parameters,screening_fit_parameters,resistance_limits,den_res]=initialization(n_sim)

%%%%%%%%%%%%%%%%%%%%% Initilize the grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distance between grid points - Steps in x and y axises (meters)
%%%%%%%%%%%%%%%%%%%%% Lattice constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Mortazavi, M., Wang, C., Deng, J., Shenoy, V. B., & Medhekar, N. V. (2014). 
%Ab initio characterization of layered MoS2 as anode for sodium-ion batteries. 
%Journal of Power Sources, 268, 279-286.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2H-MoS2 - Semiconductor a=b=0.639(A)
% Steps in x
%hx=0.2767 (nm)
hx=0.277E-9;
% Steps in y
%hy=0.3195 (nm)
hy=0.32E-9;
% Thickness of a MoS2 monolayer
hz=0.65E-9;
%% Angle in the hexagonal between sulfur positions
theta=atan(hy/hx);

%% Dimension of the square to calculate the density
% Careful! --> Use 3,4,6,7,8
% 2 is too small
% We need x axis even. 
square_size=6;

% It is important that the number of squares to calculate the density fit
% well with the simulation domain
%% Dimension of the grid
% x axis --> Wide of the device
d(1)=50E-9;
% y axis --> Channel lenght
d(2)=50E-9;

if square_size==0
d(1)=round(d(1)/hx);
% x axis should be even, because make easier to simulate the borders
if mod(d(1),2)~=0
d(1)=d(1)+1;    
end
% y axis --> Channel lenght
d(2)=round(d(2)/hy);

else
d(1)=square_size*round(d(1)/(square_size*hx));
% x axis should be even, because make easier to simulate the borders
if mod(d(1),2)~=0
d(1)=d(1)+1;    
end
% y axis --> Channel lenght
d(2)=square_size*round(d(2)/(square_size*hy));
end

[Grid_S,Nparticles]=hex_grid(d);
ni=size(Grid_S,1);
nj=size(Grid_S,2);
%% --------------------------------------------------------------------%%



%%%%%%%%%%%%%%%%%%%% Initialize the defects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
distribution=2;

%% Kind of distribution %% 

% distribution=1 --> Uniform rows
if distribution==1
%% Dose of He ions in the device: number of defects
dose=0.75;
%% dose_width = 0 nm --> the last row
% This is the half of the width
%dose_width=round(10E-9/hy);
dose_width=(2*n_sim)*1E-9/hy;

parameter_seed(1)=dose;
parameter_seed(2)=dose_width;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% defect density: %%%%%%%%%%%%%%%%%%%%%%%%% 
% Jadwiszczak, Jakub, et al. "MoS2 memtransistors fabricated by localized helium ion beam irradiation." ACS nano 13.12 (2019): 14262-14273.
% At most: ∼4.2 × 1E14 VS/cm2
% Defect density between LRS and HRS: ∼1.6 × 1E12 VS/cm2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distribution=2 --> Gaussian distribution
if distribution==2

% This is the probability we are going to use (like dose in distribution 1)
% We are going to make a distribution of probabilities
mean=0.3;
standard_deviation=0.2;
%half_damaged_region=(2+0.5*n_sim)*1E-9;
half_damaged_region=(2+0.5*4)*1E-9;

%Skewing parameter and kurtosis
skew=1.2;

% Kurt must be greater than the square of skew + 1
kurt=4.4;
%kurt=4+0.2*n_sim;

% Probability peak at the mean of the distribution
%max_probability=normpdf(mean,mean,standard_deviation)/100;
x=-2:0.01:2;
y=pearspdf(mean+x,mean,standard_deviation,skew,kurt)/100;
max_probability=max(y);
max_dose=1/max_probability;
irradiation_position=round(nj/2)+1;
parameter_seed(1)=mean;
parameter_seed(2)=standard_deviation;
parameter_seed(3)=hy;
parameter_seed(4)=half_damaged_region;
parameter_seed(5)=irradiation_position;
parameter_seed(6)=max_dose;
parameter_seed(7)=skew;
parameter_seed(8)=kurt;

% I need this position for calculating the immobile vacancies
parameters(11)=irradiation_position;
end

%% distribution=3 --> test: interaction between particles in the same row
if distribution==3
%Coordinates of the first particle
parameter_seed(1)=round(ni/2)-1;
parameter_seed(2)=round(nj/2);
end

%% distribution=4 --> test: interaction between particles in the same column
if distribution==4
%Coordinates of the first particle
parameter_seed(1)=round(ni/2);
parameter_seed(2)=round(nj/2);
end

%% distribution=5 --> test: particles at borders: top, left and right
if distribution==5
%Coordinates of the first particle
parameter_seed(1)=round(ni/2);
parameter_seed(2)=round(nj/2);
end

% distribution=6 --> Triangle
if distribution==6
%% Probability of defects in the channel 
%dose_height=0.3+n_sim*0.1;
dose_height=0.9;
%% Base of the triangle
% The half of the triangle base
triangle_base=(2*n_sim)*1E-9/hy;
%triangle_base=(2*5)*1E-9/hy;
% Triangle pointing upward --> 1
% Triangle pointing downward --> -1

orentation_triangle=1;

parameter_seed(1)=dose_height;
parameter_seed(2)=triangle_base;
parameter_seed(3)=orentation_triangle; 
end

[Grid_S,Nparticles,Vs_ij]=seed_defects(Grid_S,Nparticles,distribution,parameter_seed);
%% --------------------------------------------------------------------%%


%%%%%%%%%%%%%%%%%%%% Physical constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nu0 (s^-1) bond vibration frequency
nu0=7E13;
% Charge (e-)
q=2;
%% Screening of charge --> screening=1 --> no screening
%screnning=0.0035;
screnning=0.0015;
% Screening of the electric field due to the defects concentration
%screnning_max=0.04+0.01*n_sim;
screnning_max=0.073;
screening_min=1;
%screening_min=0.68;


%% Relative permittivity of MoS2 --> %%%%%%%%%%%%%%%%%%%%
% Davelou, D., Kopidakis, G., Kioseoglou, G., & Remediakis, I. N. (2014). 
% MoS2 nanostructures: Semiconductors with metallic edges. Solid state communications, 192, 42-46.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
er=3.7;

%% Thermal conductivity MoS2
%Yan, Rusen, et al. "Thermal conductivity of monolayer molybdenum disulfide obtained from temperature-dependent Raman spectroscopy." 
% ACS nano 8.1 (2014): 986-993.

% Zhu, Gaohua, Jun Liu, Qiye Zheng, Ruigang Zhang, Dongyao Li, Debasish Banerjee, and David G. Cahill. 
% "Tuning thermal conductivity in molybdenum disulfide by electrochemical intercalation." Nature communications 7, no. 1 (2016): 1-9.
% Different conductivity in plane and through plane

% Defects have a large impact on the phonon scattering --> %% Thermal
% diffusivity MoS2 and conductivity 
% Dolleman, Robin J., et al. "Transient thermal characterization of suspended monolayer MoS 2." 
% Physical Review Materials 2.11 (2018): 114008.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W/mK (watt/(meter*kelvin))) 
% Between 1.35 and 84 --> Some cites in the introduction of the paper
%KMoS2=34.5;
KMoS2=84;
%% Density (kg/m^3)
ro_MoS2=5060;
ro_Mo=10200;
%% Capacidad calorifica especifica (J/(kg*K))
C_esp_MoS2=373;
C_esp_Mo=251;
%% Thermal conductivity Molybdenum
% https://www.imoa.info/molybdenum/molybdenum-properties.php
% W/mK (watt/(meter*kelvin)))
%KMo=138; % Thermal conductivity for bulk Mo
KMo=138; % Thermal conductivity for bulk Mo

%% Thermal diffusivity MoS2
% Dolleman, Robin J., et al. "Transient thermal characterization of suspended monolayer MoS 2." 
% Physical Review Materials 2.11 (2018): 114008.
% Thermal diffusivity (m**2/s) 
%alpha_MoS2=1.14E-5;
%alpha_MoS2=4.45E-5;
%alpha_MoS2=0.2E-5; % Within the experimental range
alpha_MoS2=KMoS2/(ro_MoS2*C_esp_MoS2);
%% Thermal diffusivity Mo
% https://www.engineersedge.com/heat_transfer/thermal_diffusivity_table_13953.htm
% Thermal diffusivity (m**2/s) 
%alpha_Mo=5.43E-5; ---> Thermal diffusivity for bulk Mo
alpha_Mo=KMo/(ro_Mo*C_esp_Mo);
%% Factor for decreasing the heat source
quenching_heat=1/10; % The substrate and the electrodes are huge, it might be higher
%% Activate heat equation
% 2 Time-dependent heat equation: https://www.youtube.com/watch?v=V00p-TgQML0https://www.youtube.com/watch?v=V00p-TgQML0
%% Time-dependent heat equation --> It doesn't work! --> Because it doesn't fulfill the condition: Courant–Friedrichs–Lewy
% Link 1: https://www.simscale.com/blog/cfl-condition/
% Link 2: https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition
% 1 Steady equation
% 0 desactivate
heat_equation=1;


phy_const(1)=nu0;
phy_const(2)=q;
phy_const(3)=theta;
phy_const(4)=hx;
phy_const(5)=hy;
phy_const(6)=er;
phy_const(7)=screnning_max;
phy_const(8)=screening_min;
phy_const(9)=screnning;
phy_const(10)=d(1);
phy_const(11)=d(2);
phy_const(12)=KMoS2;
phy_const(13)=KMo;
phy_const(14)=hz;
phy_const(15)=alpha_MoS2;
phy_const(16)=alpha_Mo;

%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tol Poisson equation
tol=0.1;

% Applied voltage (V)
V_initial=0;

% Temperature (K) T = 358.15 K = 85 celsius
T_ambient=300;
T_experiment = 358.15;
T_top_electrode=T_experiment;
T_bottom_electrode=T_experiment;

% Degradation level of the resistance --> 0.1 means a degradation of 90%
degradation_level = 0.95;
% Type of experiment: Thermal annealing: thermal - experiment = 1
% Voltage ramp: ramp - experiment = 2
experiment = 2;

% Time (s)
time=0;
% Delta time (s)
delta_t=1;
% Delta V (V)
delta_V=0.714;
%delta_V=2.1;
% Vmax during RESET (V)
Vmax=35;
% Vmin during SET (V)
Vmin=-35;
% Number of Resistive Switching Cycles
n_RS=1;
% Polarity
% Polarity = 1 --> y=0 --> V
%                  y=L --> 0
% Polarity = -1 --> y=0 --> 0
%                  y=L --> V
polarity=1;

% Number of rows to calculate the density
n_rows_density=3;





% Screening mechanisms: 1 --> Density
%                       2 --> Resistivity
Screening_mechanisms=2;

% Limits of resistivity regimes like in:
%% [Fox2015] Fox, Daniel S., et al. "Nanopatterning and electrical tuning of MoS2 layers with a subnanometer helium ion beam." Nano letters 15.8 (2015): 5307-5313.

% A row completely full of defects --> (ni/2)/(ni*hx)
if square_size==0
max_dens=(1E-9)/(2*hx);
else   
% Maximum number of particles in a square
max_dens=max(sum(sum(Grid_S(1:square_size,1:square_size))),sum(sum(Grid_S(1:square_size,square_size+1:2*square_size))));
% Maximum density in nanometers^2
max_dens=(max_dens*1E-18)/(hx*hy*square_size^2);
min_dens=1E-18/(hx*hy*square_size^2);
end

%% Resistance of the pristine channel, maximum resistance, minimum resistance
%pristine_state_resistance=3E6;
pristine_state_resistance=3E2;
resistance_limits=zeros(6,1);
resistance_limits(1)=pristine_state_resistance;


% We are going to carry out 6 different tests with different resistivity configurations, to assess how it behaves the resistive switching in each case, namely:
resistance_configuration=4;

%% 1) Linear 
% resistance_configuration==1

%% 2) Power
% resistance_configuration==2

%% 3) Semiconductor-Insulator-Metal --> /\
if resistance_configuration==3
section_1_lim=0.6*max_dens;
section_2_lim=1*max_dens;

% 3 points 
second_point_resistance=6E5;
max_resistance=6E8;
min_resistance=2.5E5;

resistance_limits(3)=second_point_resistance;
resistance_limits(4)=max_resistance;
resistance_limits(5)=min_resistance;
end 

%% 4) Insulator-Metal-Insulator --> \/
if resistance_configuration==4
%section_1_lim=0.5*max_dens;
section_1_lim=min_dens;
%section_1_lim=(0.1*n_sim)*max_dens;
section_2_lim=1*max_dens;

% 3 points
max_resistance=1.7E9;
min_resistance=5E2;
last_point_resistance=1.7E9;
%last_point_resistance=3E7;

resistance_limits(4)=max_resistance;
resistance_limits(5)=min_resistance;
resistance_limits(6)=last_point_resistance;
end

%% 5) Semiconductor-Insulator-Metal-Insulator (Complete) [Fox2015] --> _/\/
if resistance_configuration==5
section_1_lim=0.2*max_dens;
% Between insulator and metal
section_2_lim=0.45*max_dens;

% 5 points 
first_point_resistance=2E6;
second_point_resistance=6E5;
max_resistance=3E8;
min_resistance=2.5E5;
last_point_resistance=3E8;

resistance_limits(2)=first_point_resistance;
resistance_limits(3)=second_point_resistance;
resistance_limits(4)=max_resistance;
resistance_limits(5)=min_resistance;
resistance_limits(6)=last_point_resistance;
end

%% 6) Semiconductor-Insulator-Metal (uncompleted) --> _/\
if resistance_configuration==6
section_1_lim=0.2*max_dens;
% Between insulator and metal
section_2_lim=0.45*max_dens;

% 4 points 
first_point_resistance=2E6;
second_point_resistance=6E5;
max_resistance=3E8;
min_resistance=2.5E5;

resistance_limits(2)=first_point_resistance;
resistance_limits(3)=second_point_resistance;
resistance_limits(4)=max_resistance;
resistance_limits(5)=min_resistance;
end





parameters(1)=tol;
parameters(2)=V_initial;
parameters(3)=T_ambient;
parameters(4)=time;
parameters(5)=delta_t;
parameters(6)=delta_V;
parameters(7)=Vmax;
parameters(8)=Vmin;
parameters(9)=n_rows_density;
parameters(10)=n_RS;
%parameters(11)=irradiation_position;
parameters(12)=section_1_lim;
parameters(13)=section_2_lim;
parameters(14)=Screening_mechanisms;
parameters(15)=resistance_configuration;
parameters(16)=polarity;
parameters(17)=square_size;
parameters(18)=max_dens;
parameters(19)=min_dens;
parameters(20)=T_top_electrode;
parameters(21)=T_bottom_electrode;
parameters(22)=quenching_heat;
parameters(23)=heat_equation;
parameters(24)=experiment;
parameters(25) = T_experiment;
parameters(26) = degradation_level; 

%%%%%%%%%%%%%%%%% Saving --> properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
destiny_direction='C:\Users\aldanads\OneDrive - TCDUD.onmicrosoft.com\2D device simulator project\Publications\Failure mechanism - thermal\Failure mechanism\KMC2_test\';
folder_name=strcat(num2str(parameters(10)),'RS_Sim_',num2str(n_sim));
[direction_vac,direction_phy_pl,direction_data]=save_files(destiny_direction,folder_name);





%% --- Activation energies (eV) ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wang, L., Liao, W., Wong, S. L., Yu, Z. G., Li, S., Lim, Y. F., ... & Ang, K. W. (2019). 
% Artificial synapses based on multiterminal memtransistors for neuromorphic application. 
% Advanced Functional Materials, 29(25), 1901106. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Movements
mov_actE=2.297;
% Double vacancy
doub_mov=4.149;




% Movement limitation
mov_lim=0.8;

% Toward up 
ActE(1)=mov_actE;
% Toward down
ActE(2)=mov_actE;
% Horizontal movement --> Positive x
ActE(3)=mov_actE;
% Activation energy for movement (reference)
ActE(4)=mov_actE;
% Activation energy inferior limitation
ActE(5)=mov_lim;
% Activation energy double vacancy
ActE(6)=doub_mov;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resistivity
[res_fit_parameters,den_res]=linear_fit_resistance(parameters,phy_const,resistance_limits);
[screening_fit_parameters]=linear_fit_screening(parameters,phy_const,resistance_limits);


    function[res_fit_parameters,den_res]=linear_fit_resistance(parameters,phy_const,resistance_limits)
    % Step distance between each grid point in a row
    hx=phy_const(4)*1E9;
    % Grid points in a row
    ni=phy_const(10);
    % ni is very small, hence the change in density for every new particle 
    % is quite big
    
    % Only one defect in a row --> We must work in nm, as the width of
    % the simulation domain is only 10 nm. If we work in cm, this means
    % that when there is only one particle in the row, we assume that we
    % find 1 particle periodically. So, the density would be overrated.
        
    % Only one defect in the row
    min_dens=parameters(19);
    % A row completely full of defects --> (ni/2)/(ni*hx)
    max_dens=parameters(18);


    den_res=zeros(2,round(max_dens/min_dens));
    den_res(1,:)=min_dens:min_dens:max_dens;
    

    % Separation of sections depending on the defect density
    section_1_lim=parameters(12);
    section_2_lim=parameters(13);

    % resistivity configuration
    resistance_configuration=parameters(15);

    if resistance_configuration==3
    [res_fit_parameters]=linear_fit_config_3(section_1_lim,section_2_lim,min_dens,resistance_limits);
    end

    if resistance_configuration==4
    [res_fit_parameters]=linear_fit_config_4(section_1_lim,section_2_lim,min_dens,resistance_limits);

        slope_1=res_fit_parameters(1,1);
        b_1=res_fit_parameters(1,2);
        slope_2=res_fit_parameters(2,1);
        b_2=res_fit_parameters(2,2);


    % logical indices
    idx=den_res(1,:)<=section_1_lim;
    den_res(2,idx)=den_res(1,idx).^slope_1*10^b_1;
    den_res(2,~idx)=den_res(1,~idx).^slope_2*10^b_2;

    end

    if resistance_configuration==6
    [res_fit_parameters]=linear_fit_config_6(section_1_lim,section_2_lim,max_dens,min_dens,resistance_limits);
    end

        function [res_fit_parameters]=linear_fit_config_3(section_1_lim,section_2_lim,min_dens,resistance_limits)
        
        % /\
        second_point_res=resistance_limits(3); 
        max_res=resistance_limits(4);
        min_res=resistance_limits(5);
            
        % The resistivity has two regimes
        %% First regime --> Semiconductor to insulator
        slope_1=log10(max_res/second_point_res)/log10((section_1_lim)/(min_dens));
        b_1=log10(max_res)-slope_1*log10(section_1_lim);
        %% Second regime --> Insulator to metal
        slope_2=log10(min_res/max_res)/log10(section_2_lim/(section_1_lim));
        b_2=log10(min_res)-slope_2*log10(section_2_lim);

        res_fit_parameters(1,1)=slope_1;
        res_fit_parameters(1,2)=b_1;
        res_fit_parameters(2,1)=slope_2;
        res_fit_parameters(2,2)=b_2;
        end

        function [res_fit_parameters]=linear_fit_config_4(section_1_lim,section_2_lim,min_dens,resistance_limits)
        
        % \/
        max_res=resistance_limits(4);
        min_res=resistance_limits(5);
        last_point_res=resistance_limits(6);
            
        %% First regime --> Insulator to metal
        slope_1=log10(min_res/max_res)/log10(section_1_lim/(min_dens));
        b_1=log10(min_res)-slope_1*log10(section_1_lim);

        %% Second regime --> Metal to insulator
        slope_2=log10(last_point_res/min_res)/log10(section_2_lim/section_1_lim);
        b_2=log10(last_point_res)-slope_2*log10(section_2_lim);

        % If slope_1=inf --> Only one branch
        if isinf(slope_1)
            slope_1=slope_2;
            b_1=b_2;
        end

        res_fit_parameters(1,1)=slope_1;
        res_fit_parameters(1,2)=b_1;
        res_fit_parameters(2,1)=slope_2;
        res_fit_parameters(2,2)=b_2;
        end
        
        function [res_fit_parameters]=linear_fit_config_6(section_1_lim,section_2_lim,max_dens,min_dens,resistance_limits)

        %_/\
        first_point_res=resistance_limits(2);
        second_point_res=resistance_limits(3);
        max_res=resistance_limits(4);
        min_res=resistance_limits(5);
            
        % Mapping of the resistivity in function of the dose.
        % Using this:
        %% [Fox2015] Fox, Daniel S., et al. "Nanopatterning and electrical tuning of MoS2 layers with a subnanometer helium ion beam." Nano letters 15.8 (2015): 5307-5313.
    
        % The resistivity has three regimes
        %% First regime --> Semiconductor due to the dopping effect: 6E5/2E6
        slope_1=log10(second_point_res/first_point_res)/log10(section_1_lim/min_dens);
        b_1=log10(second_point_res)-slope_1*log10(section_1_lim);

        %% Second regime --> Insulator --> Max in Fox2015: 3E10/6E5
        slope_2=log10(max_res/second_point_res)/log10((section_2_lim)/(section_1_lim));
        b_2=log10(max_res)-slope_2*log10(section_2_lim);

        %% Third regime --> Metal 2.5E5/3E8
        slope_3=log10(min_res/max_res)/log10(max_dens/(section_2_lim));
        b_3=log10(min_res)-slope_3*log10(max_dens);

        res_fit_parameters(1,1)=slope_1;
        res_fit_parameters(1,2)=b_1;
        res_fit_parameters(2,1)=slope_2;
        res_fit_parameters(2,2)=b_2;
        res_fit_parameters(3,1)=slope_3;
        res_fit_parameters(3,2)=b_3;
        end
        
    end

    function[screening_fit_parameters]=linear_fit_screening(parameters,phy_const,resistance_limits)
    
    % Be careful with the notation --> screening_min > screening_max
    % Screening_min is closer to 1, so the screening effect over the
    % electric field is lower.
    screnning_max=phy_const(7);
    screening_min=phy_const(8);

    screening_mechanism=parameters(14);

    if screening_mechanism==1


    % Only one defect in a row
    min_dens=1E7/(ni*hx);
    % A row completely full of defects --> (ni/2)*1E7/(ni*hx)
    max_dens=1E7/(2*hx);
    

    % We assume the relation between concentration of defects and screening
    % is linear. y=slope*x+b
    %% Slope
    slope=(screnning_max-screening_min)/(max_dens-min_dens);
    b=screening_min-slope*min_dens;
    screening_fit_parameters(1)=b;
    screening_fit_parameters(2)=slope;
    screening_fit_parameters=screening_fit_parameters';
    else
        
    % Max and min in resistivity (Ohmios)
    min_res=min(resistance_limits(resistance_limits>resistance_limits(1)));
    max_res=max(resistance_limits);


    % We assume the relation between resistivity and screening
    % is linear. y=slope*x+b
    % When the resistance is mininum, the screening effect is higher -->
    % the electric field is lower --> Difficult to drag particles
    % When the resistance is maximum, the screening effect is lower --> the
    % electric field is higher --> Drag particles
    
    %% Slope
    slope=(screening_min-screnning_max)/(max_res-min_res);
    b=screening_min-slope*max_res;

    screening_fit_parameters(1)=b;
    screening_fit_parameters(2)=slope;
    screening_fit_parameters=screening_fit_parameters';
    end

    end

%% Function for data
    function [direction_vac,direction_phy_pl,direction_data]=save_files(destiny_direction,folder_name)
        
        % Create a few folders
        mkdir(destiny_direction,folder_name);
        destiny_folder=strcat(destiny_direction,'\',folder_name,'\');
        mkdir(destiny_folder,'Data');
        mkdir(destiny_folder,'Vac_mov');
        direction_vac=strcat(destiny_folder,'\Vac_mov\');
        mkdir(destiny_folder,'Phy_plots');
        direction_phy_pl=strcat(destiny_folder,'\Phy_plots\');
        direction_data=strcat(destiny_direction,'\',folder_name,'\Data\');
        mkdir(direction_data,'Program');
        

        % Save the code
        path = strcat(direction_data,'Program\initialization.m');
        copyfile('initialization.m', path)

        path = strcat(direction_data,'Program\hex_grid.m');
        copyfile('hex_grid.m', path)

        path = strcat(direction_data,'Program\density.m');
        copyfile('density.m', path)

        path = strcat(direction_data,'Program\KMC.m');
        copyfile('KMC.m', path)

        path = strcat(direction_data,'Program\plot_graphs.m');
        copyfile('plot_graphs.m', path)

        path = strcat(direction_data,'Program\seed_defects.m');
        copyfile('seed_defects.m', path)

        path = strcat(direction_data,'Program\simulator_core.m');
        copyfile('simulator_core.m', path)

        path = strcat(direction_data,'Program\SolvePotentialAndField.m');
        copyfile('SolvePotentialAndField.m', path)

        path = strcat(direction_data,'Program\pearspdf.m');
        copyfile('pearspdf.m', path)
          
    end


end
