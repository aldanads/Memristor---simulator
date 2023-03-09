%-----------------------------------------------------------------%
% Simulator: Lateral memtransistors from monolayer single-crystal %
% molybdenum disulfide (MoS2)                                     %
% Starting Date: 2021/11/16                                       %
% Samuel Aldana Delgado                                           %
%-----------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Purpose: This is the core of the simulator. From here, we call all the functions

%This is necessary to obtain different random numbers for each session of MATLAB
rng('shuffle','twister');

max_sim=2;
for n_sim=2:max_sim

[Grid_S,phy_const,ActE,Nparticles,Vs_ij,parameters,direction_vac,direction_phy_pl,direction_data,res_fit_parameters,screening_fit_parameters,resistance_limits,den_res]=initialization(n_sim);

[V_initial,time,delta_t,delta_V,Vmax,Vmin,n_RS,square_size,v_resistance,v_density_Vs,v_total_res,v_current,v_time,voltage,v_temperature,R_ratio,V_set_reset,cont_SET_RESET,j,i,tmax,v_curr_map,v_den_map,v_ey,v_res_map,u,T]=variables_counters(parameters,Grid_S);

% Number of RS cycles
for cycle=1:n_RS

%V_set_reset(3,cycle)=1;
%V_set_reset(6,cycle)=1;

% s=1 --> Right branch
% s=2 --> Left branch
for s=1:2

carry_on=true;

while (carry_on)

    %% ACTIVATE AND DESACTIVE ANNEALING AND VOLTAGE RAMP EXPERIMENT
    %% STOP ANNEALING WHEN THE RESISTANCE DECREASE A SPECIFIC AMOUNT 
    % parameters(24)=experiment(1);

    if parameters(24) == 1 % Heat experiment
        delta_V = 0;
        parameters(24) = 11; % Heat experiment on
        heat_equation = 0;
    elseif parameters(24) == 2 % Ramp experiment
        delta_V = parameters(6);
        T_ambient = parameters(3);
        T = ones(size(Grid_S,1),size(Grid_S,2)) * T_ambient;
        parameters(20) = T_ambient; % T top electrode
        parameters(21) = T_ambient; % T bottom electrode
        heat_equation = parameters(23);
        parameters(24) = 22; % Ramp experiment on
        i = 1; % Start voltage ramp from 0
    end

V=V_initial+(i-1)*delta_V;
tmax=tmax+delta_t;
V, j
while time(1)<tmax

voltage(j)=V;

%% Defect density
[density_Vs,v_density_Vs]=density(Grid_S,Vs_ij,phy_const,parameters,j,v_density_Vs);
v_den_map(:,:,j)=density_Vs;

%% Resistivity
[resistance_step,v_resistance,v_total_res]=resistance_calculation(density_Vs,parameters,res_fit_parameters,resistance_limits,j,v_resistance,v_total_res);
v_res_map(:,:,j)=resistance_step;
%% Current
[v_current,current_mapping,f]=current(voltage,v_total_res,j,v_current,v_resistance,resistance_step,Grid_S,parameters,phy_const);
v_curr_map(:,:,j)=current_mapping;
%% Solve Poisson equation and electric field
[u,ex,ey,screening_j]=SolvePotentialAndField(Grid_S,phy_const,Vs_ij,parameters,V,density_Vs,resistance_step,screening_fit_parameters,u);
v_ey(:,:,j)=ey;
%% Heat equation
if heat_equation>=1
[T]=SolveHeat(f,phy_const,Grid_S,parameters,time,T);
end
v_temperature(j,1)=mean(mean(T));
v_temperature(j,2)=max(max(T));
v_temperature(j,2)

%% KMC algorithm
[Grid_S,time,Vs_ij,prob]=KMC(Grid_S,phy_const,Vs_ij,ex,ey,ActE,T,time,parameters,delta_V);


end

switch s
    
    % Ramp voltage increases until max --> Then decreases until min
    case 1

        % When the resistance degradate that amount, we stop the annealing
        % experiment and start the voltage ramp exp
        % parameters(26) = degradation level
        if (v_total_res(j) <= v_total_res(1)*parameters(26)) && (parameters(24) == 11 || parameters(24) == 1)
        parameters(24) = 2; % experiment(1)
        end

        if V>=Vmax
        i=1;
        delta_V=-delta_V;
        V_initial=V;
        carry_on=false;

        % FLag for reading the resistance ratio in the RESET curve
        %R_ratio(4,1)=0;
        %R_ratio(4,2)=0;

        end
        
    % When V min is reached, go until 0
    case 2
        if V<=Vmin
        i=1;
        delta_V=-delta_V;
        V_initial=V;
        end
        
        if (delta_V>0)&&(V==0)
        carry_on=false;
        V_initial=V;
        end
end



%% Vector of time
v_time(j)=time(1);
%% Vector of voltage

% Only one SET and RESET voltage for each branch
if (V==0)
    cont_SET_RESET=cont_SET_RESET+1;
    R_ratio(4,1)=0;
    R_ratio(4,2)=0;
end
[R_ratio,V_set_reset]=measure(R_ratio,voltage,delta_V,v_resistance,j,V_set_reset,cont_SET_RESET);


%% Plot the movements of the vacancies
plot_graphs(Grid_S,Vs_ij,voltage,direction_vac,direction_phy_pl,v_time,Nparticles,phy_const,j,parameters,v_density_Vs,ey,v_resistance,screening_j,den_res,resistance_limits,density_Vs,resistance_step,v_current,v_total_res,current_mapping)
j=j+1;
i=i+1;
end
i=2;
end

end

[R_ratio]=statistics_R_ratio(R_ratio);
% Cost function useful for I-V fitting
%[J]=exp_sim_fit(n_RS,v_current,voltage,parameters);

%Save the variables in workspace
save(strcat(direction_data,'variables.mat'))

beep;
clear;
java.lang.System.gc();
end

function [R_ratio]=statistics_R_ratio(R_ratio)
m=size(R_ratio,2);

R_ratio(5,1:m/2)=R_ratio(3,1:2:m);
R_ratio(6,1)=sum(R_ratio(3,1:2:m)>1);
R_ratio(6,2)=m/2-R_ratio(6,1);
R_ratio(6,3)=mean(R_ratio(3,1:2:m));
R_ratio(6,4)=max(R_ratio(3,1:2:m));
R_ratio(6,5)=min(R_ratio(3,1:2:m));
R_ratio(6,6)=std(R_ratio(5,1:m/2));
R_ratio(6,7)=std(abs(R_ratio(5,1:m/2-1)-R_ratio(5,2:m/2)));
R_ratio(7,1:m/2)=R_ratio(3,2:2:m);
R_ratio(8,1)=sum(R_ratio(3,2:2:m)>1);
R_ratio(8,2)=m/2-R_ratio(8,1);
R_ratio(8,3)=mean(R_ratio(3,2:2:m));
R_ratio(8,4)=max(R_ratio(3,2:2:m));
R_ratio(8,5)=min(R_ratio(3,2:2:m));
R_ratio(8,6)=std(R_ratio(7,1:m/2));
R_ratio(8,7)=std(abs(R_ratio(7,1:m/2-1)-R_ratio(7,2:m/2)));

end

function [V_initial,time,delta_t,delta_V,Vmax,Vmin,n_RS,square_size,v_resistance,v_density_Vs,v_total_res,v_current,v_time,voltage,v_temperature,R_ratio,V_set_reset,cont_SET_RESET,j,i,tmax,v_curr_map,v_den_map,v_ey,v_res_map,u,T]=variables_counters(parameters,Grid_S)

ni=size(Grid_S,1);
nj=size(Grid_S,2);

V_initial=parameters(2);
u=zeros(ni,nj); % Electric potential
T_experiment=parameters(25);
T = ones(ni,nj) * T_experiment;
time(1)=parameters(4);
time(2)=0;
delta_t=parameters(5);
delta_V=parameters(6);
Vmax=parameters(7);
Vmin=parameters(8);
n_RS=parameters(10);
square_size=parameters(17);
%size_vectors=(Vmax-Vmin)/delta_V;
if delta_V > 0
size_vectors=ceil(2*n_RS*(Vmax-Vmin)/delta_V);
else
size_vectors=100;
end
% Important variables
if square_size==0
v_resistance=zeros(nj,size_vectors);
v_density_Vs=zeros(nj,size_vectors);
else
%v_resistance=zeros(size(Grid_S,2)/square_size,size_vectors);
%v_density_Vs=zeros(size(Grid_S,2)/square_size,size_vectors);
%v_curr_map=zeros(size(Grid_S,1)/square_size,size(Grid_S,2)/square_size,size_vectors);
%v_den_map=zeros(size(Grid_S,1)/square_size,size(Grid_S,2)/square_size,size_vectors);

% Continuous density
v_resistance=zeros(nj,size_vectors);
v_density_Vs=zeros(nj,size_vectors);
v_curr_map=zeros(ni,nj,size_vectors);
v_den_map=zeros(ni,nj,size_vectors);
v_res_map=zeros(ni,nj,size_vectors);

v_ey=zeros(ni,nj,size_vectors);
end

v_total_res=zeros(size_vectors,1);
v_current=zeros(size_vectors,1);
v_time=zeros(size_vectors,1);
voltage=zeros(size_vectors,1);
v_temperature=zeros(size_vectors,2);

% Resistance ratio between LRS and HRS 
% R_ratio(1,:) --> ON RESISTANCE
% R_ratio(2,:) --> OFF RESISTANCE
% R_ratio(3,:) --> ON/OFF Ratio
% R_ratio(4,:) --> Flag to read the values only once
R_ratio=zeros(8,2*n_RS);

% V SET and V RESET
% V_set_reset(1,:) --> Vset
% V_set_reset(2,:) --> Step of the simulation j for Vset
% V_set_reset(3,:) --> Flag to save Vset only once
% V_set_reset(4,:) --> Vreset
% V_set_reset(5,:) --> Step of the simulation j for Vreset
% V_set_reset(6,:) --> Flag to save Vreset only once
V_set_reset=zeros(6,2*n_RS);
cont_SET_RESET=0;

j=1;
i=1;
tmax=0;

end
