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

max_sim=6;
for n_sim=1:max_sim

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
        delta_t = parameters(26);
    elseif parameters(24) == 2 % Ramp experiment
        delta_V = parameters(6);
        delta_t = parameters(5);
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
%v_temperature(j,2)
%% KMC algorithm
[Grid_S,time,Vs_ij,prob]=KMC(Grid_S,phy_const,Vs_ij,ex,ey,ActE,T,time,parameters,tmax);


end

switch s
    
    % Ramp voltage increases until max --> Then decreases until min
    case 1

        % When the resistance degradate that amount, we stop the annealing
        % experiment and start the voltage ramp exp
        % parameters(26) = degradation level
        %if (v_total_res(j) <= v_total_res(1)*parameters(26)) && (parameters(24) == 11 || parameters(24) == 1)
        if (time(1) >= parameters(26)) && (parameters(24) == 11 || parameters(24) == 1)

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
if (V==0) && (parameters(24) == 2 || parameters(24) == 22)
    cont_SET_RESET=cont_SET_RESET+1;
    R_ratio(4,1)=0;
    R_ratio(4,2)=0;
end

if parameters(24) == 2
[R_ratio,V_set_reset]=measure(R_ratio,voltage,delta_V,v_resistance,j,V_set_reset,cont_SET_RESET);
end


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

exit;







