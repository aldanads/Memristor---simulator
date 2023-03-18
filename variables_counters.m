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

    if parameters(24) == 1
    size_vectors = size_vectors + 30;
    end

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