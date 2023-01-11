% Function:
% - Calculate total current
% - Calculate map of current for every square
% 

function [v_current,current_mapping,f]=current(voltage,v_total_res,j,v_current,v_resistance,resistance_step,Grid_S,parameters,phy_const)

    [v_current]=total_current(voltage,v_total_res,j,v_current);

    [current_mapping]=current_map(v_resistance,v_current,j,resistance_step);

    [f]=joule_heating(current_mapping,resistance_step, Grid_S,parameters,phy_const);


end

function [v_current]=total_current(voltage,v_total_res,j,v_current)

v_current(j)=voltage(j)/v_total_res(j);

end

function [current_mapping]=current_map(v_resistance,v_current,j,resistance_step)

current_mapping=zeros(size(resistance_step));
voltage_row=v_current(j)*v_resistance(:,j);
for i=1:size(current_mapping,2)
current_mapping(:,i)=voltage_row(i)./resistance_step(:,i);
end

end

function [f]=joule_heating(current_mapping,resistance_step,Grid_S,parameters,phy_const)
    
square_size=parameters(17);
hx=phy_const(4);
hy=phy_const(5);
hz=phy_const(14);
quenching_heat=parameters(22);

%% We are averaging the density through a square with an area = hx*square_size*hy*square_size --> We assume that the current go through the cross section of this 
% Current density --> Current through the cross section (hx*square_size*hz)
% --> 2D: cross section is hx*square_size
current_mapping=current_mapping/(hx*square_size);
% ------- Resistivity=R*Area/length ----------------------------------------------
% A=hz*hx*square_size --> Cross section area
% l=hy*square_size --> length
% f is the energy per unit volumen

% --> Attention: We use this if resistance_step is resistance
%f=(hx*hz*square_size/(hy*square_size))*resistance_step.*current_mapping.^2;

% --> We use this if resistance step is resistivity --> Fox use resistivity
% Fox, Daniel S., Yangbo Zhou, Pierce Maguire, Arlene O’Neill, Cormac Ó’Coileáin, Riley Gatensby, Alexey M. Glushenkov et al. 
% "Nanopatterning and electrical tuning of MoS2 layers with a subnanometer helium ion beam." Nano letters 15, no. 8 (2015): 5307-5313.
f=(hy*square_size/(hx*hz*square_size))*resistance_step.*current_mapping.^2;

% Number of particles for each square --> We consider that the same heat is
% produced in each particle, no matter if it is a defect or not
%n_particles=sum(sum(Grid_S(1:square_size,1:square_size)>0));


% Turn f matrix in a matrix of Grid_S dimension
%    for i=1:square_size:size(Grid_S,1)
%        for j=1:square_size:size(Grid_S,2)
%    f(i:(i+square_size-1),j:(j+square_size-1))=f_aux(1+(i-1)/square_size,1+(j-1)/square_size);
%        end
%    end
%f=quenching_heat*f/n_particles;
f=quenching_heat*f;
f(Grid_S==0)=0;

end

