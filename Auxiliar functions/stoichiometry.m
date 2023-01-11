% Stoichiometry by row

%stoich=zeros(size(Grid_S,2),1);
% The MoS2 stoichiometry is 1:2 (Mo:S)
% We do not implement Mo, so we don't change the number.
% Sulfur is implemented with 1
% Vs is implemente with 2
% So the number of Mo is the half of the number of sulfur or Vs
%num_Mo=sum(Grid_S(1,:)>0)/2;
% Averaged in the y-axis, along x-axis
%num_S=sum(Grid_S(:,:)==1,1);
%stoich=num_S./(num_S+num_Mo);

% Stoichiometry by square
% Max defect density in a nm^2 = 5.64
num_Mo=5.64/2;
v_density_S=5.64-v_density_Vs;
stoich_S=v_density_S./(v_density_S+num_Mo);
stoich_Mo=num_Mo./(v_density_S+num_Mo);

h2=phy_const(5)*1E9;
y=h2:h2:size(Grid_S,2)*h2;

figure;plot(y,stoich_S(:,1))



