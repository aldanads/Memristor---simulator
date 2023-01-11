

function [J]=exp_sim_fit(n_RS,v_current,voltage,parameters)

load I_V_exp_1400cycles.mat;
m=length(I_V_exp);
J=zeros(n_RS,2);

Vmax=parameters(7);
Vmin=parameters(8);
delta_V=parameters(6);

cycle_steps=2*(Vmax-Vmin)/delta_V;

figure;semilogy(voltage,abs(v_current),I_V_exp(1:cycle_steps,1),abs(I_V_exp(1:cycle_steps,2)))
if m==size(v_current,1)
for i=1:n_RS
J(i,1)=sum((v_current(((i-1)*m+1):i*m)-I_V_exp(1:m,2)).^2)/(2*m);
end
end
J(1,2)=min(J(:,1));

end