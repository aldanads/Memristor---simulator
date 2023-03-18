function [R_ratio,V_set_reset]=measure(R_ratio,voltage,delta_V,v_resistance,j,V_set_reset,cont_SET_RESET)

sum_res_j=sum(v_resistance(:,j));

[R_ratio]=R_ratio_linear(voltage,delta_V,cont_SET_RESET,R_ratio,sum_res_j,j);

if j>1
[V_set_reset]=Looking_for_V_set_reset(v_resistance,voltage,sum_res_j,V_set_reset,cont_SET_RESET,j);
end


end

function [R_ratio]=R_ratio_linear(voltage,delta_V,cont_SET_RESET,R_ratio,sum_res_j,j)

%% Ratio between ON and OFF state
% We calculate the ratio with the values in the RESET curve during the
% the sweep 0 --> Vmin (ON Resistance) and the sweep Vmin --> 0 (OFF Resistance)
% V=-4 is low enough to not drag the defects
if (voltage(j)<0)&&(voltage(j)<=-4)&&(delta_V<0)&&(R_ratio(4,1)==0)
R_ratio(1,cont_SET_RESET)=sum_res_j;
R_ratio(4,1)=1;
end

if (voltage(j)<0)&&(voltage(j)>=-4)&&(delta_V>0)&&(R_ratio(4,2)==0)
R_ratio(2,cont_SET_RESET)=sum_res_j;
% Ratio OFF/ON should be greater than 1
R_ratio(3,cont_SET_RESET)=R_ratio(2,cont_SET_RESET)/R_ratio(1,cont_SET_RESET);
R_ratio(4,2)=1;
end

if (voltage(j)>0)&&(voltage(j)>=4)&&(delta_V>0)&&(R_ratio(4,1)==0)
R_ratio(1,cont_SET_RESET)=sum_res_j;
R_ratio(4,1)=1;
end

if (voltage(j)>0)&&(voltage(j)<=4)&&(delta_V<0)&&(R_ratio(4,2)==0)
R_ratio(2,cont_SET_RESET)=sum_res_j;
R_ratio(3,cont_SET_RESET)=R_ratio(2,cont_SET_RESET)/R_ratio(1,cont_SET_RESET);
R_ratio(4,2)=1;
end
end

function [V_set_reset]=Looking_for_V_set_reset(v_resistance,voltage,sum_res_j,V_set_reset,cont_SET_RESET,j)

% Percentage of resistance difference
dif_res=(sum(v_resistance(:,j-1))-sum_res_j)/sum_res_j;

% It mays change between ON and OFF state many times in each branch, so if
% a SET occurs, I allow to save a new RESET. If a RESET occurs, I allow to
% save a new SET. If the RESET and SET are stable, only the first V is
% saved.
% Vset
if (cont_SET_RESET<=size(V_set_reset,2))&&(dif_res>0.09)&&(V_set_reset(3,cont_SET_RESET)==0)
V_set_reset(1,cont_SET_RESET)=voltage(j-1);
V_set_reset(2,cont_SET_RESET)=j-1;

% If the device is in ON state, it might be a RESET
V_set_reset(3,cont_SET_RESET)=1;
V_set_reset(6,cont_SET_RESET)=0;
end
% Vreset
if (cont_SET_RESET<=size(V_set_reset,2))&&(dif_res<-0.1)&&(V_set_reset(6,cont_SET_RESET)==0)
V_set_reset(4,cont_SET_RESET)=voltage(j-1);
V_set_reset(5,cont_SET_RESET)=j-1;

% If the device is in OFF state, it might be a SET
V_set_reset(6,cont_SET_RESET)=1;
V_set_reset(3,cont_SET_RESET)=0;

end
end

