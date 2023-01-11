% Choose one curve each 15

%figure;
%k=1;
%for i=1:100:1000
%volt((1+(k-1)*197):k*197)=I_V_exp(1+197*(i-1):197*i,1);
%curr((1+(k-1)*197):k*197)=I_V_exp(1+197*(i-1):197*i,2);
%k=k+1;
%end
%semilogy(volt,abs(curr))


%I_V_exp_2(:,1)=I_V_exp((I_V_exp(:,4)<1E9) & (I_V_exp(:,4)>1E6),1);
%I_V_exp_2(:,2)=I_V_exp((I_V_exp(:,4)<1E9) & (I_V_exp(:,4)>1E6),2);
%I_V_exp_2(:,3)=I_V_exp((I_V_exp(:,4)<1E9) & (I_V_exp(:,4)>1E6),3);
%I_V_exp_2(:,4)=I_V_exp((I_V_exp(:,4)<1E9) & (I_V_exp(:,4)>1E6),4);

%figure;semilogy(I_V_exp_2(1:100000,3),I_V_exp_2(1:100000,4))

load('I_V_exp_1400cycles.mat')

res_sup_limit=1E9;
res_inf_limit=1E6;
n_points=5E3;
%[I_V_exp_2]=screening_data(I_V_exp,res_sup_limit,res_inf_limit,n_points);

[R_ratio]=calculate_res_ratio(I_V_exp);

[p,delta,S,x2,y,y_fit_2,delta_2]=linear_fitting(R_ratio);

%[SET_power_consumption,RESET_power_consumption]=power_consumption(I_V_exp);

function [R_ratio]=calculate_res_ratio(I_V_exp)

R_ratio=zeros(5,2400);
cont_SET_RESET=0;
voltage=I_V_exp(:,1);
Res=I_V_exp(:,4);
delta_V=0.7;
R_ratio(4,1)=1;
R_ratio(4,2)=1;
flag = 1;
%for j=1:size(I_V_exp)
for j=1:236100

% When it is starting the negative curve    
if (voltage(j) < 0 && flag == 1 && delta_V < 0)
    cont_SET_RESET=cont_SET_RESET+1;
    R_ratio(4,1)=0;
    R_ratio(4,2)=0;
    flag = 2;
end

% Only allow to enter again in the previous if after the positive curve
if voltage(j) > 0
    flag = 1;
end

if j>1
   delta_V=voltage(j)-voltage(j-1);
end


if (voltage(j)<0)&&(voltage(j)<=-10)&&(delta_V<0)&&(R_ratio(4,1)==0)
R_ratio(1,cont_SET_RESET)=Res(j);
R_ratio(4,1)=1;
end

if (voltage(j)<0)&&(voltage(j)>=-10)&&(delta_V>0)&&(R_ratio(4,2)==0)
R_ratio(2,cont_SET_RESET)=Res(j);
% Ratio OFF/ON should be greater than 1
R_ratio(3,cont_SET_RESET)=R_ratio(1,cont_SET_RESET)/R_ratio(2,cont_SET_RESET);
R_ratio(4,2)=1;
end
end


m=size(R_ratio,2);
cycles=sum(R_ratio(3,:)>0);

R_ratio(5,1)=sum(R_ratio(3,:)>1);
R_ratio(5,2)=m/2-R_ratio(5,1);
R_ratio(5,3)=mean(R_ratio(3,1:cycles));
R_ratio(5,4)=max(R_ratio(3,:));
R_ratio(5,5)=min(R_ratio(3,1:cycles));
R_ratio(5,6)=std(R_ratio(3,1:cycles));
R_ratio(5,7)=std(abs(R_ratio(3,1:cycles-1)-R_ratio(3,2:cycles)));

%figure;plot(R_ratio(3,1:cycles))

end

function [p,delta,S,x2,y,y_fit_2,delta_2]=linear_fitting(R_ratio)

initial_point=223;
x=initial_point:R_ratio(5,1);
cycles=R_ratio(5,1);
res_ratio=R_ratio(3,initial_point:cycles);
[p, S]=polyfit(x,res_ratio,1);
[y_fit, delta]=polyval(p,x,S);


x2=initial_point:cycles;
y=p(1)*x2+p(2);
figure;
plot(x,R_ratio(3,initial_point:cycles),x2,y,x2,y+2*mean(delta),'m--',x2,y-2*mean(delta),'m--')

log_x=log(x);
log_y=log(res_ratio);
[p, S]=polyfit(log_x,log_y,1);
[ylog_fit, delta_2]=polyval(p,log_x,S);

% Log representation 
log_x2=log(x2);
a=exp(p(2));
y_fit_2=a*exp(p(1)*log_x2);

figure;semilogx(x,R_ratio(3,initial_point:cycles),x2,y_fit_2, x2,y_fit_2*exp(2*mean(delta_2)),'m--',x2,y_fit_2/exp(2*mean(delta_2)),'m--')
figure;plot(x,R_ratio(3,initial_point:cycles),x2,y_fit_2, x2,y_fit_2*exp(2*mean(delta_2)),'m--',x2,y_fit_2/exp(2*mean(delta_2)),'m--')

end

function [I_V_exp_2]=screening_data(I_V_exp,res_sup_limit,res_inf_limit,n_points)
I_V_exp_2(:,1)=I_V_exp((I_V_exp(:,4)<res_sup_limit) & (I_V_exp(:,4)>res_inf_limit),1);
I_V_exp_2(:,2)=I_V_exp((I_V_exp(:,4)<res_sup_limit) & (I_V_exp(:,4)>res_inf_limit),2);
I_V_exp_2(:,3)=I_V_exp((I_V_exp(:,4)<res_sup_limit) & (I_V_exp(:,4)>res_inf_limit),3);
I_V_exp_2(:,4)=I_V_exp((I_V_exp(:,4)<res_sup_limit) & (I_V_exp(:,4)>res_inf_limit),4);

figure;semilogy(I_V_exp_2(1:n_points,3),I_V_exp_2(1:n_points,4))
xlabel('Time (s)')
ylabel('Resistance (Ohms)')
figure;plot(I_V_exp_2(1:n_points,3),I_V_exp_2(1:n_points,1))

end

function [SET_power_consumption,RESET_power_consumption]=power_consumption(I_V_exp)

SET_power_consumption = zeros(1201,1);
RESET_power_consumption = zeros(1201,1);
j=1;

for i=1:size(I_V_exp(:,1))

if I_V_exp(i,1) >= 0
SET_power_consumption(j)=SET_power_consumption(j)+abs(I_V_exp(i,1)*I_V_exp(i,2));
flag=1;
else 
RESET_power_consumption(j)=RESET_power_consumption(j)+abs(I_V_exp(i,1)*I_V_exp(i,2));
delta_V=I_V_exp(i,1)-I_V_exp(i-1,1);
    
    if (I_V_exp(i,1) < 2*abs(delta_V)) && (delta_V > 0) && (flag == 1) 
    j=j+1;
    flag=2;
    end
end

end

end
