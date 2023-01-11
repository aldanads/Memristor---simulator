

function [T]=SolveHeat(f,phy_const,Grid_S,parameters,time)

heat_equation_act=parameters(23);


% Simulation domain
nx=phy_const(10);
ny=phy_const(11);

% Grid step
hx=phy_const(4);
hy=phy_const(5);
hz=phy_const(14);

% Thermal conductivity
KMoS2=phy_const(12);
KMo=phy_const(13);
% Thermal diffusivity
alpha_MoS2=phy_const(15);
alpha_Mo=phy_const(16);

%Initial temperature
T_ambient=parameters(3);
T_top_electrode=parameters(20);
T_bottom_electrode=parameters(21);

% Convergence parameter
nunk=(nx-2)*(ny-2);
% Tolerance parameter of convergence
eps=0.1;

% Thermal conductivity matrix
[k]=Kmatrix(nx,ny,Grid_S,KMoS2,KMo);
% Thermal diffusivity matrix
[alpha]=alpha_matrix(nx,ny,Grid_S,alpha_MoS2,alpha_Mo);

% Temperature matrix
[T]=Temp_matrix(nx,ny,T_ambient);

% Dirichtlet boundary conditions
[T]=dirichlet(T,T_top_electrode,T_bottom_electrode);
% Loop iterations
max_iter=200;

if heat_equation_act == 1
%Constant initialization
Const(2)=hy^2;
Const(3)=hx^2;
Const(1)=1/(2*(Const(2)+Const(3)));
% f is the energy per unit volumen
idx=k>0;
const_heat_source=zeros(nx,ny);
const_heat_source(idx)=(f(idx)./k(idx))*(Const(2))*(Const(3));

for iter=1:max_iter
    % Neumann boundary conditions
    [T]=neumann(T,nx,iter,T_ambient);
    % Solving the heat equation
    [T,numi]=heat_equation(nx,ny,const_heat_source,Const,T,eps);
    if numi==nunk
    break
    end

end

else 
    if heat_equation_act == 2
    %Const(1)=(hx*1E9)^2;
    %Const(2)=(hy*1E9)^2;
    % Temperature per unit surface
    square_size=parameters(17);
    Const(1)=hx^2/(hx*square_size*hy*square_size);
    Const(2)=hy^2/(hx*square_size*hy*square_size);
    Const(3)=hz^2/(hx*square_size*hy*square_size);
    alpha_time=time(2)*alpha;

   
    %% TAKE CARE --> h in nm and alpha in (m**2/s)

    idx=k>0;
    const_heat_source=zeros(nx,ny);
    % f is the energy per unit volume
    % const_heat_source is heat flux= f * volume / cross_section
    const_heat_source(idx)=(hx*hy/(hx))*f(idx)./k(idx); 
    end
for iter=1:max_iter
    % Neumann boundary conditions
    [T]=neumann(T,nx,iter,T_ambient);
    % Solving the heat equation
    [T,numi]=heat_equation_t(nx,ny,T,eps,Const,const_heat_source,alpha_time);
    if numi==nunk
    break
    end
end


end
end

function [k]=Kmatrix(nx,ny,Grid_S,KMoS2,KMo)

k=ones(nx,ny)*KMoS2;
%k=zeros(nx,ny);
%k(Grid_S==1)=KMoS2;
k(Grid_S==2)=KMo;

end

function [alpha]=alpha_matrix(nx,ny,Grid_S,alpha_MoS2,alpha_Mo)

alpha=ones(nx,ny)*alpha_MoS2;
%alpha(Grid_S==1)=alpha_MoS2;
alpha(Grid_S==2)=alpha_Mo;

end

function [T]=Temp_matrix(nx,ny,T_ambient)

T=ones(nx,ny)*T_ambient;

end

function [T]=dirichlet(T,T_top_electrode,T_bottom_electrode)

T(:,1)=T_bottom_electrode;
T(:,size(T,2))=T_top_electrode;

end

function [T]=neumann(T,nx,iter,T_ambient)

if iter==1
T(1,:)=T_ambient;
T(nx,:)=T_ambient;
else
T(1,:)=(4*T(2,:)-T(3,:))/3;
T(nx,:)=(4*T(nx-1,:)-T(nx-2,:))/3;
end

end

function [T,numi]=heat_equation(nx,ny,const_heat_source,Const,T,eps)
%Nodes converged
numi=0;

    for i=2:nx-1
        for j=2:ny-1
        T_aux=Const(1)*(Const(2)*(T(i+1,j)+T(i-1,j))+Const(3)*(T(i,j+1)+T(i,j-1))+const_heat_source(i,j));
    
        % error
        error=abs(T_aux-T(i,j));
        T(i,j)=T_aux;

            if error<eps
            numi=numi+1;
            end

        end
    end


end

function [T,numi]=heat_equation_t(nx,ny,T,eps,Const,const_heat_source,alpha_time)

Tn=T;
numi=0;

for i=2:nx-1
    for j=2:ny-1
        % Last component (600-2Tn) is to implement a heat sink under the layer --> Due to the substrate
        % I assume that the diffusivity of Si is 200 times the diffusivity of MoS2
        T(i,j)=Tn(i,j)+alpha_time(i,j)*((Tn(i+1,j)+Tn(i-1,j)-2*Tn(i,j))/Const(1)+(Tn(i,j+1)+Tn(i,j-1)-2*Tn(i,j))/Const(2)+200*(600-2*Tn(i,j))/Const(3)+const_heat_source(i,j));
        %T(i,j)=Tn(i,j)+alpha_time(i,j)*((Tn(i+1,j)+Tn(i-1,j)-2*Tn(i,j))/Const(1)+(Tn(i,j+1)+Tn(i,j-1)-2*Tn(i,j))/Const(2)+const_heat_source(i,j));
       
        % error
        error=abs(Tn(i,j)-T(i,j));
        %T(i,j)
        if error < eps
            numi=numi+1;
        end
    end
end


end
