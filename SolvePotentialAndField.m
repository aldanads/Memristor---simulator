%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This algorithm solve electric potential and Poisson equation through the
% iterative process of Gauss-Seidel with first order central differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% Grid_S -> Grid with sulfur and sulfur vacancies
% er -> dielectric permitivity
% tol -> Tolerance, convergence parameter
% V -> Applied voltage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
% u -> Electric potential
% ex -> Electric field - X
% ey -> Electric field - Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [u,ex,ey,screening_j]=SolvePotentialAndField(Grid_S,phy_const,Vs_ij,parameters,V,density_Vs,resistance_step,screening_fit_parameters)

ni=size(Grid_S,1);
nj=size(Grid_S,2);
%Electril potential
u=zeros(ni,nj);

num=round((ni*nj)^(1/2));

% Charge distribution
[ro_eps]=charge(Grid_S,phy_const,Vs_ij);
% Dirichtlet conditions
[u]=dirichtlet(u,V,parameters);
% Auxiliar variable
u1=u;

% count
cont=0;

% Tolerance
tol=parameters(1);

%Constant initialization --> Poisson
h1=phy_const(4)*1E9;
h2=phy_const(5)*1E9;
er=phy_const(6);
cte=(1/(2*(h1^2+h2^2)));
cte_ro=(h1^2)*(h2^2)*ro_eps/er;
h1_2=h1^2;
h2_2=h2^2;
Const(1)=cte;
Const(2)=h1_2;
Const(3)=h2_2;

%% Algorithm

while lt(cont,100)
     u2=u;
        if lt(cont,40)
            NN = 60*num;
        end
        if and(ge(cont,40),lt(cont,70))
            NN = 30*num;
        end
        if and(ge(cont,70),lt(cont,90))
            NN = 15*num;
        end
        if and(ge(cont,90),lt(cont,100))
            NN = 5*num;
        end
        for i=1:NN
            [u]=EcPoisson(u,u1,cte_ro,Const,ni,nj);
            % Neuman boundary conditions
            [u]=neuman(u,u1);
             u1=u;
        end
        % Check convergence
        [cont]=checkf(u,u2,tol);
    
end

% The screening factor depends on the concentration of defects. So
% we calculate it for every row.
[screening_j]=spatial_screening(density_Vs,resistance_step,parameters,screening_fit_parameters,Grid_S); 
%for j=1:size(u,2)
%    u(:,j)=u(:,j)*screening_j(j);
%end
%Calculate electric field
[ex,ey]=Electric_field(u,phy_const,screening_j);

    



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------- Functions for Poisson equation ------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Function to distribute the charge in the grid
    function [ro_eps]=charge(Grid_S,phy_const,Vs_ij)
        ni=size(Grid_S,1);
        nj=size(Grid_S,2);
        h1=phy_const(4);
        h2=phy_const(5);
        er=phy_const(6);
        ro_eps=zeros(ni,nj);
        q=phy_const(2);
        screening=phy_const(9);

        % The screening factor depends on the concentration of defects. So
        % we calculate it for every row.
        %[screening_j]=spatial_screening(phy_const,density_Vs); 

        % Charge of an electron (C)
        qe=1.602176565e-19;
        % Vacuum permittivity (F/m)
        e0=8.8541878176e-12;
        % Volume --> phy_const(4)=h(1); phy_const(5)=h(2)
        vol=h1*h2;
        % Charge density phy_const(2)=q;
        % (V/m) --> For some reason, it's not dimensionally correct
        density=q*qe/(vol*e0);
        % nm better
        density=density*1E-9;

        for i=1:size(Vs_ij,1)
            %ro_eps(Vs_ij(i,1),Vs_ij(i,2))=screening_j(Vs_ij(i,2))*density/er;
            ro_eps(Vs_ij(i,1),Vs_ij(i,2))=screening*density/er;
        end

    end

    % Dirichtlet boundary conditions at both extreme of the channel.
    % Constant voltage is applied.
    function [u]=dirichtlet(u,V,parameters)
        polarity=parameters(16);
        %% External voltage in the drain
        if polarity==1
        u(:,1)=V;
        
        u(:,size(u,2))=0;
        else
        u(:,1)=0;
        
        u(:,size(u,2))=V;
        end
    end

    % Solving Poisson equation using finite differences technique
    function [u]=EcPoisson(u,u1,cte_ro,Const,ni,nj)

    %% h1~=h2 --> This change the finite difference technique for Poisson equation
    %% Only interior points, the boundary points with boundary conditions

    for i=2:ni-1
        for j=2:nj-1
            u(i,j)=Const(1)*(cte_ro(i,j)+(Const(3))*(u1(i+1,j)+u1(i-1,j))+(Const(2))*(u1(i,j+1)+u1(i,j-1)));
        end
    end

        
    end

    % Neuman conditions
    % (ver: https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference)
    function [u]=neuman(u,u1)
    ni=size(u,1);
    nj=size(u,2);

    %for j=2:size(u,2)-1
        % Forward finite difference --> error ~hx
     %   u(1,j)=u1(2,j);
     %   u(ni,j)=u1(ni-1,j);
        
        % See Lecture4 in Bibliography Ghost points --> Polynomial fitting and Taylor serie:
        % --> error ~hx^3 --> WARNING! NOT WORKING!The potential diverges!
        %u(1,j)=(4*u1(2,j)-u1(3,j))/3;
        %u(ni,j)=(4*u1(ni-1,j)-u1(ni-2,j))/3;
        
        %error_1=(aux_1-u(1,j))/aux_1
        %error_ni=(aux_ni-u(ni,j)/aux_ni)   
    %end

     u(1,2:nj-1)=u1(2,2:nj-1);
     u(ni,2:nj-1)=u1(ni-1,2:nj-1);

    end

    %% Function to check how many points have converged
    function [cont]=checkf(u,u2,tol)
    ni=size(u,1);
    nj=size(u,2);
    
    %% Points to converge
    total=ni*(nj-2);
    
    cont=0;
    for j=2:nj-1
        for i=1:ni
            if ne(u(i,j),0)
               diff=abs((u(i,j)-u2(i,j))/u(i,j)); 
            else
               diff=abs(u(i,j)-u2(i,j)); 
            end
            
            if gt(tol,diff)
               cont=cont+1; 
            end
            
            
        end
    end
    
    cont=100*cont/total;
    
    end

    %% Calculate electric field
    function [ex,ey]=Electric_field(u,phy_const,screening_j)

    h1=phy_const(4);
    h2=phy_const(5);

    [ey,ex]=gradient(u);
    ey=-ey.*screening_j/h2;
    ex=-ex.*screening_j/h1;


    %ey=zeros(ni,nj-1);
    %ex=zeros(ni-1,nj);
    %for i=1:ni-1
    %    for j=1:nj
    %            ex(i,j) = -screening_j(i,j)*(u(i+1,j)-u(i,j))/h1;
    %    end
    %end

    %for i=1:ni
    %    for j=1:nj-1
    %            ey(i,j) = -screening_j(i,j)*(u(i,j+1)-u(i,j))/h2;
    %    end
    %end

    %% ex and ey are  ex=zeros(ni-1,nj) and ey=zeros(ni,nj-1); 
    %% We consider that at the border ni and nj the electric field is the same.
    %exx=zeros(size(Grid_S,1),size(Grid_S,2));
    %eyy=zeros(size(Grid_S,1),size(Grid_S,2));

    %exx(1:(size(Grid_S,1)-1),:)=ex;
    %exx(size(Grid_S,1),:)=ex(size(ex,1),:);
    %ex=exx;

    %eyy(:,1:(size(Grid_S,2)-1))=ey;
    %eyy(:,size(Grid_S,2))=ey(:,size(ey,2));
    %ey=eyy;





        
    end

    function [screening_j]=spatial_screening(density_Vs,resistance_step,parameters,screening_fit_parameters,Grid_S)
    
    screening_mechanism=parameters(14);
    

    %% Screening based on defect density
    if screening_mechanism==1
    
    density=ones(2,size(density_Vs,1));
    density(2,:)=density_Vs;
    density=density';

    
    % Matrix and vector multiplication is more efficient in MATLAB than
    % "for" bucle.
    screening_j=density*screening_fit_parameters;


    %% Screening based on resistivity
    else
    square_size=parameters(17);

        if square_size==0
    resistance=ones(2,size(resistance_step,1));
    resistance(2,:)=resistance_step;
    resistance=resistance';
   
    % Matrix and vector multiplication is more efficient in MATLAB than for
    % bucle.
    screening_j=resistance*screening_fit_parameters;
        else
    b=screening_fit_parameters(1);
    slope=screening_fit_parameters(2);
    screening_j=resistance_step.*slope+b;
    
    %% Discrete density
    %screening_j_aux=resistance_step.*slope+b;
    %screening_j=zeros(size(Grid_S));
    % Turn screening_j into a matrix with the grid dimension
    %for i=1:square_size:size(Grid_S,1)
    %    for j=1:square_size:size(Grid_S,2)
    %screening_j(i:(i+square_size-1),j:(j+square_size-1))=screening_j_aux(1+(i-1)/square_size,1+(j-1)/square_size);
    %    end
    %end
    
        end

    end


    end
