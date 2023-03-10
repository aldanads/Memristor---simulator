%-----------------------------------------------------------------%
% Simulator: Lateral memtransistors from monolayer single-crystal %
% molybdenum disulfide (MoS2)                                     %
% Starting Date: 2021/11/16                                       %
% Samuel Aldana Delgado                                           %
%-----------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Purpose: This function is the kinetic Monte Carlo algorithm. 
%% This algorithm is used to select the occurring events.

function [Grid_S,time,Vs_ij,prob]=KMC(Grid_S,phy_const,Vs_ij,ex,ey,ActE,T,time,parameters,delta_V)

%Time max for step -> The event with bigger time step
time_max=0;
%-> delta_t=parameters(5);
delta_t=parameters(5);
annealing_time = parameters(26);

% Probability of each process
prob=zeros(1,4);

for i=1:size(Vs_ij,1)

    %% Initialization of transisiton rates --> Set to 0 for every particles
    TR(1)=0;
    TR(2)=0;
    TR(3)=0;

    %% We select the allowed movements
    [mov_banner]=mov_screening(Grid_S,Vs_ij(i,1),Vs_ij(i,2));
    
    % Electric field due to a particle --> There are fields in every direction, 
    % but in different grid points, whilst the particle is in a single grid point. 
    % This is a problem because make the particle interact with itself
    if Vs_ij(i,1)>1
       E_field_x=ex(Vs_ij(i,1)-1,Vs_ij(i,2))+ex(Vs_ij(i,1),Vs_ij(i,2));
    else
       E_field_x=ex(Vs_ij(i,1),Vs_ij(i,2));
    end
    
    if Vs_ij(i,2)>1
       E_field_y=ey(Vs_ij(i,1),Vs_ij(i,2)-1)+ey(Vs_ij(i,1),Vs_ij(i,2));
    else
       E_field_y=ey(Vs_ij(i,1),Vs_ij(i,2));
    end
    
    [ActE]=Act_energy(ActE,mov_banner,phy_const,E_field_x,E_field_y);
 
    if (mov_banner(1)==1)
        TR(1)=Trans_Rate(ActE(1),T(Vs_ij(i,1),Vs_ij(i,2)),phy_const(1));
    end
    
    if (mov_banner(2)==1)
        TR(2)=Trans_Rate(ActE(2),T(Vs_ij(i,1),Vs_ij(i,2)),phy_const(1));
    end
    
    if (mov_banner(3)==1)
        TR(3)=Trans_Rate(ActE(3),T(Vs_ij(i,1),Vs_ij(i,2)),phy_const(1));
    end
    %% We sum all the transition rates and we multiply by a random number. 
    %% This is, we select a point in a bar whose length is the sum of all TRs.
    sumTR=sum(TR)*rand;
    %% Continue to the next iteration of the loop 
    if (sumTR==0)
        continue;
    end
    
    %% We select the event through kMC algorithm
    aux=TR(1);
    s=1;
    while (aux<=sumTR)
         s=s+1;
         aux=aux+TR(s);
    end
    %% Time for at least one occurring event
    time(2)=-log(rand)/sum(TR);
    
    if delta_V ~= 0
    if time(2)>delta_t
        time(2)=-log(rand)/TR(s);
        
        if time(2)>delta_t
           time(2)=delta_t; 
        end
    end
    %The event with bigger time step
    if time(2)>time_max
       time_max=time(2); 
    end

    else 
        if time(1)+time(2) > annealing_time
            time(2) = annealing_time - time(1);
        end
    
    end
    
    P=1-exp(-TR(s)*time(2));
    if (rand<P)
    prob(s)=prob(s)+1;
    prob(4)=prob(4)+1;
    [Grid_S,Vs_ij]=processes(Grid_S,i,Vs_ij,s,mov_banner(4));
    end

end

% Update time
time(1)=time(1)+time(2);
if prob(4)>0
prob(1)=prob(1)/prob(4);
prob(2)=prob(2)/prob(4);
prob(3)=prob(3)/prob(4);
end

%% Processes involved in the system --> We apply the event selected by kMC
function [Grid_S,Vs_ij]=processes(Grid_S,ptr,Vs_ij,s,mov_banner)
    
    %% We use k instead of i for not mixing with previous notation
    k=Vs_ij(ptr,1);
    j=Vs_ij(ptr,2);
    
switch mov_banner
    
    %% Right side particle
    case 2
        
    switch s
        %% Movement to the right-up side
        case 1
        Grid_S(k+1,j+1)=2;
        Grid_S(k,j)=1;
        % Update coordinates
        Vs_ij(ptr,1)=k+1;
        Vs_ij(ptr,2)=j+1;
        
        %% Movement to the right-down side
        case 2
        Grid_S(k+1,j-1)=2;
        Grid_S(k,j)=1;
        % Update coordinates
        Vs_ij(ptr,1)=k+1;
        Vs_ij(ptr,2)=j-1;
        
        %% Movement to the left side
        case 3    
        Grid_S(k-1,j)=2;
        Grid_S(k,j)=1;
        % Update coordinates
        Vs_ij(ptr,1)=k-1;
 

    end
    
    %% Left side particle
    case 3
        
    switch s
        %% Movement to the left-up side
        case 1
        Grid_S(k-1,j+1)=2;   
        Grid_S(k,j)=1;
        % Update coordinates
        Vs_ij(ptr,1)=k-1;
        Vs_ij(ptr,2)=j+1;
                
        %% Movement to the left-down side
        case 2
        Grid_S(k-1,j-1)=2;   
        Grid_S(k,j)=1;
        % Update coordinates
        Vs_ij(ptr,1)=k-1;
        Vs_ij(ptr,2)=j-1;
                
        %% Movement to the right side
        case 3
        Grid_S(k+1,j)=2;
        Grid_S(k,j)=1;
        % Update coordinates
        Vs_ij(ptr,1)=k+1;
 
    end
end
    
end

%%  Movement processes screening 
function [mov_banner]=mov_screening(Grid_S,i,j)

%% Movement_banner --> Each particle can only take three paths 
%% 1 --> movement allowed
%% 0 --> movement restricted
%% 2 --> right particle
%% 3 --> left particle

ni=size(Grid_S,1);
nj=size(Grid_S,2);

%% mov toward up
mov_banner(1)=1;
%% mov toward down
mov_banner(2)=1;
%% horizontal mov
mov_banner(3)=1;
    %if ((i==1)&&(Grid_S(i+1,j))==0)||(i>1)&&(Grid_S(i-1,j)==1)||(i>1)&&(Grid_S(i-1,j)==2)
    if ((i<ni)&&(Grid_S(i+1,j))==0)||((i==ni)&&(Grid_S(i-1,j))~=0)
%%                  O X X O  
%%                  X O O X 
%%                  O X X O 
%%                     (?)
%%                  X O O X 
    %% The particle is in the right side
    %% mov_banner(4) = 2--> right side; 3 --> left side
    mov_banner(4)=2;

        %% If it is the upper face or vacancy at the up-right
        if (j==nj)||(i==ni)||(Grid_S(i+1,j+1)==2)
            mov_banner(1)=0;
        end
        
        %% If it is the face on below side or vacancy at the down-right
        if (j==1)||(i==ni)||(Grid_S(i+1,j-1)==2)
            mov_banner(2)=0;
        end
        
        %% If it is the left face or vacancy at the left
        if (i==1)||(Grid_S(i-1,j)==2)
            mov_banner(3)=0;
        end
        
    else
%%                  O X X O  
%%                  X O O X 
%%                  O X X O 
%%                   (?)
%%                  X O O X 
    %% The particle is in the left side
    %% mov_banner(4) = 2--> right side; 3 --> left side
    mov_banner(4)=3;
            %% If it is the upper face or vacancy at the up-left
        if (j==nj)||(i==1)||(Grid_S(i-1,j+1)==2)
            mov_banner(1)=0;
        end
        
        %% If it is the face on below side or vacancy at the down-left
        if (j==1)||(i==1)||(Grid_S(i-1,j-1)==2)
            mov_banner(2)=0;
        end
        
        %% If it is the right face or vacancy at the right
        if (i==ni)||(Grid_S(i+1,j)==2)
            mov_banner(3)=0;
        end
            
    end
end

%% Function to calculate the activation energy
% ActE --> Vector with activation energies
% Grid_S --> Grid with sulfur and sulfur vacancies
% i, j --> coordinates
% phy_const(2) --> 2 (charge in e-)
% phy_const(3) --> theta --> Angle between sulfur positions in the
% hexagonal
% angle_E --> the angle of the electric field (where is pointing out)
function [ActE]=Act_energy(ActE,mov_banner,phy_const,ex,ey)    

    q=phy_const(2);
    theta=phy_const(3);
    h1=phy_const(4);
    h2=phy_const(5);
    h=sqrt(h1^2+h2^2);
    er=phy_const(6);
    % The angle of the electric field
    [angle_E]=arctan(ey,ex);

    %if j<=(parameters(11)-4)
    E_mov=ActE(4);
    %else
    %E_mov=ActE(6);
    %end

% Module of electric field
    E_field=sqrt(ex^2+ey^2);
    %% The only movements allowed is up-right down-right or to the left
    if (mov_banner(4)==2)
        % Towards field direction --> up-right
        if (mov_banner(1)==1)
        [local_field]=local_field_dielectric(E_field,er);
        ActE(1)=E_mov-local_field*q*cos(theta-angle_E)*h/2;
        %Inferior limitation
        if ActE(1)<ActE(5)
           ActE(1)=ActE(5); 
        end
        end
        
        % Opposite to field direction --> down-right
        if (mov_banner(2)==1)
        [local_field]=local_field_dielectric(E_field,er);
        ActE(2)=E_mov-local_field*q*cos(-theta-angle_E)*h/2;
        %Inferior limitation
        if ActE(2)<ActE(5)
           ActE(2)=ActE(5); 
        end
        end
        
        %% Horizontal movement --> to the left
        if (mov_banner(3)==1)
        [local_field]=local_field_dielectric(ex,er);
        ActE(3)=E_mov+local_field*q*h2/2; 
        
        %Inferior limitation
        if ActE(3)<ActE(5)
           ActE(3)=ActE(5); 
        end
        end

    else
    %% The only movements allowed is up-left, down-left or to the right
    pi=4*atan(1);
        % Towards field direction --> up-left
        if (mov_banner(1)==1)
        [local_field]=local_field_dielectric(E_field,er);
        ActE(1)=E_mov-local_field*q*cos(pi-theta-angle_E)*h/2;
        %Inferior limitation
        if ActE(1)<ActE(5)
           ActE(1)=ActE(5); 
        end
        end
        
        % Opposite to field direction --> down-left
        if (mov_banner(2)==1)
        [local_field]=local_field_dielectric(E_field,er);
        ActE(2)=E_mov-local_field*q*cos(pi+theta-angle_E)*h/2;
        %Inferior limitation
        if ActE(2)<ActE(5)
           ActE(2)=ActE(5); 
        end
        end
        
        %% Horizontal movement --> to the right
        if (mov_banner(3)==1)
        [local_field]=local_field_dielectric(ex,er);
        ActE(3)=E_mov-local_field*q*h2/2;  
        
        %Inferior limitation
        if ActE(3)<ActE(5)
           ActE(3)=ActE(5); 
        end
        end
                   

    end
    
    function [angle_E]=arctan(ey,ex)
       pi=4*atan(1);

       if (ey>=0)&&(ex>=0)
           angle_E=atan(ey/ex);
           
         if (ex==0)&&(ey==0)
              angle_E=0; 
         end
       end
       
       if (ey>=0)&&(ex<0)
           angle_E=pi+atan(ey/ex);
       end
       
       if (ey<0)&&(ex<0)
           angle_E=pi+atan(ey/ex);
       end
       
       if (ey<0)&&(ex>=0)
           angle_E=atan(ey/ex);
           if (ex==0)
               angle_E=-pi/2;
           end
       end
       
    end
    
    function [local_field]=local_field_dielectric(E_field,er)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Padovani, A., Larcher, L., Pirrotta, O., Vandelli, L., & Bersuker, G. (2015). 
    %% Microscopic modeling of HfO x RRAM operations: From forming to switching. 
    %% IEEE Transactions on electron devices, 62(6), 1998-2006.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For the expresion %
    
    % dipole moment (eA) --> 1Debay = 0.20819 eA
     %p=2.5E-10; --> HfO2
    % enm
    %Debay = 0.020819;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Klein, J., Wierzbowski, J., Regler, A., Becker, J., Heimbach, F., Muller, K., ... & Finley, J. J. (2016). 
    % Stark effect spectroscopy of mono-and few-layer MoS2. Nano letters, 16(3), 1554-1559.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dipole moment
    %p=1.4*Debay;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Wierzbowski, J., Klein, J., Kaniber, M., Müller, K., & Finley, J. J. 
    % Polarization control in few-layer MoS2 by electric field induced symmetry breaking.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dipole moment --> (enm)    %p=40;
    %p=40;
    p=15;
    
    L=1/3;
    relative_permittivity=er-1;
    local_field=E_field*(1+L*relative_permittivity)*p;
    end
    
end


%% Function to calculate transition rates
% Ea --> Activation energy
% T --> Temperature at this point in kelvin
% phy_const(1) --> nu0 (s^-1)
function [TR]=Trans_Rate(Ea,T,phy_const)
    
%Constante de Boltzmann (eV/K)
kb=8.6173324E-5;
%Transition rate
TR=phy_const*exp(-Ea/(kb*T));
end

end