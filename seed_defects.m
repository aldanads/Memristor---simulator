%-----------------------------------------------------------------%
% Simulator: Lateral memtransistors from monolayer single-crystal %
% molybdenum disulfide (MoS2)                                     %
% Starting Date: 2021/11/16                                       %
% Samuel Aldana Delgado                                           %
%-----------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Purpose: This function initialize the defects in the grid
%%
%% -----------------------------------------------------------------------
%% -------------------- Parameters ---------------------------------------
%% dose = probability (<1) to find a defect
%% dose_width = number of rows, starting at the end of the grid

function [Grid_S,Nparticles,Vs_ij]=seed_defects(Grid_S,Nparticles,distribution,parameter_seed)


%% Coordinates for sulfur vacancies
%Vs_ij=zeros(size(Grid_S,1)*dose_width,2);


%% Distribution 1: Perfect rows in x, no dispersion of defects. Centereded at d(2) (the end of the grid)
if (distribution==1)
    [Vs_ij,Grid_S,Nparticles]=perfect_rows(Grid_S,parameter_seed,Nparticles);
end

%% Distribution 2: Gaussian distribution of defects. Distributed in rows and centered at d(2) (the end of the grid)
if (distribution==2)
    [Vs_ij,Grid_S,Nparticles]=gaussian_distribution(Grid_S,parameter_seed,Nparticles);
end

%% Distribution 3: Test for electric field and interaction between particles
if ((distribution>=3)&&(distribution<6))
    [Vs_ij,Grid_S,Nparticles]=test_interaction_part(Grid_S,Nparticles,parameter_seed,distribution);    
end

if (distribution==6)
    [Vs_ij,Grid_S,Nparticles]=triangle(Grid_S,Nparticles,parameter_seed);
end


    % Perfect rows distribution of defects
    function [Vs_ij,Grid_S,Nparticles]=perfect_rows(Grid_S,parameter_seed,Nparticles)
    dose=parameter_seed(1);
    dose_width=parameter_seed(2);

    nj=size(Grid_S,2);
    start_irradiated_region=round(nj/2-dose_width);

    finish_irradiated_region=round(nj/2+dose_width);
    
    for i=1:size(Grid_S,1)

    %% dose_width determines how many rows will have defects starting from the last one
    % backward count
    for j=finish_irradiated_region:-1:(start_irradiated_region)
    if (Grid_S(i,j)==1)
        %% dose determines the probability of finding a defect    
        if (rand<dose)
        Grid_S(i,j)=2;
        Nparticles(2)=Nparticles(2)+1;
        Nparticles(1)=Nparticles(1)-1;
        Vs_ij(Nparticles(2),1)=i;
        Vs_ij(Nparticles(2),2)=j;
        end
    else continue;

    end
    end
    end
    
    end

    % Gaussian distribution of defects centered at the end of the grid in
    % axis j
    function [Vs_ij,Grid_S,Nparticles]=gaussian_distribution(Grid_S,parameter_seed,Nparticles)

    mean = parameter_seed(1);
    standard_deviation = parameter_seed(2);
    hy=parameter_seed(3);
    damaged_region=parameter_seed(4);
    irradiated_row=parameter_seed(5);
    steps_damreg=damaged_region/hy;
    Ly=size(Grid_S,2);
    max_dose=parameter_seed(6);
    skew=parameter_seed(7);
    kurt=parameter_seed(8);

    for j=1:Ly-1
    % +- 1 standard desviation --> 68% of the data
    % +- 2 standard desviation --> 95.4% of the data
    % We divide the standard desviation between the number of step in the
    % damaged region "steps_damreg"
    x=mean+(irradiated_row-j)*standard_deviation/steps_damreg;
    
    %% Normal distribution function
    % The normpdf returns the PDF in terms of %. To normalize to fraction, divide the output of normpdf by 100.
    %dose=max_dose*normpdf(x,mean,standard_deviation)/100;

    %% Pearson distribution
    % Cite this:
    % Pierce Brady (2022). pearspdf (https://www.mathworks.com/matlabcentral/fileexchange/26516-pearspdf), MATLAB Central File Exchange. Retrieved January 24, 2022. 
    dose = max_dose*pearspdf(x,mean,standard_deviation,skew,kurt)/100;
    % Along the complete row with the same probability "dose" for finding a
    % defect
    for i=1:size(Grid_S,1)
    if (Grid_S(i,Ly-j)==1)
        %% dose determines the probability of finding a defect    
        if (rand<dose)
        Grid_S(i,Ly-j)=2;
        Nparticles(2)=Nparticles(2)+1;
        Nparticles(1)=Nparticles(1)-1;
        Vs_ij(Nparticles(2),1)=i;
        Vs_ij(Nparticles(2),2)=Ly-j;
        end
    else continue;

    end
    end
    end

    % To know the probability between two limits --> Problem! We don't know the
    % limits
    %lims = [mean Ly-hy];
    %cp = normcdf(lims, mean, standard_deviation);
    %We can calculate the probability between these two points
    %Prob = cp(1) - cp(2);
        
    end

    function [Vs_ij,Grid_S,Nparticles]=test_interaction_part(Grid_S,Nparticles,parameter_seed,distribution)
        
        % Two particles in the same row
        if (distribution==3)
        i=parameter_seed(1);    
        j=parameter_seed(2);
        
        while (Nparticles(2)<2)
        if (Grid_S(i,j)==1)&&(Grid_S(i+5,j)==1)
        Grid_S(i,j)=2;
        Nparticles(2)=Nparticles(2)+1;
        Nparticles(1)=Nparticles(1)-1;
        Vs_ij(Nparticles(2),1)=i;
        Vs_ij(Nparticles(2),2)=j; 
        Grid_S(i+5,j)=2;
        Nparticles(2)=Nparticles(2)+1;
        Nparticles(1)=Nparticles(1)-1;
        Vs_ij(Nparticles(2),1)=i+5;
        Vs_ij(Nparticles(2),2)=j; 
    
        else
            if (Grid_S(i,j)==1)&&(Grid_S(i-5,j)==1)    
            Grid_S(i,j)=2;
            Nparticles(2)=Nparticles(2)+1;
            Nparticles(1)=Nparticles(1)-1;
            Vs_ij(Nparticles(2),1)=i;
            Vs_ij(Nparticles(2),2)=j; 
            Grid_S(i-5,j)=2;
            Nparticles(2)=Nparticles(2)+1;
            Nparticles(1)=Nparticles(1)-1;
            Vs_ij(Nparticles(2),1)=i-5;
            Vs_ij(Nparticles(2),2)=j; 
        
            else 
            j=j+1;
            end    
        end
        end
        end
        
        % Two particles in the same column
        if (distribution==4)
        i=parameter_seed(1);    
        j=parameter_seed(2);
        
        while (Nparticles(2)<2)
        if (Grid_S(i,j)==1)&&(Grid_S(i,j+2)==1)
        Grid_S(i,j)=2;
        Nparticles(2)=Nparticles(2)+1;
        Nparticles(1)=Nparticles(1)-1;
        Vs_ij(Nparticles(2),1)=i;
        Vs_ij(Nparticles(2),2)=j; 
        Grid_S(i,j+2)=2;
        Nparticles(2)=Nparticles(2)+1;
        Nparticles(1)=Nparticles(1)-1;
        Vs_ij(Nparticles(2),1)=i;
        Vs_ij(Nparticles(2),2)=j+2; 
    
        else
            if (Grid_S(i,j)==1)&&(Grid_S(i,j-2)==1)    
            Grid_S(i,j)=2;
            Nparticles(2)=Nparticles(2)+1;
            Nparticles(1)=Nparticles(1)-1;
            Vs_ij(Nparticles(2),1)=i;
            Vs_ij(Nparticles(2),2)=j; 
            Grid_S(i,j-2)=2;
            Nparticles(2)=Nparticles(2)+1;
            Nparticles(1)=Nparticles(1)-1;
            Vs_ij(Nparticles(2),1)=i;
            Vs_ij(Nparticles(2),2)=j-2; 
        
            else 
            j=j+1;
            end    
        end
        end
        end
        
        % Particles in the borders (top, left and right)
        if (distribution==5)
        
        % Particle in top border    
        i=parameter_seed(1);
        j=size(Grid_S,2);
        
        while Nparticles(2)<1
        if (Grid_S(i,j)==1)
        Grid_S(i,j)=2;
        Nparticles(2)=Nparticles(2)+1;
        Nparticles(1)=Nparticles(1)-1;
        Vs_ij(Nparticles(2),1)=i;
        Vs_ij(Nparticles(2),2)=j;
        else
        i=i+1;
        end
        end
        
        % Particle in left border    
        i=1;
        j=parameter_seed(2);
        
        while Nparticles(2)<2
        if (Grid_S(i,j)==1)
        Grid_S(i,j)=2;
        Nparticles(2)=Nparticles(2)+1;
        Nparticles(1)=Nparticles(1)-1;
        Vs_ij(Nparticles(2),1)=i;
        Vs_ij(Nparticles(2),2)=j; 
        else
        j=j+1;
        end
        end
        
        % Particle in right border    
        i=size(Grid_S,1);
        j=parameter_seed(2);
        while Nparticles(2)<3
        if (Grid_S(i,j)==1)
        Grid_S(i,j)=2;
        Nparticles(2)=Nparticles(2)+1;
        Nparticles(1)=Nparticles(1)-1;
        Vs_ij(Nparticles(2),1)=i;
        Vs_ij(Nparticles(2),2)=j;
        else
        j=j+1;
        end
        end
        
        
        end
    end

    function [Vs_ij,Grid_S,Nparticles]=triangle(Grid_S,Nparticles,parameter_seed)
    dose_height=parameter_seed(1);
    triangle_base=parameter_seed(2);
    orientation_triangle=parameter_seed(3);

    nj=size(Grid_S,2);
    start_irradiated_region=round(nj/2-triangle_base);
    finish_irradiated_region=round(nj/2+triangle_base);

    % Linear fit for the dose --> Triangle shape
    % Orientation is chosen --> 1 Pointing to upward
    %                       --> -1 Pointing to downward
    if orientation_triangle==1
    slope=dose_height/(2*triangle_base);
    b=dose_height;
    else
    slope=-dose_height/(2*triangle_base);
    b=0;
    end


    %% dose_width determines how many rows will have defects starting from the last one
    % backward count
    for j=start_irradiated_region:(finish_irradiated_region)
    dose=slope*(start_irradiated_region-j)+b;

    for i=1:size(Grid_S,1)    
    if (Grid_S(i,j)==1)
        %% dose determines the probability of finding a defect    
        if (rand<dose)
        Grid_S(i,j)=2;
        Nparticles(2)=Nparticles(2)+1;
        Nparticles(1)=Nparticles(1)-1;
        Vs_ij(Nparticles(2),1)=i;
        Vs_ij(Nparticles(2),2)=j;
        end
    else continue;

    end
    end
    end

    end
end