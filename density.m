%-----------------------------------------------------------------%
% Simulator: Lateral memtransistors from monolayer single-crystal %
% molybdenum disulfide (MoS2)                                     %
% Starting Date: 2021/11/16                                       %
% Samuel Aldana Delgado                                           %
%-----------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Purpose: We calculate the linear density for each row j


function [density_Vs,v_density_Vs]=density(Grid_S,Vs_ij,phy_const,parameters,j,v_density_Vs)


square_size=parameters(17);

if square_size==0

    [density_Vs]=linear_density(Vs_ij,Grid_S,phy_const);
    v_density_Vs(:,j)=density_Vs;


else

    %[density_Vs]=square_density(Vs_ij,Grid_S,square_size,phy_const);
    [density_Vs]=continuous_density(Vs_ij,Grid_S,square_size,phy_const);
    v_density_Vs(:,j)=mean(density_Vs,1);


end
end

function [density_Vs]=linear_density(Vs_ij,Grid_S,phy_const)

% In nanometer
h1=phy_const(4)*1E9;
%% We are looking for the density for the row j in micrometers
density_Vs=zeros(size(Grid_S,2),1);

 for i=1:size(Vs_ij,1)
     % If there is a vacancy at row j, we add +1
     density_Vs(Vs_ij(i,2))=density_Vs(Vs_ij(i,2))+1;
 end
    % Linear density for each row in nm --> To turn it into micrometer: *1000
    % Linear density for each row in nm --> To turn it into centimeter: 1E7
    density_Vs=1E7*density_Vs/(size(Grid_S,1)*h1);

end
function [density_Vs]=square_density(Vs_ij,Grid_S,square_size,phy_const)
    
    hx=phy_const(4)*1E9;
    hy=phy_const(5)*1E9;

    ni=size(Grid_S,1);
    nj=size(Grid_S,2);

    nx_squares=ni/square_size;
    ny_squares=nj/square_size;
    density_Vs=zeros(nx_squares,ny_squares);

    % We start searching from the min row to the max_row --> skip
    % many empty rows
    min_row=min(Vs_ij(:,2));
    max_row=max(Vs_ij(:,2));
    j=1;

    while (j<=max_row)
    % The number of particles within the following square:
    if j>=min_row
    for i=1:square_size:(size(Grid_S,1))
    % We should undo the index Grid_S index into the square index.
    density_Vs(1+(i-1)/square_size,1+(j-1)/square_size)=sum((Vs_ij(:,1)>=i) & (Vs_ij(:,1)<i+square_size)& (Vs_ij(:,2)>=j) & (Vs_ij(:,2)<(j+square_size)));
    end
    end
    j=j+square_size;

    end

    density_Vs=density_Vs/(hx*hy*square_size^2);
end

function [density_Vs]=continuous_density(Vs_ij,Grid_S,square_size,phy_const)
    
    hx=phy_const(4)*1E9;
    hy=phy_const(5)*1E9;
    half_square=square_size/2;

    ni=size(Grid_S,1);
    nj=size(Grid_S,2);

    density_Vs=zeros(ni,nj);
    [virt_Vs_ij]=periodic_boundary(Vs_ij,half_square,ni,nj);

    for i=half_square+1:ni+half_square
        for j=half_square+1:nj+half_square
    density_Vs(i-half_square,j-half_square)=sum((virt_Vs_ij(:,1)>=(i-half_square)) & (virt_Vs_ij(:,1)<(i+half_square)) & (virt_Vs_ij(:,2)>=(j-half_square)) & (virt_Vs_ij(:,2)<(j+half_square)));
        end
    end

    density_Vs=density_Vs/(hx*hy*square_size^2);
    %% Periodic boundary conditions
    function [virt_Vs_ij]=periodic_boundary(Vs_ij,half_square,ni,nj)
        
        % We create some virtual particles outside the simulation domain
        % --> the extended simulation domain is ni+square_size and
        % nj+square_size

        virt_Vs_ij=Vs_ij+half_square;

        for k=1:half_square
        
        %% Borders    
        % Mirror particles in x=0
        mirror_particles=Vs_ij((Vs_ij(:,1)==k),1)+(half_square-2*k+1);
        mirror_particles(:,2)=Vs_ij((Vs_ij(:,1)==k),2)+half_square;
        len=size(mirror_particles,1);
        virt_Vs_ij(end+1:end+len,:)=mirror_particles;

        % Mirror particles in xmax
        mirror_particles=Vs_ij((Vs_ij(:,1)==ni+1-k),1)+half_square+2*k-1;
        mirror_particles(:,2)=Vs_ij((Vs_ij(:,1)==ni+1-k),2)+half_square;
        len=size(mirror_particles,1);
        virt_Vs_ij(end+1:end+len,:)=mirror_particles;

        % Mirror particles in y=0
        mirror_particles=Vs_ij((Vs_ij(:,2)==k),1)+half_square;
        mirror_particles(:,2)=Vs_ij((Vs_ij(:,2)==k),2)+(half_square-2*k+1);
        len=size(mirror_particles,1);
        virt_Vs_ij(end+1:end+len,:)=mirror_particles;

        % Mirror particles in xmax
        mirror_particles=Vs_ij((Vs_ij(:,2)==nj+1-k),1)+half_square;
        mirror_particles(:,2)=Vs_ij((Vs_ij(:,2)==nj+1-k),2)+half_square+2*k-1;
        len=size(mirror_particles,1);
        virt_Vs_ij(end+1:end+len,:)=mirror_particles;

        end

        %% Corners
        % Mirror particles in x=0 and y=0
        mirror_particles=virt_Vs_ij((virt_Vs_ij(:,1)<=half_square & virt_Vs_ij(:,2)>half_square & virt_Vs_ij(:,2)<=2*half_square),1);
        mirror_particles(:,2)=virt_Vs_ij((virt_Vs_ij(:,1)<=half_square & virt_Vs_ij(:,2)>half_square & virt_Vs_ij(:,2)<=2*half_square),2)-half_square;
        len=size(mirror_particles,1);
        virt_Vs_ij(end+1:end+len,:)=mirror_particles;

        % Mirror particles in xmax and y=0
        mirror_particles=virt_Vs_ij((virt_Vs_ij(:,1)>=(ni+half_square) & virt_Vs_ij(:,2)>half_square & virt_Vs_ij(:,2)<=2*half_square),1);
        mirror_particles(:,2)=virt_Vs_ij((virt_Vs_ij(:,1)>=(ni+half_square) & virt_Vs_ij(:,2)>half_square & virt_Vs_ij(:,2)<=2*half_square),2)-half_square;
        len=size(mirror_particles,1);
        virt_Vs_ij(end+1:end+len,:)=mirror_particles;

        % Mirror particles in x=0 and ymax
        mirror_particles=virt_Vs_ij((virt_Vs_ij(:,1)<=half_square & virt_Vs_ij(:,2)>nj & virt_Vs_ij(:,2)<=nj+half_square),1);
        mirror_particles(:,2)=virt_Vs_ij((virt_Vs_ij(:,1)<=half_square & virt_Vs_ij(:,2)>nj & virt_Vs_ij(:,2)<=nj+half_square),2)+half_square;
        len=size(mirror_particles,1);
        virt_Vs_ij(end+1:end+len,:)=mirror_particles;

        % Mirror particles in xmax and ymax
        mirror_particles=virt_Vs_ij((virt_Vs_ij(:,1)>=(ni+half_square) & virt_Vs_ij(:,2)>nj & virt_Vs_ij(:,2)<=nj+half_square),1);
        mirror_particles(:,2)=virt_Vs_ij((virt_Vs_ij(:,1)>=(ni+half_square) & virt_Vs_ij(:,2)>nj & virt_Vs_ij(:,2)<=nj+half_square),2)+half_square;
        len=size(mirror_particles,1);
        virt_Vs_ij(end+1:end+len,:)=mirror_particles;
    
    end

end