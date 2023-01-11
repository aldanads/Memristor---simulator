%-----------------------------------------------------------------%
% Simulator: Lateral memtransistors from monolayer single-crystal %
% molybdenum disulfide (MoS2)                                     %
% Starting Date: 2021/11/16                                       %
% Samuel Aldana Delgado                                           %
%-----------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Purpose: This function make a hexagonal grid with a 2D matrix
%%
%%------------------------------------------------------------------------
%% -------------------- Parameters ---------------------------------------
%% d(1) = dx --> dimension of the grid in x axis
%% d(2) = dy --> dimension of the grid in y axis
%% Grid_S --> 2D hexagonal grid with d(1) and d(2) dimension and h(1) and h(2) separation between grid points
%%       within the Grid_S matrix:
%%                                - 0 forbidden position, not part of the crystal
%%                                - 1 Sulfur atom, S
%% -----------------------------------------------------------------------


function [Grid_S,Nparticles]=hex_grid(d)

Grid_S=zeros(d(1),d(2));
%% Number or particles in the grid (1-S and 2-S vacancies)
Nparticles=zeros(1,2);
part_row_1=0;

%% This is the crystal structure with O (sulfur position) and X (forbidden positions)
%%                   X O O X X O O X
%%                   O X X O O X X O 
%%                   X O O X X O O X
%%                   O X X O O X X O
%%                   X O O X X O O X
%%                   1 2 3 4 5 6 7 8 ...

%% First position for Sulfur in first row: 2 and 3
i=2;
%% First position for Sulfur in second row: 1 and 2
j=1;
Grid_S(j,2)=1;
part_row_2=1;
j=j+3;

%% Until we finish the first row
while (i+1<=d(1))    
%%  Row=1
    Grid_S(i,1)=1;
    Grid_S(i+1,1)=1;
    part_row_1=part_row_1+2;
    
%%  Row=2
    Grid_S(j,2)=1;
    if(j+1<d(1))
    Grid_S(j+1,2)=1;
    part_row_2=part_row_2+2;
    else
    part_row_2=part_row_2+1;    
    end

%% Jumps of 4
i=i+4;
j=j+4;
end

%% Particles in the first two rows
Nparticles(1)=Nparticles(1)+part_row_1+part_row_2;


%% Next rows are like the first one or the second one
%% Going along all the rows
i=3;
while (i+1<=d(2))
    Grid_S(:,i)=Grid_S(:,1);
    Grid_S(:,i+1)=Grid_S(:,2);
    Nparticles(1)=Nparticles(1)+part_row_1+part_row_2;

    i=i+2;
end

if mod(d(2),2)~=0
   Grid_S(:,i)=Grid_S(:,1);
   Nparticles(1)=Nparticles(1)+part_row_1;

end
end







