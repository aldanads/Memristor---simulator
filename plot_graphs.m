%-----------------------------------------------------------------%
% Simulator: Lateral memtransistors from monolayer single-crystal %
% molybdenum disulfide (MoS2)                                     %
% Starting Date: 2021/11/16                                       %
% Samuel Aldana Delgado                                           %
%-----------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function []=plot_graphs(Grid_S,Vs_ij,voltage,direction_vac,direction_phy_pl,v_time,Nparticles,phy_const,j_total,parameters,v_density_Vs,ey,v_resistance,screening_j,den_res,resistance_limits,density_Vs,resistance_step,v_current,v_total_res,current_mapping)

h1=phy_const(4)*1E9;
h2=phy_const(5)*1E9;
square_size=parameters(17);

vac_mov=figure('visible','off');

ax1=subplot(1,1,1);
plot_vacancies(Grid_S,Vs_ij,voltage(j_total),Nparticles,phy_const,ax1);
% Direction and name of the files
W= strcat(direction_vac,'\',num2str(j_total),'_Graph','_',num2str(voltage(j_total)),'V','_',num2str(v_time(j_total)),'s','.jpg');
saveas(vac_mov,W)

phy_variables=figure('visible','off');
% Graph for superficial density
ax2=subplot(2,3,1);

% Channel resistivity and screening factor plot
ax3=subplot(2,3,2);
ax5=subplot(2,3,4);
% I-V plot
ax6=subplot(2,3,5);
% Density - Resistivity
ax7=subplot(2,3,6);
if square_size==0
plot_density_resistivity(v_density_Vs(:,j_total),phy_const,parameters,v_density_Vs,j_total,voltage,ax2,ax3,ax5,v_resistance(:,j_total),v_resistance,v_time,screening_j,ax6,den_res,ax7,resistance_limits,v_current,v_total_res)
else
plot_density_resistivity_squares(ax2,ax3,ax5,ax6,ax7,phy_const,Grid_S,screening_j,v_resistance,voltage,v_time,v_density_Vs,resistance_limits,den_res,density_Vs,square_size,resistance_step,j_total,v_current,v_total_res,current_mapping)
end
% Electric field --> Component y
ax4=subplot(2,3,3);
x_coord=h1:h1:size(Grid_S,1)*h1;
y_coord=h2:h2:size(Grid_S,2)*h2;
imagesc(ax4,y_coord,x_coord,ey);
colorbar
c = colorbar;
c.Label.String = 'Electric field (V/m)';
ylabel(ax4, 'X axis (nm)', 'FontSize',6)
xlabel(ax4,'Y axis (nm)','fontsize',6)
xlim([h2 size(Grid_S,2)*h2]);
ylim([h1 size(Grid_S,1)*h1]);





% Direction and name of the files
W= strcat(direction_phy_pl,'\',num2str(j_total),'_Phy_v','_',num2str(voltage(j_total)),'V','_',num2str(v_time(j_total)),'s','.jpg');
saveas(phy_variables,W)

end


    function []=plot_vacancies(Grid_S,Vs_ij,V,Nparticles,phy_const,ax1)

    % In nanometers
    h1=phy_const(4)*1E9;
    h2=phy_const(5)*1E9;
        
    % Red
    color='r';
    % Blue
    color2='b';

    %% Sulfur atoms coordinates
    [S_ij]=S_coordinates(Grid_S,Nparticles);

    % Stand the graph and plot both S and Vs
    hold on
    m=max(size(Grid_S));
    %scatter(Vs_ij(:,1)*h1,Vs_ij(:,2)*h2,2500/m,'MarkerEdgeColor','k','MarkerFaceColor',color, 'LineWidth', 0.5)
    scatter(ax1,Vs_ij(:,1)*h1,Vs_ij(:,2)*h2,1500/m,'MarkerEdgeColor',color, 'LineWidth', 0.5)
    scatter(ax1,S_ij(:,1)*h1,S_ij(:,2)*h2,500/m,'MarkerEdgeColor','k','MarkerFaceColor',color2, 'LineWidth', 0.5)
    
    %grid on
    %grid minor

    title([' V(V)= ',num2str(V)],'fontsize',12)
    xlim([0 size(Grid_S,1)*h1]);
    ylim([0 size(Grid_S,2)*h2]);
    
    % Proportion between axis x, y and z
    %pbaspect([1 3 1])

    xlabel('X axis (nm)','fontsize',12)
    ylabel('Y axis (nm)','fontsize',12)
    
    end

    function [S_ij]=S_coordinates(Grid_S,Nparticles)
        ni=size(Grid_S,1);
        nj=size(Grid_S,2);
        S_ij=zeros(Nparticles(1),2);
        count=1;
        for i=1:ni
            for j=1:nj
                if (Grid_S(i,j)==1)
                    S_ij(count,1)=i;
                    S_ij(count,2)=j;
                    count=count+1;
                end
            end
        end
    end

    function []=plot_density_resistivity(density_Vs,phy_const,parameters,v_density_Vs,j_total,voltage,ax2,ax3,ax5,resistivity,v_resistance,v_time,screening_j,ax6,den_res,ax7,resistance_limits,v_current,v_total_res)
        
        % In nanometers
        h2=phy_const(5)*1E9;
        
        % The number of columns with select to calculate the superficial
        % density
        n_rows=parameters(9);
        %Vmax=parameters(7);
        %Vmin=parameters(8);
        [density,y]=superficial_density(h2,n_rows,density_Vs);
        

        %y=h2:h2:size(Grid_S,2)*h2;
        plot(ax2,density,y)
        title(ax2,' Density of vacancies','fontsize',12)
        xlabel(ax2,'V_{S}/cm^2','fontsize',6)
        ylabel(ax2,'Y axis (nm)','fontsize',12)
        xlim(ax2,[0 6E14]);
        ylim(ax2,[0 size(Grid_S,2)*h2]);

        % Screening plot
        y=h2:h2:size(Grid_S,2)*h2;
        plot(ax3,screening_j,y)
        title(ax3, 'Screening')
        xlabel(ax3, 'Screening', 'FontSize',6)
        ylabel(ax3,'Y axis (nm)','fontsize',12)
        ylim(ax3,[0 size(Grid_S,2)*h2]);
        xlim(ax3,[0 1]);

        % Channel resistivity
        semilogy(ax5,v_time(1:j_total,1),v_total_res(1:j_total))
        title(ax5, 'Channel resistivity')
        xlabel(ax5,'Time (s)', 'FontSize',6)
        ylabel(ax5,'Resistivity (Ohmios)','FontSize',12)
        %ylim(ax5,[5E8 7E9]);

        % Current vs voltage
        semilogy(ax6,voltage(1:j_total,1),abs(v_current(1:j_total)))
        title(ax6, 'I-V')
        xlabel(ax6,'Voltage (V)', 'FontSize',6)
        ylabel(ax6,'Current (A)','FontSize',12)
        
        % Set properties of x label
        %xlabh=get(gca,'xlabel');
        %xlabh.Position(1)=-2000;
        %xlabh.Position(2)=-1000;
        %set(xlabh,'Position',get(xlabh,'Position') - [0 -6 0])
        
        %[peak,peak_position_y]=density_peak(v_density_Vs,j_total);
        %hold on
        %plot(ax3,voltage(1:j_total),peak)
        %title(ax3, 'Max density')
        %xlabel(ax3, 'Voltage (V)', 'FontSize',6)
        %ylabel(ax3,'Max vacancy density (V_{S}/cm^2)','fontsize',12)

        % Resistivity plot
        %semilogx(ax3,resistivity,y)
        %title(ax3, 'Resistivity')
        %xlabel(ax3, 'Resistivity (Ohmios)', 'FontSize',6)
        %ylabel(ax3,'Y axis (nm)','fontsize',12)
        %ylim(ax3,[0 size(Grid_S,2)*h2]);
        %xlim(ax3,[2.5E5 7E9]);

        %plot(ax5,voltage(1:j_total),peak_position_y*h2)
        %title(ax5, 'Max density position')
        %xlabel(ax5,'Voltage (V)', 'FontSize',6)
        %ylabel(ax5,'Max density position (nm)','FontSize',12)

        %[hAx]=plotyy(ax3,voltage(1:j_total),peak,voltage(1:j_total),peak_position_y*h2);
        %xlabel(ax3,'Voltage (V)','fontsize',10)
        %xlabh2=get(gca,'xlabel');
        %set(xlabh2,'Position',get(xlabh2,'Position') - [0 -6 0])
        %set(hAx,{'ycolor'},{'r';'b'});
        %ylabel(hAx(1),'Max vacancy density (V_{S}/cm^2)','fontsize',10);
        %ylabel(hAx(2),'Max density position (nm)','fontsize',10);
        %xlim(hAx,[Vmin Vmax]);

        [max_min_res(1), coord(1)]=max(v_resistance(:,j_total));

        % We have to avoid the minimum resistance where there is no defect
        [max_min_res(2), coord(2)]=min(v_resistance(v_resistance(:,j_total)>resistance_limits(1),j_total));
        
        % Skip 0 density points
        coord_den=find(v_density_Vs(:,j_total)>0);
        coord(2)=coord_den(coord(2));

        %cord(1)=find(v_resistivity(:,j_total)==max_min_res(1),1);
        %cord(2)=find(v_resistivity(:,j_total)==max_min_res(2),1);

        semilogy(ax7,den_res(1,:)*1E7,den_res(2,:))
        hold on;
        grid on;
        scatter(ax7,v_density_Vs(coord,j_total),max_min_res)
        title(ax7,'Resistance - defect density')
        xlabel(ax7,'V_{S}/cm', 'FontSize',6)
        ylabel(ax7,'Resistance (Ohmios)','FontSize',12)
        
        
        function [density,y]=superficial_density(h2,n_rows,density_Vs)
        % The number of devisions of the simulation domain for calculate
        % the superficial density in ribbons
        divisions=round(size(Grid_S,2)/n_rows);
        density=zeros(divisions,1);
        y=zeros(divisions,1);
        
        % The superficial density ---> area: length is x-axis and width n_rows
        j=1;
        y(j)=n_rows*h2/2;
        while j<=divisions
        for i=1:n_rows
            if n_rows*j<=size(Grid_S,2)
            density(divisions+1-j)=density(divisions+1-j)+density_Vs(size(density_Vs,1)+1-i-(j-1)*n_rows);
            end
            y(j)=y(1)+(j-1)*n_rows*h2;
        end
        j=j+1;
        end
        % The density calculated in "density_row" is a linear density
        % To turn it into micrometer: *1000
        % To turn it into centimeter: *1E7

        density=1E7*density/(n_rows*h2);
        end

        %% AVISO NAVEGANTES --> Máximo en densidad lineal, no superficial
        function [peak,peak_position_y]=density_peak(v_density_Vs,j_total)

            [peak,peak_position_y]=max(v_density_Vs(:,j_total));

            
        end




    end

    function []=plot_density_resistivity_squares(ax2,ax3,ax5,ax6,ax7,phy_const,Grid_S,screening_j,v_resistance,voltage,v_time,v_density_Vs,resistance_limits,den_res,density_Vs,square_size,resistance_step,j_total,v_current,v_total_res,current_mapping)

    % In nanometers
    h2=phy_const(5)*1E9;
    h1=phy_const(4)*1E9;
    %y=h2:h2*square_size:size(Grid_S,2)*h2;
    %x=h1:h1*square_size:size(Grid_S,1)*h1;
    x=h1:h1:size(Grid_S,1)*h1;
    y=h2:h2:size(Grid_S,2)*h2;
    x=x';
    y=y';

    % Density plot    
    imagesc(ax2,y,x,density_Vs)
    colorbar
    c = colorbar;
    c.Label.String = 'V_{S}/nm^2';
    %surf(ax2,y,x,density_Vs);
    title(ax2,' Density of vacancies','fontsize',6)
    %zlabel(ax2,'V_{S}/nm^2','fontsize',6)
    ylabel(ax2, 'X axis (nm)', 'FontSize',6)
    xlabel(ax2,'Y axis (nm)','fontsize',6)
    ylim(ax2,[0 size(Grid_S,1)*h1]);
    xlim(ax2,[0 size(Grid_S,2)*h2]);
    %zlim(ax2,[0 6E14])

    max_resistance=resistance_limits(4);
    min_resistance=resistance_limits(5);

    % Screening plot / Resistivity / Current map
    %surf(ax3,y,x,screening_j)
    imagesc(ax3,y,x,resistance_step)
    colorbar
    c = colorbar;
    c.Label.String = 'Resistivity (Ohms)';
    %bar3(ax3,abs(current_mapping))
    %bar3(ax3,resistance_step)
    %s=pcolor(ax3,y,x,resistivity_step);
    %s.FaceColor = 'interp';
    %colorbar(ax3,'northoutside')
    %caxis([1E4 1E8])
    %Label.String = 'Resistivity (Ohmios)';
    %title(ax3, 'Screening')
    %zlabel(ax3, 'Screening', 'FontSize',6)
    ylabel(ax3, 'X axis (nm)', 'FontSize',6)
    xlabel(ax3,'Y axis (nm)','fontsize',6)
    xlim(ax3,[0 size(Grid_S,2)*h2]);
    ylim(ax3,[0 size(Grid_S,1)*h1]);
    %zlim(ax3,[0 2E-5]);

    % Channel resistivity
    semilogy(ax5,v_time(1:j_total,1),v_total_res(1:j_total))
    title(ax5, 'Channel', 'FontSize',6)
    xlabel(ax5,'Time (s)', 'FontSize',6)
    ylabel(ax5,'Resistivity (Ohmios)','FontSize',12)

    % Current vs voltage
    semilogy(ax6,voltage(1:j_total,1),abs(v_current(1:j_total)))
    title(ax6, 'I-V')
    xlabel(ax6,'Voltage (V)', 'FontSize',6)
    ylabel(ax6,'Current (A)','FontSize',12)


    % Resistance vs defect density
    %[max_min_res(1), coord(1)]=max(max(resistivity_j));
    [max_res, idx_max]=max(resistance_step);
    [max_min_res(1), idx]=max(max_res);
    coord_max(1)=idx_max(idx);
    coord_max(2)=idx;
    % We have to avoid the minimum resistance where there is no defect
    idx=find(resistance_step>resistance_limits(1));
    [max_min_res(2), idx_min]=min(resistance_step(idx));
    density(1)=density_Vs(coord_max(1),coord_max(2));
    density(2)=density_Vs(idx(idx_min));

    %[max_min_res(2), coord(2)]=min(resistivity_j(resistivity_j>resistance(1)));
     
    % Skip 0 density points
    %coord_den=find(v_density_Vs(:,j_total)>0);
    %coord(2)=coord_den(coord(2));

    semilogy(ax7,den_res(1,:),den_res(2,:))
    hold on;
    grid on;
    scatter(ax7,density,max_min_res)
    %scatter(ax7,density_Vs(coord_max(1),coord_max(2)),max_min_res(1),density_Vs(idx(idx_min)),max_min_res(2),40,'kx')
    title(ax7,'Single square')
    xlabel(ax7,'V_{S}/cm', 'FontSize',6)
    ylabel(ax7,'Resistance (Ohmios)','FontSize',12)

    end