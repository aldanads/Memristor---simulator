
function [resistance_step,v_resistance,v_total_res]=resistance_calculation(density_Vs,parameters,res_fit_parameters,resistance_limits,j_total,v_resistance,v_total_res)
        
    %density_Vs=density_Vs*1E-7;

    % Separation of sections depending on the defect density
    section_1_lim=parameters(12);
    section_2_lim=parameters(13);

    % resistivity configuration
    resistance_configuration=parameters(15);



    % 3) Semiconductor-Insulator-Metal --> /\
    if resistance_configuration==3
    [resistance_step]=resistance_configuration_3(density_Vs,section_1_lim,res_fit_parameters,resistance_limits);
    end
        
    % 4) Semiconductor-Insulator-Metal --> \/
    if resistance_configuration==4
    [resistance_step,v_resistance,v_total_res]=resistance_configuration_4(density_Vs,section_1_lim,res_fit_parameters,resistance_limits,j_total,v_resistance,v_total_res,parameters);
    end

    % 6) Semiconductor-Insulator-Metal (uncompleted) --> _/\
    if resistance_configuration==6
    [resistance_step]=resistance_configuration_6(density_Vs,section_1_lim,section_2_lim,res_fit_parameters,resistance_limits);
    end

   

end

    % 3) Semiconductor-Insulator-Metal --> /\
    function [resistance_step]=resistance_configuration_3(density_Vs,section_1_lim,res_fit_parameters,resistance_limits)
    % We calculate the resistivity for every row (ohmios)
    resistance_step=zeros(size(density_Vs));

    slope_1=res_fit_parameters(1,1);
    b_1=res_fit_parameters(1,2);
    slope_2=res_fit_parameters(2,1);
    b_2=res_fit_parameters(2,2);

    for j=1:length(density_Vs)
    
    % No particles in the row    
    if density_Vs(j)==0
    %resistance(1)=pristine_state_resistance;    
    resistance_step(j)=resistance_limits(1);
    end

    % Particle concentration corresponds to the first region
    if (density_Vs(j)<=section_1_lim)&&(density_Vs(j)>0)
    resistance_step(j)=density_Vs(j)^slope_1*10^b_1;
    % Particle concentration corresponds to the second region 
    else
        % skip 0
        if (density_Vs(j)>section_1_lim)    
        resistance_step(j)=density_Vs(j)^slope_2*10^b_2;
        end
    end

    end

    end

    % 4) Semiconductor-Insulator-Metal --> \/
    function [resistance_step,v_resistance,v_total_res]=resistance_configuration_4(density_Vs,section_1_lim,res_fit_parameters,resistance_limits,j_total,v_resistance,v_total_res,parameters)
    

    square_size=parameters(17);

    slope_1=res_fit_parameters(1,1);
    b_1=res_fit_parameters(1,2);
    slope_2=res_fit_parameters(2,1);
    b_2=res_fit_parameters(2,2);

    if square_size==0
    %% THIS SHOULD BE VECTORIZED
    % We calculate the resistivity for every row (ohmios)
    resistance_step=zeros(size(density_Vs));

    for j=1:length(density_Vs)
    
    % No particles in the row    
    if density_Vs(j)==0
    %resistance(1)=pristine_state_resistance;    
    resistance_step(j)=resistance_limits(1);
    end

    % Particle concentration corresponds to the first region
    if (density_Vs(j)<=section_1_lim)&&(density_Vs(j)>0)
    resistance_step(j)=density_Vs(j)^slope_1*10^b_1;
    % Particle concentration corresponds to the second region 
    else
        % skip 0
        if (density_Vs(j)>section_1_lim)    
        resistance_step(j)=density_Vs(j)^slope_2*10^b_2;
        end
    end

    end

    % Save resistance along the simulation
    v_resistance(:,j_total)=resistance_step;
    v_total_res(j_total)=sum(v_resistance(:,j_total));
    
    else
        
        resistance_step=zeros(size(density_Vs));
        %idx=(density_Vs<=section_1_lim & density_Vs~=0);
        idx=(density_Vs<section_1_lim & density_Vs~=0);
        resistance_step(idx)=(density_Vs(idx).^slope_1)*10^b_1;
        resistance_step(~idx)=(density_Vs(~idx).^slope_2)*10^b_2;
        resistance_step(density_Vs==0)=resistance_limits(1);

        % Save resistance along the simulation
        v_resistance(:,j_total)=1./sum(1./resistance_step,1);
        v_total_res(j_total)=sum(v_resistance(:,j_total));

    end

    end

    % 6) Semiconductor-Insulator-Metal (uncompleted) --> _/\
    function [resistance_step]=resistance_configuration_6(density_Vs,section_1_lim,section_2_lim,res_fit_parameters,resistance_limits)
    
    % We calculate the resistivity for every row (ohmios)
    resistance_step=zeros(size(density_Vs));
    
    slope_1=res_fit_parameters(1,1);
    b_1=res_fit_parameters(1,2);
    slope_2=res_fit_parameters(2,1);
    b_2=res_fit_parameters(2,2);
    slope_3=res_fit_parameters(3,1);
    b_3=res_fit_parameters(3,2);

    for j=1:length(density_Vs)
    
    % We assume the first point in [Fox2015] is the pristine state and
    % there is no defects (before irradiation). This is not exactly true,
    % as usually defects are found in the material depending on the
    % fabrication process --> No particle = 3E6 ohmios
    if density_Vs(j)==0
    %resistance(1)=pristine_state_resistance;    
    resistance_step(j)=resistance_limits(1);
    end

    if (density_Vs(j)<=section_1_lim)&&(density_Vs(j)>0)
    resistance_step(j)=density_Vs(j)^slope_1*10^b_1;
    end

    if (density_Vs(j)<section_2_lim)&&(density_Vs(j)>section_1_lim)
    resistance_step(j)=density_Vs(j)^slope_2*10^b_2;
    end

    if (density_Vs(j)>=section_2_lim)
    resistance_step(j)=density_Vs(j)^slope_3*10^b_3;
    end

    end
    end