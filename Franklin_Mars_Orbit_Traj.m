%% Assessed Exercise 1 
clear all; close all; clc;
%% Script to compute total DeltaV for all the possible of the transfers of following grids of departure and arrival time

%% Initialise Variables

%Spacecraft features:
Orbiter_Isp = 310;

%Global
mu_S = getAstroConstants('Sun','mu');

mu_E = getAstroConstants('Earth','mu');
rad_E = getAstroConstants('Earth','Radius');
g_E = 9.81;                                                                 % ms^-2

mu_M = getAstroConstants('Mars','mu');
rad_M = getAstroConstants('Mars','Radius');



%%%%%%%%% Chosen BaselineDeparture Date Option
%Refined arrival and departure dates                                       for efficiency of running program
Departure_date = [2028 12 14 00 00 00];                                     % tinitial
Arrival_date = [2029 07 31 00 00 00];
%%%%%%%%%




Departure_date_mjd2 = date2mjd2000(Departure_date);                         % convert tinitial to date format 
Arrival_date_mjd2 = date2mjd2000(Arrival_date);                             % convert tfinal to date format
Delta_t_target = (Arrival_date_mjd2 - Departure_date_mjd2)*86400;           % Transfer Shoot Duration

% Generate Grid of conditions to analyse
Range = 30*5;                                                               % Set range to plot either side of departure date
Departure_Grid = Departure_date_mjd2-Range:10:Departure_date_mjd2+Range;     % Select dates from min-max departure dates in steps of 5
Arrival_Grid = Arrival_date_mjd2-Range:10:Arrival_date_mjd2+2*Range;           % Select dates from min-max arrival dates in steps of 5

%Target final Parking Orbit
R_CircPark_E = rad_E+200;                                                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%Check this value!
R_CircPark_M = rad_M+400;




%Falcon 9 Heavy (With Recovery) constraints
C3_vector = [6.25 60.14];
WetMass_vector = [6900 20000];
C3_vector = [0 5 10 15 20 25 30 40 60];
WetMass_vector = [5750 6500 4250 3600 3000 1900 1000];
C3_vector = [0.06 6 ];
WetMass_vector = [20000 6500 ];

%Arianne 5 Example
C3_vector =      [0     1    2.25  4     6.25  9     12.25 16   20.25];
WetMass_vector = [6191  6059 5844  5549  5180  4744  4251  3714 2155];

C3_vector =      [0      2.25   6.25    12.5    16    20.25     24.5];
WetMass_vector = [8000   7500   6900    6000    5800  5500      5250];

max_v_inf = sqrt(13);


%Initialise matrices where solutions will be stored
Delta_V_Jumbo = zeros(length(Departure_Grid));                              %Delta_V Initialised   
Final_Mass_Solutions = zeros(length(Departure_Grid));                       %Mass Solutions Initialised
ToF_Jumbo = zeros(length(Departure_Grid));


%Preallocate array sizes
unique_departure_dates = zeros(length(Departure_Grid));
unique_arrival_dates = zeros(length(Arrival_Grid));
C3_Dep_Jumbo = zeros(length(Departure_Grid));
Arrival_V_infl_Jumbo = zeros(length(Departure_Grid));
Arrival_V_infs_Jumbo = zeros(length(Departure_Grid));
C3_calc_Jumbo = zeros(length(Departure_Grid));
CalculatedFinalMass_Jumbo = zeros(length(Departure_Grid));
%% Software featrues
h = waitbar(0, 'Calculating...');                                           % Initialise Waitbar
tic;                                                                        % Start the timer

%% Nested loop to construct Matrix
for m = 1:length(Departure_Grid)
    for n = 1:length(Arrival_Grid)

        [r_E, v_E] = EphSS_car(3, Departure_Grid(m));                       % Generate Earth initial Position and Velocity 
        [r_M, v_M] = EphSS_car(4, Arrival_Grid(n));                         % Generate Mars initial Position and Velocity

        ToF = (Arrival_Grid(n)-Departure_Grid(m))*86400;
        % disp(['Current step in n: ',n,' in nested loop'])


        %Calculate cost to arrive at Mars parking orbit on SHORT transfer
        [rdot1_Lamberts, rdot2_Lamberts] = LambertArc_ATD2024_PARALLEL(r_E,r_M,ToF,+1, mu_S);
            % Delta_V1s = norm(rdot1_Lambertl - v_E);                             % Earth departure v+_inf long - use for direct orbit insertion
        Delta_V1s_inf = norm(v_E - rdot1_Lamberts);                             % Earth departure v+_inf           is this vinfinity?                                 
        V_parksE = sqrt(2*(mu_E/R_CircPark_E + norm(Delta_V1s_inf)^2/2));    % DeltaV cost for parking orbit
        V_circsE = sqrt(mu_E/R_CircPark_E);                                  %
        Delta_V1s = V_parksE - V_circsE;

        
        Delta_V2s_inf = norm(rdot2_Lamberts - v_M);                         % Mars Arrival v-_inf
        V_parks = sqrt(2*(mu_M/R_CircPark_M + norm(Delta_V2s_inf)^2/2));    % DeltaV cost for parking orbit
        V_circs = sqrt(mu_M/R_CircPark_M);                                  %
        Delta_V2s = V_parks - V_circs;

        Delta_V_Short = norm(Delta_V2s + Delta_V1s);   


        %Calculate cost to arrive at Mars parking orbit on LONG transfer
        [rdot1_Lambertl, rdot2_Lambertl] = LambertArc_ATD2024_PARALLEL(r_E,r_M,ToF,-1, mu_S);
            % Delta_V1l = norm(rdot1_Lambertl - v_E);                             % Earth departure v+_inf long - use for direct orbit insertion
        Delta_V1l_inf = norm(v_E - rdot1_Lambertl);                             % Earth departure v+_inf           is this vinfinity?                                 
        V_parksEl = sqrt(2*(mu_E/R_CircPark_E + norm(Delta_V1l_inf)^2/2));    % DeltaV cost for parking orbit
        V_circsEl = sqrt(mu_E/R_CircPark_E);                                  %
        Delta_V1l = V_parksEl - V_circsEl;


        Delta_V2l_inf = norm(rdot2_Lambertl - v_M);                         % Mars Arrival v-_inf long
        V_parkl = sqrt(2*(mu_M/R_CircPark_M + norm(Delta_V2l_inf)^2/2));    % DeltaV cost for parking orbit
        V_circl = sqrt(mu_M/R_CircPark_M);
        Delta_V2l = V_parkl - V_circl;
        
        Delta_V_Long = norm(Delta_V2l + Delta_V1l);                         % Total DeltaV for long transfer


        
        Delta_V_Jumbo(m,n) = min([Delta_V_Short Delta_V_Long]);             % Create matrix of DeltaV Values:           m = rows, n = columns     
        ToF_Jumbo(m,n) = (Arrival_Grid(n) - Departure_Grid(m));             % Create matrix of ToF in Days
        Arrival_V_inf_Jumbo(m,n) = min([Delta_V2s_inf Delta_V2l_inf]);      % Create matrix of Arrival V-_inf
        C3_calc_Jumbo(m,n) = min([Delta_V1l Delta_V1s]);                    % Create matrix of Departure C3
        Delta_V2_Capture(m,n) = min([Delta_V2s Delta_V2l]);                 % Create matrix of Capture Delta V
        Delta_V1_Capture(m,n) = min([Delta_V1s Delta_V1l]);                 % Create matrix of Capture Delta V
        
        %Identify mass related to long and short Mars transfers
        CalculatedFinalMass_L = interp1(C3_vector, WetMass_vector, Delta_V1l^2, 'linear');
        CalculatedFinalMass_S = interp1(C3_vector, WetMass_vector, Delta_V1s^2, 'linear');
        Worst_Case_WetMass = min([CalculatedFinalMass_L CalculatedFinalMass_S]);
        Delta_V_Capture_calc = min([Delta_V2s Delta_V2l]);                  %Calculate the DeltaV of chosen method

        Final_Payload_Mass_calc = Worst_Case_WetMass*exp(-Delta_V_Capture_calc/Orbiter_Isp*g_E);    % Calculate final Payload Mass

        CalculatedFinalMass_Jumbo(m,n) = Final_Payload_Mass_calc;

        % choose smaller DeltaV for that transfer,
        if Delta_V_Short <= Delta_V_Long                                    
            rdot1_Lambert = rdot1_Lamberts;                                 % If short has less Delta_V, store these values
            rdot2_Lambert = rdot2_Lamberts;
            rdot1_Lambert_norm = norm(rdot1_Lambert);
            rdot2_Lambert_norm = norm(rdot2_Lambert);
        else
            rdot1_Lambert = rdot1_Lambertl;                                 % If long has less Delta_V, store these values
            rdot2_Lambert = rdot2_Lambertl;
            rdot1_Lambert_norm = norm(rdot1_Lambert);
            rdot2_Lambert_norm = norm(rdot2_Lambert);
        end
        
        if rdot1_Lambert_norm <= max_v_inf                                  % Select relevant excess velocity v_inf
            C3_Ex4 = rdot1_Lambert_norm^2;
            mass_Wet = interp1(C3_vector,WetMass_vector,C3_Ex4,'linear', 'extrap');             % Find required initial wet mass m0       
        else 
            rdot2_Lambert_norm = max_v_inf;
            C3_Ex4 = rdot2_Lambert_norm^2;
            mass_Wet0 = interp1(C3_vector, WetMass_vector, C3_Ex4,'linear', 'extrap');          % Find required initial wet mass m0

            v_exhaust = 310*g_E;                                              %g_E = 9.81
            mass_Wet = mass_Wet0 * exp(-(Delta_V_Jumbo(m,n)*1000)/v_exhaust);
        end
        
        % rdot
        C3_Dep_Jumbo(m,n) = C3_Ex4/1e6;                                     % Convert from m2 to km2
        Mass_final = mass_Wet;
        Final_Mass_Solutions(m,n) = Mass_final;                             % Store mass Solution
        
        % Update waitbar
        progress = (m - 1) * length(Arrival_Grid) + n;
        total_iterations = length(Departure_Grid) * length(Arrival_Grid);
        waitbar(progress / total_iterations, h, sprintf('Calculating... %d%%', round(progress / total_iterations * 100)));
    end
end

close(h);                                                                   %Close Waitbar

% Stop the timer
elapsed_time = toc;
disp(['Elapsed time for identifying plot points: ' num2str(elapsed_time) ' seconds']);
% profile off;
% profile viewer;

%% Create and store contour plotting values
Smallest_V = min(min(Delta_V_Jumbo));
Smallest_Arrival_V_inf = min(min(Arrival_V_inf_Jumbo));
Smallest_Departure_C3 = min(min(C3_calc_Jumbo.^2));
Smallest_CalculatedFinalMass = min(min(CalculatedFinalMass_Jumbo));

% Initialize arrays to store unique departure and arrival dates
unique_departure_dates = [];
unique_arrival_dates = [];

% Iterate through Departure_Grid
for i = 1:length(Departure_Grid)
    % Convert MJD2000 date to MATLAB date value format
    departure_date = mjd20002date(Departure_Grid(i));
    arrival_date = mjd20002date(Arrival_Grid(i));
    
    % Convert date values to strings
    departure_date_str = datestr(departure_date, 'dd-mmm-yyyy');
    arrival_date_str = datestr(arrival_date, 'dd-mmm-yyyy');
    
    % Add the date to unique_departure_dates if not already present
    if ~any(strcmp(unique_departure_dates, departure_date_str))'
        unique_departure_dates = [unique_departure_dates; departure_date_str];
    end
        % Add the date to unique_arrival_dates if not already present
    if ~any(strcmp(unique_arrival_dates, arrival_date_str))'
        unique_arrival_dates = [unique_arrival_dates; arrival_date_str];
    end
end

%% Plot Support                                                             -
% Calculate the range of visible dates in the x and y directions
x_range = [Departure_Grid(1) Departure_Grid(end)];
y_range = [Arrival_Grid(1) Arrival_Grid(end)];

% Determine the number of ticks to display based on the range
num_ticks_x = sum(Departure_Grid >= x_range(1) & Departure_Grid <= x_range(2));
num_ticks_y = sum(Arrival_Grid >= y_range(1) & Arrival_Grid <= y_range(2));

%% C3 vs Wet Mass Info
C32WetMass(C3_vector, WetMass_vector, 13.34);

%% Coloured visual plot                                                     1
% Calculate the step size for displaying ticks
%%%%%%%%%%%%%%%%% Adjust this value to choose number of dates displayed on each axis %%%%%%%%%%%%%%%%%%%%%% 
tick_density = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (adjust tick density) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_fontsize = 20;
x_tick_step = max([1, ceil(num_ticks_x / tick_density)]);                   % Adjust the divisor to control the density of ticks
y_tick_step = max([1, ceil(num_ticks_y / tick_density)]);                   % Adjust the divisor to control the density of ticks

% Convert MJD2000 departure dates to date strings
departure_date_str = cell(length(Departure_Grid), 1);
for i = 1:x_tick_step:length(Departure_Grid)                                %-1 here to remove last entry for nicer visuals
    departure_date_strarr{i} = datestr(mjd20002date(Departure_Grid(i)), 'dd-mm-yyyy');
end

% Convert MJD2000 arrival dates to date strings
arrival_date_str = cell(length(Arrival_Grid), 1);
for j = 1:y_tick_step:length(Arrival_Grid)                                  %-1 here to remove last entry for nice visuals
    arrival_date_strarr{j} = datestr(mjd20002date(Arrival_Grid(j)), 'dd-mm-yyyy');
end

figure;
hold on
    axis equal tight
    contourf(Departure_Grid, Arrival_Grid, Delta_V_Jumbo', [Smallest_V:0.2:10], 'LineColor', 'none');%[0.5 0.5 0.5]);
   
    %Axis and plot Options
    %Following three commented lines for inversed colour scheme:            % swap between these and the below two
    % cmap = flipud(colormap(bone));
    % cb = colorbar;
    % colormap(cmap)
    
    %Following two commented lines for standard colour schemes:             %swap between these and the above three
    cb = colorbar;
    cmap = colormap("bone");
    
    
    %General Avix and Plot Options
    ax = gca;
    ax.XAxis.MinorTick = 'on';
    ax.YAxis.MinorTick = 'on';
    ylabel(cb, '\DeltaV (km/s)', 'Fontsize', plot_fontsize);            % colourbar text options
    set(gca, 'XTick', Departure_Grid(1:x_tick_step:end), 'XTickLabel', departure_date_strarr(1:x_tick_step:end), 'Fontsize', plot_fontsize);
    set(gca, 'YTick', Arrival_Grid(1:y_tick_step:end), 'YTickLabel', arrival_date_strarr(1:y_tick_step:end), 'Fontsize', plot_fontsize);
    xlabel('Departure Date', 'Fontsize', plot_fontsize);
    ylabel('Arrival Date', 'Fontsize', plot_fontsize);
    
    %Title Plotting:
    Departure_date_dt = datetime(mjd20002date(Departure_date_mjd2), 'ConvertFrom', 'datenum');
    Arrival_date_dt = datetime(mjd20002date(Arrival_date_mjd2), 'ConvertFrom', 'datenum');
    Departure_date_str = datestr(Departure_date_dt, 'dd-mm-yyyy');
    Arrival_date_str = datestr(Arrival_date_dt, 'dd-mm-yyyy');
    title_str = sprintf('Gradient Porkchop Plot, Selected Baseline: %s to %s', Departure_date_str, Arrival_date_str);
    title(title_str, 'Fontsize', plot_fontsize);

hold off
 %% Technical grayscale plot with annotated values on graph                 2 DeltaV w/ tof     
% Calculate the step size for displaying ticks
%%%%%%%%%%%%%%%%% Adjust this value to choose number of dates displayed on each axis %%%%%%%%%%%%%%%%%%%%%%
tick_density = 30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (adjust tick density) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_fontsize = 20;
x_tick_step = max([1, ceil(num_ticks_x / tick_density)]);                   % Adjust the divisor to control the density of ticks
y_tick_step = max([1, ceil(num_ticks_y / tick_density)]);                   % Adjust the divisor to control the density of ticks

% Convert MJD2000 departure dates to date strings
departure_date_str = cell(length(Departure_Grid), 1);
for i = 1:x_tick_step:length(Departure_Grid)                                % -1 here to remove last entry for nice visuals
    departure_date_strarr{i} = datestr(mjd20002date(Departure_Grid(i)), 'dd-mm-yyyy');
end

% Convert MJD2000 arrival dates to date strings
arrival_date_str = cell(length(Arrival_Grid), 1);
for j = 1:y_tick_step:length(Arrival_Grid)               
    arrival_date_strarr{j} = datestr(mjd20002date(Arrival_Grid(j)), 'dd-mm-yyyy');
end

figure;
hold on
    axis equal
    grid on
    Delta_V_Display = [1:0.1:7.5];
    [P,i] = contour(Departure_Grid, Arrival_Grid, Delta_V_Jumbo', [Smallest_V:0.1:7.5], 'LevelList', Delta_V_Display, 'LineColor', 'k','LabelSpacing', 500);
    [C,h] = contour(Departure_Grid, Arrival_Grid, ToF_Jumbo', 'LineColor', 'k', 'LineWidth', 0.5, 'LevelList', [180, 225, 270, 315, 360, 405, 450, 495], 'LineColor', 'r'); % Adjust the LevelList as needed
    sp = plot(Departure_date_mjd2, Arrival_date_mjd2, 'MarkerFaceColor','m','MarkerEdgeColor','m', 'MarkerSize',10, 'Marker','*', 'LineStyle','none');
    
    
    clabel(P, i, 'Fontsize', plot_fontsize-6, 'Color','k')                      % Add labels for the Porkchop contour lines
    clabel(C, h, 'FontSize', plot_fontsize-6, 'Color','r');                     % Add labels for the ToF contour lines
    
    % Plot Surroundings
    %Axis Option
    set(gca, 'XTick', Departure_Grid(1:x_tick_step:end), 'XTickLabel', departure_date_strarr(1:x_tick_step:end), 'FontSize',plot_fontsize);
    set(gca, 'YTick', Arrival_Grid(1:y_tick_step:end), 'YTickLabel', arrival_date_strarr(1:y_tick_step:end), 'FontSize',plot_fontsize);
    xlabel('Departure Date', 'FontSize',plot_fontsize);
    ylabel('Arrival Date', 'FontSize',plot_fontsize);
    
    legend([i, h, sp], {'ΔV (km/s)','Time of ToF (days)','Baseline Date'}, 'FontSize', plot_fontsize, 'Location', 'northwest');  
                                                                   
    %Title Plotting:
    Departure_date_dt = datetime(mjd20002date(Departure_date_mjd2), 'ConvertFrom', 'datenum');
    Arrival_date_dt = datetime(mjd20002date(Arrival_date_mjd2), 'ConvertFrom', 'datenum');
    
    Departure_date_str = datestr(Departure_date_dt, 'dd-mm-yyyy');
    Arrival_date_str = datestr(Arrival_date_dt, 'dd-mm-yyyy');
    
    % title_str = sprintf('ΔV with ToF Porkchop Plot, Selected Baseline: %s to %s: Initial Parking: %i km', Departure_date_str, Arrival_date_str, R_CircPark_E-rad_E);
    title_str = sprintf('ExoMars 2022 - ΔV with ToF Porkchop Plot, Selected Baseline: %s to %s', Departure_date_str, Arrival_date_str);

    title(title_str, 'FontSize',plot_fontsize);

hold off

%% Technical                                                                2.5 Arrival DeltaV Burn     
% Calculate the step size for displaying ticks
%%%%%%%%%%%%%%%%% Adjust this value to choose number of dates displayed on each axis %%%%%%%%%%%%%%%%%%%%%%
tick_density = 35;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (adjust tick density) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_fontsize = 20;
x_tick_step = max([1, ceil(num_ticks_x / tick_density)]);                   % Adjust the divisor to control the density of ticks
y_tick_step = max([1, ceil(num_ticks_y / tick_density)]);                   % Adjust the divisor to control the density of ticks

% Convert MJD2000 departure dates to date strings
departure_date_str = cell(length(Departure_Grid), 1);
for i = 1:x_tick_step:length(Departure_Grid)                                % -1 here to remove last entry for nice visuals
    departure_date_strarr{i} = datestr(mjd20002date(Departure_Grid(i)), 'dd-mm-yyyy');
end

% Convert MJD2000 arrival dates to date strings
arrival_date_str = cell(length(Arrival_Grid), 1);
for j = 1:y_tick_step:length(Arrival_Grid)               
    arrival_date_strarr{j} = datestr(mjd20002date(Arrival_Grid(j)), 'dd-mm-yyyy');
end

figure;
hold on
    axis equal
    grid on
    Delta_V_Display = [1:0.2:4];
    [P,i] = contour(Departure_Grid, Arrival_Grid, Delta_V2_Capture', [Smallest_V:1:4], 'LevelList', Delta_V_Display, 'LineColor', [0, 0.5, 0], 'LabelSpacing', 500);
    [C,h] = contour(Departure_Grid, Arrival_Grid, ToF_Jumbo', 'LineColor', 'r', 'LineWidth', 0.5, 'LevelList', [90, 180, 270, 360, 450, 540]); % Adjust the LevelList as needed
    sp = plot(Departure_date_mjd2, Arrival_date_mjd2, 'MarkerFaceColor','m','MarkerEdgeColor','m', 'MarkerSize',10, 'Marker','*', 'LineStyle','none');
    
    
    clabel(P, i, 'Fontsize', plot_fontsize-6, 'Color',[0, 0.5, 0])                      % Add labels for the Porkchop contour lines
    clabel(C, h, 'FontSize', plot_fontsize-6, 'Color','r');                     % Add labels for the ToF contour lines
    
    % Plot Surroundings
    %Axis Option
    set(gca, 'XTick', Departure_Grid(1:x_tick_step:end), 'XTickLabel', departure_date_strarr(1:x_tick_step:end), 'FontSize',plot_fontsize);
    set(gca, 'YTick', Arrival_Grid(1:y_tick_step:end), 'YTickLabel', arrival_date_strarr(1:y_tick_step:end), 'FontSize',plot_fontsize);
    xlabel('Departure Date', 'FontSize',plot_fontsize);
    ylabel('Arrival Date', 'FontSize',plot_fontsize);
    
    legend([i, h, sp], {'Capture ΔV (km/s)','Time of ToF (days)','Baseline Date'}, 'FontSize', plot_fontsize, 'Location', 'northwest');  
                                                                   
    %Title Plotting:
    Departure_date_dt = datetime(mjd20002date(Departure_date_mjd2), 'ConvertFrom', 'datenum');
    Arrival_date_dt = datetime(mjd20002date(Arrival_date_mjd2), 'ConvertFrom', 'datenum');
    
    Departure_date_str = datestr(Departure_date_dt, 'dd-mm-yyyy');
    Arrival_date_str = datestr(Arrival_date_dt, 'dd-mm-yyyy');
    
    % title_str = sprintf('Capture ΔV with ToF Porkchop Plot, Selected Baseline: %s to %s: Initial Parking:%i km', Departure_date_str, Arrival_date_str, R_CircPark_M-rad_M);
    title_str = sprintf('ExoMars 2022 - ΔV with ToF Porkchop Plot, Selected Baseline: %s to %s', Departure_date_str, Arrival_date_str);

    title(title_str, 'FontSize',plot_fontsize);

hold off

%% Technical                                                                2.6 Depart DeltaV Burn     
% Calculate the step size for displaying ticks
%%%%%%%%%%%%%%%%% Adjust this value to choose number of dates displayed on each axis %%%%%%%%%%%%%%%%%%%%%%
tick_density = 35;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (adjust tick density) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_fontsize = 20;
x_tick_step = max([1, ceil(num_ticks_x / tick_density)]);                   % Adjust the divisor to control the density of ticks
y_tick_step = max([1, ceil(num_ticks_y / tick_density)]);                   % Adjust the divisor to control the density of ticks

% Convert MJD2000 departure dates to date strings
departure_date_str = cell(length(Departure_Grid), 1);
for i = 1:x_tick_step:length(Departure_Grid)                                % -1 here to remove last entry for nice visuals
    departure_date_strarr{i} = datestr(mjd20002date(Departure_Grid(i)), 'dd-mm-yyyy');
end

% Convert MJD2000 arrival dates to date strings
arrival_date_str = cell(length(Arrival_Grid), 1);
for j = 1:y_tick_step:length(Arrival_Grid)               
    arrival_date_strarr{j} = datestr(mjd20002date(Arrival_Grid(j)), 'dd-mm-yyyy');
end

figure;
hold on
    axis equal
    grid on
    Delta_V_Display = [1:0.2:5];
    [P,i] = contour(Departure_Grid, Arrival_Grid, Delta_V1_Capture', [Smallest_V:1:5], 'LevelList', Delta_V_Display, 'LineColor', [0.5, 0.0, 0.5], 'LabelSpacing', 500);
    [C,h] = contour(Departure_Grid, Arrival_Grid, ToF_Jumbo', 'LineColor', 'r', 'LineWidth', 0.5, 'LevelList', [90, 180, 270, 360, 450, 540]); % Adjust the LevelList as needed
    sp = plot(Departure_date_mjd2, Arrival_date_mjd2, 'MarkerFaceColor','m','MarkerEdgeColor','m', 'MarkerSize',10, 'Marker','*', 'LineStyle','none');
    
    
    clabel(P, i, 'Fontsize', plot_fontsize-6, 'Color',[0.5, 0.0, 0.5])                      % Add labels for the Porkchop contour lines
    clabel(C, h, 'FontSize', plot_fontsize-6, 'Color','r');                     % Add labels for the ToF contour lines
    
    % Plot Surroundings
    %Axis Option
    set(gca, 'XTick', Departure_Grid(1:x_tick_step:end), 'XTickLabel', departure_date_strarr(1:x_tick_step:end), 'FontSize',plot_fontsize);
    set(gca, 'YTick', Arrival_Grid(1:y_tick_step:end), 'YTickLabel', arrival_date_strarr(1:y_tick_step:end), 'FontSize',plot_fontsize);
    xlabel('Departure Date', 'FontSize',plot_fontsize);
    ylabel('Arrival Date', 'FontSize',plot_fontsize);
    
    legend([i, h, sp], {'Escape ΔV (km/s)','Time of ToF (days)','Baseline Date'}, 'FontSize', plot_fontsize, 'Location', 'northwest');  
                                                                   
    %Title Plotting:
    Departure_date_dt = datetime(mjd20002date(Departure_date_mjd2), 'ConvertFrom', 'datenum');
    Arrival_date_dt = datetime(mjd20002date(Arrival_date_mjd2), 'ConvertFrom', 'datenum');
    
    Departure_date_str = datestr(Departure_date_dt, 'dd-mm-yyyy');
    Arrival_date_str = datestr(Arrival_date_dt, 'dd-mm-yyyy');
    
    % title_str = sprintf('Departure ΔV with ToF Porkchop Plot, Selected Baseline: %s to %s : Initial Parking: %i km', Departure_date_str, Arrival_date_str, R_CircPark_E-rad_E);
    title_str = sprintf('ExoMars 2022 - ΔV with ToF Porkchop Plot, Selected Baseline: %s to %s', Departure_date_str, Arrival_date_str);

    title(title_str, 'FontSize',plot_fontsize);

hold off
 %% Technical grayscale plot with annotated values on graph                  C3 w/ tof                
% Calculate the step size for displaying ticks
%%%%%%%%%%%%%%%%% Adjust this value to choose number of dates displayed on each axis %%%%%%%%%%%%%%%%%%%%%% (adjust tick density)
tick_density = 40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (adjust tick density) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_fontsize = 20;
x_tick_step = max([1, ceil(num_ticks_x / tick_density)]);                   % Adjust the divisor to control the density of ticks
y_tick_step = max([1, ceil(num_ticks_y / tick_density)]);                   % Adjust the divisor to control the density of ticks

% Convert MJD2000 departure dates to date strings
departure_date_str = cell(length(Departure_Grid), 1);
for i = 1:x_tick_step:length(Departure_Grid)           
    departure_date_strarr{i} = datestr(mjd20002date(Departure_Grid(i)), 'dd-mm-yyyy');
end

% Convert MJD2000 arrival dates to date strings
arrival_date_str = cell(length(Arrival_Grid), 1);
for j = 1:y_tick_step:length(Arrival_Grid)               %-1 here to remove last entry for nice visuals
    arrival_date_strarr{j} = datestr(mjd20002date(Arrival_Grid(j)), 'dd-mm-yyyy');
end

figure;
hold on
    axis equal
    grid on
    Arrival_V_inf_Display = [1.0:0.2:8];
    Departure_C3_Display = [1:0.5:40];
    c3_ACTUAL_Jumbo = C3_calc_Jumbo'.^2;
    
    %contour matrix       %contour handles  
    [Avinf, h1] = contour(Departure_Grid, Arrival_Grid, Arrival_V_inf_Jumbo', [Smallest_Arrival_V_inf:0.2:10], 'LevelList', Arrival_V_inf_Display, 'LineStyle',':', 'LineWidth', 0.1, 'LineColor', 'k', 'LabelSpacing',400);
    [C3_cont, h2] = contour(Departure_Grid, Arrival_Grid, c3_ACTUAL_Jumbo, [Smallest_Departure_C3:0.5:10], 'LevelList', Departure_C3_Display, 'LineStyle','-', 'LineWidth', 0.1, 'LineColor', 'b','LabelSpacing',400);
    [ToF,ht] = contour(Departure_Grid, Arrival_Grid, ToF_Jumbo',  'LevelList', [90, 180, 270, 360, 450, 540], 'LineWidth', 0.2,'LineColor', 'r','LabelSpacing',300); % Adjust the LevelList as needed
    plot(Departure_date_mjd2, Arrival_date_mjd2, 'MarkerFaceColor','m','MarkerEdgeColor','m', 'MarkerSize',10, 'Marker','*', 'LineStyle','none')
    clabel(Avinf, h1, 'Fontsize', plot_fontsize-6, 'Color','k');                            % Add labels for the v_inf Arrival contour lines
    clabel(C3_cont, h2, 'Fontsize', plot_fontsize-6, 'Color','b');                          % Add labels for the C3 Departure contour lines
    clabel(ToF, ht, 'Fontsize', plot_fontsize-6, 'Color', 'r');                             % Add labels for the ToF contour lines
    
    % Plot Surroundings
    %Axis Option
    set(gca, 'XTick', Departure_Grid(1:x_tick_step:end), 'XTickLabel', departure_date_strarr(1:x_tick_step:end), 'Fontsize', plot_fontsize);
    set(gca, 'YTick', Arrival_Grid(1:y_tick_step:end), 'YTickLabel', arrival_date_strarr(1:y_tick_step:end), 'Fontsize', plot_fontsize);
    xlabel('Departure Date', 'Fontsize', plot_fontsize);
    ylabel('Arrival Date', 'Fontsize', plot_fontsize);
    legend([h1, h2, ht], {'Arrival V_∞ (km/s)', 'Departure C3 (km^2/s^2)', 'Transfer duration (days)'}, 'FontSize', plot_fontsize, 'Location', 'northwest');  
    
    %Title Plotting:
    Departure_date_dt = datetime(mjd20002date(Departure_date_mjd2), 'ConvertFrom', 'datenum');
    Arrival_date_dt = datetime(mjd20002date(Arrival_date_mjd2), 'ConvertFrom', 'datenum');
    
    Departure_date_str = datestr(Departure_date_dt, 'dd-mm-yyyy');
    Arrival_date_str = datestr(Arrival_date_dt, 'dd-mm-yyyy');
    
    title_str = sprintf('C3, v_∞ and TOF, Selected Baseline: %s to %s', Departure_date_str, Arrival_date_str);
    title(title_str, 'Fontsize', plot_fontsize);

hold off

 %% Mass Limit Plotting graphics                                          4 Mass Contour              
% Calculate the step size for displaying ticks
%%%%%%%%%%%%%%%%% Adjust this value to choose number of dates displayed on each axis %%%%%%%%%%%%%%%%%%%%%% (adjust tick density)
tick_density = 60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (adjust tick density) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_fontsize = 20;
x_tick_step = max([1, ceil(num_ticks_x / tick_density)]);                   % Adjust the divisor to control the density of ticks
y_tick_step = max([1, ceil(num_ticks_y / tick_density)]);                   % Adjust the divisor to control the density of ticks

% Convert MJD2000 departure dates to date strings
departure_date_str = cell(length(Departure_Grid), 1);                       %-1 here to remove last entry for nice visuals
for i = 1:x_tick_step:length(Departure_Grid)           
    departure_date_strarr{i} = datestr(mjd20002date(Departure_Grid(i)), 'dd-mm-yyyy');
end

% Convert MJD2000 arrival dates to date strings
arrival_date_str = cell(length(Arrival_Grid), 1);
for j = 1:y_tick_step:length(Arrival_Grid)                                  %-1 here to remove last entry for nice visuals
    arrival_date_strarr{j} = datestr(mjd20002date(Arrival_Grid(j)), 'dd-mm-yyyy');
end

figure;
hold on
    % set(figure, 'Color', 'black');                                             % Set background color to black
    axis equal
    grid on
    
    CalculatedFinalMass_Display = [5150:25:6500];
    
    cb = colorbar;
    cmap = colormap("pink");
    % cmap_limited = cmap(round(linspace(1, size(cmap, 1), 64)), :); % Limit colormap range
    % colormap(cmap_limited)
    % cb.Color = 'pink'; % Set colorbar text color to white
    ax = gca;
    ax.XAxis.MinorTick = 'on';
    ax.YAxis.MinorTick = 'on';
    ylabel(cb, 'Final Spacecraft Total Mass [kg]', 'Fontsize', plot_fontsize);
    %contour matrix                                                         % contour handles  
    [Avinf, h1] = contourf(Departure_Grid, Arrival_Grid, CalculatedFinalMass_Jumbo', [5100:25:6500], 'LevelList', CalculatedFinalMass_Display, 'LineStyle','-', 'LineColor', 'k', 'LabelSpacing',600);
    % plot(Departure_date_mjd2, Arrival_date_mjd2, 'MarkerFaceColor','m','MarkerEdgeColor','m', 'MarkerSize',10, 'Marker','*', 'LineStyle','none')                                                  % Marker for selected launch and arrival dayes
    clabel(Avinf, h1, 'Fontsize', plot_fontsize-4, 'Color',[0.3 0.3 0.3]);                            % Add labels for the v_inf Arrival contour lines
    
    % Plot Surroundings
    %Axis Option
    set(gca, 'XTick', Departure_Grid(1:x_tick_step:end), 'XTickLabel', departure_date_strarr(1:x_tick_step:end), 'Fontsize', plot_fontsize);
    set(gca, 'YTick', Arrival_Grid(1:y_tick_step:end), 'YTickLabel', arrival_date_strarr(1:y_tick_step:end), 'Fontsize', plot_fontsize);
    xlabel('Departure Date', 'Fontsize', plot_fontsize);
    ylabel('Arrival Date', 'Fontsize', plot_fontsize);
    
    %Title Plotting:
    Departure_date_dt = datetime(mjd20002date(Departure_date_mjd2), 'ConvertFrom', 'datenum');
    Arrival_date_dt = datetime(mjd20002date(Arrival_date_mjd2), 'ConvertFrom', 'datenum');
    
    Departure_date_str = datestr(Departure_date_dt, 'dd-mm-yyyy');
    Arrival_date_str = datestr(Arrival_date_dt, 'dd-mm-yyyy');
    
    title_str = sprintf('Contour Plot of Final Spacecraft Mass', Departure_date_str, Arrival_date_str);
    title(title_str, 'Fontsize', plot_fontsize);

hold off




% %% Mass Limit Plotting graphics                                          4 Mass Contour              
% % Calculate the step size for displaying ticks
% %%%%%%%%%%%%%%%%% Adjust this value to choose number of dates displayed on each axis %%%%%%%%%%%%%%%%%%%%%% (adjust tick density)
% tick_density = 70;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (adjust tick density) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% plot_fontsize = 20;
% x_tick_step = max([1, ceil(num_ticks_x / tick_density)]);                   % Adjust the divisor to control the density of ticks
% y_tick_step = max([1, ceil(num_ticks_y / tick_density)]);                   % Adjust the divisor to control the density of ticks
% 
% % Convert MJD2000 departure dates to date strings
% departure_date_str = cell(length(Departure_Grid), 1);                       %-1 here to remove last entry for nice visuals
% for i = 1:x_tick_step:length(Departure_Grid)           
%     departure_date_strarr{i} = datestr(mjd20002date(Departure_Grid(i)), 'dd-mm-yyyy');
% end
% 
% % Convert MJD2000 arrival dates to date strings
% arrival_date_str = cell(length(Arrival_Grid), 1);
% for j = 1:y_tick_step:length(Arrival_Grid)                                  %-1 here to remove last entry for nice visuals
%     arrival_date_strarr{j} = datestr(mjd20002date(Arrival_Grid(j)), 'dd-mm-yyyy');
% end
% 
% % Create a dark-themed plot
% fig = figure;
% set(fig, 'Color', 'black'); % Set background color to black
% hold on
% axis equal
% grid on
% 
% CalculatedFinalMass_Display = [4500:50:6500];
% 
% % Define custom colormap with limited range
% cmap = colormap("pink");
% cmap_limited = cmap(round(linspace(1, size(cmap, 1), 64)), :); % Limit colormap range
% 
% % Plot using limited colormap
% cb = colorbar;
% colormap(cmap_limited);
% cb.Color = 'white'; % Set colorbar text color to white
% ylabel(cb, 'Final Spacecraft Total Mass (kg)', 'Fontsize', plot_fontsize, 'Color', 'white');
% 
% % contour matrix                                                         % contour handles  
% [Avinf, h1] = contourf(Departure_Grid, Arrival_Grid, CalculatedFinalMass_Jumbo', [5000:100:6500], 'LevelList', CalculatedFinalMass_Display, 'LineStyle','-', 'LineColor', 'k', 'LabelSpacing',300);
% plot(Departure_date_mjd2, Arrival_date_mjd2, 'MarkerFaceColor','m','MarkerEdgeColor','m', 'MarkerSize',10, 'Marker','*', 'LineStyle','none')
% clabel(Avinf, h1, 'Fontsize', plot_fontsize-4, 'Color',[0.0 0.0 0.0]);                            % Add labels for the v_inf Arrival contour lines
% 
% % Plot Surroundings
% %Axis Option
% ax = gca; % Get current axes handle
% set(ax, 'XTick', Departure_Grid(1:x_tick_step:end), 'XTickLabel', departure_date_strarr(1:x_tick_step:end), 'Fontsize', plot_fontsize, 'Color', 'white');
% set(ax, 'YTick', Arrival_Grid(1:y_tick_step:end), 'YTickLabel', arrival_date_strarr(1:y_tick_step:end), 'Fontsize', plot_fontsize);
% xlabel('Departure Date', 'Fontsize', plot_fontsize, 'Color', 'white');
% ylabel('Arrival Date', 'Fontsize', plot_fontsize, 'Color', 'white');
% 
% % Set tick label colors to white
% % ax.XAxis.TickLabels.Color = 'white';
% % ax.YAxis.TickLabels.Color = 'white';
% 
% % ax.YAxis.TickValuesMode.Color = 'white';
% 
% %Title Plotting:
% Departure_date_dt = datetime(mjd20002date(Departure_date_mjd2), 'ConvertFrom', 'datenum');
% Arrival_date_dt = datetime(mjd20002date(Arrival_date_mjd2), 'ConvertFrom', 'datenum');
% 
% Departure_date_str = datestr(Departure_date_dt, 'dd-mm-yyyy');
% Arrival_date_str = datestr(Arrival_date_dt, 'dd-mm-yyyy');
% 
% title_str = sprintf('Contour Plot of Final Spacecraft Mass', Departure_date_str, Arrival_date_str);
% title(title_str, 'Fontsize', plot_fontsize, 'Color', 'white');
% 
% hold off