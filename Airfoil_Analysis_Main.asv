%%ASEN 2502 Wind Tunnel Lab
% Author:
% Lab Description: Employ basic experimental wind tunnel testing procedures, 
% conduct flow measurements, and develop an awareness for sources of 
% error/error analysis as part of a team.  Validate the AES PILOT 
% Low-SpeedWind Tunnel through the evaluation of a Clark Y-14 airfoil 
% and comparingthe results to the published National Advisory Committee 
% for Aeronautics (NACA) results from 1938.
% 
% Inputs:
% - Raw Wind Tunnel CSV data files for tests
% - Airfoil Static Port Locations
% - Comparison NACA Clark Y Data
% 
% Outputs:
% 
%
%% Intialize Workspace
clear; clc; close all; % Clear workspace, command window, and close all figures

addpath(genpath('15 mps Data Files')); %Adds 15 m/s test data files folder and subfolders
addpath(genpath('30 mps Data Files')); %Adds 30 m/s test data files folder and subfolders

%% Import Airfoil Port Locations
PortsTable = readtable('Port_Locations.xlsx','Sheet','Port_Locations'); %Read in CSV file with port locations
Ports = table2array(PortsTable);
SegmentsTable = readtable('Port_Locations.xlsx','Sheet','Segments'); %Read in segment information from CSV file
Segments = table2array(SegmentsTable(:,2:3));
%% Search Data Folders, Pull File Names & Count Data Files
% Get filenames for test data files
fileLoc = '15 mps Data Files/';
list = dir([fileLoc, '*AoA*']); % This lists all files in fileLoc with 'WTData' in the file name
numFiles = length(list);

for i = 1:numFiles
    fileNames15{i} = [fileLoc, list(i).name]; % This makes a string of the complete file name with the path in front of it
    AoA15(i) = str2num(char(extractBetween(list(i).name, 'AoA_', '.csv'))); % Finds the angle of attack value in the name and make it a usable number
end
fileNames15 = string(fileNames15); % Convert to char array from cell array
AoA_Count15 = length(AoA15); %Counts number of AoAs tested
%{
fileLoc = '30 mps Data Files/';
list = dir([fileLoc, '*AoA*']); % This lists all files in fileLoc with 'WTData' in the file name
numFiles = length(list);

for i = 1:numFiles
    fileNames30{i} = [fileLoc, list(i).name]; % This makes a string of the complete file name with the path in front of it
    AoA30(i) = str2num(char(extractBetween(list(i).name, 'AoA_', '.csv'))); % Finds the angle of attack value in the name and make it a usable number
end
fileNames30 = string(fileNames30); % Convert to char array from cell array
AoA_Count30 = length(AoA30); %Counts number of AoAs tested
%}
%% Ingest Data Files and Data Conditioning
% Averaging Raw Data Samples for each Velocity & AOA Tested
% Derived Values for each Velocity & AoA tested
    % Average Test Section Static Pressure (done by student code)
    % Average Airfoil Port Local Velocity (done by student code)
    % Average Airfoil Port Local Static Pressure (done by student code)

% Initialize Storage Arrays
Data15 = zeros(numFiles,25); %Conditioned wind tunnel data file 15 m/s
%%%%Data30 = zeros(numFiles,25); %Conditioned wind tunnel data file 30 m/s

%% Ingest and Condition Data
% Averaging Raw Data Samples for each Velocity & AOA Tested
% Derived Values for each Velocity & AoA tested
    % Average Test Section Static Pressure (done by student code)
    % Average Airfoil Port Local Velocity (done by student code)
    % Average Airfoil Port Local Static Pressure (done by student code)

for j = 1:numFiles
    RawData15 = readmatrix(fileNames15(j),'NumHeaderLines',1); % load the data
    %Condition data
    Data15(j,1) = AoA15(j); %Sets AoA (deg) as 1st column
    Data15(j,2) = mean(RawData15(:,4)); %Sets mean of velocity measurements as 2nd column
    Data15(j,3) = mean(RawData15(:,2)); %Sets mean of total pressure (atmosphere) as 3rd column
    Data15(j,4) = mean(RawData15(:,1)); %Sets mean of temperature (atmosphere) as 4th column
    Data15(j,5) = mean(RawData15(:,3)); %Sets mean of calculaed density (atmosphere) as 5th column
    Data15(j,6) = mean(RawData15(:,5)); %Sets mean of test section dynamic pressure as 6th column
    Data15(j,7) = Data15(j,3) - Data15(j,6); %Calculates test section static pressure by calculating Po - q

    %Calculate airfoil port static pressures from measured differential
    %pressues (P_port - P_static_test) & sets columns 8 - 24 as airfoil static port pressure values (calculated) 

    %To get the port static pressure, we have to subtract out the static
    %pressure in  the atmosphere that was calculated in data15(j,7)

    %Since scanivalve ports are differential pressure from the 
    %test section, we dont need to subtract out anything and the TE should
    %be 0
    Data15(j,8) = mean(RawData15(:,15)); %Airfoil scanivalve port 1
    Data15(j,9) = mean(RawData15(:,16)); %Airfoil scanivalve port 2
    Data15(j,10) = mean(RawData15(:,17)); %Airfoil scanivalve port 3
    Data15(j,11) = mean(RawData15(:,18)); %Airfoil scanivalve port 4
    Data15(j,12) = mean(RawData15(:,19)); %Airfoil scanivalve port 5
    Data15(j,13) = mean(RawData15(:,20)); %Airfoil scanivalve port 6
    Data15(j,14) = mean(RawData15(:,21)); %Airfoil scanivalve port 7
    Data15(j,15) = mean(RawData15(:,22)); %Airfoil scanivalve port 8
    Data15(j,16) = mean(RawData15(:,23)); %Airfoil scanivalve port 9
    Data15(j,17) = 0; %TE - Assume Airfoil trailing edge pressure = test section static pressure
    Data15(j,18) = mean(RawData15(:,24)); %Airfoil scanivalve port 10
    Data15(j,19) = mean(RawData15(:,25)); %Airfoil scanivalve port 11
    Data15(j,20) = mean(RawData15(:,26)); %Airfoil scanivalve port 12
    Data15(j,21) = mean(RawData15(:,27)); %Airfoil scanivalve port 13
    Data15(j,22) = mean(RawData15(:,28)); %Airfoil scanivalve port 14
    Data15(j,23) = mean(RawData15(:,29)); %Airfoil scanivalve port 15
    Data15(j,24) = mean(RawData15(:,30)); %Airfoil scanivalve port 16
    Data15(j,25) = Data15(j,8); %Repeats port 1 (leading edge)
end
Data15 = sortrows(Data15,1); %Sorts data by increasing AoA
%{
for j = 1:numFiles
    RawData30 = readmatrix(fileNames30(j),'NumHeaderLines',1); % load the data
    %Condition data
    Data30(j,1) = AoA30(j); %Sets AoA (deg) as 1st column
    Data30(j,2) = ; %Sets mean of velocity measurements as 2nd column
    Data30(j,3) = ; %Sets mean of total pressure (atmosphere) as 3rd column
    Data30(j,4) = ; %Sets mean of temperature (atmosphere) as 4th column
    Data30(j,5) = ; %Sets mean of calculaed density (atmosphere) as 5th column
    Data30(j,6) = ; %Sets mean of test section dynamic pressure as 6th column
    Data30(j,7) = ; %Calculates test section static pressure by calculating Po - q
    %Calculate airfoil port static pressures from measured differential
    %pressues (P_port - P_static_test) & sets columns 8 - 24 as airfoil static port pressure values calculated 
    Data30(j,8) = ; %Airfoil scanivalve port 1
    Data30(j,9) = ; %Airfoil scanivalve port 2
    Data30(j,10) = ; %Airfoil scanivalve port 3
    Data30(j,11) = ; %Airfoil scanivalve port 4
    Data30(j,12) = ; %Airfoil scanivalve port 5
    Data30(j,13) = ; %Airfoil scanivalve port 6
    Data30(j,14) = ; %Airfoil scanivalve port 7
    Data30(j,15) = ; %Airfoil scanivalve port 8
    Data30(j,16) = ; %Airfoil scanivalve port 9
    Data30(j,17) = ; %TE  - Assume Airfoil trailing edge pressure = test section static pressure
    Data30(j,18) = ; %Airfoil scanivalve port 10
    Data30(j,19) = ; %Airfoil scanivalve port 11
    Data30(j,20) = ; %Airfoil scanivalve port 12
    Data30(j,21) = ; %Airfoil scanivalve port 13
    Data30(j,22) = ; %Airfoil scanivalve port 14
    Data30(j,23) = ; %Airfoil scanivalve port 15
    Data30(j,24) = ; %Airfoil scanivalve port 16
    Data30(j,25) = ; %Repeats port 1 (leading edge)

end
Data30 = sortrows(Data30,1); %Sorts data by increasing AoA
%}


%% Determine Forces & Analyze Results
% Pressure Distribution for each Velocity & AoA tested (done by student code)
Calculated_Values = zeros(30,5);

dX = Ports(1,2) - Ports(10,2);
dY = Ports(1,3) - Ports(10,3);
ChordLength = sqrt(dX^2 + dY^2);

for i = 1:30

    N_upper_sum = 0;
    A_upper_sum = 0;

    N_lower_sum = 0;
    A_lower_sum = 0;

    %Loops over upper surface, from leading edge to trailing edge
    for j = 8:16
        Index = j-7; % 1...9
        N_upper_sum = N_upper_sum  - 0.5*(Data15(i,j) + Data15(i,j+1)) * Segments(Index,1);
        A_upper_sum = A_upper_sum  + 0.5*(Data15(i,j) + Data15(i,j+1)) * Segments(Index,2);
    end

    %Loops over lower surface from leading edge to trailing edge
    for j = 25:-1:18
        Index = -1*j + 35;
        N_lower_sum = N_lower_sum + 0.5*(Data15(i,j) + Data15(i,j-1)) * Segments(Index,1);
        A_lower_sum = A_lower_sum - 0.5*(Data15(i,j) + Data15(i,j-1)) * Segments(Index,2);
    end
    
    NormalForce = N_upper_sum + N_lower_sum;
    AxialForce  = A_upper_sum + A_lower_sum;

    %%Storing all calculated values into array
    %Order of storage: AoA, Fa, Fn, Fl, Cl
    Calculated_Values(i,1) = Data15(i,1);
    Calculated_Values(i,2) = AxialForce;
    Calculated_Values(i,3) = NormalForce;

    %Lift force calculations
    Lift_Force = NormalForce * cosd(Data15(i,1)) - AxialForce * sind(Data15(i,1));
    Calculated_Values(i,4) = Lift_Force;

    %Coefficient of Lift Calculation
    FreeStreamDynPressure = 0.5 * Data15(i,5) * Data15(i,2)^2;
    Coefficient_lift = Lift_Force / (FreeStreamDynPressure * ChordLength);
    Calculated_Values(i,5) = Coefficient_lift;

end





% Normal and Axial Force components Velocity & AoA tested (done by student code)
% Lift & Coefficient of Lift Velocity & AoA tested (done by student code)

%% Plots
% Velocity vs normalized chord (x/c)
figure;
subplot(3,1,1);


% Coefficient of Pressure vs normalized chord (x/c)
subplot(3,1,2);

% Coefficient of Lift vs Angle of Attack
subplot(3,1,3);
plot(Calculated_Values(:,1) , Calculated_Values(:,5));

