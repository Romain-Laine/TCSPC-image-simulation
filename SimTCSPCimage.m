%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TCSPC data simulator
% This code generates TCSPC data with Gaussian IRF and multiexcitation
% peaks, background (afterpulsing). Choose between photon scan, tau scan, uniform, or freestyle simulation options. 
% The data is saved as a TIFF stack via various methods. 
% 
% Laser Analytics Group: http://laser.ceb.cam.ac.uk/ 
% Dr Romain Laine rfl30@cam.ac.uk
% 2017-05-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

% Set the size of image you want to simulate. This also relates to the number of repeats for each condition -------
n_repeats = 1024;     % number of pixels in the vertical direction
n_conditions = 1024;  % number of pixels in the horizontal direction

% Change acquisition parameters --------------------
R = 80;      % Repetition rate of the laser in MHz units
T = 12.5;    % Acquisition window 0-T in ns units
n = 8;       % number of bits coding the TAC n = 8 --> 256 bins
Ap = 5;      % Afterpulsing in % --> background in TCSPC

% Change parameters for IRF-----------------------------
t0 = 1.5;    % offset in ns units
s = 0.15;    % standard deviation of Gaussian function in ns units
N_irf = 1e8; % MAXIMUM 10^8! (150,000 for 256 bins to get close to 16 bit histograms)

% 4 different types of simulations can be generated below

% 1. For Photon number scan ----------------------
% Here the photon counts are linearly spaced between N_min and N_max across
% a number of n_conditions
N_min = 100;    % minimum number of photons
N_max = 5000;   % maximum number of photons
Tau = 5;        % in ns units

% 2. For Tau scan -------------------------------
% Here the lifetimes are linearly spaced between Tau_min and Tau_max.
% All pixels have the same number of photons N. 
Tau_min = 0.05;  % minimum lifetime
Tau_max = 6;     % minimum lifetime
N = 5000;        % number of photons

% 3. For Uniform --------------------------------
% Here a uniform lifetime and number of photons is generated using the
% variables Tau and N set above

% 4. For Freestyle -----------------------------
% Here the list of lifetimes and photon counts are defined below in the
% code
FolderNameFreestyle = '_2 lines no BG';


%Save files to folder
Save_ON = 1;
FileName_append = '';
Folder_for_Save = 'D:\F3-CMM\TCSPC-image-simulation-master\TCSPC-image-simulation-master\';

% Simulation type -----------------------------
str = cell(1,4);
str{1} = 'Photon # scan';
str{2} = 'Tau scan';
str{3} = 'Uniform';
str{4} = 'Freestyle';

[Selection,ok] = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',str);

Sim_type = str{Selection};
if ~ok 
    return
end


%% -------------------------------------------------------------------------------------------------------
FolderName = [Sim_type, ...
    ' n=',num2str(n),'_T=',num2str(T),'ns_R=',num2str(R),'MHz_Ap=',num2str(Ap),'%',...
    '_IRF t0=',num2str(t0),'ns_s=',num2str(s),'ns_Nirf=',num2str(N_irf),'phot',...
    '_',num2str(n_repeats),'x',num2str(n_conditions)];

if strcmp(Sim_type, 'Tau scan')
    FolderName = [FolderName, '_Tau=',num2str(Tau_min),'-', num2str(Tau_max),'ns_N=',num2str(N),'phot'];
elseif strcmp(Sim_type, 'Photon # scan')
    FolderName = [FolderName, '_Phot=',num2str(N_min),'-', num2str(N_max),'phot_Tau=',num2str(Tau),'ns'];
elseif strcmp(Sim_type, 'Freestyle')
    FolderName = [FolderName, FolderNameFreestyle];
end

Folder_for_Save = [Folder_for_Save, FolderName];

% -------------------------------------------------------------------------

Sim_param = [t0, s, n, T, R, Ap];
h_wait = waitbar(0,'Wait for the data to be simulated...') ;

disp('Simulating data...');

tic
if strcmp(Sim_type,'Photon # scan')
    N_list = round(linspace(N_min,N_max,n_conditions));
    Tau_list = Tau*ones(1,n_conditions);
    
elseif strcmp(Sim_type,'Tau scan')
    N_list = N*ones(1,n_conditions);
    Tau_list = linspace(Tau_min,Tau_max,n_conditions);
    
elseif strcmp(Sim_type,'Uniform')
    N_list = N*ones(1,n_conditions);
    Tau_list = Tau*ones(1,n_conditions);
    
elseif strcmp(Sim_type,'Freestyle')
    % Here the list of photons and lifetime can be freely created. The vertical
    % line #i will be composed of repeats of a single conditions set by
    % Tau_list(i) and N_list(i)
    
    %     N_list = 50*ones(1,n_conditions);
    %     N_list(20:40) = N;
    %     Tau_list = Tau*ones(1,n_conditions);
    
    %     N_list = 50*ones(1,n_conditions);
    %     X = -15:15;
    %     CosIntensity = cos((pi/2)*X/15);
    %     N_list(66:(66+30)) = N_list(66:(66+30)) - round(34*CosIntensity);
    %     N_list(161:(161+30)) = N_list(161:(161+30)) + round(67*CosIntensity);
    %     Tau_list = 1.5*ones(1,n_conditions);
    %     Tau_list(66:(66+30)) = 0.5;
    %     Tau_list(161:(161+30)) = 3.5;
    
    N_list = 50*ones(1,n_conditions);
    N_list(20) = 5000;
    N_list(40) = 50000;
    Tau_list = 2.5*ones(1,n_conditions);
    
    
end

% Generate the IRF curve -----------------------------------------------------
[ t, Phot_number_IRF ] = SimDecayPhotonCount_Multiple( Sim_param, 0, N_irf, 1);

TCSPC_image = zeros(2^n, n_repeats, n_conditions);
for i = 1:n_conditions
    waitbar(i / n_conditions);
    [ t, Phot_number ] = SimDecayPhotonCount_Multiple( Sim_param, Tau_list(i), N_list(i), n_repeats );
    TCSPC_image(:,:,i) = Phot_number;
end
toc

if max(Phot_number_IRF(:)) >= 2^16
    disp('IRF DATA RESCALED!!');
    Phot_number_IRF = (2^16-1)*Phot_number_IRF/max(Phot_number_IRF(:));
end

if max(TCSPC_image(:)) >= 2^16
    disp('IMAGE DATA RESCALED!!');
    TCSPC_image = (2^16-1)*TCSPC_image/max(TCSPC_image(:));
end

Phot_number_IRF = uint16(Phot_number_IRF);
TCSPC_image = uint16(TCSPC_image);
close(h_wait);

% Save the data as tif stack
dt = 1000*T/2^n; % in ps
if Save_ON == 1
    if (exist(Folder_for_Save, 'dir') == 0)
        mkdir(Folder_for_Save);
    end
    SaveAsOMETIFF( permute(Phot_number_IRF,[2,3,1]), [Folder_for_Save, '\', 'IRF Stack', FileName_append], dt);
    SaveAsOMETIFF( permute(TCSPC_image,[2,3,1]), [Folder_for_Save, '\', 'Data Stack', FileName_append], dt);

%     SaveTIFFStack( Phot_number_IRF, ['IRF Stack',FileName_append], Folder_for_Save, SaveMethod, dt);
%     SaveTIFFStack( TCSPC_image, ['Data Stack', FileName_append], Folder_for_Save, SaveMethod, dt);
end


% ------------------------------------------------------------------------

Decay_display = reshape(TCSPC_image(:,1,:),2^n,n_conditions);

figure('Color','white','Units','normalized','position',[0.2 0.1 0.5 0.8],'name','IRF decay plot');
subplot(2,1,1)
plot(t,Phot_number_IRF); % display the first of all the repeats
xlabel 'Time (ns)'
ylabel 'Photon counts'
subplot(2,1,2)
semilogy(t,Phot_number_IRF); % display the first of all the repeats
xlabel 'Time (ns)'
ylabel 'Photon counts'

figure('Color','white','Units','normalized','position',[0.2 0.1 0.5 0.8],'name','Decay plots');
subplot(2,1,1)
plot(t,Decay_display); % display the first of all the repeats
xlabel 'Time (ns)'
ylabel 'Photon counts'

subplot(2,1,2)
semilogy(t,Decay_display); % display the first of all the repeats
xlabel 'Time (ns)'
ylabel 'Photon counts'

disp('------------------------');
disp('All done.')

