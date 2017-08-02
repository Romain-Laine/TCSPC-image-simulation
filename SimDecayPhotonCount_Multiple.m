function [ t, Phot_number ] = SimDecayPhotonCount_Multiple( Sim_param, Tau, N, n_decay )
%SIMDECAYPHOTONCOUNT Summary of this function goes here
% n_decay is the number of decay to simulate

Display = 0;

% IRF parameters
t0 = Sim_param(1);    % offset in ns
s = Sim_param(2);    % standard deviation of Gaussian function in ns

% Acquisition parameters
n = Sim_param(3);       % number of bits coding the TAC n = 8 --> 256 bins
T = Sim_param(4);      % Acquisition window 0-T in ns
R = Sim_param(5);      % Repetition rate of the laser in MHz
Ap = Sim_param(6);      % Afterpulsing in % --> background in TCSPC


% Generating the arrival times
if Tau > 0
    u_f = rand(N,n_decay);
    t_f = Tau*log(1./(1-u_f));
else
    t_f = zeros(N,n_decay);
end


if t0 <= 0 || s <= 0
    t_irf = zeros(N,n_decay);
else
    u_irf = rand(N,n_decay);
    t_irf = t0 - s*erfinv(erf(t0/s) - 2*u_irf);
end


t_tot = t_f + t_irf;

% Incomplete decays
T_rep = 1000/R;
t_tot = mod(t_tot,T_rep);

% Afterpulsing background
N_bg = ceil(Ap/100*N);
t_Ap = T*rand(N_bg,n_decay);
t_tot = [t_tot; t_Ap];

% Histogramming
t = linspace(0,T,2^n+1);
t(end) = [];
Phot_number = histc(t_tot,[t t(end)+t(2)-t(1)]); % histc works in columns
Phot_number(end,:) = [];


if Display == 1
    % Displaying
    figure('Color','white');
    semilogy(t,Phot_number);
    xlabel 'Time (ns)'
    ylabel 'Photon counts'
    title(['N = ',num2str(N),' photons, Tau = ',num2str(Tau),' ns'])
    
    figure('Color','white');
    plot(t,Phot_number);
    xlabel 'Time (ns)'
    ylabel 'Photon counts'
    title(['N = ',num2str(N),' photons, Tau = ',num2str(Tau),' ns'])
end

end

