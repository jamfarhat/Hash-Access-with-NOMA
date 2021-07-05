close all; clear; clc;
%--------------------------------------------------------------------------------
% Parameters used to plot the figures
FigureWidth = 16;           
FigureHeight = 10;           
LegendFontSize = 12;    % 
AxisFontSize = 12;
AxisLabelFontSize = 12;
LineWidth = 1.5;
MarkerSize = 8;
colorHA = [0.7410, 0.4470, 0];
colorNOHA = [0, 0.4470, 0.7410];

% 
% ThroughputvsSNRvsDifficulty = 1;
% ThroughputvsDifficulty = 1;
% ThroughputvsN = 1;
% ThroughputvsSNRvsR = 1;
% GainvsR = 1;

%% PART 1: Throughput vs SNR (Various kappa)

SNRdBmVec = [-60:5:20];             % SNR (dBm)
SNRVec = 10.^(SNRdBmVec./10)*1e-3;  % SNR (linear scale)
sigma = (10.^(-104./10))*1e-3;      % noise variance (W)
DifficultyVec = [0.85 0.9 0.95];    % difficulty to access the channel - kappa
% alpha = 4;                        % path-loss exponent
R = 1;                              % transmission rate (bps/Hz)
D0 = 100;                           % radius of the circular area (m)
N = 20;                             % number of iot devices (active/inactive) in the area

ThroughputHA = zeros(length(DifficultyVec),length(SNRVec));
ThroughputHATheo = zeros(length(DifficultyVec),length(SNRVec));
ThroughputNOHA = zeros(length(DifficultyVec),length(SNRVec));
ThroughputNOHATheo = zeros(length(DifficultyVec),length(SNRVec));
for diff_id = 1:length(DifficultyVec)
    disp([ 'Part 1/5 : [' num2str(diff_id) '/' num2str(length(DifficultyVec)) '] ']);       
    Difficulty = DifficultyVec(diff_id);
    for snr_id = 1:length(SNRVec)
        SNR = SNRVec(snr_id);
        [ThroughputHATheo(diff_id,snr_id),...
         ThroughputHA(diff_id,snr_id)] = calcThroughput(Difficulty,sigma/SNR,D0,N,R,'HA');
        [ThroughputNOHATheo(diff_id,snr_id),...
         ThroughputNOHA(diff_id,snr_id)] = calcThroughput(Difficulty,sigma/SNR,D0,N,R,'NOHA');
    end
    
end

figure1 = figure;
set(gcf, 'Units', 'centimeters');
afFigurePosition = [2 2 FigureWidth FigureHeight]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
step = 2;
plot(SNRdBmVec, ThroughputNOHATheo(1,:), '-', 'LineWidth', LineWidth, 'Color', colorNOHA);
hold on;
plot(SNRdBmVec(1:step:end), ThroughputNOHA(1,1:step:end), 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize , 'Color', colorNOHA);
plot(SNRdBmVec, ThroughputNOHATheo(2,:), '-', 'LineWidth', LineWidth, 'Color', colorNOHA);
plot(SNRdBmVec(1:step:end), ThroughputNOHA(2,1:step:end), 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize , 'Color', colorNOHA);
plot(SNRdBmVec, ThroughputNOHATheo(3,:), '-', 'LineWidth', LineWidth, 'Color', colorNOHA);
plot(SNRdBmVec(1:step:end), ThroughputNOHA(3,1:step:end), 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize , 'Color', colorNOHA);
grid on;
xlabel('Transmit Power $P$ (dBm)','Interpreter','LaTeX','Fontsize',AxisFontSize);
ylabel('Throughput $\rho^{\mathrm{(NOHA)}}$','Interpreter','LaTeX','Fontsize',AxisFontSize);
legend({'Analytical from (18)','Simulation'},...
       'fontsize',LegendFontSize,'location','northwest');
%xlim([0.65 1]);
%annotation(figure1,'textarrow',[0.759933774834437 0.759933774834437],...
%    [0.674603174603172 0.727513227513227],'String',{'\kappa = 0.95'},'fontsize',AxisLabelFontSize);
%annotation(figure1,'textarrow',[0.849337748344371 0.849337748344371],...
%    [0.537037037037036 0.6005291005291],'String',{'\kappa = 0.85'},'fontsize',AxisLabelFontSize);
%annotation(figure1,'textarrow',[0.844370860927152 0.844370860927152],...
%    [0.875661375661371 0.812169312169312],'String',{'\kappa = 0.9'},'fontsize',AxisLabelFontSize);
set(gca,'fontsize',AxisLabelFontSize);

%% PART 2: Throughput vs Difficulty
clc;
DifficultyVec = [0.65:0.01:1]; 
sigma = (10.^(-104./10))*1e-3;
D0 = 100;
N = 20; 

SNRdBm1 = 0; 
SNR1 = 10.^(SNRdBm1./10)*1e-3;
R1 = 1;

SNRdBm2 = 0; 
SNR2 = 10.^(SNRdBm2./10)*1e-3;
R2 = 0.5;

ThroughputHA1 = zeros(1,length(DifficultyVec));
ThroughputHATheo1 = zeros(1,length(DifficultyVec));
ThroughputHA2 = zeros(1,length(DifficultyVec));
ThroughputHATheo2 = zeros(1,length(DifficultyVec));
ThroughputNOHA1 = zeros(1,length(DifficultyVec));
ThroughputNOHATheo1 = zeros(1,length(DifficultyVec));
ThroughputNOHA2 = zeros(1,length(DifficultyVec));
ThroughputNOHATheo2 = zeros(1,length(DifficultyVec));

for diff_id = 1:length(DifficultyVec)
        disp([ 'Part 2/5 : [' num2str(diff_id) '/' num2str(length(DifficultyVec)) '] ']);
        Difficulty = DifficultyVec(diff_id);
        [ThroughputHATheo1(diff_id),...
         ThroughputHA1(diff_id),...
         OptimalKappaHA1,...
         MaximumThroughputHA1] = calcThroughput(Difficulty,sigma/SNR1,D0,N,R1,'HA');
     
         [ThroughputHATheo2(diff_id),...
         ThroughputHA2(diff_id),...
         OptimalKappaHA2,...
         MaximumThroughputHA2] = calcThroughput(Difficulty,sigma/SNR2,D0,N,R2,'HA');
        
        [ThroughputNOHATheo1(diff_id),...
         ThroughputNOHA1(diff_id),...
         OptimalKappaNOHA1,...
         MaximumThroughputNOHA1] = calcThroughput(Difficulty,sigma/SNR1,D0,N,R1,'NOHA');
     
     [ThroughputNOHATheo2(diff_id),...
         ThroughputNOHA2(diff_id),...
         OptimalKappaNOHA2,...
         MaximumThroughputNOHA2] = calcThroughput(Difficulty,sigma/SNR2,D0,N,R2,'NOHA');
end

figure2 = figure;
set(gcf, 'Units', 'centimeters');
afFigurePosition = [2 2 FigureWidth FigureHeight]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
step = 2;
plot(DifficultyVec, ThroughputHATheo1, '--', 'LineWidth', LineWidth, 'Color', colorHA);
hold on;
plot(DifficultyVec(1:step:end), ThroughputHA1(1:step:end), 's', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize ,'Color', colorHA);
plot(OptimalKappaHA1, MaximumThroughputHA1, 's', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize ,'MarkerFaceColor', colorHA, 'Color', colorHA);

plot(DifficultyVec, ThroughputNOHATheo1, '-', 'LineWidth', LineWidth, 'Color', colorNOHA);
plot(DifficultyVec(1:step:end), ThroughputNOHA1(1:step:end), 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize , 'Color', colorNOHA);
plot(OptimalKappaNOHA1, MaximumThroughputNOHA1, 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize ,'MarkerFaceColor', colorNOHA, 'Color', colorNOHA);

plot(DifficultyVec, ThroughputHATheo2, '--', 'LineWidth', LineWidth, 'Color', colorHA);
plot(DifficultyVec(1:step:end), ThroughputHA2(1:step:end), 's', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize ,'Color', colorHA);
plot(OptimalKappaHA2, MaximumThroughputHA2, 's', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize ,'MarkerFaceColor', colorHA, 'Color', colorHA);

plot(DifficultyVec, ThroughputNOHATheo2, '-', 'LineWidth', LineWidth, 'Color', colorNOHA);
plot(DifficultyVec(1:step:end), ThroughputNOHA2(1:step:end), 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize , 'Color', colorNOHA);
plot(OptimalKappaNOHA2, MaximumThroughputNOHA2, 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize ,'MarkerFaceColor', colorNOHA, 'Color', colorNOHA);

grid on;
xlabel('Difficulty $\kappa$','Interpreter','LaTeX','Fontsize',AxisFontSize);
ylabel('Throughput $\rho^{\mathrm{(sch)}}$','Interpreter','LaTeX','Fontsize',AxisFontSize);
legend({'HA (Analytical)','HA (Simulation)', '\kappa_{\ast}^{(HA)} from (4)',...
       'NOHA (Analytical)','NOHA (Simulation)', '\kappa_{\ast}^{(NOHA)} from (20)'},...
       'fontsize',LegendFontSize,'location','northwest');
xlim([0.65 1]);
set(gca,'fontsize',AxisLabelFontSize);
% Create arrow
annotation(figure2,'arrow',[0.417395833333334 0.501731770833333],...
    [0.365208333333334 0.170416666666667]);

% Create arrow
annotation(figure2,'arrow',[0.733443708609272 0.743377483443708],...
    [0.724867724867725 0.806878306878306],'Color',[0.8 0.8 0.8]);

% Create arrow
annotation(figure2,'arrow',[0.735099337748344 0.791390728476821],...
    [0.62962962962963 0.304232804232804],'Color',[0.8 0.8 0.8]);

% Create arrow
annotation(figure2,'arrow',[0.748344370860927 0.781456953642384],[0.62962962962963 0.58994708994709],...
    'Color',[0.8 0.8 0.8]);

% Create arrow
annotation(figure2,'arrow',[0.721854304635762 0.725165562913907],...
    [0.626984126984127 0.497354497354497],'Color',[0.8 0.8 0.8]);

% Create arrow
annotation(figure2,'arrow',[0.519854442604857 0.55794701986755],...
    [0.678424933862433 0.288359788359788]);

% Create arrow
annotation(figure2,'arrow',[0.52151007174393 0.591059602649007],[0.675779431216929 0.6005291005291]);

% Create textbox
annotation(figure2,'textbox',...
    [0.648886416942607 0.64021164021164 0.170650006898454 0.066300595238094],'String',{'Optimal values'},...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1],...
    'BackgroundColor',[0.941176470588235 0.941176470588235 0.941176470588235]);

% Create textbox
annotation(figure2,'textbox',...
    [0.458711196192055 0.708994708994709 0.177050393211919 0.0478105158730178],...
    'String',{'$\xi=3$ bps/Hz'},...
    'Interpreter','latex',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure2,'textbox',...
    [0.24541666666667 0.386243386243386 0.188515624999998 0.0451107804232834],...
    'String',{'$\xi=1$ bps/Hz'},...
    'Interpreter','latex',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create arrow
annotation(figure2,'arrow',[0.417218543046358 0.503311258278146],...
    [0.365079365079365 0.306878306878307]);

%% PART 3: Throughput vs Number of IoT devices
clc;
SNRdBm = 0; 
SNR = 10.^(SNRdBm1./10)*1e-3;
R = 1;
D0 = 100;
Nvec = [2:10 20:10:100];                   
ThroughputHA = zeros(1,length(Nvec));
ThroughputNOHA = zeros(1,length(Nvec));
MaximumThroughputHATheo = zeros(1,length(Nvec));
MaximumThroughputNOHATheo = zeros(1,length(Nvec));
MaximumThroughputHA = zeros(1,length(Nvec));
MaximumThroughputNOHA = zeros(1,length(Nvec));

for n_id = 1:length(Nvec)
    disp([ 'Part 3/5 : [' num2str(n_id) '/' num2str(length(Nvec)) '] ']);
    N = Nvec(n_id);
    % HA: fixed difficulty
    [Unused,...
     ThroughputHA(n_id),...
     OptimalKappaHA,...
     Unused] = calcThroughput(0.9,sigma/SNR,D0,N,R,'HA');
    % HA: optimal difficulty
    [MaximumThroughputHATheo(n_id),...
     MaximumThroughputHA(n_id),...
     Unused,Unused] = calcThroughput(OptimalKappaHA,sigma/SNR,D0,N,R,'HA');
     
    % NOHA: fixed difficulty
    [Unused,...
     ThroughputNOHA(n_id),...
     OptimalKappaNOHA,...
     Unused] = calcThroughput(0.9,sigma/SNR,D0,N,R,'NOHA');
    % NOHA: optimal difficulty
     [MaximumThroughputNOHATheo(n_id),...
      MaximumThroughputNOHA(n_id),...
      Unused,Unused] = calcThroughput(OptimalKappaNOHA,sigma/SNR,D0,N,R,'NOHA');
  
end

figure
set(gcf, 'Units', 'centimeters');
afFigurePosition = [2 2 FigureWidth FigureHeight]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
step = 1;
semilogx(Nvec, MaximumThroughputHATheo, '--', 'LineWidth', LineWidth, 'Color', colorHA);
hold on;
semilogx(Nvec(1:step:end), MaximumThroughputHA(1:step:end), 's', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize ,'Color', colorHA,'MarkerFaceColor', colorHA);
semilogx(Nvec, ThroughputHA, 's', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize ,'Color', colorHA);

semilogx(Nvec, MaximumThroughputNOHATheo, '-', 'LineWidth', LineWidth, 'Color', colorNOHA);
semilogx(Nvec(1:step:end), MaximumThroughputNOHA(1:step:end), 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize , 'Color', colorNOHA, 'MarkerFaceColor', colorNOHA);
semilogx(Nvec, ThroughputNOHA, 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize ,'Color', colorNOHA);

grid on;
xlabel('Number of IoT Devices $N$','Interpreter','LaTeX','Fontsize',AxisFontSize);
ylabel('Throughput $\rho^{\mathrm{(sch)}}$','Interpreter','LaTeX','Fontsize',AxisFontSize);
legend({'HA (Analytical) - Optimal \kappa','HA (Simulation) - Optimal \kappa',...
       'HA (Simulation) - \kappa = 0.9', 'NOHA (Analytical) - Optimal \kappa',...
       'NOHA (Simulation) - Optimal \kappa','NOHA (Simulation) - \kappa = 0.9'},...
       'fontsize',LegendFontSize,'location','northeast');
set(gca,'fontsize',AxisLabelFontSize);
xlim([2 Nvec(end)])
ylim([0 2])

%% PART 4: Throughput vs SNR (Various R)
clc;
SNRdBmVec = [-40:5:30];             % SNR (dBm)
SNRVec = 10.^(SNRdBmVec./10)*1e-3;  % SNR (linear scale)
sigma = (10.^(-104./10))*1e-3;      % noise variance (W)
RVec = [1 3 5];                     % spectral efficiency (bps/Hz)
Difficulty = 0.95;                  % difficulty to access the channel - kappa
% alpha = 4;                        % path-loss exponent
D0 = 100;                           % radius of the circular area (m)
N = 20;                             % number of iot devices (active/inactive) in the area

ThroughputHA = zeros(length(RVec),length(SNRVec));
ThroughputHATheo = zeros(length(RVec),length(SNRVec));
ThroughputNOHA = zeros(length(RVec),length(SNRVec));
ThroughputNOHATheo = zeros(length(RVec),length(SNRVec));

for r_id = 1:length(RVec)
    disp([ 'Part 4/5 : [' num2str(r_id) '/' num2str(length(RVec)) '] ']);       
    R = RVec(r_id);
    for snr_id = 1:length(SNRVec)
        SNR = SNRVec(snr_id);
        [ThroughputHATheo(r_id,snr_id),...
         ThroughputHA(r_id,snr_id)] = calcThroughput(Difficulty,sigma/SNR,D0,N,R,'HA');
        [ThroughputNOHATheo(r_id,snr_id),...
         ThroughputNOHA(r_id,snr_id)] = calcThroughput(Difficulty,sigma/SNR,D0,N,R,'NOHA');
    end
end

figure;
set(gcf, 'Units', 'centimeters');
afFigurePosition = [2 2 FigureWidth FigureHeight]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
step = 1;
plot(SNRdBmVec, ThroughputNOHATheo(1,:), '-', 'LineWidth', LineWidth, 'Color', colorNOHA);
hold on;
plot(SNRdBmVec, ThroughputNOHATheo(2,:), '--', 'LineWidth', LineWidth, 'Color', colorNOHA);
plot(SNRdBmVec, ThroughputNOHATheo(3,:), '-.', 'LineWidth', LineWidth, 'Color', colorNOHA);
plot(SNRdBmVec(1:step:end), ThroughputNOHA(1,1:step:end), 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize , 'Color', colorNOHA);
plot(SNRdBmVec(1:step:end), ThroughputNOHA(2,1:step:end), 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize , 'Color', colorNOHA);
plot(SNRdBmVec(1:step:end), ThroughputNOHA(3,1:step:end), 'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize , 'Color', colorNOHA);
grid on;
xlabel('Transmit Power $P$ (dBm)','Interpreter','LaTeX','Fontsize',AxisFontSize);
ylabel('Throughput $\rho^{\mathrm{(NOHA)}}$','Interpreter','LaTeX','Fontsize',AxisFontSize);
legend({'Analytical, \xi=1 bps/Hz','Analytical, \xi=3 bps/Hz','Analytical, \xi=5 bps/Hz','Simulation'},...
       'fontsize',LegendFontSize,'location','northwest');

set(gca,'fontsize',AxisLabelFontSize);

%% PART 5: Relative Gain vs R
clc;
SNRdBm = 0;             % SNR (dBm)
SNR = 10.^(SNRdBm./10)*1e-3;  % SNR (linear scale)
sigma = (10.^(-104./10))*1e-3;      % noise variance (W)
RVec = [0.1 0.5 1 3 5];             % spectral efficiency (bps/Hz)
% alpha = 4;                        % path-loss exponent
D0 = 100;                           % radius of the circular area (m)
N = 20;                             % number of iot devices (active/inactive) in the area

ThroughputHA = zeros(1,length(RVec));
ThroughputHATheo = zeros(1,length(RVec));
ThroughputNOHA = zeros(1,length(RVec));
ThroughputNOHATheo = zeros(1,length(RVec));

for r_id = 1:length(RVec)
    disp([ 'Part 5/5 : [' num2str(r_id) '/' num2str(length(RVec)) '] ']);       
    R = RVec(r_id);
    [Unused,...
     Unused,...
     OptimalKappaHA,...
     Unused] = calcThroughput(0.9,sigma/SNR,D0,N,R,'HA');
 
    [ThroughputHATheo(1,r_id),...
     ThroughputHA(1,r_id)] = calcThroughput(OptimalKappaHA,sigma/SNR,D0,N,R,'HA');
     
    [Unused,...
     Unused,...
     OptimalKappaNOHA,...
     Unused] = calcThroughput(0.9,sigma/SNR,D0,N,R,'NOHA');
 
    [ThroughputNOHATheo(1,r_id),...
     ThroughputNOHA(1,r_id)] = calcThroughput(OptimalKappaNOHA,sigma/SNR,D0,N,R,'NOHA');
end
clc;
GainSimul = ThroughputNOHA./ThroughputHA
GainTheo = ThroughputNOHATheo./ThroughputHATheo
