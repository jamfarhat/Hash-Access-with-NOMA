function [ThroughputTheo,...
          ThroughputSimul,...
          OptimalKappa,...
          MaximumThroughput] = calcThroughput(diff,sigmaP,D0,N,R,sch)
%==========================================================================
% Inputs:
%           diff: Difficulty in accesing the channel
%           sigmaP: sigma^2/P 
%           D0: Radius of the circular area (m)
%           N: Number of IoT devices in the area
%           R: Spectral efficiency (bps/Hz)
%           sch: scheme 'HA' (default) or 'NOHA'.
% Outputs:
%           ThroughputTheo: Analytical throughput 
%           ThroughputSimul: Simulated throughput (Monte Carlo)
%           OptimalKappa: Optimal difficulty that maximizes the throughput
%           MaximumThroughput: Throughput achieved by OptimalKappa
%==========================================================================
monteCarloSetups = 1e4;
alpha = 4;                  % path-loss exponent
%R = 1;                      % transmission rate (bps/Hz)

%% MONTE CARLO
    hashIoT = hash_sim(N, monteCarloSetups);
    r = checkhash(hashIoT,diff);
    HA = sum(r, 1);

    Probability_NHA_1 = sum(HA==1)/monteCarloSetups;    % Pr{NHA = 1}, Eq. (13)
    Probability_NHA_2 = sum(HA==2)/monteCarloSetups;    % Pr{NHA = 2}, Eq. (13)
   
    SNR_IoT=[];
    for ii = 1:monteCarloSetups
        Raio = D0*sqrt(rand(1, 2));
        Theta = 2*pi*rand(1, 2);
        positionsX = Raio.*cos(Theta);
        positionsY = Raio.*sin(Theta);

        Distances=(sqrt(positionsX.^2+positionsY.^2));
        h = (abs(sqrt(0.5)*(randn(1, 2) +1i*randn(1, 2)))).^2;  
        SNR_MonteCarlo=(h(1))./(sigmaP.*(Distances(1).^(alpha)));
        SNR_MonteCarlo(2)=(h(2))./(sigmaP.*(Distances(2).^(alpha)));

        SNR_IoT = [SNR_IoT SNR_MonteCarlo'];
    end
    OutageMC=(sum((log2(1+SNR_IoT(1, :))<=R)~=0)/monteCarloSetups);
    SNR_NOMA = sort(SNR_IoT, 'ascend');
    
    C_SINR = log2(1+SNR_NOMA(2, :)./(1+SNR_NOMA(1, :)));
    C_SNR = log2(1+SNR_NOMA(1, :));
    
    %OutageSIC1=(sum(C_SINR<R~=0)/monteCarloSetups);
    %OutageSIC2=(sum(C_SNR<R~=0)/monteCarloSetups);
    OutageSIC2not1=(sum((C_SINR>=R).*(C_SNR<R))/monteCarloSetups);
    OutageSICNone=(sum((C_SINR>=R).*(C_SNR>=R))/monteCarloSetups);
   
     ThroughputSimul = R*(1-OutageMC)*Probability_NHA_1;  % Throughput HA, Eq. (10)
     if(strcmp(sch,'NOHA') == 1)
        ThroughputSimul= R*(1-OutageMC)*Probability_NHA_1+... % Throughput NOHA, Eq. (18)
                     R*OutageSIC2not1*Probability_NHA_2+...
                     2*R*OutageSICNone*Probability_NHA_2;
     end
       
    
%% ANALYSIS
    Prob = @(n,k) nchoosek(N,n)*((1-k).^(n)).*(k).^(N-n);

    Probability_NHA_1_Theo = Prob(1,diff);  % Pr{NHA = 1}, Eq. (13)
    Probability_NHA_2_Theo = Prob(2,diff);  % Pr{NHA = 2}, Eq. (13)
    
    gamma0 = 2^R-1;
    Outage = 1 - (sqrt(pi)*erf(D0^2*sqrt(gamma0*sigmaP)))/(2*D0^2*sqrt(gamma0*sigmaP));
    ThroughputTheo = R.*(1-Outage).*Probability_NHA_1_Theo;  
    OptimalKappa=(N-1)/(N);
    MaximumThroughput = R*Prob(1,OptimalKappa)*(1-Outage);

    if(strcmp(sch,'NOHA') == 1)   
        OutageSIC2Theo = 2*Outage-(Outage)^2;
        L = 5;
        Somatorio = 0; 
        for(nn = 1:L)
            Thetan = cos(((2*nn-1)/(2*L))*pi);
            cn = 1+((D0/2)*Thetan+(D0/2))^alpha;
            for(kk = 1:L)
                Thetak = cos(((2*kk-1)/(2*L))*pi);
                ck = 1+((D0/2)*Thetak+(D0/2))^alpha;

                Somatorio=Somatorio+(cn*(1/(cn + ck) - exp(-ck*gamma0*sigmaP)/(cn+gamma0*ck))*(pi^2)*(1 + Thetan)*...
                                    sqrt(1 - Thetan^2)*(1 + Thetak)*sqrt(1 - Thetak^2))/(2*L^2);

          end
        end  
        OutageSIC1Theo = Somatorio;

        ThroughputTheo = R.*(1-Outage).*Probability_NHA_1_Theo+...
                         R.*(1-OutageSIC1Theo).*(OutageSIC2Theo).*Probability_NHA_2_Theo+...
                         2.*R.*(1-OutageSIC1Theo).*(1-OutageSIC2Theo).*Probability_NHA_2_Theo;
        A = R*N*(1-Outage);                                             % Eq. (21)
        B = R*nchoosek(N,2)*(1-OutageSIC1Theo)*(2-OutageSIC2Theo);      % Eq. (22)
        OptimalKappa = (1/(2*N*(A-B)))*((N-1)*(A-2*B)...                % Eq. (20)
                    +sqrt(A^2*(N-1)^2 + 4*B*(B-A)));
        MaximumThroughput = (A/N)*Prob(1,OptimalKappa) + (B/nchoosek(N,2))*Prob(2,OptimalKappa);
    
    end
    
end