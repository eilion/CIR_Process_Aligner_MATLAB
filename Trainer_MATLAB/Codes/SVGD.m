function [data,SAMPLES] = SVGD(data,param,SAMPLES,cal,setting,MAX_ITERS,nSamples)

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
eta = 1e-3;

A1 = setting.a;
B1 = setting.b;

% eps = 1e-2;

THSD = 1e-4;

YC = data.C14;

D1 = data.depth.*data.R;
N1 = size(D1,1);
D1 = sort(D1,'ascend');

N0 = size(SAMPLES.D,1);


AP = 0:setting.interval/nSamples:(setting.interval-setting.interval/nSamples);

CNT = zeros(N1,1);
HH = zeros(N1,N0);
for n = 1:N1
    CNT(n) = floor(D1(n)/setting.interval);
    HH(n,1:CNT(n)) = setting.interval;
end
HH0 = - AP;
HH1 = min(D1-CNT.*setting.interval+AP,setting.interval);
HH2 = max(D1-CNT.*setting.interval+AP-setting.interval,0);

UQ = unique(CNT);
QQ = zeros(length(UQ),N1);
for n = 1:N1
    QQ(UQ==CNT(n),n) = QQ(UQ==CNT(n),n) + 1;
end


BETA = param.BETA;
RHO = param.RHO;

RHO = sqrt(RHO^setting.interval);


SAM = SAMPLES.F;
SAM_BIAS = SAMPLES.BIAS;

if size(SAM,2) ~= nSamples
    RAND_SEED = ceil(size(SAM,2).*rand(nSamples,1));

    SAM = SAM(:,RAND_SEED,:);
    SAM_BIAS = SAM_BIAS(:,RAND_SEED);
end


Mw = 0;
Vw = 0;

Mw_BIAS = 0;
Vw_BIAS = 0;


for rr = 1:MAX_ITERS


    SAM_T = pagetranspose(SAM);

    KK = sum(SAM_T.^2,2) + pagetranspose(sum(SAM_T.^2,2)) - 2.*pagemtimes(SAM_T,SAM);
    KK_BIAS = (SAM_BIAS'-SAM_BIAS).^2;

    hh = median(sqrt(sum(KK,3)+KK_BIAS),"all").^2./log(nSamples);

    KK = exp(-KK./hh);
    PK = (2./hh.*SAM).*sum(KK,1) - pagemtimes(2./hh.*SAM,KK);

    KK_BIAS = exp(-KK_BIAS./hh);
    PK_BIAS = (2./hh.*SAM_BIAS).*sum(KK_BIAS,1) - (2./hh.*SAM_BIAS)*KK_BIAS;


    % Initialization:
    PDEV = zeros(N0,nSamples,setting.K);
    PDEV_BIAS = zeros(1,nSamples);

    % Transition Model:
    SAM_G = sum(SAM.^2,3);

    PDEV(1:end-1,:,:) = PDEV(1:end-1,:,:) - 2.*BETA./(1-RHO.^2).*RHO.*(RHO.*SAM(1:end-1,:,:)-SAM(2:end,:,:));
    PDEV(2:end,:,:) = PDEV(2:end,:,:) + 2.*BETA./(1-RHO.^2).*(RHO.*SAM(1:end-1,:,:)-SAM(2:end,:,:));
    PDEV(1,:,:) = PDEV(1,:,:) - 2.*BETA.*SAM(1,:,:);

    % Emission Model:
    AA_C = HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:) + SAM_BIAS.^2 - 0.1; % ages

    % 14C:
    MU = interp1(cal.A,cal.MU,AA_C,'linear','extrap');
    SIG = interp1(cal.A,cal.SIG,AA_C,'linear',1./THSD);

    PMU = interp1(cal.A,cal.PMU,AA_C,'linear','extrap');
    PSIG = interp1(cal.A,cal.PSIG,AA_C,'linear','extrap');

    WW = sqrt(SIG.^2+YC(:,2).^2+YC(:,4).^2);
    ZZ = (MU+YC(:,3)-YC(:,1))./WW;

    if isinf(A1)&&isinf(B1)
        DELTA = - ZZ./WW.*PMU + SIG.*(ZZ./WW).^2.*PSIG - SIG./(WW.^2).*PSIG;
    else
        DELTA = - (0.5+A1)./(B1+0.5.*ZZ.^2).*(ZZ./WW.*PMU-SIG.*(ZZ./WW).^2.*PSIG) - SIG./(WW.^2).*PSIG;
    end

    
    PDEV = PDEV + (HH'*DELTA).*(2.*SAM);
    PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(DELTA,1).*(2.*SAM(1,:,:));
    PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*DELTA)).*(2.*SAM(UQ+1,:,:));
    PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*DELTA)).*(2.*SAM(UQ+2,:,:));
  
    PDEV_BIAS = PDEV_BIAS + sum(DELTA,1).*(2.*SAM_BIAS);


    PHI = pagemtimes(PDEV,KK) + PK;
    PHI_BIAS = PDEV_BIAS*KK_BIAS + PK_BIAS;

    PHI = PHI./nSamples;
    PHI_BIAS = PHI_BIAS./nSamples;


    % Update:
    Mw = beta1.*Mw - (1-beta1).*PHI;
    Vw = beta2.*Vw + (1-beta2).*PHI.*PHI;

    Mw_BIAS = beta1.*Mw_BIAS - (1-beta1).*PHI_BIAS;
    Vw_BIAS = beta2.*Vw_BIAS + (1-beta2).*PHI_BIAS.*PHI_BIAS;

    SAM = SAM - (eta*sqrt(1-beta2^rr)/(1-beta1^rr)).*Mw./(sqrt(Vw)+epsilon);
    SAM_BIAS = SAM_BIAS - (eta*sqrt(1-beta2^rr)/(1-beta1^rr)).*Mw_BIAS./(sqrt(Vw_BIAS)+epsilon);
end


SAMPLES.F = SAM;
SAM_G =  sum(SAM.^2,3);
SAMPLES.AGE = SAM_BIAS.^2 - 0.1 + HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:);
SAMPLES.BIAS = SAM_BIAS;


SAMPLES.AGE_D = zeros(N0-1,nSamples);
D0 = SAMPLES.D0;



CNT = zeros(N0-1,1);
HH = zeros(N0-1,N0);
for n = 1:N0-1
    CNT(n) = floor(D0(n)/setting.interval);
    HH(n,1:CNT(n)) = setting.interval;
end
HH0 = - AP;
HH1 = min(D0-CNT.*setting.interval+AP,setting.interval);
HH2 = max(D0-CNT.*setting.interval+AP-setting.interval,0);

SAMPLES.AGE_D = SAM_BIAS.^2 - 0.1 + HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:);

SAMPLES.ACC_RATE = SAM_G;


end