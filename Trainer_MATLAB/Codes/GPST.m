function [data,SAMPLES] = GPST(data,param,SAMPLES,cal,setting,MAX_ITERS,nSamples)

A1 = setting.a;
B1 = setting.b;

eps = 1e-3;
M = 30;

THSD = 1e-2;

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
SAM_G =  sum(SAM.^2,3);
SAM_BIAS = SAMPLES.BIAS;

if size(SAM,2) ~= nSamples
    RAND_SEED = ceil(size(SAM,2).*rand(nSamples,1));

    SAM = SAM(:,RAND_SEED,:);
    SAM_G =  sum(SAM.^2,3);
    SAM_BIAS = SAM_BIAS(:,RAND_SEED);
end

AA_C = SAM_BIAS.^2 - 0.1 + HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:);

PHI = normrnd(0,1,[N0,nSamples,setting.K]);
PHI_BIAS = normrnd(0,1,[1,nSamples]);

LOGLIK_old = sum(sum(-0.5*PHI.^2-0.5*log(2*pi),3),1);
LOGLIK_old = LOGLIK_old - 0.5*PHI_BIAS.^2 - 0.5*log(2*pi);

LOGLIK_old = LOGLIK_old + sum(-BETA.*SAM(end,:,:).^2-BETA./(1-RHO.^2).*sum((SAM(1:end-1,:,:)-RHO.*SAM(2:end,:,:)).^2,1),3);


% 14C:
MU = interp1(cal.A,cal.MU,AA_C,'linear','extrap');
SIG = interp1(cal.A,cal.SIG,AA_C,'linear',1./THSD);

WW = sqrt(SIG.^2+YC(:,2).^2+YC(:,4).^2);
ZZ = (MU+YC(:,3)-YC(:,1))./WW;

if isinf(A1)&&isinf(B1)
    LOGLIK_old = LOGLIK_old + sum(-0.5.*ZZ.^2-log(WW),1);
else
    LOGLIK_old = LOGLIK_old + sum(-(0.5+A1).*log(1+ZZ.^2./(2.*B1))-log(WW),1);
end

rand_seed = log(rand(MAX_ITERS,nSamples));


APP_RATE = zeros(MAX_ITERS,1);


for rr = 1:MAX_ITERS

    SAM_G =  sum(SAM.^2,3);

    AA_C = SAM_BIAS.^2 - 0.1 + HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:);
   
    PDEV = zeros(N0,nSamples,setting.K);
    PDEV_BIAS = zeros(1,nSamples);
    
    PDEV(1:end-1,:,:) = PDEV(1:end-1,:,:) - 2.*BETA./(1-RHO.^2).*RHO.*(RHO.*SAM(1:end-1,:,:)-SAM(2:end,:,:));
    PDEV(2:end,:,:) = PDEV(2:end,:,:) + 2.*BETA./(1-RHO.^2).*(RHO.*SAM(1:end-1,:,:)-SAM(2:end,:,:));
    PDEV(1,:,:) = PDEV(1,:,:) - 2.*BETA.*SAM(1,:,:);


    % 14C:
    MU = interp1(cal.A,cal.MU,AA_C,'linear','extrap');
    SIG = interp1(cal.A,cal.SIG,AA_C,'linear',1./THSD);

    PMU = interp1(cal.A,cal.PMU,AA_C,'linear','extrap');
    PSIG = interp1(cal.A,cal.PSIG,AA_C,'linear','extrap');

    WW = sqrt(SIG.^2+YC(:,2).^2+YC(:,4).^2);
    ZZ = (MU+YC(:,3)-YC(:,1))./WW;

    if isinf(A1)&&isinf(B1)
        ALPHA = - ZZ./WW.*PMU + SIG.*(ZZ./WW).^2.*PSIG - SIG./(WW.^2).*PSIG;
    else
        ALPHA = - (0.5+A1)./(B1+0.5.*ZZ.^2).*(ZZ./WW.*PMU-SIG.*(ZZ./WW).^2.*PSIG) - SIG./(WW.^2).*PSIG;
    end

    PDEV = PDEV + (HH'*ALPHA).*(2.*SAM);
    PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(ALPHA,1).*(2.*SAM(1,:,:));
    PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*ALPHA)).*(2.*SAM(UQ+1,:,:));
    PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*ALPHA)).*(2.*SAM(UQ+2,:,:));

    PDEV_BIAS = PDEV_BIAS + sum(ALPHA,1).*(2.*SAM_BIAS);

    PHI_new = PHI + 0.5*eps*PDEV;
    PHI_BIAS_new = PHI_BIAS + 0.5*eps*PDEV_BIAS;
    SAM_new = SAM;
    SAM_BIAS_new = SAM_BIAS;

    for m = 1:M-1
        SAM_new = SAM_new + eps*PHI_new;
        SAM_BIAS_new = SAM_BIAS_new + eps*PHI_BIAS_new;

        SAM_G_new =  sum(SAM_new.^2,3);

        AA_C = SAM_BIAS_new.^2 - 0.1 + HH*SAM_G_new;
        AA_C = AA_C + HH0.*SAM_G_new(1,:) + HH1.*SAM_G_new(CNT+1,:) + HH2.*SAM_G_new(CNT+2,:);

        PDEV = zeros(N0,nSamples,setting.K);
        PDEV_BIAS = zeros(1,nSamples);
        
        PDEV(1:end-1,:,:) = PDEV(1:end-1,:,:) - 2.*BETA./(1-RHO.^2).*RHO.*(RHO.*SAM_new(1:end-1,:,:)-SAM_new(2:end,:,:));
        PDEV(2:end,:,:) = PDEV(2:end,:,:) + 2.*BETA./(1-RHO.^2).*(RHO.*SAM_new(1:end-1,:,:)-SAM_new(2:end,:,:));
        PDEV(1,:,:) = PDEV(1,:,:) - 2.*BETA.*SAM_new(1,:,:);

        % 14C:
        MU = interp1(cal.A,cal.MU,AA_C,'linear','extrap');
        SIG = interp1(cal.A,cal.SIG,AA_C,'linear',1./THSD);

        PMU = interp1(cal.A,cal.PMU,AA_C,'linear','extrap');
        PSIG = interp1(cal.A,cal.PSIG,AA_C,'linear','extrap');

        WW = sqrt(SIG.^2+YC(:,2).^2+YC(:,4).^2);
        ZZ = (MU+YC(:,3)-YC(:,1))./WW;

        if isinf(A1)&&isinf(B1)
            ALPHA = - ZZ./WW.*PMU + SIG.*(ZZ./WW).^2.*PSIG - SIG./(WW.^2).*PSIG;
        else
            ALPHA = - (0.5+A1)./(B1+0.5.*ZZ.^2).*(ZZ./WW.*PMU-SIG.*(ZZ./WW).^2.*PSIG) - SIG./(WW.^2).*PSIG;
        end

        PDEV = PDEV + (HH'*ALPHA).*(2.*SAM_new);
        PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(ALPHA,1).*(2.*SAM_new(1,:,:));
        PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*ALPHA)).*(2.*SAM_new(UQ+1,:,:));
        PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*ALPHA)).*(2.*SAM_new(UQ+2,:,:));
        
        PDEV_BIAS = PDEV_BIAS + sum(ALPHA,1).*(2.*SAM_BIAS_new);

        PHI_new = PHI_new + eps*PDEV;
        PHI_BIAS_new = PHI_BIAS_new + eps*PDEV_BIAS;
    end

    SAM_new = SAM_new + eps*PHI_new;
    SAM_BIAS_new = SAM_BIAS_new + eps*PHI_BIAS_new;

    SAM_G_new =  sum(SAM_new.^2,3);

    AA_C = SAM_BIAS_new.^2 - 0.1 + HH*SAM_G_new + HH0.*SAM_G_new(1,:) + HH1.*SAM_G_new(CNT+1,:) + HH2.*SAM_G_new(CNT+2,:);

    PDEV = zeros(N0,nSamples,setting.K);
    PDEV_BIAS = zeros(1,nSamples);

    PDEV(1:end-1,:,:) = PDEV(1:end-1,:,:) - 2.*BETA./(1-RHO.^2).*RHO.*(RHO.*SAM_new(1:end-1,:,:)-SAM_new(2:end,:,:));
    PDEV(2:end,:,:) = PDEV(2:end,:,:) + 2.*BETA./(1-RHO.^2).*(RHO.*SAM_new(1:end-1,:,:)-SAM_new(2:end,:,:));
    PDEV(1,:,:) = PDEV(1,:,:) - 2.*BETA.*SAM_new(1,:,:);

    % 14C:
    MU = interp1(cal.A,cal.MU,AA_C,'linear','extrap');
    SIG = interp1(cal.A,cal.SIG,AA_C,'linear',1./THSD);

    PMU = interp1(cal.A,cal.PMU,AA_C,'linear','extrap');
    PSIG = interp1(cal.A,cal.PSIG,AA_C,'linear','extrap');

    WW = sqrt(SIG.^2+YC(:,2).^2+YC(:,4).^2);
    ZZ = (MU+YC(:,3)-YC(:,1))./WW;

    if isinf(A1)&&isinf(B1)
        ALPHA = - ZZ./WW.*PMU + SIG.*(ZZ./WW).^2.*PSIG - SIG./(WW.^2).*PSIG;
    else
        ALPHA = - (0.5+A1)./(B1+0.5.*ZZ.^2).*(ZZ./WW.*PMU-SIG.*(ZZ./WW).^2.*PSIG) - SIG./(WW.^2).*PSIG;
    end

    PDEV = PDEV + (HH'*ALPHA).*(2.*SAM_new);
    PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(ALPHA,1).*(2.*SAM_new(1,:,:));
    PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*ALPHA)).*(2.*SAM_new(UQ+1,:,:));
    PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*ALPHA)).*(2.*SAM_new(UQ+2,:,:));

    PDEV_BIAS = PDEV_BIAS + sum(ALPHA,1).*(2.*SAM_BIAS_new);

    PHI_new = PHI_new + 0.5*eps*PDEV;
    PHI_BIAS_new = PHI_BIAS_new + 0.5*eps*PDEV_BIAS;

    AA_C = SAM_BIAS_new.^2 - 0.1 + HH*SAM_G_new + HH0.*SAM_G_new(1,:) + HH1.*SAM_G_new(CNT+1,:) + HH2.*SAM_G_new(CNT+2,:);

    LOGLIK_new = sum(sum(-0.5*PHI_new.^2-0.5*log(2*pi),3),1);
    LOGLIK_new = LOGLIK_new - 0.5*PHI_BIAS_new.^2 - 0.5*log(2*pi);

    LOGLIK_new = LOGLIK_new + sum(-BETA.*SAM_new(end,:,:).^2-BETA./(1-RHO.^2).*sum((SAM_new(1:end-1,:,:)-RHO.*SAM_new(2:end,:,:)).^2,1),3);

    % 14C:
    MU = interp1(cal.A,cal.MU,AA_C,'linear','extrap');
    SIG = interp1(cal.A,cal.SIG,AA_C,'linear',1./THSD);

    WW = sqrt(SIG.^2+YC(:,2).^2+YC(:,4).^2);
    ZZ = (MU+YC(:,3)-YC(:,1))./WW;

    if isinf(A1)&&isinf(B1)
        LOGLIK_new = LOGLIK_new + sum(-0.5.*ZZ.^2-log(WW),1);
    else
        LOGLIK_new = LOGLIK_new + sum(-(0.5+A1).*log(1+ZZ.^2./(2.*B1))-log(WW),1);
    end

    
    ID = (~isnan(LOGLIK_new))&(LOGLIK_new-LOGLIK_old>rand_seed(rr,:));

    SAM(:,ID,:) = SAM_new(:,ID,:);
    SAM_BIAS(ID) = SAM_BIAS_new(ID);

    PHI(:,ID,:) = PHI_new(:,ID,:);
    PHI_BIAS(:,ID) = PHI_BIAS_new(:,ID);
    LOGLIK_old(ID) = LOGLIK_new(ID);

    LOGLIK_old = LOGLIK_old - sum(sum(-0.5*PHI.^2-0.5*log(2*pi),3),1);
    LOGLIK_old = LOGLIK_old + 0.5*PHI_BIAS.^2 + 0.5*log(2*pi);
    PHI = normrnd(0,1,[N0,nSamples,setting.K]);
    PHI_BIAS = normrnd(0,1,[1,nSamples]);
    LOGLIK_old = LOGLIK_old + sum(sum(-0.5*PHI.^2-0.5*log(2*pi),3),1);
    LOGLIK_old = LOGLIK_old - 0.5*PHI_BIAS.^2 - 0.5*log(2*pi);

    APP_RATE(rr) = 100.*sum(ID)./nSamples;

    if APP_RATE(rr) < 90
        eps = eps./1.01;
    elseif APP_RATE(rr) > 95
        eps = eps.*1.01;
    end
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


SAMPLES.APP_RATE = APP_RATE;
SAMPLES.ACC_RATE = SAM_G;


end