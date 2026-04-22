function [data,SAMPLES] = SVGD(data,SAMPLES,target,setting,MAX_ITERS,nSamples)

stack = target.stack;
cal = target.calibration;

STACK_SIG = stack.SIG.^(-1);

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
eta = 1e-3;

A1 = setting.a_14C;
B1 = setting.b_14C;

A2 = setting.a_d18O;
B2 = setting.b_d18O;

% eps = 1e-2;

THSD = 1e-4;

YY = data.d18O;
YC = data.C14;
YA = data.ETC;

if ~isempty(YY)
    ID_Y = YY(:,end);
    YY = YY(:,1:end-1);
else
    ID_Y = [];
end
if ~isempty(YC)
    ID_C = YC(:,end);
    YC = YC(:,1:end-1);
else
    ID_C = [];
end
if ~isempty(YA)
    ID_A = YA(:,end);
    YA = YA(:,1:end-1);
else
    ID_A = [];
end


D1 = data.depth.*data.R;
N1 = size(D1,1);
D1 = sort(D1,'ascend');

N0 = size(SAMPLES.F,1);

INT = setting.interval*data.R;


AP = 0:INT/nSamples:(INT-INT/nSamples);

CNT = zeros(N1,1);
HH = zeros(N1,N0);
for n = 1:N1
    CNT(n) = floor((D1(n)-setting.start_depth)/INT);
    HH(n,1:CNT(n)) = INT;
end
HH0 = - AP;
HH1 = min(D1-setting.start_depth-CNT.*INT+AP,INT);
HH2 = max(D1-setting.start_depth-CNT.*INT+AP-INT,0);

UQ = unique(CNT);
QQ = zeros(length(UQ),N1);
for n = 1:N1
    QQ(UQ==CNT(n),n) = QQ(UQ==CNT(n),n) + 1;
end

UQ_Y = unique(ID_Y);
QQ_Y = zeros(size(UQ_Y,1),size(ID_Y,1));
for n = 1:size(ID_Y,1)
    QQ_Y(UQ_Y==ID_Y(n),n) = QQ_Y(UQ_Y==ID_Y(n),n) + 1;
end

UQ_C = unique(ID_C);
QQ_C = zeros(size(UQ_C,1),size(ID_C,1));
for n = 1:size(ID_C,1)
    QQ_C(UQ_C==ID_C(n),n) = QQ_C(UQ_C==ID_C(n),n) + 1;
end

UQ_A = unique(ID_A);
QQ_A = zeros(size(UQ_A,1),size(ID_A,1));
for n = 1:size(ID_A,1)
    QQ_A(UQ_A==ID_A(n),n) = QQ_A(UQ_A==ID_A(n),n) + 1;
end


BETA = setting.beta;
RHO = setting.rho;

RHO = sqrt(RHO^INT);


N0 = size(SAMPLES.F,1);
SAM = ones(N0,nSamples,2.*setting.alpha)./sqrt(2.*setting.alpha);
SAM = SAM + normrnd(0,0.01,[N0,nSamples,2.*setting.alpha]);

SAM_BIAS = ones(1,nSamples) + normrnd(0,0.1,[1,nSamples]);


Mw = 0;
Vw = 0;

Mw_BIAS = 0;
Vw_BIAS = 0;


for rr = 1:MAX_ITERS


    SAM_T = pagetranspose(SAM);

    KK = sum(SAM_T.^2,2) + pagetranspose(sum(SAM_T.^2,2)) - 2.*pagemtimes(SAM_T,SAM);
    KK_BIAS = (SAM_BIAS'-SAM_BIAS).^2;

    hh = median(sum(KK,3)+KK_BIAS,"all")./log(nSamples);

    KK = exp(-KK./hh);
    PK = (2./hh.*SAM).*sum(KK,1) - pagemtimes(2./hh.*SAM,KK);

    KK_BIAS = exp(-KK_BIAS./hh);
    PK_BIAS = (2./hh.*SAM_BIAS).*sum(KK_BIAS,1) - (2./hh.*SAM_BIAS)*KK_BIAS;


    % Initialization:
    PDEV = zeros(N0,nSamples,2*setting.alpha);
    PDEV_BIAS = zeros(1,nSamples);

    % Transition Model:
    SAM_G = sum(SAM.^2,3);

    PDEV(1:end-1,:,:) = PDEV(1:end-1,:,:) - 2.*BETA./(1-RHO.^2).*RHO.*(RHO.*SAM(1:end-1,:,:)-SAM(2:end,:,:));
    PDEV(2:end,:,:) = PDEV(2:end,:,:) + 2.*BETA./(1-RHO.^2).*(RHO.*SAM(1:end-1,:,:)-SAM(2:end,:,:));
    PDEV(1,:,:) = PDEV(1,:,:) - 2.*BETA.*SAM(1,:,:);

    % Emission Model:
    AA = HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:) + SAM_BIAS.^2 + setting.min_start_age; % ages

    AA_Y = AA(ID_Y,:);
    AA_C = AA(ID_C,:);
    AA_A = AA(ID_A,:);

    % d18O:
    if ~isempty(AA_Y)
        MU = interp1(stack.A,stack.MU,AA_Y,'linear','extrap');
        SIG = interp1(stack.A,STACK_SIG,AA_Y,'linear',THSD);

        PMU = interp1(stack.A,stack.PMU,AA_Y,'linear','extrap');
        PSIG = interp1(stack.A,stack.PSIG,AA_Y,'linear','extrap');

        ZZ = (YY-data.scale.*MU-data.shift)./data.scale;
        if isinf(A2)&&isinf(B2)
            ALPHA = ZZ.*SIG.^2.*PMU + ZZ.^2.*SIG.^3.*PSIG - SIG.*PSIG;
        else
            ALPHA = (2.*A2+1).*(ZZ.*SIG.^2.*PMU+ZZ.^2.*SIG.^3.*PSIG)./(2.*B2+ZZ.^2.*SIG.^2) - SIG.*PSIG;
        end

        ALPHA = QQ_Y*ALPHA;

        DELTA = zeros(N1,nSamples);
        DELTA(UQ_Y,:) = ALPHA;

        PDEV = PDEV + (HH'*DELTA).*(2.*SAM);
        PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(DELTA,1).*(2.*SAM(1,:,:));
        PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*DELTA)).*(2.*SAM(UQ+1,:,:));
        PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*DELTA)).*(2.*SAM(UQ+2,:,:));

        PDEV_BIAS = PDEV_BIAS + sum(DELTA,1).*(2.*SAM_BIAS);
    end


    % 14C:
    if ~isempty(AA_C)

        MU = zeros(size(AA_C));
        SIG = zeros(size(AA_C));

        PMU = zeros(size(AA_C));
        PSIG = zeros(size(AA_C));

        for k = 1:length(cal)
            ID = (YC(:,end)==k);

            if sum(ID) > 0
                MU(ID,:) = interp1(cal(k).A,cal(k).MU,AA_C(ID,:),'linear','extrap');
                SIG(ID,:) = interp1(cal(k).A,cal(k).SIG,AA_C(ID,:),'linear',1./THSD);

                PMU(ID,:) = interp1(cal(k).A,cal(k).PMU,AA_C(ID,:),'linear','extrap');
                PSIG(ID,:) = interp1(cal(k).A,cal(k).PSIG,AA_C(ID,:),'linear','extrap');
            end
        end

        WW = sqrt(SIG.^2+YC(:,2).^2+YC(:,4).^2);
        ZZ = (MU+YC(:,3)-YC(:,1))./WW;

        if isinf(A1)&&isinf(B1)
            ALPHA = - ZZ./WW.*PMU + SIG.*(ZZ./WW).^2.*PSIG - SIG./(WW.^2).*PSIG;
        else
            ALPHA = - (0.5+A1)./(B1+0.5.*ZZ.^2).*(ZZ./WW.*PMU-SIG.*(ZZ./WW).^2.*PSIG) - SIG./(WW.^2).*PSIG;
        end

        ALPHA = QQ_C*ALPHA;

        DELTA = zeros(N1,nSamples);
        DELTA(UQ_C,:) = ALPHA;

        PDEV = PDEV + (HH'*DELTA).*(2.*SAM);
        PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(DELTA,1).*(2.*SAM(1,:,:));
        PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*DELTA)).*(2.*SAM(UQ+1,:,:));
        PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*DELTA)).*(2.*SAM(UQ+2,:,:));

        PDEV_BIAS = PDEV_BIAS + sum(DELTA,1).*(2.*SAM_BIAS);
    end


    % ADD_AGE:
    if ~isempty(AA_A)
        ALPHA = - (AA_A-YA(:,1)).*YA(:,2).^(-2);

        ALPHA = QQ_A*ALPHA;

        DELTA = zeros(N1,nSamples);
        DELTA(UQ_A,:) = ALPHA;

        PDEV = PDEV + (HH'*DELTA).*(2.*SAM);
        PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(DELTA,1).*(2.*SAM(1,:,:));
        PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*DELTA)).*(2.*SAM(UQ+1,:,:));
        PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*DELTA)).*(2.*SAM(UQ+2,:,:));

        PDEV_BIAS = PDEV_BIAS + sum(DELTA,1).*(2.*SAM_BIAS);
    end


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


SAMPLES.F = SAM.*sqrt(data.R);
SAM_G =  sum(SAM.^2,3);
SAMPLES.AGE = SAM_BIAS.^2 + setting.min_start_age + HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:);
SAMPLES.BIAS = SAM_BIAS;


SAMPLES.AGE_D = zeros(N0-1,nSamples);
D0 = SAMPLES.D0.*data.R;



CNT = zeros(N0-1,1);
HH = zeros(N0-1,N0);
for n = 1:N0-1
    CNT(n) = floor(D0(n)/INT);
    HH(n,1:CNT(n)) = INT;
end
HH0 = - AP;
HH1 = min(D0-CNT.*INT+AP,INT);
HH2 = max(D0-CNT.*INT+AP-INT,0);

SAMPLES.AGE_D = SAM_BIAS.^2 + setting.min_start_age + HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:);
SAMPLES.ACC_RATE = SAM_G.*data.R;


end