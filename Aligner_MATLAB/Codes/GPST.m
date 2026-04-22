function [data,SAMPLES] = GPST(data,SAMPLES,target,setting,MAX_ITERS,nSamples)

stack = target.stack;
cal = target.calibration;

STACK_SIG = stack.SIG.^(-1);

A1 = setting.a_14C;
B1 = setting.b_14C;

A2 = setting.a_d18O;
B2 = setting.b_d18O;

eps = 1e-3;
M = 30;

THSD = 1e-2;

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



SAM = SAMPLES.F./sqrt(data.R);
SAM_G =  sum(SAM.^2,3);
SAM_BIAS = SAMPLES.BIAS;

if size(SAM,2) ~= nSamples
    RAND_SEED = ceil(size(SAM,2).*rand(nSamples,1));

    SAM = SAM(:,RAND_SEED,:);
    SAM_G =  sum(SAM.^2,3);
    SAM_BIAS = SAM_BIAS(:,RAND_SEED);
end

AA = SAM_BIAS.^2 + setting.min_start_age + HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:);
AA_Y = AA(ID_Y,:);
AA_C = AA(ID_C,:);
AA_A = AA(ID_A,:);

PHI = normrnd(0,1,[N0,nSamples,2.*setting.alpha]);
PHI_BIAS = normrnd(0,1,[1,nSamples]);

LOGLIK_old = sum(sum(-0.5*PHI.^2-0.5*log(2*pi),3),1);
LOGLIK_old = LOGLIK_old - 0.5*PHI_BIAS.^2 - 0.5*log(2*pi);

LOGLIK_old = LOGLIK_old + sum(-BETA.*SAM(end,:,:).^2-BETA./(1-RHO.^2).*sum((SAM(1:end-1,:,:)-RHO.*SAM(2:end,:,:)).^2,1),3);

% d18O:
if ~isempty(AA_Y)
    MU = interp1(stack.A,stack.MU,AA_Y,'linear','extrap');
    SIG = interp1(stack.A,STACK_SIG,AA_Y,'linear',THSD);

    ZZ = (YY-data.scale.*MU-data.shift)./data.scale;

    if isinf(A2)&&isinf(B2)
        LOGLIK_old = LOGLIK_old + sum(-0.5.*ZZ.^2.*SIG.^2+log*(SIG)-log(data.scale),1);
    else
        LOGLIK_old = LOGLIK_old + sum(-(A2+1./2).*log(1+1./(2.*B2).*ZZ.^2.*SIG.^2)+log(SIG)-log(data.scale),1);
    end
end

% 14C:
if ~isempty(AA_C)
    MU = zeros(size(AA_C));
    SIG = zeros(size(AA_C));

    for k = 1:length(cal)
        ID = (YC(:,end)==k);

        if sum(ID) > 0
            MU(ID,:) = interp1(cal(k).A,cal(k).MU,AA_C(ID,:),'linear','extrap');
            SIG(ID,:) = interp1(cal(k).A,cal(k).SIG,AA_C(ID,:),'linear',1./THSD);
        end
    end

    WW = sqrt(SIG.^2+YC(:,2).^2+YC(:,4).^2);
    ZZ = (MU+YC(:,3)-YC(:,1))./WW;

    if isinf(A1)&&isinf(B1)
        LOGLIK_old = LOGLIK_old + sum(-0.5.*ZZ.^2-log(WW),1);
    else
        LOGLIK_old = LOGLIK_old + sum(-(0.5+A1).*log(1+ZZ.^2./(2.*B1))-log(WW),1);
    end
end

% ADD_AGE:
if ~isempty(AA_A)
    LOGLIK_old = LOGLIK_old + sum(-0.5.*(AA_A-YA(:,1)).^2.*YA(:,2).^(-2),1);
end


rand_seed = log(rand(MAX_ITERS,nSamples));


APP_RATE = zeros(MAX_ITERS,1);


for rr = 1:MAX_ITERS

    SAM_G =  sum(SAM.^2,3);

    AA = SAM_BIAS.^2 + setting.min_start_age + HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:);
    AA_Y = AA(ID_Y,:);
    AA_C = AA(ID_C,:);
    AA_A = AA(ID_A,:);

    PDEV = zeros(N0,nSamples,2.*setting.alpha);
    PDEV_BIAS = zeros(1,nSamples);

    PDEV(1:end-1,:,:) = PDEV(1:end-1,:,:) - 2.*BETA./(1-RHO.^2).*RHO.*(RHO.*SAM(1:end-1,:,:)-SAM(2:end,:,:));
    PDEV(2:end,:,:) = PDEV(2:end,:,:) + 2.*BETA./(1-RHO.^2).*(RHO.*SAM(1:end-1,:,:)-SAM(2:end,:,:));
    PDEV(1,:,:) = PDEV(1,:,:) - 2.*BETA.*SAM(1,:,:);


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


    PHI_new = PHI + 0.5*eps*PDEV;
    PHI_BIAS_new = PHI_BIAS + 0.5*eps*PDEV_BIAS;
    SAM_new = SAM;
    SAM_BIAS_new = SAM_BIAS;

    for m = 1:M-1
        SAM_new = SAM_new + eps*PHI_new;
        SAM_BIAS_new = SAM_BIAS_new + eps*PHI_BIAS_new;

        SAM_G_new =  sum(SAM_new.^2,3);

        AA = SAM_BIAS_new.^2 + setting.min_start_age + HH*SAM_G_new + HH0.*SAM_G_new(1,:) + HH1.*SAM_G_new(CNT+1,:) + HH2.*SAM_G_new(CNT+2,:);
        AA_Y = AA(ID_Y,:);
        AA_C = AA(ID_C,:);
        AA_A = AA(ID_A,:);

        PDEV = zeros(N0,nSamples,2.*setting.alpha);
        PDEV_BIAS = zeros(1,nSamples);

        PDEV(1:end-1,:,:) = PDEV(1:end-1,:,:) - 2.*BETA./(1-RHO.^2).*RHO.*(RHO.*SAM_new(1:end-1,:,:)-SAM_new(2:end,:,:));
        PDEV(2:end,:,:) = PDEV(2:end,:,:) + 2.*BETA./(1-RHO.^2).*(RHO.*SAM_new(1:end-1,:,:)-SAM_new(2:end,:,:));
        PDEV(1,:,:) = PDEV(1,:,:) - 2.*BETA.*SAM_new(1,:,:);

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

            PDEV = PDEV + (HH'*DELTA).*(2.*SAM_new);
            PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(DELTA,1).*(2.*SAM_new(1,:,:));
            PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*DELTA)).*(2.*SAM_new(UQ+1,:,:));
            PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*DELTA)).*(2.*SAM_new(UQ+2,:,:));

            PDEV_BIAS = PDEV_BIAS + sum(DELTA,1).*(2.*SAM_BIAS_new);
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

            PDEV = PDEV + (HH'*DELTA).*(2.*SAM_new);
            PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(DELTA,1).*(2.*SAM_new(1,:,:));
            PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*DELTA)).*(2.*SAM_new(UQ+1,:,:));
            PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*DELTA)).*(2.*SAM_new(UQ+2,:,:));

            PDEV_BIAS = PDEV_BIAS + sum(DELTA,1).*(2.*SAM_BIAS_new);
        end

        % ADD_AGE:
        if ~isempty(AA_A)
            ALPHA = - (AA_A-YA(:,1)).*YA(:,2).^(-2);

            ALPHA = QQ_A*ALPHA;

            DELTA = zeros(N1,nSamples);
            DELTA(UQ_A,:) = ALPHA;

            PDEV = PDEV + (HH'*DELTA).*(2.*SAM_new);
            PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(DELTA,1).*(2.*SAM_new(1,:,:));
            PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*DELTA)).*(2.*SAM_new(UQ+1,:,:));
            PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*DELTA)).*(2.*SAM_new(UQ+2,:,:));

            PDEV_BIAS = PDEV_BIAS + sum(DELTA,1).*(2.*SAM_BIAS_new);
        end

        PHI_new = PHI_new + eps*PDEV;
        PHI_BIAS_new = PHI_BIAS_new + eps*PDEV_BIAS;
    end

    SAM_new = SAM_new + eps*PHI_new;
    SAM_BIAS_new = SAM_BIAS_new + eps*PHI_BIAS_new;

    SAM_G_new =  sum(SAM_new.^2,3);

    AA = SAM_BIAS_new.^2 + setting.min_start_age + HH*SAM_G_new + HH0.*SAM_G_new(1,:) + HH1.*SAM_G_new(CNT+1,:) + HH2.*SAM_G_new(CNT+2,:);
    AA_Y = AA(ID_Y,:);
    AA_C = AA(ID_C,:);
    AA_A = AA(ID_A,:);

    PDEV = zeros(N0,nSamples,2.*setting.alpha);
    PDEV_BIAS = zeros(1,nSamples);

    PDEV(1:end-1,:,:) = PDEV(1:end-1,:,:) - 2.*BETA./(1-RHO.^2).*RHO.*(RHO.*SAM_new(1:end-1,:,:)-SAM_new(2:end,:,:));
    PDEV(2:end,:,:) = PDEV(2:end,:,:) + 2.*BETA./(1-RHO.^2).*(RHO.*SAM_new(1:end-1,:,:)-SAM_new(2:end,:,:));
    PDEV(1,:,:) = PDEV(1,:,:) - 2.*BETA.*SAM_new(1,:,:);

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

        PDEV = PDEV + (HH'*DELTA).*(2.*SAM_new);
        PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(DELTA,1).*(2.*SAM_new(1,:,:));
        PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*DELTA)).*(2.*SAM_new(UQ+1,:,:));
        PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*DELTA)).*(2.*SAM_new(UQ+2,:,:));

        PDEV_BIAS = PDEV_BIAS + sum(DELTA,1).*(2.*SAM_BIAS_new);
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

        PDEV = PDEV + (HH'*DELTA).*(2.*SAM_new);
        PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(DELTA,1).*(2.*SAM_new(1,:,:));
        PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*DELTA)).*(2.*SAM_new(UQ+1,:,:));
        PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*DELTA)).*(2.*SAM_new(UQ+2,:,:));

        PDEV_BIAS = PDEV_BIAS + sum(DELTA,1).*(2.*SAM_BIAS_new);
    end

    % ADD_AGE:
    if ~isempty(AA_A)
        ALPHA = - (AA_A-YA(:,1)).*YA(:,2).^(-2);

        ALPHA = QQ_A*ALPHA;

        DELTA = zeros(N1,nSamples);
        DELTA(UQ_A,:) = ALPHA;

        PDEV = PDEV + (HH'*DELTA).*(2.*SAM_new);
        PDEV(1,:,:) = PDEV(1,:,:) + HH0.*sum(DELTA,1).*(2.*SAM_new(1,:,:));
        PDEV(UQ+1,:,:) = PDEV(UQ+1,:,:) + (QQ*(HH1.*DELTA)).*(2.*SAM_new(UQ+1,:,:));
        PDEV(UQ+2,:,:) = PDEV(UQ+2,:,:) + (QQ*(HH2.*DELTA)).*(2.*SAM_new(UQ+2,:,:));

        PDEV_BIAS = PDEV_BIAS + sum(DELTA,1).*(2.*SAM_BIAS_new);
    end

    PHI_new = PHI_new + 0.5*eps*PDEV;
    PHI_BIAS_new = PHI_BIAS_new + 0.5*eps*PDEV_BIAS;

    AA = SAM_BIAS_new.^2 + setting.min_start_age + HH*SAM_G_new + HH0.*SAM_G_new(1,:) + HH1.*SAM_G_new(CNT+1,:) + HH2.*SAM_G_new(CNT+2,:);
    AA_Y = AA(ID_Y,:);
    AA_C = AA(ID_C,:);
    AA_A = AA(ID_A,:);

    LOGLIK_new = sum(sum(-0.5*PHI_new.^2-0.5*log(2*pi),3),1);
    LOGLIK_new = LOGLIK_new - 0.5*PHI_BIAS_new.^2 - 0.5*log(2*pi);

    LOGLIK_new = LOGLIK_new + sum(-BETA.*SAM_new(end,:,:).^2-BETA./(1-RHO.^2).*sum((SAM_new(1:end-1,:,:)-RHO.*SAM_new(2:end,:,:)).^2,1),3);

    % d18O:
    if ~isempty(AA_Y)
        MU = interp1(stack.A,stack.MU,AA_Y,'linear','extrap');
        SIG = interp1(stack.A,STACK_SIG,AA_Y,'linear',THSD);

        ZZ = (YY-data.scale.*MU-data.shift)./data.scale;

        if isinf(A2)&&isinf(B2)
            LOGLIK_new = LOGLIK_new + sum(-0.5.*ZZ.^2.*SIG.^2+log*(SIG)-log(data.scale),1);
        else
            LOGLIK_new = LOGLIK_new + sum(-(A2+1./2).*log(1+1./(2.*B2).*ZZ.^2.*SIG.^2)+log(SIG)-log(data.scale),1);
        end
    end

    % 14C:
    if ~isempty(AA_C)
        MU = zeros(size(AA_C));
        SIG = zeros(size(AA_C));

        for k = 1:length(cal)
            ID = (YC(:,end)==k);

            if sum(ID) > 0
                MU(ID,:) = interp1(cal(k).A,cal(k).MU,AA_C(ID,:),'linear','extrap');
                SIG(ID,:) = interp1(cal(k).A,cal(k).SIG,AA_C(ID,:),'linear',1./THSD);
            end
        end

        WW = sqrt(SIG.^2+YC(:,2).^2+YC(:,4).^2);
        ZZ = (MU+YC(:,3)-YC(:,1))./WW;

        if isinf(A1)&&isinf(B1)
            LOGLIK_new = LOGLIK_new + sum(-0.5.*ZZ.^2-log(WW),1);
        else
            LOGLIK_new = LOGLIK_new + sum(-(0.5+A1).*log(1+ZZ.^2./(2.*B1))-log(WW),1);
        end
    end

    % ADD_AGE:
    if ~isempty(AA_A)
        LOGLIK_new = LOGLIK_new + sum(-0.5.*(AA_A-YA(:,1)).^2.*YA(:,2).^(-2),1);
    end


    ID = (~isnan(LOGLIK_new))&(LOGLIK_new-LOGLIK_old>rand_seed(rr,:));

    SAM(:,ID,:) = SAM_new(:,ID,:);
    SAM_BIAS(ID) = SAM_BIAS_new(ID);

    PHI(:,ID,:) = PHI_new(:,ID,:);
    PHI_BIAS(:,ID) = PHI_BIAS_new(:,ID);
    LOGLIK_old(ID) = LOGLIK_new(ID);

    LOGLIK_old = LOGLIK_old - sum(sum(-0.5*PHI.^2-0.5*log(2*pi),3),1);
    LOGLIK_old = LOGLIK_old + 0.5*PHI_BIAS.^2 + 0.5*log(2*pi);
    PHI = normrnd(0,1,[N0,nSamples,2.*setting.alpha]);
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

SAMPLES.F = SAM.*sqrt(data.R);
SAM_G =  sum(SAM.^2,3);

SAMPLES.AGE = SAM_BIAS.^2 + setting.min_start_age + HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:);
SAMPLES.BIAS = SAM_BIAS;


SAMPLES.AGE_D = zeros(N0-1,nSamples);
D0 = SAMPLES.D0.*data.R;

CNT = zeros(N0-1,1);
HH = zeros(N0-1,N0);
for n = 1:N0-1
    CNT(n) = floor((D0(n)-setting.start_depth)/INT);
    HH(n,1:CNT(n)) = INT;
end
HH0 = - AP;
HH1 = min(D0-setting.start_depth-CNT.*INT+AP,INT);
HH2 = max(D0-setting.start_depth-CNT.*INT+AP-INT,0);

SAMPLES.AGE_D = SAM_BIAS.^2 + setting.min_start_age + HH*SAM_G + HH0.*SAM_G(1,:) + HH1.*SAM_G(CNT+1,:) + HH2.*SAM_G(CNT+2,:);

SAMPLES.APP_RATE = APP_RATE;
SAMPLES.ACC_RATE = SAM_G.*data.R;


end