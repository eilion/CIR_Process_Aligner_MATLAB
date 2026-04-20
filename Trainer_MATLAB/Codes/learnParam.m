function param = learnParam(data,param,SAMPLES,setting)

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma = 1e-3;

L = length(data);

max_iters = 10000;

INT = setting.interval;
BETA = param.BETA;
RHO = param.RHO;

Mw = 0;
Vw = 0;

MM_0 = 0;
MM_1 = 0;
F_00 = 0;
F_11 = 0;
F_01 = 0;
F_22 = 0;

for ll = 1:L
    MM_0 = MM_0 + 0.5.*size(SAMPLES(ll).F,1).*setting.K;
    MM_1 = MM_1 + 0.5.*(size(SAMPLES(ll).F,1)-1).*setting.K;
    F_00 = F_00 + sum(SAMPLES(ll).F(1:end-1,:,:).^2,"all")./size(SAMPLES(ll).F,2);
    F_11 = F_11 + sum(SAMPLES(ll).F(2:end,:,:).^2,"all")./size(SAMPLES(ll).F,2);
    F_01 = F_01 + sum(SAMPLES(ll).F(1:end-1,:,:).*SAMPLES(ll).F(2:end,:,:),"all")./size(SAMPLES(ll).F,2);
    F_22 = F_22 + sum(SAMPLES(ll).F(end,:,:).^2,"all")./size(SAMPLES(ll).F,2);
end


for r = 1:max_iters
    PDEVs = zeros(2,1);

    RHO_INT = sqrt(RHO^INT);

    % BETA:
    PDEVs(1) = MM_0./BETA - 1./(1-RHO_INT.^2).*F_00 - RHO_INT.^2./(1-RHO_INT.^2).*F_11 + 2./(1-RHO_INT.^2).*RHO_INT.*F_01 - F_22;
    
    % RHO:
    PDEVs(2) = - BETA.*INT.*RHO_INT.^2./RHO./((1-RHO_INT.^2).^2).*(F_00+F_11) + BETA.*INT.*RHO_INT./RHO.*(1+RHO_INT.^2)./((1-RHO_INT.^2).^2).*F_01 + MM_1.*INT.*RHO_INT.^2./RHO./(1-RHO_INT.^2);


    Mw = beta1*Mw - (1-beta1)*PDEVs;
    Vw = beta2*Vw + (1-beta2)*PDEVs.*PDEVs;

    BETA = BETA - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1)./(sqrt(Vw(1))+epsilon);
    RHO = RHO - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(2)./(sqrt(Vw(2))+epsilon);
end

param.BETA = BETA;
param.RHO = RHO;


end