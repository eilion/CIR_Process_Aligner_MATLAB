function getFigures(outputFile)

PATH = ['Outputs/',outputFile,'/results.mat'];
AA = load(PATH);
data = AA.data;
setting = AA.setting;
SAMPLES = AA.SAMPLES;
stack = AA.target.stack;
C14_CI = AA.C14_CI;

L = length(data);


for ll = 1:L
    % CIR:
    fig = figure;

    if isempty(data(ll).d18O)
        set(fig,'Position',[20 20 800 450]);
    else
        set(fig,'Position',[20 20 800 900]);
    end
    movegui(fig,'center');
    drawnow;

    % age model
    if isempty(data(ll).d18O)
        subplot(1,3,1:2);
    else
        subplot(2,3,1:2);
    end
    hold on;
    tt = data(ll).name;
    tt(tt=='_') = '-';
    title(tt,'FontSize',14);

    DD = [SAMPLES(ll).D0;data(ll).depth];
    LB = [quantile(SAMPLES(ll).AGE_D,0.025,2);quantile(SAMPLES(ll).AGE,0.025,2)];
    MB = [quantile(SAMPLES(ll).AGE_D,0.5,2);quantile(SAMPLES(ll).AGE,0.5,2)];
    UB = [quantile(SAMPLES(ll).AGE_D,0.975,2);quantile(SAMPLES(ll).AGE,0.975,2)];

    [~,ORDER] = sort(DD,'ascend');
    DD = DD(ORDER);
    LB = LB(ORDER);
    MB = MB(ORDER);
    UB = UB(ORDER);

    [DD,ID,~] = unique(DD);
    LB = LB(ID);
    MB = MB(ID);
    UB = UB(ID);

    TT = linspace(min(DD),max(DD),1001)';

    LB = interp1(DD,LB,TT);
    MB = interp1(DD,MB,TT);
    UB = interp1(DD,UB,TT);

    plot(SAMPLES(ll).AGE_D,SAMPLES(ll).D0,'color',[0 0 0 0.015]);

    plot(LB,TT,'--k');
    plot(MB,TT,'r');
    plot(UB,TT,'--k');

    if ~isempty(C14_CI(ll).depth)
        for m = 1:size(C14_CI(ll).depth,1)
            plot([C14_CI(ll).lower(m,1),C14_CI(ll).upper(m,2)],C14_CI(ll).depth(m).*ones(1,2),'b','LineWidth',2);
        end
    end

    if ~isempty(data(ll).ETC)
        for n = 1:size(data(ll).ETC,1)
            plot([data(ll).ETC(n,1)-1.96.*data(ll).ETC(n,2),data(ll).ETC(n,1)+1.96.*data(ll).ETC(n,2)],data(ll).depth(data(ll).ETC(n,3))*ones(1,2),'g','LineWidth',2);
        end
    end

    MIN_Y = min(DD);
    MAX_Y = max(SAMPLES(ll).D0);

    MIN_X = setting(ll).min_start_age;
    MAX_X = ceil(max(UB)/5)*5;

    xlim([MIN_X MAX_X]);
    ylim([MIN_Y MAX_Y]);

    xlabel('Age (kyr)','FontSize',12);
    ylabel('Depth (m)','FontSize',12);


    % log sedimentation rates
    if isempty(data(ll).d18O)
        subplot(1,3,3);
    else
        subplot(2,3,3);
    end
    hold on;
    title('log-sedimentation rates','FontSize',14);

    DD = SAMPLES(ll).D;
    TT = linspace(0,data(ll).depth(end),1001)';

    ACC_RATE = zeros(length(TT),size(DD,2));
    for k = 1:size(DD,2)
        ACC_RATE(:,k) = interp1(DD(:,k),-log(SAMPLES(ll).ACC_RATE(:,k)),TT);
    end

    LB = quantile(ACC_RATE,0.025,2);
    MB = quantile(ACC_RATE,0.5,2);
    UB = quantile(ACC_RATE,0.975,2);

    plot(ACC_RATE,TT,'color',[0 0 0 0.015]);

    plot(LB,TT,'--k');
    plot(MB,TT,'r');
    plot(UB,TT,'--k');

    MIN = min(LB(~isnan(LB))) - 0.5;
    MAX = max(UB(~isnan(UB))) + 0.5;

    if ~isempty(data(ll).C14)
        plot(MIN.*ones(length(sum(~isnan(data(ll).C14(:,1)))),1),data(ll).depth(data(ll).C14(:,end)),'>b','LineWidth',1);
    end

    if ~isempty(data(ll).ETC)
        plot(MIN.*ones(length(sum(~isnan(data(ll).ETC(:,1)))),1),data(ll).depth(data(ll).ETC(:,end)),'square','Color','g','LineWidth',1);
    end

    xlim([MIN MAX]);
    ylim([MIN_Y MAX_Y]);

    xlabel('log-sed.rate','FontSize',12);



    % alignment to stack
    if ~isempty(data(ll).d18O)
        subplot(2,3,4:6);
        hold on;
        title('Alignment to the Stack','FontSize',14);

        xx = [stack.A;flipud(stack.A)];
        yy = [stack.MU-1.96.*stack.SIG;flipud(stack.MU+1.96.*stack.SIG)];
        patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.10,'EdgeColor','none');
        plot(stack.A,stack.MU,':k','LineWidth',2);

        YY = (data(ll).d18O(:,1)-data(ll).shift)/data(ll).scale;

        AGE = SAMPLES(ll).AGE(data(ll).d18O(:,2),:);
        AGE_lower = quantile(AGE,0.025,2);
        AGE_median = quantile(AGE,0.5,2);
        AGE_upper = quantile(AGE,0.975,2);

        for n = 1:size(AGE,1)
            plot([AGE_lower(n),AGE_upper(n)],YY(n)*ones(1,2),'r','LineWidth',1);
        end

        [~,order] = sort(AGE_median,'ascend');

        plot(AGE_median(order),YY(order),'-*r','LineWidth',1);

        if ~isempty(C14_CI(ll).depth)
            plot(C14_CI(ll).median,5.5.*ones(length(C14_CI(ll).median),1),'^b','LineWidth',1);
        end

        if ~isempty(data(ll).ETC)
            plot(data(ll).ETC(:,1),5.5.*ones(size(data(ll).ETC,1),1),'square','Color','g','LineWidth',1);
        end

        xlim([MIN_X MAX_X]);
        ylim([2.5 5.5]);

        set(gca,'ydir','reverse');

        xlabel('Age (kyr)','FontSize',12);
        ylabel(['\delta^1^8O (',char(8240),')'],'FontSize',12);
    end
end

end