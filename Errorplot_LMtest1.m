clear all
close all
%% test details
odir    = 'modelfiles_testing';
arun    = 'hef_LMtest1';
modstr  = 'JF10';
freq    = '1.0';
%fine
acvols  = 2:2:4;
lgtauMRs = 7:0.2:8;
% custom
% acvols = [10];
% lgtauMRs = [7.5];

if1hz   = 0;
ifwt    = 1;
ifplot  = 1;
redo    = 1;


%% parms
Rearth = 6371e3;
Zcore = 2891e3;
cutoff_Q = 650;
mod_or_trans = 'trans'; % 'mod' or 'trans' for whether the comparison is the translated points or the smooth model

%% starting
cd('~/Dropbox/CIDER_grain_size_project/Seismowork/TPd_to_VQ/LMtesting');
addpath('../matlab_guts');
if if1hz; hzstr = '_1hz'; else hzstr = ''; end
% get PREM and QL6
prm = prem('depths',[0:5:2891]'); 
qrm = ql6('depths',[0:5:2891]'); 

%% weight
wt.z = [0:1e3:Rearth]';
W = ones(size(wt.z));
if ifwt
W(wt.z > 1000e3) = 0.7;
W(wt.z > 2000e3) = 0.1;
W(wt.z > 500e3 & wt.z < 800e3) = 2;
W(wt.z < 35e3) = 0; % none in crust
W(wt.z > Zcore) = 0; % none in core
wtstr = '';
else
wtstr = 'Un';
end
wt.wt = W;
figure(9), clf, 
plot(wt.wt,wt.z/1000,'Linewidth',2)
set(gca,'ydir','reverse','xlim',[0 2.5],'ylim',[0 Zcore]./1000,'Fontsize',14)
xlabel('weight'), ylabel('depth (km)')

%% loop
if redo
mft = nan(length(acvols),length(lgtauMRs));
misfit = struct('Vs',mft,'Vp',mft,'Qm',mft);

for ia = 1:length(acvols);
for it = 1:length(lgtauMRs);
acvstr = num2str(acvols(ia),'%.0f');
taustr = num2str(lgtauMRs(it),'%.1f');
if strcmp(mod_or_trans,'trans')
    tfile = ['trans_' arun '_V' acvstr '_lgTau' taustr '_' modstr '_' freq];

    fprintf('Reading translated aspect output file\n')
    fprintf('for V = %.f and lgTauMR = %.2f\n',acvols(ia),lgtauMRs(it))
    [ ASPPOST ] = read_asppost( [odir '/' tfile] ); 

    %% Mask crazy Qm
    indq_lm = abs(ASPPOST.Qm)>cutoff_Q & ASPPOST.P > 24.3e9;
    ASPPOST.Qm(indq_lm) = cutoff_Q;
    indq_um = abs(ASPPOST.Qm)>cutoff_Q & ASPPOST.P <= 24.3e9;
    ASPPOST.Qm(indq_um) = cutoff_Q;
    ASPPOST.Qm(ASPPOST.Qm==cutoff_Q)=nan;

    %% calc differences to models
    dVs = ASPPOST.Vs_1hz - interp1(prm.depth*1e3,prm.vs*1e3,ASPPOST.Z);
    dVp = ASPPOST.Vp_1hz - interp1(prm.depth*1e3,prm.vp*1e3,ASPPOST.Z);
    dQm = ASPPOST.Qm - interp1(qrm.depth*1e3,qrm.qu,ASPPOST.Z);
    Z = ASPPOST.Z;

    %% plot
    if ifplot
    figure(1), clf, set(gcf,'pos',[20 600 600 600]), hold on
    plot(ASPPOST.Vs_1hz,ASPPOST.Z,'bo')
    plot(1e3*prm.vs,1e3*prm.depth,'r','Linewidth',2.5)
    plot(interp1(prm.depth*1e3,prm.vs*1e3,ASPPOST.Z),ASPPOST.Z,'g.')
    plot(abs(dVs),ASPPOST.Z,'ko')
    set(gca,'ydir','reverse')
    title(sprintf('for V = %.f and lgTauMR = %.2f\n',acvols(ia),lgtauMRs(it)),'FontSize',18)
    figure(2), clf, set(gcf,'pos',[620 600 600 600]), hold on
    plot(ASPPOST.Qm,ASPPOST.Z,'bo')
    plot(qrm.qu,1e3*qrm.depth,'r','Linewidth',2.5)
    plot(interp1(qrm.depth*1e3,qrm.qu,ASPPOST.Z),ASPPOST.Z,'g.')
    plot(abs(dQm),ASPPOST.Z,'ko')
    set(gca,'ydir','reverse')
    title(sprintf('for V = %.f and lgTauMR = %.2f\n',acvols(ia),lgtauMRs(it)),'FontSize',18)
    pause
    end

elseif strcmp(mod_or_trans,'mod')
    fprintf('Reading translated aspect output file\n')
    fprintf('for V = %.f and lgTauMR = %.2f\n',acvols(ia),lgtauMRs(it))
    
    mfile = ['1D_' arun '_V' acvstr '_lgTau' taustr '_' modstr '_' freq hzstr '.model'];
    try
    MODEL = read_modelfile([odir,'/',mfile]);
    catch
    MODEL.Vs = nan; MODEL.Vp=nan;MODEL.Qm=nan;MODEL.Z=nan;
    end
    
    %% calc differences to models
    dVs = MODEL.Vs - interp1(prm.depth*1e3,prm.vs*1e3,MODEL.Z);
    dVp = MODEL.Vp - interp1(prm.depth*1e3,prm.vp*1e3,MODEL.Z);
    dQm = MODEL.Qm - interp1(qrm.depth*1e3,qrm.qu,MODEL.Z);
    Z = MODEL.Z;

    % plot
    if ifplot
    figure(1), clf, set(gcf,'pos',[20 600 600 600]), hold on
    plot(MODEL.Vs,MODEL.Z,'bo')
    plot(1e3*prm.vs,1e3*prm.depth,'r','Linewidth',2.5)
    plot(interp1(prm.depth*1e3,prm.vs*1e3,MODEL.Z),MODEL.Z,'g.')
    plot(abs(dVs),MODEL.Z,'ko')
    set(gca,'ydir','reverse','ylim',[0 Zcore],'xlim',[-1 10]*1e3,...
        'XTick',[-1e3:1e3:1e4],'XTickLabel',num2str([-1:1:10]'))
    title(sprintf('for V = %.f and lgTauMR = %.2f\n',acvols(ia),lgtauMRs(it)),'FontSize',18)
    figure(2), clf, set(gcf,'pos',[620 600 600 600]), hold on
    plot(MODEL.Qm,MODEL.Z,'bo')
    plot(qrm.qu,1e3*qrm.depth,'r','Linewidth',2.5)
    plot(interp1(qrm.depth*1e3,qrm.qu,MODEL.Z),MODEL.Z,'g.')
    plot(abs(dQm),MODEL.Z,'ko')
    set(gca,'ydir','reverse','ylim',[0 Zcore],'xlim',[00 1000])
    title(sprintf('for V = %.f and lgTauMR = %.2f\n',acvols(ia),lgtauMRs(it)),'FontSize',18)
    pause
    end
end % if mod/trans

% add in weights
if ~exist('weight','var'), weight = linterp(wt.z,wt.wt,Z); end
dVs = weight.*dVs;
dVp = weight.*dVp;
dQm = weight.*dQm;

%% calc misfits
misfit.Vs(ia,it) = norm(dVs(~isnan(dVs)))/sum(~isnan(dVs));
misfit.Vp(ia,it) = norm(dVp(~isnan(dVp)))/sum(~isnan(dVp));
misfit.Qm(ia,it) = norm(dQm(~isnan(dQm)))/sum(~isnan(dQm));

end %loop on TauMRs
end % loop on acvos

end % if redo


figure(10),clf
contourf(misfit.Vs,30)
shading flat
% contourf((misfit.Vs-mean(misfit.Vs,2)*ones(1,length(TauMRs)))./(mean(misfit.Vs,2)*ones(1,length(TauMRs))))
set(gca,'YTick',[1:length(acvols)],'YTickLabel',num2str(acvols'),...
        'XTick',[1:length(lgtauMRs)],'XTickLabel',num2str(lgtauMRs'),...
        'FontSize',12)
ylabel('Activation Volume'), xlabel('lgTauMR')
title(sprintf('%sWeighted Vs misfit',wtstr),'FontSize',18)

figure(11),clf
contourf(misfit.Qm,30)
shading flat
set(gca,'YTick',[1:length(acvols)],'YTickLabel',num2str(acvols'),...
        'XTick',[1:length(lgtauMRs)],'XTickLabel',num2str(lgtauMRs'),...
        'FontSize',12)
ylabel('Activation Volume'), xlabel('lgTauMR')
title(sprintf('%sWeighted Vs misfit',wtstr),'FontSize',18)
colorbar
