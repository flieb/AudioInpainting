%% plotting results consecutive...
%close all;

load('results_Consecutive_JanssenAlgo.mat');
load('Experiments_ConsecutiveDoerflerData_GABERB.mat');


tauvec = [2.5e-2 5e-2 7.5e-2 1e-1 2.5e-1 5e-1 7.5e-1 1 5 10];
gaplengthvec = [2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30];
transvec = [{'gab'},{'erb'},{'wav'}];

res2 = squeeze(max(res,[],6));
%res2 = squeeze(res(:,:,:,:,8));

%res = zeros(length(snum),length(gaplengthvec),length(transvec),length(threshvec),length(tauvec));

col = get(groot,'DefaultAxesColorOrder');

col2 = [col(1,:); col(4,:); col(3,:); col(2,:)];

for kkk=1:4
    figure( kkk ), cla; hold on;
    
    semilogx(gaplengthvec,squeeze(res2(kkk,:,1,1)),'*-','Color',col(1,:),'MarkerSize',4);
    semilogx(gaplengthvec,squeeze(res2(kkk,:,1,2)),'o--','Color',col(1,:),'MarkerSize',4);
    semilogx(gaplengthvec,squeeze(res2(kkk,:,2,1)),'*-','Color',col(2,:),'MarkerSize',4);
    semilogx(gaplengthvec,squeeze(res2(kkk,:,2,2)),'o--','Color',col(2,:),'MarkerSize',4);
    plot(gaplengthvec,res_janssen(kkk,:),'s-','Color',col(5,:),'MarkerSize',4);
    hold off;
    ax20
    %legend('gab-lasso','gab-ew','gab-pew','erb-lasso','erb-ew','erb-pew','Location','best');
   
    
    if kkk==2
        ax = gca;
        ax.YTick= [6,8,10,12,14,16];
%         hh = legend('GAB-LASSO','GAB-PEW','ERB-LASSO','ERB-PEW','Janssen','Location','northeast');
%         hh.Interpreter = 'latex';
%         tmp2 = hh.Position;
%         hh.Position = [tmp2(1)-0.03 tmp2(2)-.06 tmp2(3) tmp2(4)];
    end
    if kkk==1
        hh = legend('GAB-LASSO','GAB-PEW','ERB-LASSO','ERB-PEW','Janssen','Location','northeast');
        hh.Interpreter = 'latex';
        tmp2 = hh.Position;
        hh.Position = [tmp2(1)-0.03 tmp2(2)-.06 tmp2(3) tmp2(4)];
    end
    %title(['Signal: ' num2str(kkk)]);
    
    xlabel('Gap Length (ms)');
    ylabel('$\textrm{SNR}_{\textrm{M}}$ (dB)');
    
    filename = ['C:\Users\lieb\Documents\Promotion\__Dissertation__\figures\audioinpainting\SNR_gaps_' num2str(kkk)];
    %fig2pdf(gcf,filename,55);
    
end


%% plotting results consecutive SOLUTIONS
if exist('fig','var')
    try 
        close(fig);
    catch
    end
end

load('Example_Solution_Testsignal4_22p5msGaplength_GabL_GabP_ErbL_ErbP.mat');

alpha = .7;
lw = 0.7;
t = linspace(0,100000/44100,100000)*1000;

col = get(groot,'DefaultAxesColorOrder');


fig = figure(1);
plot(t,s,'-','LineWidth',lw,'Color',[0.3,0.3,0.3]); hold on;
%plot(t,sol1,'LineWidth',lw); %gab-lasso
plot(t,sol2,'-','LineWidth',lw,'Color',[col(1,:), alpha]); %gab-pew
%plot(t,sol3,'LineWidth',lw); %erb-lasso
plot(t,sol4,'-','LineWidth',lw,'Color',[col(2,:),alpha]); %erb-pew
%plot(t,x,'-','LineWidth',lw,'Color',[col(5,:),alpha]);
hold off;
xlim([1.748 1.772]*1000);
ylim([-0.5, 0.7]);

%hh = legend('Original','GAB-PEW','ERB-PEW','Janssen','Location','northeast');
hh = legend('Original','GAB-PEW','ERB-PEW','Location','northeast');
hh.Interpreter = 'latex';
tmp2 = hh.Position;
hh.Position = [tmp2(1)+0.02 tmp2(2)+.01 tmp2(3) tmp2(4)];


xlabel('Time (ms)');
ylabel('Amplitude');


filename = 'C:\Users\lieb\Documents\Promotion\__Dissertation__\figures\audioinpainting\ExampleSol4';
fig2pdf(gcf,filename,1);

%plotting timing results consecutive
% load('results_Consecutive_JanssenAlgo_11000.mat');
% 
% res = squeeze(res_timing);
% M = zeros(4,12,2,2);
% S = M;
% 
% for jj=1:4
%     for kk=1:12
%         for ll=1:2
%             for ii=1:2
%                 tmp = squeeze(res(jj,kk,ll,ii,:));
%                 tmp2 = tmp(tmp>0);
%                 M(jj,kk,ll,ii) = mean(tmp);
%                 S(jj,kk,ll,ii) = std(tmp);
%             end
%         end
%     end
%     
%     figure(jj)
%     plot(squeeze(M(jj,:,1,1))); hold on
%     plot(squeeze(M(jj,:,1,2)));
%     plot(squeeze(M(jj,:,2,1)));
%     plot(squeeze(M(jj,:,2,2))); hold off;
% end

%% results random



load('Experiments_randomMask_New.mat');

%res = zeros(length(snum),length(transvec),length(threshvec),length(tauvec));
res2 = squeeze(max(res,[],4));

transvec = [{'gab'},{'erb'},{'wav'}];

threshvec = [{'lasso'},{'ew'},{'wgl'},{'pew'}];
    
fprintf('           1          2          3          4\n');
for k = [1,3,2]
    fprintf('%3s \n',transvec{k});
    for j = 1:4
        fprintf('%5s',threshvec{j});
        for l=1:4
            fprintf('    %2.2f  %2.2f  %2.2f  %2.2f\n',res2(l,k,j));
        end
        fprintf('\n');
    end
end


%% results syn vs ana vs dr vs fista

load('Experiments_DRvsFistavsSynvsAna.mat');


res2 = squeeze(max(res,[],5));
transvec = [{'gab'},{'erb'},{'wav'}];
synthesisvec = [{'synthsis'}, {'analysis'}];



for m=1:4
    fprintf('------------------------ %d: --------------------------\n',m);
    fprintf('                         FISTA |  DR   \n');
    fprintf('            GAB    WAV    ERB  |  GAB    WAV    ERB    \n');
    for k = 1:2
        fprintf([synthesisvec{k} ': ']);
        for j = 1:2
            for l=[1 3 2]
                fprintf(' %5.2f ',res2(m,j,k,l));
            end
            fprintf('|');
        end
        fprintf('\n');
    end
end


%% results convergence dr vs fista

load('Experiments_DRvsFistavsSynvsAna_relnorm.mat');

sign = 4; %1-tauk=11
tauk = 11;
maxit = 200;


figure(1), cla; hold off;

semilogy(squeeze(res_relnorm{sign,1,1,1,tauk}),'LineWidth', 1);
hold on;
semilogy(squeeze(res_relnorm{sign,1,2,1,tauk}),'LineWidth', 1);
semilogy(squeeze(res_relnorm{sign,2,1,1,tauk}),'LineWidth', 1);
semilogy(squeeze(res_relnorm{sign,2,2,1,tauk}),'LineWidth', 1);

ylim([2e-4 1]);

hh = legend('FISTA-Synthesis','FISTA-Analysis','DR-Synthesis','DR-Analysis','Location','best');
hh.Interpreter = 'latex';
xlabel('Number of Iterations');
ylabel('$\|z\|_{\mathrm{rel}}$');


filename = 'C:\Users\lieb\Documents\Promotion\__Dissertation__\figures\audioinpainting\convergenceRate';
fig2pdf(gcf,filename,1);









