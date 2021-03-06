%% audio_inpainting tests D�rfler Datasets
% this script reproduces the results presented in the paper. 
% Please note, that the LTFAT toolbox (http://ltfat.sourceforge.net/) has 
% to be included and added to the path.

clear;
close all;

addpath('Methods');
if (exist('dgtreal') ~= 2)
    error('Please include LTFAT toolbox');
else
    warning off; %some methods in the ltfat toolbox might produce warnings
end

%%

% 1: Strings
% 2: Piano
% 3: Percussion
% 4: Jazz
snum = 4;

filename = ['./wav/sig_' num2str(snum) '.wav'];
[s,fs] = audioread(filename);

if snum==4
    s = s(1:200000);
else
    s = s(1:100000);
end
L = length(s);


soundsc(s,fs)




%%   automatic testrun, RANDOM



p = 0.8; %percentage of coefficients to delete




snum = 1:4;
tauvec = [5e-3 1e-2 2.5e-2 5e-2 7.5e-2 1e-1 2.5e-1 5e-1 7.5e-1 1 5 10];


transvec = [{'gab'},{'erb'},{'wav'}];

threshvec = [{'lasso'},{'ew'},{'wgl'},{'pew'}];


%result matrix
res = zeros(length(snum),length(transvec),length(threshvec),length(tauvec));


% Gabor
settings.a = 160;
settings.M = 3125;
settings.width = 23.22;

% ERBLET
settings.ebins = 18;
settings.qvar = 0.08;

% Wavelet
settings.fmin = 100;
settings.wbins = 120;
settings.bw = 3;


settings.plotting = 0;


zzz = 1;
zzztotal = numel(res(:));
fix(clock)
for iii=1:length(snum)
    %get signal
    filename = ['./wav/sig_' num2str(snum(iii)) '.wav'];
    [s,fs] = audioread(filename);
    if snum(iii)==4
        s = s(1:200000);
    else
        s = s(1:100000);
    end
    
    if snum(iii) == 1
        settings.perslength = 72;
        settings.percussionpersistence = 0;
    elseif snum(iii) == 2
        settings.perslength = 48;
        settings.percussionpersistence = 0;
    elseif snum(iii) == 3
        settings.perslength = 24;
        settings.percussionpersistence = 1;
    else
        settings.perslength = 48;
        settings.percussionpersistence = 0;
    end
    L = length(s);
    %get mask
    rng(10);
    Mask = rand(size(s)) > p;
    depl_s = Mask.*s;
    M = logical(1-Mask);
    snr_m = @(sol) 20 *log10(std(s(M))/std(sol(M)-s(M)));
    
    for jjj=1:length(transvec)
        settings.trans = transvec{jjj};
        
        for kkk=1:length(threshvec)
            settings.thres = threshvec{kkk};
            
            for lll=1:length(tauvec)
                tau = tauvec(lll);
                
                sol = audio_inpainting(s,depl_s,Mask,fs,tau,settings);
                res(iii,jjj,kkk,lll) = snr_m(sol);
                fprintf('\n -- %d / % d -- \n',zzz,zzztotal);
                fprintf('%d: %s (%5s) - SNR: %2.2f',snum(iii),settings.trans,settings.thres, snr_m(sol));
                zzz = zzz + 1;
            end
        end 
    end
end
fix(clock)

%squeeze(res)
    
save('Experiments_randomMask_New','res');


%%   automatic testrun, CONSECUTIVE

gaplengthvec = [2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30];

snum = [1,2,3,4];
tauvec = [7.5e-3 1e-2 2.5e-2 5e-2 7.5e-2 1e-1 2.5e-1 5e-1 7.5e-1 1 2.5 5 7.5 10 25 50 75 100];

transvec = [{'gab'},{'erb'}];

threshvec = [{'lasso'},{'pew'}];


%result matrix
res = zeros(length(snum),length(gaplengthvec),length(transvec),length(threshvec),length(tauvec));
res_timing = res;

% Gabor
settings.a = 250;
settings.M = 3125;
settings.width = 150; %%23.22

% ERBLET
settings.ebins = 18;
settings.qvar = 0.08;


settings.plotting = 0;


zzz = 1;
zzztotal = numel(res(:));
fix(clock)

for iii=1:length(snum)
    %get signal
    filename = ['./wav/sig_' num2str(snum(iii)) '.wav'];
    [s,fs] = audioread(filename);
    s = s(1:100000);
    if snum(iii) == 1
        settings.perslength = 72;
        settings.percussionpersistence = 0;
    elseif snum(iii) == 2
        settings.perslength = 48;
        settings.percussionpersistence = 0;
    elseif snum(iii) == 3
        settings.perslength = 24;
        settings.percussionpersistence = 0;
        settings.ebins = 8;
        settings.qvar = 0.2;
    else
        settings.perslength = 48;
        settings.percussionpersistence = 0;
    end
    L = length(s);
    
    
    for mmm = 1:length(gaplengthvec)
        %get mask
        pos = 11000;
        startpos = pos;
        alle = 300; %ms
        gaplength = gaplengthvec(mmm);
        numc = round(gaplength*1e-3*fs);
        Mask = ones(size(s));
        eras = zeros(1,numc);
        while (pos + numc < L - startpos)
            Mask(pos:pos+numc-1) = eras;
            pos = pos + round(alle*1e-3*fs);
        end
        %figure(1), plot(Mask);
        depl_s = Mask.*s;
        M = logical(1-Mask);
        snr_m = @(sol) 20 *log10(std(s(M))/std(sol(M)-s(M)));

        for jjj=1:length(transvec)
            settings.trans = transvec{jjj};

            for kkk=1:length(threshvec)
                settings.thres = threshvec{kkk};
                
                tmpmax = 0;
                for lll=1:length(tauvec)
                    tau = tauvec(lll);
                    
                    tic
                    sol = audio_inpainting(s,depl_s,Mask,fs,tau,settings);
                    res_timing(iii,mmm,jjj,kkk,lll) = toc;
                    tmp = snr_m(sol);
                    res(iii,mmm,jjj,kkk,lll) = tmp;
                    
                    fprintf('\n -- %d / % d -- \n',zzz,zzztotal);
                    fprintf('%d: %s (%5s) - SNR: %2.2f',snum(iii),settings.trans,settings.thres, snr_m(sol));
                    zzz = zzz + 1;
                    
                    if tmp>tmpmax
                        tmpmax = tmp;
                    elseif abs(tmp)< 0.8*tmpmax
                        fprintf('\ndecreasing value, cont.\n');
                        zzz = zzz + (length(tauvec)-lll);
                        break;
                    end
                    
                    
                end
            end 
        end
    end
    %fix(clock)
end
fix(clock)


%squeeze(res)

save('Experiments_ConsecutiveDoerflerData_newERB','res');


%% Janssen algorithm...

%note that the inpainting toolbox from Adler et al. has to be included

snum = [1,2,3,4];
gaplengthvec = [2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30];

res_janssen = zeros(length(snum),length(gaplengthvec));

fix(clock)
for iii=1:length(snum)
    %get signal
    filename = ['./wav/sig_' num2str(snum(iii)) '.wav'];
    [s,fs] = audioread(filename);
    s = s(1:100000);
    L = length(s);
    for mmm = 1:length(gaplengthvec)
        %get mask
        pos = 11000;
        startpos = pos;
        alle = 300; %ms
        gaplength = gaplengthvec(mmm);
        numc = round(gaplength*1e-3*fs);
        Mask = ones(size(s));
        eras = zeros(1,numc);
        while (pos + numc < L - startpos)
            Mask(pos:pos+numc-1) = eras;
            pos = pos + round(alle*1e-3*fs);
        end
        %figure(1), plot(Mask);
        depl_s = Mask.*s;
        M = logical(1-Mask);
        snr_m = @(sol) 20 *log10(std(s(M))/std(sol(M)-s(M)));
        
        
        problemData.IMiss = M;
        problemData.x = depl_s;
        
        try
            sol = inpaintFrame_janssenInterpolation(problemData);
        catch
            error('Please include inpainting toolbox from A.Adler, Audio Inpainting');
        end
        
        tmp = snr_m(sol);
        res_janssen(iii,mmm) = tmp;
        
        fprintf('%d: Gaplength - %2.1f -> %2.2f\n',iii,gaplengthvec(mmm),tmp);
    end
end
fix(clock)



%% synthesis vs. analysis vs DR vs FISTA

snum = [1,2,3,4];
p = 0.8; %percentage of coefficients to delete

tauvec = [7.5e-5 1e-4 2.5e-4 5e-4 7.5e-4 1e-3 2.5e-3 5e-3 1e-2 2.5e-2 5e-2 7.5e-2 1e-1 2.5e-1 5e-1 7.5e-1 1 5 10];

transvec = [{'gab'},{'erb'},{'wav'}];

algovec = [{'forward-backward'},{'douglas-rachford'}];
synthesisvec = [{'synthesis'}, {'analysis'}];

% Gabor
settings.a = 160;
settings.M = 3125;
settings.width = 23.22;

% ERBLET
settings.ebins = 18;
settings.qvar = 0.08;

% Wavelet
settings.fmin = 100;
settings.wbins = 120;
settings.bw = 3;

settings.thres = 'lasso';

settings.plotting = 0;

res = zeros(length(snum),length(algovec),length(synthesisvec),length(transvec),length(tauvec));
res_relnorm = cell(length(snum),length(algovec),length(synthesisvec),length(transvec),length(tauvec));

zzz = 1;
zzztotal = numel(res(:));
fix(clock)
for iii=1:length(snum)
    %get signal
    filename = ['./wav/sig_' num2str(snum(iii)) '.wav'];
    [s,fs] = audioread(filename);
    s = s(1:100000);
    L = length(s);
    %get mask
    rng(10);
    Mask = rand(size(s)) > p;
    depl_s = Mask.*s;
    M = logical(1-Mask);
    snr_m = @(sol) 20 *log10(std(s(M))/std(sol(M)-s(M)));
    
    for jjj=1:2
        settings.algo = algovec{jjj};
        for kkk=1:2
            settings.synthesis = mod(kkk,2);
            for lll=1:length(transvec)
                settings.trans = transvec{lll};
                tmpmax = 0;
                for mmm=1:length(tauvec)
                    tau = tauvec(mmm);
                    [sol,relnormvec] = audio_inpainting(s,depl_s,Mask,fs,tau,settings);
                    tmp = snr_m(sol);
                    res (iii,jjj,kkk,lll,mmm) = tmp;
                    res_relnorm{iii,jjj,kkk,lll,mmm} = relnormvec;
                    fprintf('\n -- %d / % d -- \n',zzz,zzztotal);
                    fprintf('%d: %s (%5s) - SNR: %2.2f',snum(iii),settings.trans,settings.thres, snr_m(sol));
                    zzz = zzz + 1;
                    
                    if tmp>tmpmax
                        tmpmax = tmp;
                    elseif abs(tmp)< 0.66*tmpmax
                        fprintf('\ndecreasing value, cont.\n');
                        zzz = zzz + (length(tauvec)-mmm);
                        break;
                    end
                end
            end
        end
    end
end


save('Experiments_DRvsFistavsSynvsAna_relnorm','res_relnorm');


%%   automatic testrun, RANDOM with vary percentage p

pvec = [0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95];

snum = 1:4;
tauvec = [1e-3 2.5e-3 5e-3 7.5e-3 1e-2 2.5e-2 5e-2 7.5e-2 1e-1 2.5e-1 5e-1];
tauvec = [5e-3 7.5e-3 1e-2 2.5e-2 5e-2 7.5e-2 1e-1 2.5e-1 5e-1 1 2.5];

transvec = [{'gab'},{'erb'},{'wav'}];
transvec = [{'gab'},{'erb'}];

threshvec = [{'lasso'},{'ew'},{'wgl'},{'pew'}];

threshvec_ns = [{'lasso'},{'pew'}];
threshvec_gab= [{'lasso'},{'pew'}];


%result matrix
res = zeros(length(pvec), length(snum),length(transvec),length(threshvec),length(tauvec));


% Gabor
settings.a = 160;
settings.M = 3125;
settings.width = 23.22;

% ERBLET
settings.ebins = 18;
settings.qvar = 0.08;

% Wavelet
settings.fmin = 100;
settings.wbins = 120;
settings.bw = 3;



%settings.perslength = .20;
settings.plotting = 0;


zzz = 1;
zzztotal = length(pvec)*length(snum)*length(transvec)*2;


fix(clock)
for ggg = 1:length(pvec)
    p = pvec(ggg);
    for iii=1:length(snum)
        %get signal
        filename = ['./wav/sig_' num2str(snum(iii)) '.wav'];
        [s,fs] = audioread(filename);
        if snum(iii)==6
            s = s(1:200000); %260000
        else
            s = s(1:100000);
        end

        if snum(iii) == 1
            settings.perslength = 72;
            settings.percussionpersistence = 0;
        elseif snum(iii) == 2
            settings.perslength = 48;
            settings.percussionpersistence = 0;
        elseif snum(iii) == 3
            settings.perslength = 24;
            settings.percussionpersistence = 1;
        else
            settings.perslength = 48;
            settings.percussionpersistence = 0;
        end
        L = length(s);
        %get mask
        rng(10);
        Mask = rand(size(s)) > p;
        depl_s = Mask.*s;
        M = logical(1-Mask);
        snr_m = @(sol) 20 *log10(std(s(M))/std(sol(M)-s(M)));

        for jjj=1:length(transvec)
            settings.trans = transvec{jjj};

            if (strcmp(settings.trans,'gab'))
                threshvec = threshvec_gab;
            else
                threshvec = threshvec_ns;
            end

            for kkk=1:length(threshvec)
                settings.thres = threshvec{kkk};

                lll = 6;
                tau = tauvec(lll); %middle element
                
                sol = audio_inpainting(s,depl_s,Mask,fs,tau,settings);
                tempsnr_prev = snr_m(sol);
                res(ggg,iii,jjj,kkk,lll) = tempsnr_prev;
                
                fprintf('\n -- %d / % d -- \n',zzz,zzztotal);
                fprintf('%2.2f -> %d: %s (%5s) @tau=%f - SNR: %2.2f',pvec(ggg),snum(iii),settings.trans,settings.thres, tau, tempsnr_prev);
                
                
                lll = lll-1;
                tau = tauvec(lll); %next element
                
                sol = audio_inpainting(s,depl_s,Mask,fs,tau,settings);
                tempsnr = snr_m(sol);
                res(ggg,iii,jjj,kkk,lll) = tempsnr;
                
                fprintf('\n -- %d / % d -- \n',zzz,zzztotal);
                fprintf('%2.2f -> %d: %s (%5s) @tau=%f - SNR: %2.2f',pvec(ggg), snum(iii),settings.trans,settings.thres, tau, tempsnr);
                
                if (tempsnr_prev < tempsnr)
                    while (tempsnr > tempsnr_prev)
                        tempsnr_prev = tempsnr;
                        lll = lll - 1;
                        tau = tauvec(lll);
                        sol = audio_inpainting(s,depl_s,Mask,fs,tau,settings);
                        tempsnr = snr_m(sol);
                        res(ggg,iii,jjj,kkk,lll) = tempsnr;

                        fprintf('\n -- %d / % d -- \n',zzz,zzztotal);
                        fprintf('%2.2f -> %d: %s (%5s) @tau=%f - SNR: %2.2f',pvec(ggg), snum(iii),settings.trans,settings.thres, tau, tempsnr);
                    end
                else
                    lll = lll + 1;
                    tmp = tempsnr;
                    tempsnr = tempsnr_prev;
                    tempsnr_prev = tmp;
                    while(tempsnr > tempsnr_prev)
                        tempsnr_prev = tempsnr;
                        lll = lll + 1;
                        tau = tauvec(lll);
                        sol = audio_inpainting(s,depl_s,Mask,fs,tau,settings);
                        tempsnr = snr_m(sol);
                        res(ggg,iii,jjj,kkk,lll) = tempsnr;

                        fprintf('\n -- %d / % d -- \n',zzz,zzztotal);
                        fprintf('%2.2f -> %d: %s (%5s) @tau=%f - SNR: %2.2f',pvec(ggg), snum(iii),settings.trans,settings.thres, tau, tempsnr);
                    end
                end
              
                
                zzz = zzz + 1;
%                 for lll=1:length(tauvec)
%                     tau = tauvec(lll);
% 
%                     sol = audio_inpainting(s,depl_s,Mask,fs,tau,settings);
%                     res(ggg,iii,jjj,kkk,lll) = snr_m(sol);
%                     fprintf('\n -- %d / % d -- \n',zzz,zzztotal);
%                     fprintf('%d: %s (%5s) - SNR: %2.2f',snum(iii),settings.trans,settings.thres, snr_m(sol));
%                     zzz = zzz + 1;
%                 end
            end 
        end
    end
end
fix(clock)

save('Experiments_30-95percent_JazzQuintet','res');


