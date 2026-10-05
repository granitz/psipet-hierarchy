%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_hopf_subject_refinement.m
%
% 1. Load optimal parameters (gBest, aBest) from step 1 (or specify manually)
% 2. Load Ceff_best
% 3. For each subject, run iterative fitting with the same (gBest, aBest)
% 4. Save final subject-level connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rng(42)

maxNumCompThreads(12);

if ~exist('GRIDON','var')
    GRIDON = false;
end

% ------------- load best (g,a) from step 1 --------------
basepath = getenv('PSIPET_MODEL_ROOT');
if isempty(basepath)
    error('Set the environment variable PSIPET_MODEL_ROOT to the run folder (see README.md).');
end
if basepath(end) ~= '/', basepath = [basepath '/']; end
addpath(genpath(pwd));   % this repository's matlab/ folder, incl. external/
assert(exist('HopfModel_LNA','file')==2, 'HopfModel_LNA.m not found -- see matlab/external/README.md');
assert(exist('gaussfilt','file')==2,     'gaussfilt.m not found -- see matlab/external/README.md');

data_dir = sprintf('%sdata/',basepath);
res_dir = sprintf('%sresults/',basepath);
tspath = sprintf('%stseries/',data_dir); % followed by rest or restvisual and / ses-01 or ses-02 and then filename

if GRIDON
    grpPrefix = 'gridGroupResults';
    subjPrefix = 'subjectRefinement';
else
    grpPrefix = 'no_grid_GroupResults';
    subjPrefix = 'no_grid_subjectRefinement';
end

% load group‐level optimum
resdir = fullfile(res_dir, 'results_group_only');
bestfile = fullfile(resdir, ...
    sprintf('%s_group-%s_task-%s_ses-%s_tau-%d.mat', ...
        grpPrefix, GROUP, TASK, SES, Tau));

bestdat  = load(bestfile);
gBest    = bestdat.gBest;
aBest    = bestdat.aBest;
Ceff_best= bestdat.Ceff_best; % group-level connectivity for optimal solution

%%%
% --------------- Define parameters ------------
N        = str2double(ATLASNAME);
TR       = 2.205;  % Repetition time
sigma    = 0.01;   % noise amplitude
epsFC    = 0.0004;
epsFCtau = 0.0001;
maxC     = 0.2;
max_iters= 5000;

% Load SC
load(sprintf('%s/atlas_SchaeferTian116v3_sift_invnodevol_radius2_count_connectivity_ses-%s.mat',data_dir,SES));
% load(sprintf('%s/atlas_SchaeferTian116v3_sift_radius2_count_connectivity_ses-%s.mat',data_dir,SES));

SC = squeeze(mean(SC,1));
SC = SC / max(SC(:)) * maxC;

IgnoreNodes = [20,47,48,49,95,96,107];

indexN = setdiff(1:N, IgnoreNodes);
SC     = SC(indexN,indexN);

disp('Loaded structural connectivity')

% Load subjects name
load(sprintf('%sPSIPET-ID_group-%s_task-%s.mat',data_dir, GROUP, TASK),'subjects');

% Load timeseries
NSUB = length(subjects);
TS = [];  % Will be NSUB x N x T

for iSub = 1:NSUB
    subID = subjects{iSub};

    fname = sprintf('%s/%s/ses-%s/%s_ses-%s_task-%s_feature-%s_atlas-SchaeferTian116_timeseries.tsv', tspath, TASK, SES, subID, SES, TASK,FEATURE);

    if isfile(fname)
        % data = load(fname);
        % Get the timeseries variable (node x time)
        ts_sub = readmatrix(fname, ...
              'FileType','text', ...
              'Delimiter','\t')';

        % Keep only valid nodes: indexN x time
        ts_sub = ts_sub(indexN, :);
        % Reshape to 1 x N x T and concatenate
        ts_sub = reshape(ts_sub, 1, length(indexN), size(ts_sub, 2));
        TS = cat(1, TS, ts_sub);
    else
        fprintf('Warning: File not found for %s\n', subID);
    end
end

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                  % lowpass frequency of filter (Hz)
fhi = 0.08;                   % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

% Compute frequencies from subject timeseries
Nvalid = length(indexN);
f_diff_sub = zeros(NSUB, Nvalid);
FCemp = zeros(NSUB, Nvalid, Nvalid);

for sub_i = 1:NSUB
    sub_i
    clear signal_filt Power_Areas;
    % Extract timeseries for this subject (N x T)
    ts = squeeze(TS(sub_i, :, :));

    for seed=1:Nvalid
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) = filtfilt(bfilt,afilt,ts(seed,:));
    end
    signal_filt=signal_filt(:,10:end-10);

    [Ns, Tmaxred] = size(signal_filt);
    TT = Tmaxred;
    Ts = TT * TR;  % Total sampling time
    freq = (0:TT/2-1) / Ts;  % Frequency range
    nfreqs = length(freq);

    % Perform spectral analysis on filtered data
    for seed = 1:Nvalid
        pw = abs(fft(zscore(signal_filt(seed, :))));
        PowSpect = pw(1:floor(TT/2)).^2 / (TT/TR);
        Power_Areas = gaussfilt(freq, PowSpect, 0.005);
        [~, index] = max(Power_Areas);
        index = squeeze(index);
        f_diff_sub(sub_i, seed) = freq(index);
    end
end

% Average frequencies across subjects
f_diff = mean(f_diff_sub, 1);
wo = 2*pi * f_diff(:);

disp('Computed frequencies and FC from subject timeseries')

% aVec
aVec = aBest * ones(Nvalid,1);

% lags
lags_model = [TR*Tau];

% ---------- subject-level refinement -----------
Ceff_Subj      = zeros(NSUB, Nvalid, Nvalid);
fittFC_subj    = zeros(NSUB,1);
fittCOVtau_subj= zeros(NSUB,1);

for nsub = 1:NSUB
    disp(['Refining subject ' num2str(nsub) ' ...']);

    ts = TS(nsub,:,:);
    ts = squeeze(ts);
    ts = ts(:, 10:end-10);

    % Empirical FC
    FCemp_sub   = corrcoef(ts');
    Cov_sub     = cov(ts');
    CovTau_emp  = zeros(Nvalid,Nvalid);
    sigratio    = zeros(Nvalid,Nvalid);
    for i=1:Nvalid
        for j=1:Nvalid
            sigratio(i,j) = 1/sqrt(Cov_sub(i,i))/sqrt(Cov_sub(j,j));
            [clag, lags]  = xcov(ts(i,:), ts(j,:), Tau);
            idx = find(lags==Tau);
            CovTau_emp(i,j) = clag(idx)/size(ts,2);
        end
    end
    CovTau_emp = CovTau_emp .* sigratio;

    % FCemp_sub(isnan(FCemp_sub)) = 0;
    % CovTau_emp(isnan(CovTau_emp)) = 0;

    % Start from optimal group-level Ceff
    Cnew_sub = Ceff_best;
    olderror_sub = Inf;

    for iterS=1:max_iters
        [FCsim_sub, Cov0_sub, Ct_sub] = HopfModel_LNA(...
            Cnew_sub, aVec, gBest, wo, sigma, lags_model);

        % Cov at TR
        CovTau_sim0 = squeeze(Ct_sub);

        sigratiosim_sub = zeros(Nvalid,Nvalid);
        for i=1:Nvalid
            for j=1:Nvalid
                sigratiosim_sub(i,j) = 1/sqrt(Cov0_sub(i,i))/sqrt(Cov0_sub(j,j));
            end
        end
        CovTau_sim = CovTau_sim0 .* sigratiosim_sub;

        % errors
        errFC_s  = mean( (FCemp_sub(:)     - FCsim_sub(:)).^2 );
        errCOV_s = mean( (CovTau_emp(:)    - CovTau_sim(:)).^2 );
        totErr_s = errFC_s + errCOV_s;

        if mod(iterS,100)<0.1
            iterS
            if (olderror_sub - totErr_s)/totErr_s < 1e-3
                disp('No improvement => stop.');
                break;
            end
            if totErr_s > olderror_sub
                disp('Error up => stop.');
                break;
            end
            olderror_sub = totErr_s;
        end

        % gradient update
        for i=1:Nvalid
            for j=1:Nvalid
                if SC(i,j)>0 % only update existing connections
                    Cnew_sub(i,j) = Cnew_sub(i,j) ...
                      + epsFC   * (FCemp_sub(i,j)  - FCsim_sub(i,j)) ...
                      + epsFCtau* (CovTau_emp(i,j) - CovTau_sim(i,j));
                    if Cnew_sub(i,j)<0
                        Cnew_sub(i,j)=0;
                    end
                end
            end
        end
        % rescale
        Cnew_sub = Cnew_sub / max(Cnew_sub(:)) * maxC;
    end

    Ceff_Subj(nsub,:,:) = Cnew_sub;

    % final correlations
    fc_ = corr2(FCemp_sub(tril(true(Nvalid), -1)), ...
                FCsim_sub(tril(true(Nvalid), -1)));
    fittFC_subj(nsub) = fc_;
    cv_ = corr2(CovTau_emp(tril(true(Nvalid), -1)), ...
                CovTau_sim(tril(true(Nvalid), -1)));
    fittCOVtau_subj(nsub) = cv_;
end

% Expand Ceff_Subj to full size (NSUB x 416 x 416) with NaN for ignored nodes
Ceff_Subj_full = nan(NSUB, N, N);
for nsub = 1:NSUB
    Ceff_Subj_full(nsub, indexN, indexN) = squeeze(Ceff_Subj(nsub, :, :));
end

% save subject‐level results with matching prefix
resdir2 = fullfile(res_dir, 'results_subject_refinement');
if ~exist(resdir2, 'dir')
    mkdir(resdir2);
end

%% Hierarchy Ceff

%% PB
hierarchicallevels_full = nan(NSUB, N);  % Full N=416 array with NaN

for nsub=1:NSUB
    Ceff=squeeze(Ceff_Subj(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels(nsub,:)=gamma';

    % Expand gamma to full size with NaN for ignored nodes
    hierarchicallevels_full(nsub, indexN) = gamma';

    H=(meshgrid(gamma)-meshgrid(gamma)'-1).^2;
    F0=sum(sum((A.*H)))/sum(sum(A));
    trophiccoherence(nsub)=1-F0;
end

outname = sprintf('%s_group-%s_task-%s_ses-%s_tau-%d_atlas-%s.mat', ...
    subjPrefix, GROUP, TASK, SES, Tau, ATLASNAME);

save(fullfile(resdir2, outname), ...
     'Ceff_Subj','Ceff_Subj_full','hierarchicallevels','hierarchicallevels_full','trophiccoherence', 'fittFC_subj', 'fittCOVtau_subj', ...
     'gBest', 'aBest', 'Tau', 'sigma','subjects','indexN');
