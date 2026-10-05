%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_hopf_grid_group.m
%
% 1. For each (g,a), run the group-level iterative procedure
% 2. Save final group-level Ceff & error metrics
% 3. Find optimal (g,a) by min error or max correlation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rng(42)

maxNumCompThreads(12);

if ~exist('GRIDON','var')
    GRIDON = false;
end

% --------------- Define ranges ----------------
if GRIDON
   % First,  coarse, logarithmic, wide enough that the optimum can't sit at an edge:
   % Gvals = logspace(-2, 1, 50);
   % Then on pass 2: fine, linear, if it comes out near 0.75 for example:
   Gvals = 0.20:0.005:0.80;
   Avals = -0.02;
else
   Gvals = 0.2;
   Avals = -0.02;
end

% (The published selected values, for reference:
%  active_01 g=0.435, active_02 g=0.420, placebo_01 g=0.405, placebo_02 g=0.400.
%  To refit at a single g rather than search, set Gvals to that value here.)

% --------------- Define parameters ------------
N        = str2double(ATLASNAME);
TR       = 2.205;  % Repetition time
sigma    = 0.01;   % noise amplitude
epsFC    = 0.0004;
epsFCtau = 0.0002;
maxC     = 0.2;
max_iters= 5000;

% --------------- Load data (SC, TS, etc.) -----

basepath = getenv('PSIPET_MODEL_ROOT');
if isempty(basepath)
    error('Set the environment variable PSIPET_MODEL_ROOT to the run folder (see README.md).');
end
if basepath(end) ~= '/', basepath = [basepath '/']; end
addpath(genpath(pwd));   % this repository's matlab/ folder, incl. external/
assert(exist('HopfModel_LNA','file')==2, 'HopfModel_LNA.m not found -- see matlab/external/README.md');
assert(exist('gaussfilt','file')==2,     'gaussfilt.m not found -- see matlab/external/README.md');

data_dir = sprintf('%sdata/',basepath);
res_dir  = sprintf('%sresults/',basepath);

tspath = sprintf('%stseries/',data_dir); % followed by rest or restvisual and / ses-01 or ses-02 and then filename

%%%%%%%%%%%
% Load SC %
%%%%%%%%%%%
load(sprintf('%s/atlas_SchaeferTian116v3_sift_invnodevol_radius2_count_connectivity_ses-%s.mat',data_dir,SES));
% load(sprintf('%s/atlas_SchaeferTian116v3_sift_radius2_count_connectivity_ses-%s.mat',data_dir,SES));

SC = squeeze(mean(SC,1));
SC     = SC / max(SC(:))*maxC;

IgnoreNodes = [20,47,48,49,95,96,107];

indexN = setdiff(1:N, IgnoreNodes);
SC     = SC(indexN,indexN);

disp('Loaded structural connectivity')

% Load subjects name
% group: active, placebo; task: rest, restvisual; ses: 01, 02
load(sprintf('%sPSIPET-ID_group-%s_task-%s.mat',data_dir, GROUP, TASK),'subjects');

% Load timeseries
NSUB = length(subjects);
TS = [];  % NSUB x N x T

for iSub = 1:NSUB
    subID = subjects{iSub};

    fname = sprintf('%s/%s/ses-%s/%s_ses-%s_task-%s_feature-%s_atlas-SchaeferTian116_timeseries.tsv', tspath, TASK, SES, subID, SES, TASK, FEATURE);

    if isfile(fname)
        % Get the timeseries variable (node x time)

        ts_sub = readmatrix(fname, 'FileType','text', 'Delimiter','\t')';  % node x time

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

if ~exist(fullfile(data_dir,'derived_data'),'dir'), mkdir(fullfile(data_dir,'derived_data')); end
save(fullfile(data_dir, '/derived_data/', sprintf('/wo_group-%s_task-%s_ses-%s_tau-%d.mat', GROUP, TASK, SES, Tau)), 'f_diff');

disp('Computed frequencies from subject timeseries')

[FCemp_group, COVtauemp_group, FCPB, COVtauPB] = compute_group_empirical(TS, Tau, indexN);

v = [data_dir sprintf('/derived_data/FC_group-%s_task-%s_ses-%s_tau-%d.mat',GROUP,TASK,SES,Tau)];
save (v,'FCPB','COVtauPB')

disp('NAN values for FCemp')
[r, c] = find(isnan(FCemp_group));
for k = 1:numel(r)
    disp(fprintf('NaN at row %d, column %d\n', r(k), c(k)));
end

disp('NAN values for COVtauemp_group')
[r, c] = find(isnan(COVtauemp_group));
for k = 1:numel(r)
    disp(fprintf('NaN at row %d, column %d\n', r(k), c(k)));
end

disp('Computed fc and covtau')

% lags for HopfModel_LNA in seconds
lags_model = [Tau*TR];

% --------------- Prepare storage --------------
nG = length(Gvals);
nA = length(Avals);
ErrorMat = zeros(nG,nA);
CeffStore = cell(nG,nA); % store group-level connectivity
errFCTraj  = cell(nG,nA); % per-cell FC error trajectory
errCOVTraj = cell(nG,nA); % per-cell COVtau error trajectory
FCsimStore  = cell(nG,nA); % simulated FC at loop exit, per cell
COVsimStore = cell(nG,nA); % simulated (normalised) COVtau at loop exit

% --------------- Loop over G, A ---------------
for iG = 1:nG
  gVal = Gvals(iG);
  for iA = 1:nA
    aVal = Avals(iA);
    aVec = aVal * ones(Nvalid,1);

    disp(sprintf('g=%.6f a=%.6f', gVal, aVal))

    % === Run group-level model ===
    [Ceff_group, errFC, errCOVtau, errTotal, FCsim_g, COVtauSim_g] = fit_group_level(...
       SC, aVec, gVal, wo, sigma, ...
       FCemp_group, COVtauemp_group, ...
       epsFC, epsFCtau, max_iters, maxC, lags_model);

    % store results
    CeffStore{iG,iA}  = Ceff_group;
    errFCTraj{iG,iA}  = errFC;
    errCOVTraj{iG,iA} = errCOVtau;
    FCsimStore{iG,iA}  = FCsim_g;
    COVsimStore{iG,iA} = COVtauSim_g;
    ErrorMat(iG,iA)  = errTotal;  % final total error
  end
end

% --------------- Find best (g*, a*) -----------
[minVal, idx] = min(ErrorMat(:));  % lower = better
[iGbest, iAbest] = ind2sub(size(ErrorMat), idx);
gBest = Gvals(iGbest);
aBest = Avals(iAbest);
Ceff_best = CeffStore{iGbest, iAbest};
errFC_best     = errFCTraj{iGbest, iAbest};
errCOVtau_best = errCOVTraj{iGbest, iAbest};
FCsim_best     = FCsimStore{iGbest, iAbest};
COVtauSim_best = COVsimStore{iGbest, iAbest};

disp('====================================================');
disp(['Best group-level parameters => g=' num2str(gBest) ...
      ', a=' num2str(aBest) ', error=' num2str(minVal)]);
disp('====================================================');

% --------------- Expand Ceff_best to full size (416x416) with NaN ------
Ceff_best_full = nan(N, N);  % N=416
Ceff_best_full(indexN, indexN) = Ceff_best;  % Fill in valid nodes

% --------------- Save all results -------------
results_dir = sprintf('%s/results_group_only/', res_dir);
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

results_file = fullfile(results_dir, ...
    sprintf('%s_group-%s_task-%s_ses-%s_tau-%d.mat', ...
        ternary(GRIDON,'gridGroupResults','no_grid_GroupResults'), ...
        GROUP, TASK, SES, Tau));

save(results_file, ...
     'Gvals','Avals','ErrorMat','CeffStore', ...
     'gBest','aBest','Ceff_best','Ceff_best_full','minVal','sigma', ...
     'errFC_best','errCOVtau_best','errFCTraj','errCOVTraj', ...
     'FCsim_best','COVtauSim_best','FCemp_group','COVtauemp_group');

%--- Helper function ---%
function out = ternary(cond,trueVal,falseVal)
    if cond, out = trueVal; else out = falseVal; end
end
