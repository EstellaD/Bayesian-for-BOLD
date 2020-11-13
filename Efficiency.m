% %%%%%%%%%%%% Efficiency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % I will explore the stimulation time interval of the Experiment Design.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Event Generation (alphas, periods, jitter)                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fab = linspace(0,3,361); %%%%%%%%%%%%%%%%%%%%%%% Build Alphas %%%%%%%%%%%%%
index=[21 41 61 81 101 121 141 161 181 201 221 241 261 281 301 321 341 361];
alphas = zeros(1,length(index));
for ii = 1:length(index)
    alphas(ii) = fab(index(ii));
end
n_alphas = length(alphas); %18

models = alphas;         %%%%%%%%%%%%%%%%%%%%%%% Build Models & Ymodels%%%%
n_models = length(models); 
% Here we want a data set for each value of alpha
n_Ymodels = n_models; % In my case they are the same, cuz only one beta. 

starting_time = 1.0;
total_time = 2400.0; %40 mins fMRI experiment. % Build Periods & Events%%%%
periods_list = 1:60; % Realistic?
num_periods = length(periods_list);% 60
total_time_vec = total_time*ones(1,num_periods);
events_list = floor(total_time_vec./periods_list); % [2400 1200 600 300 240...]  


jittermags = [0 1 2];    %%%%%%%%%%%%%%%%%%%%%%% Build Effforms %%%%%%%%%%%
n_jittermags = length(jittermags);
effforms = cell(1,n_jittermags);
for ii = 1:n_jittermags
effforms{ii} = zeros(num_periods,n_alphas);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Enhance Resolution Constants                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt= 0.1;
n_scans = 800;
TR = 3; % the repitition time is 3 seconds chosen.
n_res = n_scans*TR / dt; %24000

% Compute the HRF for a temporal resolution of 1/10 seconds.%%
hrf = spm_hrf(0.1); %by default
n_conv = floor(n_res) + length(hrf) - 1; 
   
% time vectors at high resolution and with TR onsets
t_highRes = 0:dt:(n_conv - 1)*dt; %%%%1x12320: 0,0.1,0.2,..., 1231.9
t_lowRes = 0:TR:(n_scans - 1)*TR; % 1x800: 0,3,6,...,1197(399*3) = time (s) of volume (TR = 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Many for-loops                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for v = 1:n_jittermags % Start with one efficiency matrix
    
for l = 1: num_periods % Choose a period thus a n_events./ build jitter 
    n_events = events_list(l); 
events = randi(10,1, n_events); 
rng(24);% To replicate
jitter = 0.2*periods_list(l)*abs(rand(1, n_events)); % set the max jitter to be 0.2 of the period.
t = linspace(starting_time, total_time + starting_time, n_events);
event_onset = t + jittermags(v)*jitter; % to make non-equi.dist.
length(event_onset); % = n_events

% Initialising the alpha columns in a matrix for function.
%Neurometric Power Law
func = bsxfun(@power, events(:), alphas);

%%%%%%%%%%%%%%%%%%%%%%%Matrices for high resolution result%%%%%%%%%%%%%%%%
% Initialize metrices for high res convolution:
func_highRes = zeros(floor(n_res), n_alphas);
events_highRes = zeros(floor(n_res), 1);

% jittered
for i = 1:n_events
    for j = 1: n_alphas
    func_highRes(round(event_onset(i)/dt), j) = func(i, j);
    events_highRes(round(event_onset(i)/dt)) = events(i);
    end
end

% Cut-off the surplus rows
func_highRes = func_highRes(1:floor(n_res), :);
events_highRes = events_highRes(1:floor(n_res), :);

%%%%%%%%%%%%%%%% CONVOLVING WITH HRF AND INTERPOLATED %%%%%%%%%%%%%%%%%%%%%
% Initialize matrices for storing convoluted values
func_highRes_Conv = zeros(n_conv, n_alphas); % # x 7 (here)
% For the interpolated values we need matrices of the same size as the
% number of scans
func_inter = zeros(n_scans, n_alphas);
% Convolution of delta functions and hrf, and interpolation down to low res.

for i = 1:n_alphas
% Probabilities
    func_highRes_Conv(:, i) = conv(hrf, func_highRes(:,i)); 
        %convolve each of the 3 high-resolution columns with the HRF.
        % Each column having events in the size of the choice-vals.
    func_inter(:, i) = interp1(t_highRes, func_highRes_Conv(:, i), t_lowRes); %interpolate to low resolution 
end

% Convolving the event stick functions.
events_highRes_Conv = conv(hrf, events_highRes); % Only vector for the events
events_inter = interp1(t_highRes, events_highRes_Conv, t_lowRes); % Use high resolution convolution 
% to interpolate time points of scans for events.

% Standardize the regressors %length = 800 = n_scans 
std_events_inter = zscore(events_inter); % null model  %length = 800
std_func_inter = zscore(func_inter); %standardized interpolation for 10 different alphas. ALL AT ONCE WOW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE DATA, COMPUTE LIKELIHOOD AND INVESTIGE THE DETECTION OF ALPHA  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE DESIGN (X) MATRICES, Y = X*B
% Create design matrices for each alpha.
% initialize the cell that holds the matrices.                         
matrices = cell(1, n_models); 

% Insert the null-model (as the first model)
matrices{1} = std_events_inter'; % means

% Insert all matrices that regards alpha.
for jj = 2:n_models
    matrices{jj} = [std_events_inter', std_func_inter(:, jj - 1)]; 
end

% Beta values
n_betas  = 1;
betas = [1; 1];  
% 2 * 1 matrix, one row for multiplication with mean, 
% the other row for beta value, ground truth and no scaling effect here.

Y = []; % zero(n_scans, n_Ymodel*n_beta) 800*8
 
Y(:, 1) = matrices{1}*betas(1, :); % The model with mean regressor
% matrices{1} = std_event_inter' * all the ones in betas 800*20

% The rest are constructed in a loop over alphas, such that all betas are
% done at once.
for ii = 2:n_alphas % Doing the alphas separately.
    Y = [Y, matrices{ii}*betas];
end

%%%%%%%%%%%% DERIVE AND COLLECT RESULTS (LIKELIHOOD VALUES) %%%%%%%%%%%%%%%
% Set seed for random number generator (for replication).
rng(31);
n_sims = 500;
sim_data = zeros(n_scans, n_Ymodels, n_sims); %3-D Model 800*361*50
sim_likelihoods = zeros(n_Ymodels* n_betas, n_models, n_sims); % 361*361*50 
sim_alphas_estimates = zeros(1, n_Ymodels, n_sims);% 1*361*50 
sim_evidence = zeros(n_Ymodels * n_betas, n_sims);% 361*50
    % For each simulation, we have P(alpha_MAP |Y) for 
    % n_Ymodels * n_betas = 361 data sets. -> For Bayes Factors
    
for kk = 1:n_sims
    % Adding noise to the data, this is the error in Y = XB+E
    sim_data(:, :, kk) = Y + normrnd(0, 1, n_scans, n_Ymodels);
    % error here is a plane 400*361, here we are generating a whole x-y plane for repeat round 1.

    %%% COMPUTING LIKELIHOOD
    % Compute likelihood:
    sim_likelihoods(:, :, kk) = log_likelihoodbyline(matrices, sim_data(:, :, kk));
    [~, min_ind] = min(sim_likelihoods(:, :, kk), [], 2); % It returns the min in each row in a column. % min(A, [], 1) is for column.
    min_ind = reshape(min_ind, 1, n_models); % reshape 18*1 into 1*18
    % In this case you can just invert.
    
    % Estimating alpha from here. 
    sim_alphas_estimates(:, :, kk) = models(min_ind);
end

%% COMPUTING MEAN ACROSS THE SIMULATIONS:
% Mean estimates (for all models!!)
mean_alphas_estimates = mean(sim_alphas_estimates, 3); %3 here is for the third dimension.
twodsimalphas = reshape(sim_alphas_estimates,n_sims,n_models); 
%% COMPUTING EFF ACROSS THE SIMULATIONS:
% Mean estimates (for all models!!)

for hh = 1:n_models
    sum_var = 0;
    for mm = 1:n_sims
        sum_var = sum_var + (twodsimalphas(mm,hh) - alphas(hh))^2;
    end
    expe = sum_var/n_sims;
    eff = 1 / expe;
    effforms{v}(l,hh) = eff;
end
end
end
%% PLOT OF THE ESTIMATES, Posterior Distribution of Sigma
% This plot is a heat map, where the axes specifies the combination of
% alpha and beta that was used to generate the data, while the color
% specifiees the estimated alphas.

o= effforms{1};
p = effforms{2};
q = effforms{3};
w = [o,p,q];
figure;  
cmap = colormap(parula(n_models)); 
colormap(cmap);
imagesc(alphas, periods_list, w); 

% should be alpha here, since i am ploting against alpha at any point.
cbar = colorbar('Ticks', alphas, ...
                'TickLabels', alphas, ...
                'Location', 'southoutside');
xlabel(cbar, 'Efficiency')
ax = gca;
ax.XAxisLocation = 'top';
xlabel('\alpha');
ylabel('Periods');
% Create title
title({'no jitter                                                              jitter                                                               2*jitter'});
