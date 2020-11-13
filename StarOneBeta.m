%%%%%%%%%%%%%%%%%%%%%% Event Generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
events = randi(10,1,100);
      %generate 85 stimulations, strongest with magnitude of 10
n_events = length(events);

%%% Create vector of corresponding time points 
% Equally spaces time points between 0 and 2400 seconds
t = linspace(0.0, 2400.0, n_events); 
    % length(t) = 85.

% Create jitter to make t NOT equidistant
rng(24); % To replicate
jitter = abs(4*randn(1, n_events)); 

% Decision onset
event_onset = t + jitter; % jitter is added just to make non-equi.dist.
length(event_onset); %85

%%%%%%%%%%%%%%%%%%%%%%%%% Set Alphas and  %%%%%%%%%%
alphas = linspace(0,3,361); 
n_alphas = length(alphas); %361

% Initialising the alpha columns in a matrix for function.
%Neurometric Power Law
%%func = zeros(n_events, n_alphas); % #scans x #alphas. X-Matrix
func = bsxfun(@power, events(:), alphas);
%stem(event_onset, func(:, 2)); % to see the plot

%%%%%%%%%%%%%%%%%%%%%%%Matrices for high resolution result%%%%%%%%%%%%%%%%%
dt= 0.1;
n_scans = 800;
TR = 3; % the repitition time is 3 seconds chosen.
n_res = n_scans*TR / dt; %12000

% Initialize metrices for high res convolution:
func_highRes = zeros(n_res, n_alphas);
events_highRes = zeros(floor(n_res), 1);

for i = 1:n_events
    for j = 1: n_alphas
    func_highRes(round(event_onset(i)/dt), j) = func(i, j);
    events_highRes(round(event_onset(i)/dt)) = events(i);
    end
end

% Cut-off the surplus rows
func_highRes = func_highRes(1:floor(n_res), :);
events_highRes = events_highRes(1:floor(n_res), :);

%%%%%%%%%%%%%%%% CONVOLVING WITH HRF AND INTERPOLATED %%%%%%%%%%%%
% Compute the HRF for a temporal resolution of 1/10 seconds.
hrf = spm_hrf(0.1); %by default
plot(hrf);
n_conv = floor(n_res) + length(hrf) - 1; %length of convolution = 12320.
    % The convolution of u and v has length(u)12000 + length(v)321 - 1 elements.
    % (See matlab documentation
    
% Initialize matrices for storing convoluted values
func_highRes_Conv = zeros(n_conv, n_alphas); % 12320 x 7 (here)

% time vectors at high resolution and with TR onsets
t_highRes = 0:dt:(n_conv - 1)*dt; %1x12320: 0,0.1,0.2,..., 1231.9
t_lowRes = 0:TR:(n_scans - 1)*TR; % 1x400: 0,3,6,...,1197(399*3) = time (s) of volume (TR = 3)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING THE STICK FUNCTIONS AND CONVOLUTIONS
%% Plot the indicator function for the events and their corresponding
%% convolution with the hrf. 
%Note: the plots are NOT of the standardized regressors!!!
% This is the main effect regressor
figure(1);
clf;                % Clear the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,1);     % This is just to make the plots line up prettily
stem(event_onset, events, ...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532], ...
    'MarkerSize',3, ...
    'Color',[0 0.447058826684952 0.74117648601532]);

title('Time-series of events');
axis([-10 1200 0 11]); %the end time of t_lowRes is 1197
xlabel('Time (seconds)');
ylabel('Stimulus magnitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2);
plot(t_lowRes,events_inter);
title('Stimulus time-series convolved with HRF');
axis([-10 1200 -0.025 0.25]);
xlabel('Time (seconds)');
ylabel('fMRI signal');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE DATA, COMPUTE LIKELIHOOD AND INVESTIGE THE DETECTION OF ALPHA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE DESIGN (X) MATRICES, REGRESSION COEFFICIENTS AND ERROR TERM 
% Define the models from their respective alphas.
% The 1 thus corresponds to the null-model here.
models = alphas; % 361 alphas.
n_models = length(models); % 361
k0 = find(models== 1/3); %index 41
k = find(models==0.5); %index 61
k2 = find(models==1); %index 121
k3 = find(models==2); % index 241
k4 = find(models==3); %index 361

% Number of data matrices, i.e. how many of the alphas will we use to generate data 
n_Ymodels = n_models; % Here we want a data set for each value of alpha

% Create design matrices for each alpha and the null model.
% initialize the cell that holds the matrices.                         
matrices = cell(1, n_models); 

% Insert the null-model (as the first model)
matrices{1} = std_events_inter'; % means

% Insert all matrices that regards alpha.
for jj = 2:n_models
    matrices{jj} = [std_events_inter', std_func_inter(:, jj - 1)]; 
end

% Beta values
n_betas  = 1; % 1 number.
betas = [1; 1];  
% 2 * 1 matrix, one row for multiplication with mean, 
% the other row for beta value, ground truth and no scaling effect here.

Y = []; % zero(n_scans, n_Ymodel*n_beta) 400*8
 
Y(:, 1) = matrices{1}*betas(1, :); % The model with mean regressor
% matrices{1} = std_event_inter' * all the ones in betas 400*20

% The rest are constructed in a loop over alphas, such that all betas are
% done at once.
for ii = 2:n_alphas % Doing the alphas separately.
    Y = [Y, matrices{ii}*betas];
end

%%%%%%%%%%%% DERIVE AND COLLECT RESULTS (LIKELIHOOD VALUES) %%%%%%%%%%%%%%%
% Set seed for random number generator (for replication).
rng(31);
n_sims = 1000;
sim_data = zeros(n_scans, n_Ymodels, n_sims); %3-D Model 400*361*50
sim_likelihoods = zeros(n_Ymodels* n_betas, n_models, n_sims); % 361*361*50 
sim_alphas_estimates = zeros(1, n_Ymodels, n_sims);% 1*361*50 
% sim_evidence = zeros(n_Ymodels * n_betas, n_sims);% 361*50
    % For each simulation, we have P(alpha_MAP |Y) for 
    % n_Ymodels * n_betas = 361 data sets.
    
for kk = 1:n_sims
    % Adding noise to the data, this is the error in Y = XB+E
    sim_data(:, :, kk) = Y + normrnd(0, 1, n_scans, n_Ymodels); 
    % error here is a plane 400*8, here we are generating a whole x-y plane for repeat round 1.

    %%% COMPUTING LIKELIHOOD
    % Compute likelihood:
    sim_likelihoods(:, :, kk) = log_likelihoodbyline(matrices, sim_data(:, :, kk));
    %likelihood_values = log_likelihoodbyline(matrices, Y); Output is of n_voxels (n_Ymodels*n_betas) x n_models

    % Collect the index of the minimum per each data series (column in Y), and
    % reshape it into a n_models x n_betas)
 
    % return the min element in each row, so dimension 8*1
    % We also collect the model evidence (the likelihood P(alpha_MAP |Y)) for the same estimate (i.e. the
    % probability), since this is used when computing Bayes factors later!!!
    [evidence_MAP, min_ind] = min(sim_likelihoods(:, :, kk), [], 2); % It returns the min in each row in a column. % min(A, [], 1) is for column.
    min_ind = reshape(min_ind, 1, n_models); % reshape 8*1 into 1*8 
    % In this case you can just invert.
    
    % Estimating alpha from here. 
    sim_alphas_estimates(:, :, kk) = models(min_ind);
%     sim_evidence(:, kk) = evidence_MAP;
end

%% COMPUTING MEAN ACROSS THE SIMULATIONS:
% Reshape data into 4-D matrix
sim_data = reshape(sim_data, n_scans, n_betas, n_Ymodels, n_sims);
% Mean across simulations:
mean_data = mean(sim_data, 4);

% Mean likelihood "functions"
mean_likelihoods = mean(sim_likelihoods, 3); % 1 columnwise, 2 rowwise, 3 zwise
    
% Mean estimates (for all models!!)
mean_alphas_estimates = mean(sim_alphas_estimates, 3); %3 here is for the third dimension.
    
% Mean evidence for the MAP-estimates:
% mean_evidence = mean(sim_evidence, 2); % -> to plot MAP; mean(,1) is columnwise & mean(,2) is rowwise; so get 361 number here.

%%%%%%%%%%%%%%%%%%%%%%%%Plot the errors for alpha estimated%%%%%%%%%%%%%%%
figure
% Alphas vs. All Simulation Alphas
subplot(1,3,1);
ReshapedAEsti = reshape(sim_alphas_estimates, n_models, n_sims);
Reshapedforplot = ReshapedAEsti';
plot(alphas, Reshapedforplot,'.');
axis([-0.5,3.5,-0.5,3.5]);
xlabel('Alpha values')
ylabel('Estimated Alpha values')
title('Alpha parameter recovery')


%Alphas vs. Alpha^
subplot(1,3,2);
plot(alphas, mean_alphas_estimates, 'o');
hold on
axis([-0.5,3.5,-0.5,3.5]);
lsline;
xlabel('Alpha values')
ylabel('Estimated Alpha values')
title('Alpha parameter recovery')


% Only the 5 alpha values are of interest
subplot(1,3,3);
plot([alphas(k0) alphas(k) alphas(k2) alphas(k3) alphas(k4)], [mean_alphas_estimates(:, k0) mean_alphas_estimates(k) mean_alphas_estimates(k2) mean_alphas_estimates(k3-1) mean_alphas_estimates(:, k4)], 'o');
axis([-0.5,3.5,-0.5,3.5]);
lsline;
hold on
plot([alphas(k0) alphas(k) alphas(k2) alphas(k3) alphas(k4)], [Reshapedforplot(:, k0) Reshapedforplot(:, k) Reshapedforplot(:, k2) Reshapedforplot(:, k3-1) Reshapedforplot(:, k4)],'.');
xlabel('Alpha values for 0.33 0.5 1 2 3')
ylabel('Estimated Alpha values for 0.33 0.5 1 2 3')
title('Alpha parameter recovery for 0.33 0.5 1 2 3')


% Neurometric Power Law 
K = 1;
X = vec2mat(events,1);
Y0 = X.^mean_alphas_estimates(k0);
Y1 = X.^mean_alphas_estimates(k);
Y2 = X.^mean_alphas_estimates(k2);
Y3 = X.^mean_alphas_estimates(k3);
Y4 = X.^mean_alphas_estimates(k4);

f0 = fit(X, Y0, 'power1');
f1 = fit(X, Y1, 'power1');
f2 = fit(X, Y2, 'power1');
f3 = fit(X, Y3, 'power1');  
f4 = fit(X, Y4, 'power1'); 

plot (f0, X, Y0,'*'); 
hold on
plot (f1, X, Y1,'o'); 
hold on
plot (f2, X, Y2,'d'); 
hold on
plot (f3, X, Y3); 
hold on
plot (f4, X, Y4); 
hold on
axis ([1,4,1,4]);
hold off
lgd = legend('al^pha = 0.33','al^pha = 2','al^pha = 0.5', 'al^pha = 3', 'al^pha = 1'); 
title('Neurometric Power Law')
xlabel('Events')
ylabel('f(x) = K*x^\alpha')

%  %Likelihood values for different alphas -> Not necessary
% c = categorical({'Alpha = 0.5','Alpha = 1','Alpha = 2'});
% Likeli = mean_evidence;
% Likelisubset = [Likeli(k) Likeli(122) Likeli(242)];
% bar(c,Likelisubset)
% hold on
% ylim([1200,1250]);
% xlabel('Predicted Alphas')
% title('Likelihood values for different alphas');
% error_over_repetition = std(sim_evidence, 0, 2); % calculate the standard deviation along each row.
% error_subset = [error_over_repetition(k) error_over_repetition(k2) error_over_repetition(k3)];
% errorbar(Likelisubset,error_subset)

%alpha pred with standard deviaion as error bar 
figure;
subplot(2,3,1);
c = categorical({'Alpha = 0.33','Alpha = 0.5','Alpha = 1','Alpha = 2', 'Alpha = 3'});
stdAlphasEsti = std (sim_alphas_estimates,0,3);
bar(c, [mean_alphas_estimates(k0), mean_alphas_estimates(k),mean_alphas_estimates(k2),mean_alphas_estimates(k3), mean_alphas_estimates(k4)]);
hold on
xlabel('Predicted Alphas')
title('Alpha Pred with Standard Deviaion as Error Bar')
errorbar([mean_alphas_estimates(k0), mean_alphas_estimates(k),mean_alphas_estimates(k2),mean_alphas_estimates(k3), mean_alphas_estimates(k4)],...
    [stdAlphasEsti(k0), stdAlphasEsti(k),stdAlphasEsti(k2),stdAlphasEsti(k3), stdAlphasEsti(k4)]); 

m = max(-mean_likelihoods, [], 2);
likelihood = exp(bsxfun(@plus,-mean_likelihoods, -m));
act_lik_ap = bsxfun(@times, likelihood, 1./sum(likelihood,2));

%MAP for alpha = 0.33 
subplot(2,3,2);
[pks0,locs0] = findpeaks(act_lik_ap(:, k0),alphas);
[~, max_ind0] = max(pks0);
Alpha_pred0 = locs0(max_ind0);
pks0 = max(pks0);locs0 = Alpha_pred0;
plot(alphas,act_lik_ap(:, k0));
hold on 
text(locs0, pks0, sprintf('AlphaPredict = %6.3f', locs0))
hold on
xlabel('\alpha')
ylabel('MAP: P(\alpha|Y_{\alpha=0.33})')
title('MAP distribution for alpha = 0.33');


%MAP for alpha = 0.5
subplot(2,3,3);
[pks,locs] = findpeaks(act_lik_ap(:, k),alphas);
[~, max_ind] = max(pks);
Alpha_pred = locs(max_ind);
pks = max(pks);locs = Alpha_pred;
plot(alphas,act_lik_ap(:, k));
hold on 
text(locs, pks, sprintf('AlphaPredict = %6.3f', locs))
hold on
xlabel('\alpha')
ylabel('MAP: P(\alpha|Y_{\alpha=0.5})')
title('MAP distribution for alpha = 0.5');

%MAP for alpha = 1
subplot(2,3,4);
[pks2,locs2] = findpeaks(act_lik_ap(:, k2),alphas);
[~, max_ind2] = max(pks2);
Alpha_pred2 = locs2(max_ind2);
pks2 = max(pks2);locs2 = Alpha_pred2;
plot(alphas,act_lik_ap(:, k2));
hold on 
text(locs2, pks2, sprintf('AlphaPredict = %6.3f', locs2))
hold on
xlabel('\alpha')
ylabel('MAP: P(\alpha|Y_{\alpha=1})')
title('MAP distribution for alpha = 1');

%MAP for alpha = 2
subplot(2,3,5);
[pks3,locs3] = findpeaks(act_lik_ap(:, k3),alphas);
[~, max_ind3] = max(pks3);
Alpha_pred3 = locs3(max_ind3);
pks3 = max(pks3);locs3 = Alpha_pred3;
plot(alphas,act_lik_ap(:, k3));
hold on 
text(locs3, pks3, sprintf('AlphaPredict = %6.3f', locs3))
hold on
xlabel('\alpha')
ylabel('MAP: P(\alpha|Y_{\alpha=2})')
title('MAP distribution for alpha = 2');

%MAP for alpha = 3
subplot(2,3,6);
[pks4,locs4] = findpeaks(act_lik_ap(:, k4),alphas);
[~, max_ind4] = max(pks4);
Alpha_pred4 = locs4(max_ind4);
pks4 = max(pks4);locs4 = Alpha_pred4;
plot(alphas,act_lik_ap(:, k4));
hold on 
text(locs4, pks4, sprintf('AlphaPredict = %6.3f', locs4))
hold on
xlabel('\alpha')
ylabel('MAP: P(\alpha|Y_{\alpha=3})')
title('MAP distribution for alpha = 3');