

%% 2D Problems
[A, b, x_true, ProbInfo] = PRtomo;
%[A, b, x_true, ProbInfo] = PRspherical;
%[A, b, x_true, ProbInfo] = PRseismic;


%% Adds Noise
rng(0);
bn = PRnoise(b, 10^(-2)); %adds noise
x0 = zeros(size(A,2),1); %initial guess

%% Run solvers
max_iter=50;
[x_approx_lslu_opt,output_lslu_opt] = hybrid_LSLU(A,bn,x0, x_true, max_iter,'optimal_tik');
[x_approx_lslu_wgcv,output_lslu_wgcv] = hybrid_LSLU(A,bn,x0, x_true, max_iter,'wgcv');

options_lsqr.MaxIter=max_iter;
options_lsqr.RegParam = 'wgcv';
options_lsqr.x0=x0;
options_lsqr.x_true = x_true;
options_lsqr.NoStop = 'on';
[x_approx_lqr_wgcv,output_lsqr_wgcv] = IRhybrid_lsqr(A,bn,options_lsqr);

%% Plot relative error histories
figure
plot([1;output_lslu_opt.Enrm],'LineWidth',5)
hold on
plot([1;output_lslu_wgcv.Enrm],'LineWidth',5)
plot([1;output_lsqr_wgcv.Enrm],'LineWidth',5)
legend('H-LSLU opt','H-LSLU wgcv','H-LSQR wgcv')
set(gca,'FontSize',30)
