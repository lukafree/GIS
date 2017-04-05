clear all
close all
clc
%%%%%%%%%%%% %%%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%%  %%%%%%%%%%%% %%%%%%%
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('                      DEMO 2: GIS ESTIMATORS') 
disp(' ')
disp('Reference: ')
disp('- L. Martino,V. Elvira, F. Louzada, G. Camps-Valls,') 
disp('Group Importance Sampling for particle filtering and MCMC, 2017.') 
disp(' ')
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
pause(1)
disp(' ')
cprintf('blue','       Group Importance Sampling (GIS)')
disp('   ')
disp('   ')
cprintf('blue','GIS: weigthing properly clouds of weighted samples providing consistent IS estimators.')
disp('   ')
cprintf('blue','     ')
disp('   ')
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp(' ')
cprintf('blue',    'Press a key\n' );
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Definition of the target pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig=1;
target=@(x,mu) 1/sqrt(2*pi*sig^2).*exp(-(x-mu).^2/(2*sig^2));
muT=-5; %%% mean of the target 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% number of samples per groups and number of groups
N=[2 200 20];
M=length(N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
fprintf('TARGET FUNCTION:\n')
disp(' ')
fprintf('1/sqrt(2*pi*sig^2).*exp(-(x-mu).^2/(2*sig^2)) \n')
fprintf('with mu =%d\n',muT)
fprintf('and sig =%d\n',sig)
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%  
disp(' ')
fprintf('Number of Groups =%d',length(N))
disp(' ')
disp(' ')
aux_integer=[[1:length(N)]',N'];
fprintf('Number of particles of the group no %d = %6.2f\n',aux_integer')
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
%%%%% proposal parameters
sp=[4 4 1];
mup=[0 0 2];
%%%%%
aux_integer=[[1:length(N)]',mup',sp'];
disp(' ')
fprintf('Parameters of the Gaussian proposal associated to the group no %d,  mean = %6.2f  std =%6.2f\n',aux_integer')  
disp(' ')
disp('(All these values can be changed in the code)') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
N_res_per_group=100; %%% Number of  particles resampled per group
disp(' ')
fprintf('Number of  resampled particles per group = %d \n',N_res_per_group)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp(' ')
fprintf(' Estimation of the expected value of the target')
disp(' ')
disp(' ')
fprintf('TRUE VALUE  = %d \n',muT)
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SIMU=1000;
%%%%%  
  for j=1:SIMU
 %%%%%%%%% for each group %%%
  for m=1:M
   %%%% generate samples 
    xp{m}=mup(m)+sp(m)*randn(1,N(m));
   %%%% evaluate proposal and target
    P=(1/sqrt(2*pi*sp(m).^2)).*exp(-(xp{m}- mup(m)).^2./(2*sp(m).^2));    
    T=target(xp{m},muT);
   %%%% compute importance weights of each sample  
    w{m}=T./P;    
    wn{m}=w{m}./sum(w{m});
   %%% marginal likelihood estimator 
    Z(m)=mean(w{m});
    xest(m)=sum(wn{m}.*xp{m});
    %%%% resample N_res_per_group particles per group %%%%
    xr(m,:)=randsrc(1,N_res_per_group,[xp{m}; wn{m}]);
  end %%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% approximanting moments %%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% First moment %%%%
    Est1(j)=mean(mean(xr));  %%%% Wrong estimator
    Est2(j)=1/(N_res_per_group)*sum(sum(repmat((N'.*Z')/sum(N'.*Z'),1,N_res_per_group).*xr)); %%%% Proper estimator
    %%% Second moment %%%%
    Est3(j)=mean(mean(xr.^2));  %%%% Wrong estimator
    Est4(j)=1/(N_res_per_group)*sum(sum(repmat((N'.*Z')/sum(N'.*Z'),1,N_res_per_group).*xr.^2)); %%%% Proper estimator
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  end %%%% end simulations (different indepedent runs)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% showing Results %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf('(results averaged over %d independent simulations) ',SIMU)
  disp(' ')
  disp(' ')
  disp('---------------------------------------------------------------------------')
  disp('Wrong estimator - without GIS weighting')
  disp(' ')
  disp(['Wrong estimation = ',num2str(mean(Est1))])
  disp(['MSE  = ',num2str(mean((Est1-muT).^2))])
  disp(['BIAS = ',num2str(mean(Est1)-muT)])
  disp(['VAR  = ',num2str(var(Est1,1))])
  disp(['(test: BIAS^2+VAR  = ',num2str(var(Est1,1)+(mean(Est1)-muT)^2),')'])
  disp('----------------------------------------------------------------------------') 
  disp(' ')
  disp('Proper GIS estimator using the proper weighting of the resampled particles')
  disp(' ')
  disp(['GIS estimation = ',num2str(mean(Est2))])
  disp(['MSE  = ',num2str(mean((Est2-muT).^2))])
  disp(['BIAS = ',num2str(mean(Est2)-muT)])
  disp(['VAR  = ',num2str(var(Est2,1))])
   disp(['(test: BIAS^2+VAR  = ',num2str(var(Est2,1)+(mean(Est2)-muT)^2),')']) 
  disp('----------------------------------------------------------------------------') 
cprintf('blue',    'Press a key\n' );
%%%%%%%%%
pause
%%%%%%%%%
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp(' ')
fprintf(' Estimation of the second moment of the target') 
disp(' ')
disp(' ')
fprintf('TRUE VALUE = %d \n',muT^2+sig^2)
TRUE_value_secMom=muT^2+sig^2;
disp(' ')
fprintf('(results averaged over %d independent simulations) ',SIMU)
disp(' ')
disp(' ')
disp('---------------------------------------------------------------------------')
disp('Wrong estimator - without GIS weighting')
disp(' ')
disp(['Wrong estimation = ',num2str(mean(Est3))])
disp(['MSE  = ',num2str(mean((Est3-TRUE_value_secMom).^2))])
disp(['BIAS = ',num2str(mean(Est3)-TRUE_value_secMom)])
disp(['VAR  = ',num2str(var(Est3,1))])
disp(['(test: BIAS^2+VAR  = ',num2str(var(Est3,1)+(mean(Est3)-TRUE_value_secMom)^2),')'])
disp('----------------------------------------------------------------------------') 
disp(' ')
disp('Proper GIS estimator using the proper weighting of the resampled particles')
disp(' ')
disp(['GIS estimation = ',num2str(mean(Est4))])
disp(['MSE  = ',num2str(mean((Est4-TRUE_value_secMom).^2))])
disp(['BIAS = ',num2str(mean(Est4)-TRUE_value_secMom)])
disp(['VAR  = ',num2str(var(Est4,1))])
disp(['(test: BIAS^2+VAR  = ',num2str(var(Est4,1)+(mean(Est4)-TRUE_value_secMom)^2),')']) 
disp('----------------------------------------------------------------------------')
