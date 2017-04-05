clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('       DEMO 1:  Estimators of the marginal likelihood in Sequential Importance Resampling (SIR)') 
disp(' ')
disp('Reference: ')
disp('- L. Martino,V. Elvira, F. Louzada, G. Camps-Valls,') 
disp('Group Importance Sampling for particle filtering and MCMC, 2017.') 
disp(' ')
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
pause(1)
disp(' ')
disp('   ')
cprintf('blue','Using GIS, two estimators are possible in SIR (as in SIS), and they are equivalent ')
disp('   ')
disp('   ')
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp(' ')
cprintf('blue', 'Press a key' );
pause
%%%%%%%%%%%%%
disp('   ')
disp('   ')
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('   ')
disp('Choose a Wrong or a Right (GIS) unnormalized weights after Resampling '); 
disp(' ')
disp(' ')
cprintf('blue', 'Type 1 for GIS Weighting (unnormalized weights =Z)\n'); 
cprintf('blue', 'Type 0 for Wrong Weighting (unnormalized weights =1)\n'); 
cprintf('blue', 'Type -1 for Random Wrong Weighting\n'); 
disp(' ')
disp(' ')
RightRes = input('Choose Weighting; type -1 or 0 or 1? ');
disp(' ')
%%%%%%%%%%%%
if  RightRes==1 
disp('GIS WEIGHTING')
elseif RightRes==-1
    disp('RANDOM WRONG WEIGHTING')
elseif RightRes==0
    disp('WRONG WEIGHTING') 
else
    disp('....we consider that you choose....a (right) GIS WEIGHTING...')
    pause(0.3)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('   ')
disp('   ')
disp('Choose the threshold parameter for adaptive-resampling, between 0 and 1, ');
disp('We are using the formula  ESS=1./max(wn(i,:))')
disp('but you can change for the standard formula: ESS=1./sum(wn(i,:).^2);');
disp('   ')
disp('For adaptive resampling, we suggest [0.01,0.1] with the formula ESS=1./max(wn(i,:))');
disp(' ')
ep=input('type a number between 0 and 1 (zero resampling=0; always resampling =1) ?');
%ep=1; %%%% Always Resampling 
%ep=0; %%%% NO Resampling 
%ep=0.01; %%%% adaptative resampling
disp(' ')
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% control %%%%
if ep<0 | ep>1
cprintf('err','MUST BE BETWEEN 0 AND 1\n')
disp('we set the threshold=1 (always-resampling)')
 ep=1;   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ep==1
    disp('YOU HAVE CHOSEN "ALWAYS RESAMPLING",i.e., boostrap PF ')
elseif ep==0
     disp('YOU HAVE CHOSEN "NO RESAMPLING", i.e., SIS scenario')
else
    disp('YOU HAVE CHOSEN "ADAPTIVE-RESAMPLING"')
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%PIECE OF TARGET %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig=1;
f_piece=@(x,mu) exp(-(x-mu).^2/(2*sig^2));
%%%%%%%%%%% play the role of measurements %%%%%%%%%%%%%%%%%%%%%%%%%%
%mu_d=[1 -1 3 0 -3 -2];
mu_d=[1 -1 3 0 -3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DIM=length(mu_d);
disp('-----------------------------------------------------------------------------------')
pause(0.2)
fprintf('We consider a target formed by D=%d pieces', DIM);
disp(' ')
disp(' ')
aux_integer=[[1:DIM]',mu_d'];
disp('f(x)=\prod_{d=1}^{D} exp(-(x-mu_d).^2/(2*sig^2))')
disp(' ')
fprintf(' where mu_%d = %6.2f\n',aux_integer')
disp('(these values can be changed directly in the code)')
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SIS  and SIR %%%% 
disp('-----------------------------------------------------------------------------------')
N=50000; %% number of particles 
fprintf(' we are using N=%d number of particles\n', N)
disp('-----------------------------------------------------------------------------------')
%%%%%%
x_p=zeros(1,N); 
w(1,:)=ones(1,N);
beta(1,:)=ones(1,N);
wn(1,:)=1/N*ones(1,N);
Z1(1)=1;
Z2(1)=1;
VALID_only_WITHOUT_RES(1)=1;
VALID_only_WITH_ALWAYS_RES(1)=1;
Z2alt(1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countRES=0; %%%% resampling counter  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:length(mu_d)+1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% particle generation %%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sp=2;
    x_p(i,:)=x_p(i-1,:)+sp*randn(1,N);
  %%% evaluate proposal and target  
    P(i,:)=(1/sqrt(2*pi*sp^2)).*exp(-(x_p(i,:)-x_p(i-1,:)).^2/(2*sp.^2));
    T(i,:)=f_piece(x_p(i,:),mu_d(i-1));
 %%%% Recursion for the weights
    beta(i,:)=T(i,:)./P(i,:);
    w(i,:)=w(i-1,:).*T(i,:)./P(i,:);
 %%%% normalized weight    
    wn(i,:)=w(i,:)./sum(w(i,:));
 %%%% (wrong and right) marginal likelihood estimators    
    Z1(i)=mean(w(i,:));
    VALID_only_WITHOUT_RES(i)=mean(prod(beta(1:end,:)));
    VALID_only_WITH_ALWAYS_RES(i)=prod(mean(beta(1:end,:),2));
    Z2(i)=prod(sum(wn(1:end-1,:).*beta(2:end,:),2));
    %%%  Z2alt(i)=Z2alt(i-1).*sum(wn(i-1,:).*beta(i,:)); %% alternative recursve formulation of Z2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Effective Sample Size (ESS) approx
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ESS=1./max(wn(i,:));
    %ESS=1./sum(wn(i,:).^2);
   %%%%%%%%%%%%%%%%%%
   %%% RESAMPLING
   %%%%%%%%%%%%%%%%%%  
  if ESS <=ep*N
     %%%   
      countRES=countRES+1;
     disp (['Resampling at ',[num2str(i-1)],'-th iteration, please wait....'])   
     x_p(i,:)=randsrc(1,N,[x_p(i,:); wn(i,:)]);
   %%%%%%%%%%%%%%%%%%
   switch RightRes
       case 1  
          disp ('Right GIS unnormalized weights after Resampling (w=Z)')
           w(i,:)=Z1(i).*ones(1,N);
      case 0
         disp ('Wrong unnormalized weights after Resampling (w=1)')
          w(i,:)=1.*ones(1,N);
       case  -1
         disp ('Random unnormalized weights after Resampling')
          w(i,:)=10*rand(1,1).*ones(1,N); 
       otherwise
       disp ('We consider that you choose the GIS weighting after Resampling (w=Z)')
      w(i,:)=Z1(i).*ones(1,N);
   end %%% switch case
   %%%normalized weight (since it is not a partial resampling they are always equal) 
     wn(i,:)=1/N.*ones(1,N);
  %%%%%%%%%%%%%%
  else%%% no resampling
     disp ('No - Resampling')
 end 
 %%%%% end-resampling     
end %%%% end mian loop - pieces of target
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Total number of resampling steps')
countRES
disp('TOTAL NUMBER OF ITERATIONS')
length(mu_d)
%%%%
if countRES==length(mu_d)
 disp('(THEN, we ALWAYS did RESAMPLING)')
elseif countRES==0
    disp('(THEN, no RESAMPLING)')
end
%%%
cprintf('blue',    'Press a key' );
disp(' ')
disp(' ')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('---------------------------------------------------------------------------------------------------------')
disp('Below we show all the estimations of the partial marginal likelihoods, Z_d, d=1,....,D')
disp(['with D=',num2str(DIM)])
disp('---------------------------------------------------------------------------------------------------------')
disp('VALID only WITH NO RESAMPLING (i.e., in SIS)')
disp(VALID_only_WITHOUT_RES)
disp('---------------------------------------------------------------------------------------------------------')
%%%
cprintf('blue','WITH A GIS WEIGHTING, The two estimators below are equal IN ANY CASE');
disp(' ')
disp(' ')
disp(Z1)
disp(Z2)
disp('---------------------------------------------------------------------------------------------------------')
disp('VALID_only WITH ALWAYS RESAMPLING')
disp(VALID_only_WITH_ALWAYS_RES)
disp('---------------------------------------------------------------------------------------------------------')
Ztot_true=(sqrt(2*pi*sig^2))^DIM;
disp(['True value of normalazing constant (marginal likelihood) of the Target function= ',num2str(Ztot_true)])





