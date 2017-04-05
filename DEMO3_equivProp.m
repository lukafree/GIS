clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('---------------------------------------------------------------------')
disp('                DEMO 3: EQUIVALENT PROPOSAL PDF         ')
disp('                    ')
disp('Reference: ')
disp('- L Martino, V. Elvira, G. Camps-Valls, ')
disp('Group Importance Sampling for particle filtering and MCMC, 2017 ')
disp('---------------------------------------------------------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Unnormalized TARGET DEFINITION (BIMODAL)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=2;
T=@(x) exp(-(x.^2-4).^2/(2*s^2));
%%% normalizing the target via determinist Riemann quadrature
stepT=0.1;
x=-20:stepT:20;
Zest_detQuad=sum(T(x).*stepT);
%%%%%%%%%%%%%%%%%%%%
%%% PROPOSAL PDF %%%
%%%%%%%%%%%%%%%%%%%%
sp=3;
mup=2;
P=@(x) 1/(sqrt(2*pi*sp^2)).*exp(-(x-mup).^2/(2*sp^2));
%%%%%%%%%%%%%%%
%%%% PLOT %%%%%
%%%%%%%%%%%%%%%
subplot(2,1,1)
plot(x,1/Zest_detQuad.*T(x),'b--','LineWidth',2)
hold on
plot(x,P(x),'r--','LineWidth',2)
legend('Target pdf \pi(x)','Proposal pdf q(x)')
box on
set(gca,'FontSize',15,'FontWeight','Bold')
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
subplot(2,1,2)
plot(x,1/Zest_detQuad.*T(x),'b--','LineWidth',2)
hold on
plot(x,P(x),'r--','LineWidth',2)
box on
set(gca,'FontSize',15,'FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Obtaining an approximation of the 
%%%% THEREOTICAL DEFINITION of EQUIVALENT PROPOSAL \tilde{q}(x) for N=2 (
%%%% (q_equiv in the code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:length(x)
  for j=1:length(x)
      wahora=(T(x(k))./P(x(k)));
    S=wahora+(T(x(j))./P(x(j)));
    if S==0
      S=10^(-300);
    end
    q(k,j)=(wahora./S).*P(x(k)).*P(x(j));
  end
   q_equiv(k)=2*sum(q(k,:).*stepT);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Drawing N_Approx samples from q_equiv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=2;
N_Approx=10000;
disp(['Drawing ',num2str(N_Approx),' samples from the equivalent proposal (when N=2)... '])
%%%% LOOP: drawing samples N_Approx from the equivalent proposal (when N=2)
for i=1:N_Approx
xp=mup+sp*randn(1,N);
w=T(xp)./P(xp)+10^(-300);
wn=w./(sum(w));
x_from_equiv(i)=randsrc(1,1,[xp; wn]);
end 
%%%% end drawing samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PLOTS %%%%
%%% histogram of the samples from q_equiv
[a,b]=hist(x_from_equiv,50);
aux=b(2)-b(1);
aux2=sum(a.*aux);
bar(b,1/aux2*a,'g')
%%% plotting equivalent proposal pdf
plot(x,q_equiv,'k','LineWidth',2)
%%%%
legend('Target pdf \pi(x)','Proposal pdf q(x)','Histogram by Sim.','Equiv. Prop. (N=2)')
set(gca,'FontSize',15,'FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%