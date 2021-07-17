clear all
clc;
load sampledata.dat
smin=10;
smax=200; 
N_1=30;



theta_1=0;


signal1=notch_5(1000:2000,1);

q=linspace(-5,5,101);

[n,Fq,tau,alpha,f]=F_ALPHA(signal1,smin,smax,N_1,theta_1,q);

%save n_notch_5.dat n -ascii
save alpha_output.dat alpha -ascii
save f_alpha_output.dat f -ascii
%save Fq_notch_5.dat Fq -ascii
