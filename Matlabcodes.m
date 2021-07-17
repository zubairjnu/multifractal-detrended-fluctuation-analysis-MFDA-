clear all

load sampledata.dat
scmin=16;
scmax=1024;
scres=19;
exponents=linspace(log2(scmin),log2(scmax),scres);
scale=round(2.^exponents);
q=linspace(-5,5,101);
q_1 = q';
save q_output.dat q_1 -ascii
m=1;
signal1=output_data(1:400,1); % inhA gene
%signal2=DS_2;
%signal3=DS_3;
[Hq1,tq1,hq1,Dq1,Fq1]=MFDFA1(signal1,scale,q,m,1);
% [Hq2,tq2,hq2,Dq2,Fq2]=MFDFA1_strain2(signal2,scale,q,m,1);
% [Hq3,tq3,hq3,Dq3,Fq3]=MFDFA1_strain3(signal3,scale,q,m,1);


% %Example: MFDFA2-----------------------------------------
% load fractaldata
% scale=[7,9,11,13,15,17];
% m=2;
% signal1=multifractal;
% signal2=monofractal;
% signal3=whitenoise;
% [Ht1,Htbin1,Ph1,Dh1] = MFDFA2(signal1,scale,m,1);
% [Ht2,Htbin2,Ph2,Dh2] = MFDFA2(signal2,scale,m,1);
% [Ht3,Htbin3,Ph3,Dh3] = MFDFA2(signal3,scale,m,1);
% %---------------------------------------------------------
% 







