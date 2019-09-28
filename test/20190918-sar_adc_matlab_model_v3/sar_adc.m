function [data]=sar_adc(N,std_ramp,opt_comp,Ns,fin,fs,AMP,vref)
%% This is a sar adc matlab model
% This model is based on 
% 1. Monotonic Capacitor Switching Procedure 
% 2. Top Plate sampling method
% including cap mismatch, comparator non-ideality
% data: output data 
% N: bits of ADC
% std_cap: the Standard Deviation (std) of an unit cap from foundry doc
% opt_comp: choose whether ideal comparator is used. 1: nonideal comparator
% 0: ideal comparator;
% AMP: amplitude of input signal
% fin: input frequency
% fs: sampling frequency
% vref: reference voltage amplitude
% Ns: length of output code
%% input argument preprocess
if nargin<1 % ideal adc
    error('the bits N of adc must be defined');
end
if nargin<2
    std_ramp=0.516; % default value 0: no cap mismatch
end
if nargin<3
    opt_comp=0; % default value 0: no comparator non-ideality
end
if nargin<4
    Ns=65536; 
end
if nargin<5
    fin=200; 
end
if nargin<6
    fs=2000; 
end
if nargin<7
    AMP=1; 
end
if nargin<8
    vref=1; 
end
    
%% CDAC
c_norp=zeros(1,N-1); % number of unit cap per CDAC
% CBW
for n=1:N-1
    c_norp(n)=2^(N-1-n); 
end   
c_norp=[c_norp,1]; % add the abandunt unit cap
c_norn=c_norp;

%% cap mismatch 
ncu=10; % ncu=Cu/Co=WL, the area of unit cap
std_cap=std_ramp./sqrt(ncu.*c_norp);
c_devp=std_cap.*c_norp.*randn(1,N);
c_devn=std_cap.*c_norn.*randn(1,N);
cp=c_norp+c_devp;
cn=c_norn+c_devn;
cp_tot=sum(cp); % total capacitance of the cdac 
cn_tot=sum(cn);

%% control logic
% gnd=0; 
vcm=0; 
data=zeros(Ns,N);
for n=0:Ns-1
    do=zeros(1,N);
    % Sampling
    vin=AMP*sin(2*pi*fin/fs*n); 
    vinp=vcm+0.5*vin;
    vinn=vcm-0.5*vin;
    % MSB decision 
    % CDAC switch
    bp=ones(N,1);% bit switch for vxp, 1 => vref choose; 0 => gnd choose 
    bn=ones(N,1);% bit switch for vxn
    vxp=vinp; % the bootstrap characterisitcs proposed in Liu et al.'s work
    vxn=vinn; % compare vxp & vxn instead of vinp & vinn
    % comparator decision
    do(1)=comp(vxp, vxn, opt_comp); % MSB output
    % next bit CDAC switch
    if comp(vxp, vxn, opt_comp)==1
        bp(1)=0; % switch to gnd
    else
        bn(1)=0;
    end
    for i=2:N
        vxp=cp*bp*vref/cp_tot - vref + vinp;
        vxn=cn*bn*vref/cn_tot - vref + vinn;
        vx=vxp-vxn;
        do(i)=comp(vxp, vxn, opt_comp);
        if comp(vxp, vxn, opt_comp)==1
            if i==N
                break;
            else
                bp(i)=0;
            end
        else
            if i==N
                break;
            else
                bn(i)=0;
            end
        end
    end
    for i=1:N
        data(n+1,i)=do(i);
    end
end

