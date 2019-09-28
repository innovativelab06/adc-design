function [data]=sar_adc_9b(Ns,fin,fs,AMP)
%% This is a sar adc matlab model
% 1. imcs switching method
% 2. split CDAC
% 3. +VREF and -VREF reference source
% 4. binary to thermometer code
% 5. redundancy bit 
% 6. cap mismatch
%% split CDAC
ncu=10; % WL blocks
std_ramp=0.01; % cap std slope = 0.518
Nb=9; % effective bits=8
Nr=10; % real bits=9
%% P/N array
ch=[2 2 2 1 1]; % [8 2 2 2 1 1] 
cl=[16 8 4 2 1];
ch_len=length(ch);
cl_len=length(cl);
ch_tot=sum(ch);
cl_tot=sum(cl);
cb_i=1;
ch_ui=ones(1,ch_tot);
cl_ui=ones(1,cl_tot);
%% P array
% cap mismatch
std_cap_cpb=std_ramp/sqrt(ncu*cb_i);
std_cap_cph=std_ramp./sqrt(ncu.*ch_ui);
std_cap_cpl=std_ramp./sqrt(ncu.*cl_ui);
cpb_d=std_cap_cpb.*cb_i.*randn(1,1);
cph_d=std_cap_cph.*ch_ui.*randn(1,ch_tot);
cpl_d=std_cap_cpl.*cl_ui.*randn(1,cl_tot);
% fianl cap value
cpb=cb_i+cpb_d;
cph_u=ch_ui+cph_d;
cpl_u=cl_ui+cpl_d;
cph_tot_u=sum(cph_u);
cpl_tot_u=sum(cpl_u);
cp_lsb_u=cpb*cpl_tot_u/(cpl_tot_u+cpb);
cpk_u=cph_tot_u/(cph_tot_u+cp_lsb_u);
cp_tot_u=(cpl_tot_u+cpb)*cph_tot_u;
%% N array
% cap mismatch
std_cap_cnb=std_ramp/sqrt(ncu*cb_i);
std_cap_cnh=std_ramp./sqrt(ncu.*ch_ui);
std_cap_cnl=std_ramp./sqrt(ncu.*cl_ui);
cnb_d=std_cap_cnb.*cb_i.*randn(1,1);
cnh_d=std_cap_cnh.*ch_ui.*randn(1,ch_tot);
cnl_d=std_cap_cnl.*cl_ui.*randn(1,cl_tot);
% fianl cap value
cnb=cb_i+cnb_d;
cnh_u=ch_ui+cnh_d;
cnl_u=cl_ui+cnl_d;
cnh_tot_u=sum(cnh_u);
cnl_tot_u=sum(cnl_u);
cn_lsb_u=cnb*cnl_tot_u/(cnb+cnl_tot_u);
cnk_u=cnh_tot_u/(cnh_tot_u+cn_lsb_u);
cn_tot_u=(cnl_tot_u+cnb)*cnh_tot_u;
%% IMCS switching method
vcm=0;
vref=0.5;
data=zeros(Ns,Nb);
for n=0:Ns-1
    do=zeros(1,Nr);
    % sampling
    bnl=0.5*ones(cl_len,1); 
    bpl=0.5*ones(cl_len,1);
    bnl_u=b2u_al(bnl,cl);
    bpl_u=b2u_al(bpl,cl);
    vin=AMP*sin(2*pi*fin/fs*n);
    vinn=-0.5*vin;
    vinp=0.5*vin;
    % MSB decision
    bnh=0.5*ones(ch_len,1);
    bph=0.5*ones(ch_len,1);
    bnh_u=b2u_al(bnh,ch);
    bph_u=b2u_al(bph,ch);
    % CDAC decision
    for i=1:Nr
        vxp=(1+cnk_u)*vcm+cnk_u*(-vinn);
        vxp=vxp+cnk_u*(cnh_u*(2*bnh_u-1)/cnh_tot_u*vref); 
        vxp=vxp+cnk_u*(cnl_u*(2*bnl_u-1)/cn_tot_u*vref);
        vxn=(1+cpk_u)*vcm+cpk_u*(-vinp);
        vxn=vxn+cpk_u*(cph_u*(2*bph_u-1)/cph_tot_u*vref);
        vxn=vxn+cpk_u*(cpl_u*(2*bpl_u-1)/cp_tot_u*vref);
        do(i)=comp(vxp,vxn);
        if i==Nr
            break;
        else 
            if do(i)==1 
                if (1<=i)&&(i<=ch_len-1)
                    bnh(i)=0;
                    bph(i)=1;
                else
                    bnl(i-(ch_len-1))=0;
                    bpl(i-(ch_len-1))=1;
                end
            else
                if (1<=i)&&(i<=ch_len-1)
                    bnh(i)=1;
                    bph(i)=0;
                else
                    bnl(i-(ch_len-1))=1;
                    bpl(i-(ch_len-1))=0;
                end
            end
        end
        bnl_u=b2u_al(bnl,cl);
        bpl_u=b2u_al(bpl,cl);
        bnh_u=b2u_al(bnh,ch);
        bph_u=b2u_al(bph,ch);
    end  
    % redundancy to decimal based on specific structure
    % s(1)+ 4b + 1b -1 
    count=8;
    for i=1:4
        count=count+(2*do(i)-1)*ch(i);
    end
    count=count+do(5)-1;
    % decimal to binary
    for i=1:4
        data(n+1,4-i+1)=mod(count,2);
        count=fix(count/2);
    end
    % low 4b binary
    for i=Nb-cl_len+1:Nb
        data(n+1,i)=do(i+1);
    end
end