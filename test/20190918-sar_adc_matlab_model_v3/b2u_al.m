function output = b2u_al(bits, cap_binary)
%% binary to unary code - arithmetic unit for imcs
    N=sum(cap_binary);
    output=0.5*ones(N,1);
    Ns=length(cap_binary);
    num.h=0; % number of high voltage level
    num.m=0; % number of vcm middle voltage level
    num.l=0; % number of low voltage level
    for i=1:Ns
        if bits(i)==1
            num.h=num.h+cap_binary(i);
        end
        if bits(i)==0
            num.l=num.l+cap_binary(i);
        end
    end
    count=max(num.h,num.l);
    for i=1:count
        if num.h>0   
            output(num.h)=1;
            num.h=num.h-1;
        end
        if num.l>0
            output(N-num.l+1)=0;
            num.l=num.l-1;
        end
    end
            
    