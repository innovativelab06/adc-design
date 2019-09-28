function data=ideal_dac(output,N)
%% This is an ideal dac
% data: real value of output
% output: the bit data from adc
% vref: dac reference voltage
% N: bit of dac
%% ideal dac
Ns=length(output);
data=zeros(1,Ns);
for i=1:Ns
    j=1;
    if output(i,j)==1
        data(i)=0;
    else
        data(i)=-1;
    end
    for j=2:N
        data(i)=data(i)+output(i,j)*2^(N-j+1)/2^N;
    end
end