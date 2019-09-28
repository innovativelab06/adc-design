function data=comp(vinp, vinn, opt_comp)
%% comparator module
% 1. offset and hysteresis
% 2. sampling aperture, timing resolution, uncertainty window
% 3. regeneration gain, voltage sensitivity, metastability
% 4. random decision errors, input-referred noise
%% input argument preprocess
if nargin<3
    opt_comp.offset=0;
end

%% comparator input offset voltage
comp_offset=1e-3;

%% comparator noise
comp_noise=0;

%% comparator metastability
btr=0;

%% comparator decision
if opt_comp.offset==1
    if (vinp + comp_offset)>vinn
        data = 1;
    else
        data = 0;
    end
else
    if vinp > vinn
        data=1;
    else
        data=0;
    end
end
