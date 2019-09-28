function [fsig,amp,vos,snr,thd,sndr,enob,phi]=fft_hann(data, bw_low, bw_high, fs, fig, tit, opt)
% [fsig,amp,vos,snr,thd,sndr,enob,phi]=fft_hann(data, bw_low, bw_high, fs, fig, tit, opt);
% [fsig,amp,vos,snr,thd,sndr,enob,phi]=fft_hann(data, 20, 20e3, 65536, '', 'No title', opt);
% opt=fft_hann();  % return default option.
% fsig: signal frequency
% amp: sig amplitude
% vos: offset voltage
% thd: total harmonic distortion
% sndr: signal to noise and distortion ratio
% enob: efficient number of bits
% phi: signal phase
% data: data to be fft
% bw_low: lowband of inband frequency 
% bw_high: highband of inband frequency
% fs : sampling rate
% fig = combination if 'i','f','w','I','F','W', stands for inband fig,
% fullband fig, wave fig, Upper case will generate tit.png file
% tit: title of figure.
% opt: additional options:
% opt.nhd: number of harmonic distortion force to caculate
% opt.nspur: number of spur to caculate
% opt.RL: RL to caculate power;
% opt.win = 'barthannwin','bartlett','blackman','blackmanharris','bohmanwin',
% 'chebwin', 'flattopwin', 'gausswin', 'hamming', 'hann', 'kaiser', 
% 'nuttallwin', 'parzenwin', 'rectwin', 'taylorwin', 'triang', 'tukeywin',
% opt.coherent: coherent sampling, minimize skirt
% opt.sig_leakage: when coherent sampling, still caculate as much as possible skirt
% opt.sig_add_span: when coherent sampling, add additional skirt to besides
% of window skirt, no effect when opt.sig_leakage = 1;
% Example:
% x = sin((1:65536)*2*pi*111/2048);
% x = x + 0.001*x.^2 + 0.002*x.^3 + 0.0001*sin([1:65536]*2*pi*31/2048);
% x = x + randn(size(x))*1e-5 + 0.01;
% fft_hann(x,200,20e3,65536,'fiw','title');
% opt=fft_hann();
% opt.coherent=1;
% fft_hann(x,0,20e3,65536,'i','title',opt);
% opt.win = 'rectwin';
% fft_hann(x,0,20e3,65536,'i','title',opt);
% opt.sig_leakage = 1;
% fft_hann(x,0,20e3,65536,'i','title',opt);
% opt.win = 'blackmanharris';
% fft_hann(x,0,20e3,65536,'i','title',opt);
% Change log:
% v1.1
% updated by LiuYongping, 2010/12/4;
% increase the precision of frequency caculation.
% v1.2
% update signal detection function by LiuYongping 2012/10/28;
% v1.3
% add window parameter by LiuYongping 2013/1/6;
% v2.0
% 1. add opt structure to encapsulate other parameter.
% 2. add coherent sampling option. 
% Author: Liu Yongping, liu_yongping@126.com
%% input argument preprocess.
if nargin<7
    opt.nhd = 6;  % how many harmonics to caculate?
    opt.nspur = 99;  % how many spur to caculate?
    opt.RL = 100; % 100 Ohm
    opt.win = 'hann'; % hann window
    opt.coherent = 0; % set 1 if coherent sampling;
    opt.sig_leakage = 0; % when opt.coherent = 1 and still slightly leakage in sig, set this option;
    opt.sig_add_span = 0; % add additional span to sig;
end
if nargin<1
    fsig = opt;  % return default option when no inputs.
    return;
end
if nargin<2
    bw_low = 20; % default value 20
end
if nargin<3
    bw_high = 20e3; % default value 20e3
end
if nargin<4
    fs = 65536;  % default fs = 65536
end
if nargin <5
    fig='';
end
if nargin <6
    tit='No title';
end
copyright = 'Copyright 2006-2013@LiuYongping';
inband =0;
wave = 0;
fullband = 0;
if(bw_low<0)
    bw_low=0;
end
if(bw_high>fs/2)
    bw_high=fs/2;
end
N=length(data);
N=2^(floor(log2(N)));
s=size(data);
if(N<4)
    display('data length < 4')
    return;
end
if(length(s) ~= 2)
    display('data is not a array');
end
if(s(2) == 1)
    data=data';
end
% get last N points.
data=data(end-N+1:end);
%% calculate property of window.
w = eval([opt.win '(N)'])';  % hann window;
% plot(w); % plot window 
K = 64;
tmp = eval([opt.win '(K)'])/K;
tmp = fft(tmp);
tmp = tmp.*conj(tmp);
dcspan = sum(tmp(1:K/2)>1e-6)-1; % 
win_gain = sqrt(sum(tmp)); % window gain for amp correction
%% fft process.
mag = fft(data.*w/N)/win_gain;  %fft gain correction
mag = mag(1:N/2+1); % spectrum to fs/2
mag(2:N/2) = mag(2:N/2)*sqrt(2);
magph = mag; 
mag=mag.*conj(mag);  % calculate power spectrum to fs/2;
magdb=10*log10(mag); 
magdb(magdb<-499) = -499;  % delete -inf data to prevent wrong plot.
ind_lowband=round(bw_low/fs*N)+1;
if(ind_lowband<1)
    ind_lowband = 1;
end
ind_highband=round(bw_high/fs*N)+1;
if(ind_highband>N/2+1)
    ind_highband = N/2+1;
end
f=(0:N/2)/N*fs;  % frequency
df = 1/N*fs;
magt=mag;
vos_inds = findSig(magt,1,0,inf,1,opt.coherent,dcspan);
% fprintf(1,'dcspan=%d',dcspan);
vos = sqrt(sum(mag(vos_inds)));
magt(vos_inds) = 0;
sig_ind = findMax(magt,ind_lowband,ind_highband);  % find max power in signal band, it's signal;
sig_inds = findSig(magt,sig_ind,ind_lowband,ind_highband,1,0,dcspan);
fsig = sum(((sig_inds-1)*df).*magt(sig_inds))/sum(magt(sig_inds)) ; 
phi = sum(angle(magph(sig_inds)).*magt(sig_inds))/sum(magt(sig_inds)) ; 
% angle(magph(sig_inds))
sig_ind = fsig/df+1;
% phi = angle(magph(round(sig_ind)));
% fsigt = sum(((sig_inds-1)*df).*magt(sig_inds))/sum(magt(sig_bins))
% find harmonic location, important: harmonic maybe out of signal band,
% cause it will be fold back to signal band.
hds=round((1:opt.nhd)*(sig_ind-1));
hds=mod(hds,N);  % if harmonic frequency larger than fs, then fold back;
for k=1:length(hds)
    if(hds(k)>=N/2)
        hds(k)=N-hds(k);
    end
end
hds = hds+1;     % real harmonic location and it's side band;
hds = hds(hds>=ind_lowband); %calculate inband harmonic only
hds = hds(hds<=ind_highband);
% find signal
hd_inds=[];
hd_inds_temp = findSig(magt,hds(1),ind_lowband,ind_highband,1,opt.coherent&&~opt.sig_leakage, dcspan + opt.sig_add_span);
if(length(hd_inds_temp)>0)  % signal is found;
    hd_inds{1} = hd_inds_temp;
    sig_pow=sum(magt(hd_inds{1}));
    magt(hd_inds{1}) = 0; %exclude signal
else
    sig_pow = 0;
    hds=[];
end
spur_ind = findMax(magt,ind_lowband,ind_highband);
spur_inds = findSig(magt,spur_ind,ind_lowband,ind_highband,1,opt.coherent,dcspan);  %% try to find max spur.
spur_pow = sum(mag(spur_inds)); % caculator spur power;
hds_pows=zeros(1,length(hds));
for k=2:length(hds)
    force = 1;
    if(k>6)
        force = 0;
    end
	hd_inds{k} = findSig(magt,hds(k),ind_lowband,ind_highband,force,opt.coherent,dcspan); %% try to find harmonic.
	hds_pows(k) = sum(magt(hd_inds{k}));
	magt(hd_inds{k}) = 0;
end
hds_pow = sum(hds_pows(2:end));
pow_spur=[];
ind_spur=[];
ind_spurs=[];
for i=1:opt.nspur
    ind_spur(i)=findMax(magt,ind_lowband,ind_highband);
    ind_spurs_temp = findSig(magt,ind_spur(i),ind_lowband,ind_highband,0,opt.coherent,dcspan); % try to find spur in non force mode.
    if(length(ind_spurs_temp) < 1)  %% no spur found
        break;
    end
    ind_spurs{i} = ind_spurs_temp;
    pow_spur(i) = sum(magt(ind_spurs{i}));
    magt(ind_spurs{i}) = 0;
end
all_pow=sum(mag(ind_lowband:ind_highband)); % calculate all power;
noi_pow=sum(magt(ind_lowband:ind_highband)); % calculate noise power;
noi_floor=db10(mean(magt(ind_lowband:ind_highband)))*ones(1,length(magt(ind_lowband:ind_highband)));
thd=db10(hds_pow/sig_pow);
snr=db10(sig_pow/noi_pow);
sndr=db10(sig_pow/(noi_pow+hds_pow));
sfdr=db10(sig_pow/spur_pow);
sigdb=db10(sig_pow);
sigdbm=db10(sig_pow/opt.RL*1000);
amp=sqrt(sig_pow);
enob=(sndr-1.76)/6.02;
% print message
[y,e,u]=engunits(fsig,'latex');
str=sprintf('Freq: %0.2f%sHz',y,u);
[y,e,u]=engunits(vos,'latex');
str=sprintf('%s, Vos: %0.2f%sV\n',str,y,u);
[y,e,u]=engunits(amp,'latex');
str=sprintf('%sAmp: %0.2fdBV,%0.2fdBm,%0.2f%sVrms',str,sigdb,sigdbm,y,u);
[y,e,u]=engunits(amp*2*sqrt(2),'latex');
str=sprintf('%s, %0.2f%sVpp\n',str,y,u);
str=sprintf('%sSNDR: %0.2fdB, ENOB: %0.2f\n',str,sndr,enob);
str=sprintf('%sSNR: %0.2fdB, THD: %0.2fdB, SFDR: %0.2fdB\n',str,snr,thd,sfdr);
[y,e,u]=engunits(sqrt(noi_pow),'latex');
str=sprintf('%sNoise: %0.2f%sVrms,%3.1fdBm,%3.1fdBm/Hz\n',str,y,u,db10(noi_pow/opt.RL)+30,db10(noi_pow/opt.RL/(bw_high-bw_low))+30);
str=sprintf('%sNfft: %d points,RL:%0.2fOhm\n',str,N,opt.RL);
msgstr=sprintf('%sNumber of inband point: %d',str,ind_highband-ind_lowband+1);
for k=1:length(fig)
    switch fig(k)
        case {'i','I'}
            h=figure();
            hold on;
            ind_add = ceil((ind_highband - ind_lowband + 1)*0.1);
            ind_low = max(1,ind_lowband - ind_add);
            ind_high = min(N/2,ind_highband + ind_add);
            plot(f(ind_low:ind_high),magdb(ind_low:ind_high),'b'); % plot all in blue
            plot(f(ind_lowband:ind_highband),magdb(ind_lowband:ind_highband),'Color',[1 1 1]*0.5); % plot inband in black
            magt_dbm = smooth(magt,max(floor((ind_highband-ind_lowband)/20),8));
            df = 1/N*fs;
            magt_dbm = db10(magt_dbm/opt.RL/df) + 30;
            plot(f(ind_lowband:ind_highband),magt_dbm(ind_lowband:ind_highband),'Color',[1 1 0]*0.3); % plot dbm noise in black
            
            if(length(hd_inds) < 1)
                return;
            end % no signal find, return;
%             plot(f(vos_inds),magdb(vos_inds),'-b');  % plot dc signal in blue
            plot(f(hd_inds{1}),magdb(hd_inds{1}),'-r');  % plot signal in red
            for kk=2:length(hd_inds)  % plot all harmonics in mag
                text(max(f(hd_inds{kk})),max(magdb(hd_inds{kk})),sprintf('%d',kk))
                plot(f(hd_inds{kk}),magdb(hd_inds{kk}),'-m');
            end
            if(length(spur_inds)>1)
                plot(f(spur_inds),magdb(spur_inds),'-c');
            else
                plot(f(spur_inds),magdb(spur_inds),'*c');
            end
            xcorr=f(hd_inds{1}(end));
            if(xcorr>(bw_low+bw_high)/2)
                xcorr = bw_low;
            end
            text(xcorr,max(sigdb,noi_floor(1)+35),msgstr,'Color','b','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','bold');
%             plot(f(ind_lowband:ind_highband),noi_floor);
            grid on;
            set(gca,'LineWidth',0.1);
            xlabel('Frequency(Hz)');
            ylab = sprintf('Magnitude(dB/%gHz)',fs/N);
            ylabel(ylab);
            plot(f(hds),magdb(hds),'r*');
            if(length(ind_spurs)>0)    
                for i=1:length(ind_spurs)
                    plot(f(ind_spurs{i}),db10(mag(ind_spurs{i})),'-c');
                    if(length(ind_spurs{i})==1)
                        plot(f(ind_spurs{i}),db10(mag(ind_spurs{i})),'c*');
                    end
                    str=sprintf('%0.2f',db10(pow_spur(i)));
                    text(f(ind_spur(i)),db10(mag(ind_spur(i))),str);
                end
            end
    
            str=[tit '-inband'];
            title(str,'Interpreter','none');
            set(h,'Name',[str copyright]);
            if(fig(k) == 'I') % upper case to print file.
                print(h,[str '.png'],'-dpng');
            end
        case {'w','W'}
            h=figure();
            t=[0:N-1]/fs;
            plot(t,data);
            str=[tit '-wave'];   
            title(str,'Interpreter','none');
            set(h,'Name',[str copyright]);
            if(fig(k) == 'W') % upper case to print file.
                print(h,[str '.png'],'-dpng');
            end 
        case {'f','F'}
            h=figure();
            plot(f,magdb),grid;
            str=[tit '-fullband'];
            title(str,'Interpreter','none');
            set(h,'Name',[str copyright]);
            if(fig(k) == 'F') % upper case to print file.
                print(h,[str '.png'],'-dpng');
            end 
    end
end
function help_msg
   help fft_imd
function inds=findSig(mag,ind,ind_lowband,ind_highband,force,coherent,span)
% mag: input signal; ind: center; ind_lowband:
    if(nargin < 5)
        force = 1;
    end
    hill = 3; % 3dB to climb
    ind_low = max(ind_lowband,1);
    ind_high= min(ind_highband,length(mag));
    if(ind < 1)
	    error('findSig: ind < 1');
    end
    if(ind > length(mag))
	    error('findSig: ind > length(mag)');
    end
    pow=mag(ind);% power at ind
    k = ind;
    falling = 0; % no falling at the begining;
    if(coherent)
        ind_low = max(ind-span, ind_low);
        ind_high = min(ind+span, ind_high);
    else
        while(1)
            k = k - 1;    
            if(k<ind_low)
                k = k + 1;           
                break; % out of signal band, quit;
            end
             if(mag(k) <= 0)
                break;
            end       
            if((db10(mag(k))-db10(pow) > hill))
                if(falling == 1)
            k = k + 1;
                    break; % find rising after falling, quit
                end
            else
                falling = 1;   % find falling
            end
            pow = mag(k);
        end
        ind_low = k;
        k=ind;
        pow=mag(k); % 
        falling = 0; % no falling at the begining;    
        while(1)
            k = k + 1;    
            if(k>ind_high)
                k = k - 1;         
                break; % out of signal band, quit;
            end
            if(mag(k) <= 0)
                break;
            end
            if((db10(mag(k))-db10(pow) > hill))
                if(falling ==1)
                    k = k - 1;
                    break; % find rising after falling, quit
                end
            else
                falling = 1;   % find falling
            end
            pow = mag(k);
        end
        ind_high = k;
    end
    inds=ind_low:ind_high;
    % don't consider as signal.
    maxsig=max(mag(inds));
    if(force)
        return;
    end
    %% determin if is required sig.
    %% should be improved
%     if((db10(mag(ind_low))>db10(maxsig)-20) || (db10(mag(ind_high))>db10(maxsig)-20))
    bw = floor((inds(end) - inds(1))/2);
    bw = max(bw,3);
    vbw = [(inds(1)-2*bw):inds(1)-1 inds(end)+1:(inds(end)+2*bw)];
    vbw = vbw(vbw>ind_lowband);
    vbw = vbw(vbw<ind_highband);    
    
    if(max(db10(mag(inds))) - max(db10(mag(vbw))) < 6)
        inds=[];
    end
function ind=findMax(mag,ind_lowband,ind_highband)
    if(ind_lowband < 1)
	    error('findMax: ind_lowband < 1');
    end
    if(ind_highband > length(mag))
	    error('findMax: ind_highband < length(mag)');
    end
    if(ind_highband < 1)
	    error('findMax: ind_highband < 1');
    end
    if(ind_highband > length(mag))
	    error('findMax: ind_highband < length(mag)');
    end
    if(ind_lowband > ind_highband)
	    error('findMax: ind_lowband > highband');
    end
    [amp,ind]=max(mag(ind_lowband:ind_highband));
    ind  = ind_lowband + ind -1 ;
function y = db10(x)
    y =10*log10(x);
    
function w = hann(N)
    w=0.5*(1-cos((0:N-1)/N*2*pi))';  % hann window;
        
%%%%
%% change log:
% improve FindSig to discover spurs. may caused problem if two signal are
% too close.
% dc offset may gen problem
% add dbm/Hz noise line(usefull for tx outband noise.)
