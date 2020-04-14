%% Analyzing Burst FUS PCD data 
%2.5MHz H-147 Transducer, Sonic Concepts
%Y107   PCD, Sonic Concepts

clear all;
%% Put directory which contains all mat.file here
%cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_5_Burst_25mV/UC_Carriers')
%cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_6_Burst_FUS_100mV/UC_Carriers')
cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_4_Burst_FUS_200mV/UC_Carriers')
%% isert variables
Fs=125000000; % 125MS/s
recording_window=20; % 20 ms
pulse_duration=2 % ms
t = 0:1/Fs:1-1/Fs;
window_length=Fs/(1000/recording_window)
pulse_size=Fs*(pulse_duration/1000) % number of samples for recording_window
%fprintf('Windows Length in ms :%i\n',pulse_size )
delay_after_TTL=77*125 % 77 us
psdx_size=pulse_size/2

%% AUC Inputs
%center_f=2.5*1000000
%fn=center_f*3 % 7.5 MHz
%fn_window=center_f/50 % in kHz
start_freq=12400
stop_freq=12600
%% Load .mat files in an array and take FFT
matfiles = dir('*.mat');
N = length(matfiles);
filenames = cell(N,1);
data_array=zeros(N,window_length);
psdx_array=zeros(N,psdx_size);
data_auc_array= zeros(N,1); 
%%
for i = 1:N
   %thisfig = figure();
   filenames{i} = matfiles(i).name;
   transit = load(filenames{i});
   data_array(i,:)=transpose(transit.A);
   data_tr=data_array(i,:);
   size_A=size(transit.A);
   start_pulse=round(2*size_A(1)/20)+delay_after_TTL; % pulse starts from here, based on picoscope TTL
   x=data_tr([start_pulse:start_pulse+pulse_size-1]);
   x=x.*hamming(length(x))'; % applying hamming window
   n = length(x);
   xdft = fft(x);
   xdft = xdft(1:n/2);
   psdx = (1/(Fs*n)) * abs(xdft).^2;
   freq = 1:Fs/length(x):Fs/2;
   freq=freq/1000000;
   data_auc=psdx;
   freq([start_freq+1:stop_freq+1]);
   data_auc([start_freq+1:stop_freq+1]);
   auc=trapz(1000000.*freq([start_freq+1:stop_freq+1]),data_auc([start_freq+1:stop_freq+1]));
   data_auc_array(i,1)= auc;
   psdx_array(i,:)=psdx;  % store in array
end

%% Bubbles Mean

figure();
x0=10;
y0=10;
width=1200;
height=600
set(gcf,'position',[x0,y0,width,height]);
set(gca,'fontsize',20);

mean_psdx=mean(psdx_array);
ucc_mean = mean(data_auc_array)
plot(freq([start_freq+1:stop_freq+1]),10*log10(mean_psdx([start_freq+1:stop_freq+1])),'red');
%plot(freq,10*log10(mean_psdx),'red');
xlim([6.2 6.3])
%ylim([-170 -60])
xlabel('Frequency MHz')
ylabel('Power')
set(gca,'fontsize',16)
hold on

%%
%cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_5_Burst_25mV/Baseline')
%cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_6_Burst_FUS_100mV/Baseline')
cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_4_Burst_FUS_200mV/Baseline')
%% Load .mat files in an array and take FFT
matfiles = dir('*.mat');
N = length(matfiles);
filenames = cell(N,1);
data_array=zeros(N,window_length);
psdx_array=zeros(N,psdx_size);
data_auc_array= zeros(N,1); 
%%
for i = 1:N
   %thisfig = figure();
   filenames{i} = matfiles(i).name;
   transit = load(filenames{i});
   data_array(i,:)=transpose(transit.A);
   data_tr=data_array(i,:);
   size_A=size(transit.A);
   start_pulse=round(2*size_A(1)/20)+delay_after_TTL; % pulse starts from here, based on picoscope TTL
   x=data_tr([start_pulse:start_pulse+pulse_size-1]);
   x=x.*hamming(length(x))'; % applying hamming window
   n = length(x);
   xdft = fft(x);
   xdft = xdft(1:n/2);
   psdx = (1/(Fs*n)) * abs(xdft).^2;
   freq = 1:Fs/length(x):Fs/2;
   freq=freq/1000000;
   data_auc=psdx;
   freq([start_freq+1:stop_freq+1]);
   data_auc([start_freq+1:stop_freq+1]);
   auc=trapz(1000000.*freq([start_freq+1:stop_freq+1]),data_auc([start_freq+1:stop_freq+1]));
   data_auc_array(i,1)= auc;
   psdx_array(i,:)=psdx;  % store in array
end
mean_psdx=mean(psdx_array);
plot(freq([start_freq+1:stop_freq+1]),10*log10(mean_psdx([start_freq+1:stop_freq+1])),'blue');
%plot(freq,10*log10(mean_psdx),'blue');
baseline_mean = mean(data_auc_array)

ucc_mean/baseline_mean