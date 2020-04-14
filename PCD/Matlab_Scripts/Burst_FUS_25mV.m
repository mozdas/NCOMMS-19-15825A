%% Analyzing Burst FUS PCD data 
%2.5MHz H-147 Transducer, Sonic Concepts
%Y107   PCD, Sonic Concepts

clear all;
%% Put directory which contains all mat.file here
cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_5_Burst_25mV/UC_Carriers')
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
%% Load .mat files in an array and take FFT
matfiles = dir('*.mat');
N = length(matfiles);
filenames = cell(N,1);
data_array=zeros(N,window_length);
psdx_array=zeros(N,psdx_size);


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
   psdx_array(i,:)=psdx;  % store in array
   freq = 1:Fs/length(x):Fs/2;
   freq=freq/1000000;
end

%% Bubbles Mean

figure();
x0=10;
y0=10;
width=1200;
height=600
set(gcf,'position',[x0,y0,width,height]);
set(gca,'fontsize',20);

mean_psdx_bubbles=mean(psdx_array);

plot(freq,10*log10(mean_psdx_bubbles),'red');
xlim([0 16])
ylim([-170 -60])
xlabel('Frequency MHz')
ylabel('Power')
set(gca,'fontsize',16)
hold on


%% Saline
%cd('/Users/mehmetozdas/Desktop/Choi_lab/Test_Data/NBBB179/25mV_Saline_1')
%cd('/Users/mehmetozdas/Desktop/Choi_lab/Test_Data/NBBB179/100mV_Saline')
cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_5_Burst_25mV/Baseline')

%% Load .mat files in an array and take FFT
matfiles = dir('*.mat');
N = length(matfiles);
filenames = cell(N,1);
data_array=zeros(N,window_length);
psdx_array=zeros(N,psdx_size);


for i = 1:N
   %thisfig = figure();
   filenames{i} = matfiles(i).name;
   transit = load(filenames{i});
   data_array(i,:)=transpose(transit.A);
   data_tr=data_array(i,:);
   size_A=size(transit.A);
   start_pulse=round(2*size_A(1)/20); % pulse starts from here, based on picoscope TTL
   x=data_tr([start_pulse:start_pulse+pulse_size-1])+delay_after_TTL;
   x=x.*hamming(length(x))'; % applying hamming window
   n = length(x);
   xdft = fft(x);
   xdft = xdft(1:n/2);
   psdx = (1/(Fs*n)) * abs(xdft).^2;
   psdx_array(i,:)=psdx;  % store in array
   freq = 1:Fs/length(x):Fs/2;
   freq=freq/1000000;
end

mean_psdx_baseline=mean(psdx_array);

plot(freq,10*log10(mean_psdx_baseline),'blue');


figure();

x0=10;
y0=10;
width=1200;
height=600
set(gcf,'position',[x0,y0,width,height]);
set(gca,'fontsize',20);


net_psdx=mean_psdx_baseline-mean_psdx_bubbles
plot(freq,10*log10(net_psdx),'black');
title('0.22 MPa Bubbles-Baseline')
xlim([0 16])
ylim([-170 -60])
xlabel('Frequency MHz')
ylabel('Power')
set(gca,'fontsize',16)
fname = '/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/Figures/mean_FFTs'
saveas(gca, fullfile(fname, '0.22 MPa Bubbles-Baseline'), 'png');