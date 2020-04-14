%% Analyzing AU-FUS Radiation FUS PCD data 
%2.5MHz H-147 Transducer, Sonic Concepts
%Y107   PCD, Sonic Concepts

clear all;
%% Put directory which contains all mat.file here
cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_3_AUFUS/UC_Carriers')

%% isert variables
Fs=125000000; % 125MS/s
recording_window=200; % 20 ms
pulse_duration=50 % ms
t = 0:1/Fs:1-1/Fs;
window_length=Fs/(1000/recording_window)
pulse_size=Fs*(pulse_duration/1000) % number of samples for recording_window
%fprintf('Windows Length in ms :%i\n',pulse_size )
delay_after_TTL=77*125 % 77 us
uncage_pulse=50000 % 1000 Cycles
psdx_size=uncage_pulse/2  % 1000 cycles

%% AUC Inputs
%center_f=2.5*1000000
%fn=center_f*3 % 7.5 MHz
%fn_window=center_f/50 % in kHz
start_freq=2480
stop_freq=2520
%% Load .mat files in an array and take FFT
matfiles = dir('*.mat');
N = length(matfiles);
filenames = cell(N,1);
data_array=zeros(N,window_length);
psdx_array=zeros(N,25000);
threshold=0.01
data_auc_array= zeros(N,1); 
for i = 1:N
   filenames{i} = matfiles(i).name;
   transit = load(filenames{i});
   data_array(i,:)=transpose(transit.A(1:25000000,:));
   data_tr=data_array(i,:);
   size_A=size(transit.A);
   start_pulse=round(3*size_A(1)/4); 
   x=data_tr([start_pulse:start_pulse+pulse_size-1]);
   cross=find(x > threshold);
   x=x([cross(1):round(cross(1)+uncage_pulse-1)]);
   x=x.*hamming(length(x))'; % applying hamming window
   %plot(x);
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

%% Put directory which contains all mat.file here

figure();
x0=10;
y0=10;
width=1200;
height=600
set(gcf,'position',[x0,y0,width,height]);
set(gca,'fontsize',20);

mean_psdx=mean(psdx_array);
mean_psdx_bubbles=mean(psdx_array);
ucc_mean = mean(data_auc_array)
plot(freq([start_freq+1:stop_freq+1]),10*log10(mean_psdx([start_freq+1:stop_freq+1])),'red');
%plot(freq,10*log10(mean_psdx),'red');
xlim([6.2 6.3])
%ylim([-170 -60])
xlabel('Frequency MHz')
ylabel('Power')
set(gca,'fontsize',16)
hold on

%% Load Baseline Data
cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_3_AUFUS/Baseline')
%% Load .mat files in an array and take FFT
matfiles = dir('*.mat');
N = length(matfiles);
filenames = cell(N,1);
data_array=zeros(N,window_length);
psdx_array=zeros(N,25000);
threshold=0.01
data_auc_array= zeros(N,1);

for i = 1:N
   filenames{i} = matfiles(i).name;
   transit = load(filenames{i});
   data_array(i,:)=transpose(transit.A(1:25000000,:));
   data_tr=data_array(i,:);
   size_A=size(transit.A);
   start_pulse=round(3*size_A(1)/4); 
   x=data_tr([start_pulse:start_pulse+pulse_size-1]);
   cross=find(x > threshold);
   x=x([cross(1):round(cross(1)+uncage_pulse-1)]);
   x=x.*hamming(length(x))'; % applying hamming window
   %plot(x);
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
mean_psdx_baseline=mean(psdx_array);

%plot(freq,10*log10(mean_psdx),'blue');
plot(freq([start_freq+1:stop_freq+1]),10*log10(mean_psdx([start_freq+1:stop_freq+1])),'blue');
%plot(freq,10*log10(mean_psdx),'blue');
baseline_mean = mean(data_auc_array)

ucc_mean/baseline_mean

figure();

x0=10;
y0=10;
width=1200;
height=600
set(gcf,'position',[x0,y0,width,height]);
set(gca,'fontsize',20);


net_psdx=mean_psdx_bubbles-mean_psdx_baseline
plot(freq,10*log10(net_psdx),'black');

xlim([0 16])
ylim([-170 -60])
xlabel('Frequency MHz')
ylabel('Power')
set(gca,'fontsize',16)
