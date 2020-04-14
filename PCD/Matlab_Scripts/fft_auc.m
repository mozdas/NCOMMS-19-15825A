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

