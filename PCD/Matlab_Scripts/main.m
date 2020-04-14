%% Main Function for analyzing Burst FUS PCD data 
% Author Mehmet Ozdas, zdasm@ethz.ch, ETH Zurich, Neurotechnology
% Laboratory, 2020
% 2.5MHz H-147 Transducer, Sonic Concepts
% Power amplifier: E&I 325LA for Transducer.50dB 
% Function generator: Agilent 33210A to triger Transducer
% Y107   PCD, Sonic Concepts
% 20dB, Broadband RF Amplifier 1MHz-1GHz, Ramsey Electronics
% Sampling Freq=125MS/s, 15 bits resolution 
% Picoscope: 5042 for data acquisation
% PicoScope: 3205B used for TTL to drive PRF of Picoscope: 5042

%% Insert the data path of .mat files for 

for i=1:3
    cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_6_Burst_FUS_100mV/UC_Carriers')
    sommaprodotto( 10,20)
    pause(2)
    cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_6_Burst_FUS_100mV/Baseline')
    sommaprodotto( 2,20)
    pause (5)
end
