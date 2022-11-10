close all
clear
clc
% =+=+=+=+=+=+=+=+=+=+=+=+ User data +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
% User data can be addedd here and should replace the test data below:
% =+=+=+=+=+=+=+=+=+=+ End of user data =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+

% =+=+=+=+=+=+=+=+=+=+=+=+ Test data +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
testDat = load('TestData.mat');
ECG_device  = testDat.ECG_device;
ECG_device_BP = testDat.ECG_device_BP;
rPeaks = testDat.rPeaks;
fs = testDat.fs;

% =+=+=+=+=+=+=+=+=+=+ End of test data =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
NEWC = NoiseExtractClass;
[noiseEstimated, tVec] = NEWC.ecgNoiseExtractor(ECG_device, ECG_device_BP, rPeaks, fs);
    
% Reduce noise record to same length as calculated noise record
NEWC.TestPlot(noiseEstimated, tVec, fs, ECG_device, testDat.ECG_std)