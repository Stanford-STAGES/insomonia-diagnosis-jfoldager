clear; 
close all; 
clc;
cd('C:\Users\Jonathan\Dropbox\Skole\_Master Thesis\MATLAB\');
addpath(genpath(pwd));
% Plotting paramters
plotprintflag   = false;    % True outputs .png files to figures
fsz             = 20;       % Fontsize
set(groot, 'DefaultAxesFontSize', fsz);
set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultAxesTickLabelInterpreter', 'LaTeX');
set(groot, 'DefaultAxesFontName', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');
set(groot, 'DefaultTextColor', [0, 0, 0]);
set(groot, 'DefaultAxesXColor', [0, 0, 0]);
set(groot, 'DefaultAxesYColor', [0, 0, 0]);
set(groot, 'DefaultLineLineWidth', 1.2);
set(groot, 'DefaultFigureUnits','normalized')
set(groot, 'DefaultFigurePosition',[0.25 0.25 0.35 0.5])
set(groot, 'DefaultAxesBox', 'off')  
% set(groot, 'axis', 'tight') 
% set(0,'DefaultFigureColor', '')
plotPath1 = '..\Plots\';
plotPath2 = 'C:\Users\Jonathan\Dropbox\Apps\ShareLaTeX\Thesis\figure\Step1\';
plotPath3 = 'C:\Users\Jonathan\Dropbox\Apps\ShareLaTeX\Thesis\figure\Step2\';
clc
%% Set periods of night and circadian periods
global PPeriods CPeriods
PPeriods    = [0 0.25; 0.25 0.50; 0.50 0.75; 0.75 1];
CPeriods    = [ 18 00 00, 23 59 59; ...
                00 00 00, 01 59 59; ...
                02 00 00, 03 59 59; ...
                04 00 00, 05 59 59; ...
                06 00 00, 12 59 59; ...
                ];
%% Set frequency intervals
freqs.SW    = [0.5 1];
freqs.DE    = [1 4];
freqs.TH    = [4 8];
freqs.AL    = [8 12];
freqs.SS    = [12 13.5];
freqs.FS    = [13.5 15];
freqs.SB    = [15 20];  
freqs.FB    = [20 30];    
freqs.TOTAL = [0.5 30];     


