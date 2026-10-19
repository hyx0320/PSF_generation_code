%% This script is used for generate different config file with Nnum
% MW 2021/08/16
clc,clear;
close all;
Nnum = 1;
config_pic = zeros(Nnum,Nnum);
count = 1;
for i = Nnum-1:-1:0
    
    line = (Nnum-1:-1:0) + i * Nnum;
    config_pic(Nnum-i,:) = line;
    imwrite(uint16(config_pic),strcat('./',num2str(Nnum),'x',num2str(Nnum),'.conf.sk.png'));
end