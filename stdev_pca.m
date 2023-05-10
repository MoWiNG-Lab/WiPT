clc; clear; constants; wipt = WiPT;  close all; 
DATA = DATA_5JanTest4xRandVert;
fullcsv = readmatrix(csvfiles(DATA));
csirange = csiranges(2*DATA-1):csiranges(2*DATA); fullcsv = fullcsv(csirange, :);
H = fullcsv(:, 3:66); l = height(H);


[coeff, ~, ~, ~, explained] = pca(H);
P = H * coeff;
w = 100;
sdvs = zeros(round(l, 0), 64);
for i=1:1:l-w
    sdvs(i, :) = std(P(i:i+w, :));
end




figure;
plot(sdvs);