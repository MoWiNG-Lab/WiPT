clc; clear; constants; wipt = WiPT;  close all; 
% {DATA_SeqTest4x, DATA_SeqTrain20x, Data_RandTest4x, Data_RandTrain20x}

disp('Reading SO-ET dataset ...');
[H1, tactlbls1, tlbls1, tcts1] = getData(DATA_SeqTrain20x);
disp('Reading RO-UT dataset ...');
[H2, tactlbls2, tlbls2, tcts2] = getData(Data_RandTrain20x);
disp('All dataset are read.\nGenerating the figures...');



disp('Pitch:');
plotApeSeq(H1, wipt, 50, 25, 3, 10, 1, tactlbls1, tcts1, 1, 'Pitch'); %tr-seq-pitch : 10.4758
plotApeSeq(H2, wipt, 50, 25, 2,  5, 1, tactlbls2, tcts2, 0, 'Pitch'); %tr-rand-pitch : 100-88.263986 = 11.736014

disp('Yaw:');
plotApeSeq(H1, wipt, 50, 33, 30,12, 2, tactlbls1, tcts1, 1, 'Yaw'); %tr-seq-yaw : 13.4722
plotApeSeq(H2, wipt, 50, 27, 30, 6, 2, tactlbls2, tcts2, 0, 'Yaw'); %tr-rand-yaw :  100-88.835 = 11.1647

disp('Roll:');
plotApeSeq(H1, wipt, 50, 34, 1, 11, 3, tactlbls1, tcts1, 1, 'Roll'); %tr-seq-roll : 10.8781
plotApeSeq(H2, wipt, 50, 34, 15, 5, 3, tactlbls2, tcts2, 0, 'Roll'); %tr-rand-roll : 100-88.8079906 = 11.1920094

disp('Three-fingers:');
plotApeSeq(H1, wipt, 50, 21, 30, 9, 4, tactlbls1, tcts1, 1, '3F'); %tr-seq-3f : 10.7487
plotApeSeq(H2, wipt, 50, 21, 30, 4, 4, tactlbls2, tcts2, 0, '3F'); %tr-rand-3f : 100-91.148603 = 8.851397

disp('V-sign:');
plotApeSeq(H1, wipt, 50, 15, 2, 10, 5, tactlbls1, tcts1, 1, 'V'); %tr-seq-v : 12.5202
plotApeSeq(H2, wipt, 50, 14, 2,  4, 5, tactlbls2, tcts2, 0, 'V'); %tr-rand-v : 100-89.8227462 = 10.1772538

disp('OK-sign:');
plotApeSeq(H1, wipt, 50, 34, 2, 10, 6, tactlbls1, tcts1, 1, 'OK'); %tr-seq-ok : 11.6288
plotApeSeq(H2, wipt, 50, 30, 13, 4, 6, tactlbls2, tcts2, 0, 'OK'); %tr-rand-ok : 100-91.0169136 = 8.9830864


function plotApeSeq(H, wipt, windowSize, startPC, numPC, ...
    numRepAvg, tag, true_action_labels, true_counts, toDraw, name)
    [Pw, ~] = wipt.getAveragePCASeries(H, startPC, numPC, numRepAvg, windowSize);
    [count_acc, bars]=wipt.CountReps(Pw, tag, true_action_labels, true_counts);

    diffs = abs(bars(:,2)-bars(:,1)); % |true-pred|
    ape = zeros(height(diffs), 1);
    for i=1:height(diffs)
        ape(i) = (100 * diffs(i)) / bars(i, 2);
    end
    fprintf('%f , %f\n', mean(ape), count_acc);
    ape = sort(ape);

    figName = sprintf('%d. %s', tag, name);
    marker = '-o';
    if toDraw==1
        figure('Name', figName, 'units','points', 'Position', [500,500,160,135]);
        hold on; 
        marker = '-^';
    end
    plot(ape, marker, 'LineWidth', 1.5, 'MarkerIndices', [4,8,12,16,20]);
    set(gca, 'FontSize', 11); % , 'xlim', [0, 20], 'ylim', [0,50]
    xlbl = sprintf('Sorted Segments (%s)', name);
    xlabel(xlbl); ylabel('APE (%)');
end

function [H, true_action_labels, true_labels, true_counts] = getData(DATA)
    constants; 
    fullcsv = readmatrix(csvfiles(DATA));
    csirange = csiranges(2*DATA-1):csiranges(2*DATA); fullcsv = fullcsv(csirange, :);
    H = fullcsv(:, 3:66); l = height(H);
    
    % 'PITCH': 1, 'YAW': 2, 'ROLL': 3, 'THREE_FINGERS': 4, 'V_SIGN': 5, 'OK_SIGN': 6, 'STATIC': 7, 'none': 7
    true_action_labels = fullcsv(:,2); true_labels = zeros(l, 1);
    for i=1:1:l
        l = true_action_labels(i);
        if l==0; l=5; elseif l==7; l=TAG_NONACT; else; l = TAG_ACT; end
        true_labels(i) = l;
    end
    true_counts = split(true_counts_all, '#');
    true_counts = str2num(string(true_counts(DATA)))'; %#ok<ST2NM> 
end


