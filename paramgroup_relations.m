%clc; clear; close all; 
constants; wipt = WiPT;  
DATA = Data_RandTrain20x; % DATA_5JanTest4xRandVert;
fullcsv = readmatrix(csvfiles(DATA));
csirange = csiranges(2*DATA-1):csiranges(2*DATA); fullcsv = fullcsv(csirange, :);
% fullcsv = {fullcsv(89187:213560, :), fullcsv(268684:389840, :)}; % 14-april 1st 10x & last 10x for training
% fullcsv = cat(1, fullcsv{:});
H = fullcsv(:, 3:66); l = height(H);
global TAG_ACT TAG_NONACT; %#ok<GVMIS> 

% Prepare the true-labels
true_action_labels = fullcsv(:,2); true_labels = zeros(l, 1);
for i=1:1:l
    l = true_action_labels(i);
    if l==0; l=5; elseif l==7; l=TAG_NONACT; else; l = TAG_ACT; end
    true_labels(i) = l;
end

% Prepare true-counts
true_counts = split(true_counts_all, '#');
true_counts = str2num(string(true_counts(DATA)))'; %#ok<ST2NM> 

%% Forgiveness Value for Segmentation
% fgv_plot = zeros(forgiveness_max/forgiveness_inc+1, 1); f=1;
% for fgv=0:forgiveness_inc:forgiveness_max
%     [Pw, ~] = wipt.getAveragePCASeries(H, 2, 2, 5, 100);
%     X = wipt.SegmentByLocalAvg(Pw, 100, 19, 2);
% 
%     [acc, f1] = wipt.calcAccuracyF1(X, true_labels, fgv, 0, 0);
%     disp([acc f1 fgv]);
%     fgv_plot(f)=acc;
%     f = f + 1;
% end
% figure; plot(fgv_plot);
%%

%% Merging Factor for Segmentation
% merging_acc = zeros(21, 1);
% for y=0:1:20
%     [Pw, ~] = wipt.getAveragePCASeries(H, 2, 2, 5, 100);
%     Pw = Pw(~all(isnan(Pw),2));
%     X = wipt.SegmentByLocalAvg(Pw, 100, 19, y);
%     
%     [acc, f1] = wipt.calcAccuracyF1(X, true_labels, 60, 0, 0);
%     disp([y acc f1]);
%     merging_acc(y+1) = acc;
% end
% figure; plot(merging_acc);
%%

%% Num. of Repetitive Averages for Counting (Pitch Action)
% n_count_acc = zeros(24, 1);
% for n=1:1:24
%     [Pw, ~] = wipt.getAveragePCASeries(H, 25, 32, n, 50);
%     [count_acc, ~]=wipt.CountReps(Pw, 1, true_action_labels, true_counts);
%     
%     disp([n count_acc]);
%     n_count_acc(n) = count_acc;
% end
% figure; plot(n_count_acc);
%%

%% Num. of Repetitive Averages for Segmentation
% nlist = zeros(24, 1);
% for n=1:1:24
%     [Pw, ~] = wipt.getAveragePCASeries(H, 2, 2, n, 100);
%     Pw = Pw(~all(isnan(Pw),2));
%     X = wipt.SegmentByLocalAvg(Pw, 100, 19, 2);
%     
%     [acc, f1] = wipt.calcAccuracyF1(X, true_labels, 60, 0, 0);
%     disp([n acc f1]);
%     nlist(n) = acc;
% end
% figure; plot(nlist);
%%


%% PC group for Counting
pgrid = zeros(32, 32);
for ps=1:1:32
    for p=1:1:32
        [Pw, ~] = wipt.getAveragePCASeries(H, ps, p, 7, 50);
        [count_acc, ~]=wipt.CountReps(Pw, 1, true_action_labels, true_counts);

        disp([ps p count_acc]);
        pgrid(ps, p) = count_acc;
    end
end
imagesc(pgrid);colorbar;
%%

%% PC group for Segmentation
% pgrid = zeros(32, 32);
% for ps=1:1:32
%     for p=1:1:32
%         [Pw, ~] = wipt.getAveragePCASeries(H, ps, p, 5, 100);
%         Pw = Pw(~all(isnan(Pw),2));
%         X = wipt.SegmentByLocalAvg(Pw, 100, 21, 2);
%  
%         [acc, f1] = wipt.calcAccuracyF1(X, true_labels, 60, 0, 0);
%         disp([ps p acc f1]);
%         pgrid(ps, p) = acc;
%     end
% end
% imagesc(pgrid);colorbar;
%% 


%% Window size & multiplier combo for Counting Pitch Action Rep.
wlist = zeros(32, 1);
for w=50:50:1600
    [Pw, ~] = wipt.getAveragePCASeries(H, 25, 32, 7, w);
    [count_acc, ~]=wipt.CountReps(Pw, 1, true_action_labels, true_counts);

    disp([w count_acc]);
    wlist(w/50) = count_acc;
end
figure('Position', [500,500,132,136]); plot(wlist);
%%

%% Window size & multiplier combo for Segmentation
% wgrid = zeros(32, 40);
% for w=50:50:1600
%     for wm=1:1:40
%         [Pw, ~] = wipt.getAveragePCASeries(H, 3, 3, 3, w);
%         X = wipt.SegmentByLocalAvg(Pw, w, wm, 2);
%         [accuracy, f1] = wipt.calcAccuracyF1(X, true_labels, 40, 0, 0);
%     
%         disp([w wm accuracy f1]);
%         wgrid(w/50, wm) = accuracy;
%     end
% end
% imagesc(wgrid);colorbar;
%%
