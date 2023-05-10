clc; clear; close all; 
constants; wipt = WiPT;
% Select the dataset
% DATA_SeqTest4x;   % DATA_SeqTrain20x
% Data_RandTest4x;  % Data_RandTrain20x
DATA = Data_RandTrain20x;


% Construct the CSI series from the selected experiment's CSV file
fullcsv = readmatrix(csvfiles(DATA));
csirange = csiranges(2*DATA-1):csiranges(2*DATA);%  1:251460;%   251460:315102;% 64502:315102;% 228660:540705;% 167714:228660;% 
fullcsv = fullcsv(csirange, :);
% fullcsv = {fullcsv(1:64588, :), fullcsv(79142:141484, :), ...
%     fullcsv(158891:221638, :), fullcsv(236212:298554, :)};
% fullcsv = {fullcsv(64588:79142, :), fullcsv(141484:158891, :), ...
%     fullcsv(221638:236212, :), fullcsv(298554:315102, :)};
%-fullcsv = fullcsv(213560:268684, :); % 14-april mid-4x-split for validation
%-fullcsv = {fullcsv(89187:213560, :), fullcsv(268684:389840, :)}; % 14-april 1st 10x & last 10x for training
%-fullcsv = cat(1, fullcsv{:});
H = fullcsv(:, 3:66); l = height(H);

global TAG_ACT TAG_NONACT; %#ok<GVMIS> 

WARMUP_OFFSET = 0;
COOLDOWN_OFFSET = 0;

% Prepare the true-labels
% 'PITCH': 1, 'YAW': 2, 'ROLL': 3, 'THREE_FINGERS': 4, 'V_SIGN': 5, 'OK_SIGN': 6, 'STATIC': 7, 'none': 7
true_action_labels = fullcsv(:,2); true_labels = zeros(l, 1);
for i=1:1:l
    l = true_action_labels(i);
    if l==0; l=5; elseif l==7; l=TAG_NONACT; else; l = TAG_ACT; end
    true_labels(i) = l;
end


% Prepare true-counts
true_counts = split(true_counts_all, '#');
true_counts = str2num(string(true_counts(DATA)))'; %#ok<ST2NM> 

% Use either full-param-space-search approach or hyper-param-optim. approach

% full_param_space_search;

% METHOD_TO_OPTIMIZE = { METHOD_OVERALL_AVG_THRESHOLD,
    % METHOD_LOCAL_AVG_THRESHOLD, METHOD_LOOK_AHEAD_DELTA,
    % METHOD_BASE_COUNTING}
METHOD_TO_OPTIMIZE = METHOD_BASE_COUNTING;
MAX_ITER = 300;
MAX_TIME = 20 * 60; % in seconds
if METHOD_TO_OPTIMIZE == METHOD_BASE_COUNTING
    for t=[PITCH YAW ROLL THREE_FINGERS_OR_FIST V_SIGN OK_SIGN]
        tag = t;
        hyperparam_optim;
    end
else
    for md=[METHOD_LOCAL_AVG_THRESHOLD METHOD_OVERALL_AVG_THRESHOLD]
        METHOD_TO_OPTIMIZE = md;
        hyperparam_optim;
    end
end


% plotSegmentedResult(H, wipt, 7,10,22, 100,22, 5, METHOD_LOCAL_AVG_THRESHOLD, true_labels, 0, 0); % Rand-train-LA
% plotSegmentedResult(H, wipt, 3, 3, 6, 100,19, 4, -1, METHOD_LOCAL_AVG_THRESHOLD, true_labels, 0, 0); % Seq-train-LA
% plotSegmentedResult(H, wipt, 3, 3, 3, 100,23, 2, METHOD_LOCAL_AVG_THRESHOLD, true_labels, 0, 0); % Seq-train-LA
% plotSegmentedResult(H, wipt, 5,1,15,  150,15, 4, METHOD_LOCAL_AVG_THRESHOLD, true_labels, 1800, 800); % Rand-test-LA
% plotSegmentedResult(H, wipt, 5,4,20,  50, -1,-1, -1, METHOD_BASE_COUNTING, true_labels, 300, 200); % Rand-test-Count
% plotSegmentedResult(H, wipt, 5,8,7,   50, -1,-1, -1, METHOD_BASE_COUNTING, true_labels, 300, 200); % Seq-test/train-Count

% plotCountsByGroundTruths(H, wipt, 50, 25, 3, 10, 1, true_action_labels, true_counts); %tr-seq-pitch : 10.4758
% plotCountsByGroundTruths(H, wipt, 50, 33, 30,12, 2, true_action_labels, true_counts); %tr-seq-yaw : 13.4722
% plotCountsByGroundTruths(H, wipt, 50, 34, 1, 11, 3, true_action_labels, true_counts); %tr-seq-roll : 10.8781
% plotCountsByGroundTruths(H, wipt, 50, 21, 30, 9, 4, true_action_labels, true_counts); %tr-seq-3f : 10.7487
% plotCountsByGroundTruths(H, wipt, 50, 15, 2, 10, 5, true_action_labels, true_counts); %tr-seq-v : 12.5202
% plotCountsByGroundTruths(H, wipt, 50, 34, 2, 10, 6, true_action_labels, true_counts); %tr-seq-ok : 11.6288
% plotCountsByGroundTruths(H, wipt, 50, 32, 12,  5, 1, true_action_labels, true_counts); %tr-rand-pitch : 100-91.155927 = 8.8441
% plotCountsByGroundTruths(H, wipt, 50, 27, 30, 6, 2, true_action_labels, true_counts); %tr-rand-yaw :  100-88.835 = 11.1647
% plotCountsByGroundTruths(H, wipt, 50, 34, 15, 5, 3, true_action_labels, true_counts); %tr-rand-roll : 100-88.8079906 = 11.1920094
% plotCountsByGroundTruths(H, wipt, 50, 21, 30, 4, 4, true_action_labels, true_counts); %tr-rand-3f : 100-91.148603 = 8.851397
% plotCountsByGroundTruths(H, wipt, 50, 14, 2,  4, 5, true_action_labels, true_counts); %tr-rand-v : 100-89.8227462 = 10.1772538
% plotCountsByGroundTruths(H, wipt, 50, 30, 13, 4, 6, true_action_labels, true_counts); %tr-rand-ok : 100-91.0169136 = 8.9830864

"Batch Completed" %#ok<NOPTS>


function plotDurations(X)
global TAG_ACT;
    l = height(X);
    durations = zeros(1,2);
    ct = 1;

    i=1;
    while i<=l
        if X(i)==TAG_NONACT
            last = i;
            for j=i:1:l
                if X(j)~=TAG_NONACT
                    last = j;
                    break;
                end
            end
            durations(ct, 1) = i + (last-i)/2;
            durations(ct, 2) = (last-i)/100; % (nFrames)/100Hz = seconds
            ct = ct +1;
            i=last+1;
        else 
            i = i + 1;
        end
    end

    figure; hold on;
    xaxis = durations(:,1);
    dd = durations(:,2); dd(dd<1) = 12;
    ps=plot(xaxis,dd, 'LineWidth', 2); set(gca, 'FontSize', 14);
    % ps = scatter(xaxis, ); set(gca, 'FontSize', 14);
    xlabel('Activity Segment Number'); ylabel('Action Duration (seconds)');
    pp = plot(xaxis, 12*ones(ct-1, 1), 'Color','red', 'LineWidth', 2);set(gca, 'FontSize', 14);
    legend([ps;pp], 'Predicted Duration', 'Actual Duration', ...
        'FontSize', 14, 'Orientation','horizontal');
end


function plotDurationError(X, true_labels)
global TAG_ACT;
    l = height(X);
    durations = zeros(1,2);
    ct = 1;

    i=1;
    while i<=l
        if X(i)==TAG_ACT
            last = i;
            for j=i:1:l
                if X(j)~=TAG_ACT
                    last = j;
                    break;
                end
            end
            durations(ct, 1) = i + (last-i)/2;
            durations(ct, 2) = (last-i)/100; % nFrames/100Hz = seconds
            ct = ct + 1;
            i = last + 1;
        else 
            i = i + 1;
        end
    end

    true_durations = zeros(1,2);
    ct = 1;
    i=1; l=height(true_labels);
    while i<=l
        if true_labels(i)==TAG_ACT
            last = i;
            for j=i:1:l
                if true_labels(j)~=TAG_ACT
                    last = j;
                    break;
                end
            end
            true_durations(ct, 1) = i + (last-i)/2;
            true_durations(ct, 2) = (last-i)/100; % nFrames/100Hz = seconds
            ct = ct + 1;
            i = last + 1;
        else 
            i = i + 1;
        end
    end

    figure; hold on;
    xaxis = durations(:,1);
    dd = durations(:,2); % dd(dd<1) = 12;
    tdd = true_durations(:,2); % dd(dd<1) = 12;
    ps=plot(xaxis, dd, 'LineWidth', 2); set(gca, 'FontSize', 14);
    % ps = scatter(xaxis, ); set(gca, 'FontSize', 14);
    xlabel('Activity Segment Number'); ylabel('Action Duration (seconds)');
    pp = plot(xaxis, tdd, 'Color','red', 'LineWidth', 2);set(gca, 'FontSize', 14);
    legend([ps;pp], 'Predicted Duration', 'Actual Duration', ...
        'FontSize', 14, 'Orientation','horizontal');

    diff = true_durations(:,2) - durations(:,2);
    figure; hold on; plot(diff); plot(true_durations(:,2));
end

function plotCountsByGroundTruths(H, wipt, windowSize, startPC, numPC, ...
    numRepAvg, tag, true_action_labels, true_counts)
    [Pw, ~] = wipt.getAveragePCASeries(H, startPC, numPC, numRepAvg, windowSize);
    [count_acc, bars]=wipt.CountReps(Pw, tag, true_action_labels, true_counts);
    disp(count_acc);

    diffs = abs(bars(:,2)-bars(:,1)); % |true-pred|
    ape = zeros(height(diffs), 1);
    for i=1:height(diffs)
        ape(i) = (100 * diffs(i)) / bars(i, 2);
    end
    disp(mean(ape));
    ape = sort(ape);

    figName = sprintf('Action-%d', tag);
    figure('Name', figName,'units','points', 'Position', [500,500,132,136]);
    hold on; plot(ape, '^');
    set(gca, 'FontSize', 10, 'xlim', [5, 15], 'ylim', [5,15]);
    xlabel('Sorted Segments'); ylabel('APE (%)');

%     figure('Name', figName, 'units','points', 'Position', [500,500,132,136]);
%     bar(bars);
%     set(gca, 'FontSize', 10);
%     xlabel('Segment'); ylabel('Counts');

    figure('Name', figName,'units','points', 'Position', [500,500,132,136]);
%     model = fitlm(bars(:,2), bars(:,1));
%     for x=0:20
%         h = height(bars) + 1;
%         y = predict(model, x);
%         bars(h,2) = y;
%         bars(h,1) = y;
%     end
    plot(fitlm(bars(:,2), bars(:,1)));
    set(gca, 'FontSize', 10, 'xlim', [5, 15], 'ylim', [5,20]);
    xlabel('True Counts'); ylabel('Predicted Counts');
end

function plotSegmentedResult(H, wipt, startPC, numPC, numRepAvg, ...
    windowSize, windowMultiplier, smoothingFactor, ...
    method, true_labels, warmup_offset, cooldown_offset)
constants;
    [Pw, ~] = wipt.getAveragePCASeries(H, startPC, numPC, numRepAvg, ...
        windowSize); % 2,51,460

    switch method
        case METHOD_OVERALL_AVG_THRESHOLD
            X = wipt.SegmentByOverallAvg(Pw, windowSize, smoothingFactor);
%         case METHOD_LOOK_AHEAD_DELTA
%             X = wipt.SegmentByLookAheadDelta(Pw, windowSize, ...
%                 windowMultiplier, smoothingFactor, delta);
        case METHOD_LOCAL_AVG_THRESHOLD
            X = wipt.SegmentByLocalAvg(Pw, windowSize, ...
                windowMultiplier, smoothingFactor);
        case METHOD_BASE_COUNTING
            % [count_acc, bars]=wipt.CountReps(Pw, TAG_ACT, true_labels, true_counts);
            return;
    end
    pred_segs = getSegments(X, warmup_offset, cooldown_offset);
    true_segs = getSegments(true_labels, warmup_offset, cooldown_offset);
    % diffs = true_segs - pred_segs;


    bestF1=0; bestAcc=0;
    for fgv=0:forgiveness_inc:forgiveness_max
        [acc, f1] = wipt.calcAccuracyF1(X, true_labels, fgv, ...
            warmup_offset, cooldown_offset);
        disp([acc f1 fgv]);
        if f1>bestF1; bestF1=f1; end
        if acc>bestAcc; bestAcc=acc; end
    end
    Pw = Pw(~all(isnan(Pw),2));

    fig = figure;
    hold on;
    ppw = plot((Pw-mean(Pw))/5, 'Color','blue', 'LineWidth', 1);
    set(gca, 'FontSize', 14);
    ptl = plot((true_labels-mean(true_labels)), 'color','black', 'LineWidth', 1);
    xlabel('CSI Frame Number'); % ylabel('Averaged PCA');
    px = plot((X-mean(X))/2, 'Color','red', 'LineWidth', 1);

    leg = legend([ppw;ptl;px], 'Averaged PCs', 'True Segments', 'Predicted Segments', ...
        'FontSize', 14,'Orientation','horizontal');%,'Location','south'
    % func_gen_multifig_legend(fig, leg, "seg__legend", "epsc");
end



%%% Returns 2-D array consisting the starting & ending indices of the
%%% segments as specified by per frame in the input `tags` 1-D array
function segments = getSegments(tags, warmup_offset, cooldown_offset)
global TAG_NONACT; %#ok<GVMIS> 
    l = height(tags);
    segments = zeros(1,2);

    ct = 1; i = 1+warmup_offset;
    while i<=l-cooldown_offset
        if tags(i)==TAG_NONACT
            if ct>1; segments(ct-1, 2) = i-1; end
            last = i;
            for j=i:1:l
                if tags(j)~=TAG_NONACT
                    last = j-1;
                    break;
                end
            end
            segments(ct, 1) = i;
            segments(ct, 2) = last;
            segments(ct + 1, 1) = last+1;
            % segments(ct+1, 2) = ??;
            ct = ct + 2;
            i=last+1;
        else 
            i = i + 1;
        end
    end
    segments(ct-1, 2) = l-cooldown_offset;
end


