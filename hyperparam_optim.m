% This script is supposed to be called only from/after the `main.m` script. %

[counting_acc_file, avgpca_acc_file, localavg_acc_file, look_ahead_acc_file] = ...
    wipt.initFiles(string(datetime('now', 'Format', 'dMMMyy_HHmmss')), ...
    METHOD_TO_OPTIMIZE);

w = optimizableVariable('windowSize50Multiple', [1, 3], 'Type','integer'); % window_size = w*50
wm = optimizableVariable('windowMultiplier', [1, 32], 'Type','integer');
pc_start = optimizableVariable('startPC', [1, 32], 'Type','integer');
p = optimizableVariable('numPC', [1, 32], 'Type','integer');
n = optimizableVariable('numRepAvg', [1, 24], 'Type','integer');
y = optimizableVariable('smoothingFactor', [0, 20], 'Type','integer');
d = optimizableVariable('delta', [0.05, 0.5], 'Type','real');


switch METHOD_TO_OPTIMIZE
    case METHOD_BASE_COUNTING
        obj_fun = @(vars) objective_fun_counting(vars, H, wipt, tag, ...
            true_action_labels, true_counts, counting_acc_file);
        var_arr = [pc_start, p, n, w];

    case METHOD_OVERALL_AVG_THRESHOLD
        obj_fun = @(vars) objective_fun_overall_avg(vars, ...
            H, wipt, true_labels, avgpca_acc_file);
        var_arr = [pc_start, p, n, w, y];

    case METHOD_LOOK_AHEAD_DELTA
        obj_fun = @(vars) objective_fun_lookahead(vars, ...
            H, wipt, true_labels, look_ahead_acc_file);
        var_arr = [pc_start, p, n, w, wm, y, d];

    otherwise
        obj_fun = @(vars) objective_fun_local_avg(vars, ...
            H, wipt, true_labels, localavg_acc_file);
        var_arr = [pc_start, p, n, w, wm, y];
        
end
results = bayesopt(obj_fun, var_arr, 'MaxObjectiveEvaluations', ...
    MAX_ITER, 'MaxTime',MAX_TIME, 'UseParallel',false);


% Get optimal params from results & Plot X, Pw & true-labels
x = table2array(bestPoint(results, 'Criterion', 'min-observed'));
if METHOD_TO_OPTIMIZE==METHOD_BASE_COUNTING
%     plotCountingResult(H, wipt, tag, x(1), x(2), x(3), x(4)*50, ...
%         true_action_labels, true_counts);
else
    if width(x)>5; tmp=x(5);x(5)=x(6);x(6)=tmp; end
    plotOptimalResult(H, wipt, x, METHOD_TO_OPTIMIZE, true_labels, ...
        WARMUP_OFFSET, COOLDOWN_OFFSET);
end


function plotOptimalResult(H, wipt, res, method, true_labels, ...
    warmup_offset, cooldown_offset)
constants;
    o_w = res(4) * 50;
    o_y = res(5);
    o_spc = res(1);
    o_npc = res(2);
    o_nr = res(3);
    [Pw, ~] = wipt.getAveragePCASeries(H, o_spc, o_npc, o_nr, o_w);

    switch method
        case METHOD_OVERALL_AVG_THRESHOLD
            X = wipt.SegmentByOverallAvg(Pw, o_w, o_y);
        case METHOD_LOOK_AHEAD_DELTA
            o_wm = res(6); o_d = res(7);
            X = wipt.SegmentByLookAheadDelta(Pw, o_w, o_wm, o_y, o_d);
        case METHOD_LOCAL_AVG_THRESHOLD
            o_wm = res(6);
            X = wipt.SegmentByLocalAvg(Pw, o_w, o_wm, o_y);
    end
    bestF1=0; bestAcc=0;
    for fgv=0:forgiveness_inc:forgiveness_max
        [acc, f1] = wipt.calcAccuracyF1(X, true_labels, fgv, ...
            warmup_offset, cooldown_offset);
        disp([acc f1 fgv]);
        if f1>bestF1; bestF1=f1; end
        if acc>bestAcc; bestAcc=acc; end
    end
    Pw = Pw(~all(isnan(Pw),2));

    figure;
    hold on;
    ppw = plot((Pw-mean(Pw)) * 7, 'Color','blue', 'LineWidth', 1);
    set(gca, 'FontSize', 14);
    ptl = plot(true_labels-mean(true_labels), 'color','black', 'LineWidth', 1);
    xlabel('CSI Frame Number'); % ylabel('Averaged PCA');
    px = plot(X-mean(X), 'Color','red', 'LineWidth', 1);
    legend([ppw;ptl;px], 'Averaged PCA', 'True Segments', 'Predicted Segments', ...
        'FontSize', 14,'Orientation','horizontal');%,'Location','south'
end

function plotCountingResult(H, wipt, tag, startPC, numPC, numRepAvg, ...
    windowSize, trueLabels, trueCounts)
    [Pw, ~] = wipt.getAveragePCASeries(H, startPC, numPC, ...
        numRepAvg, windowSize);
    [~, counts_for_bar] = wipt.CountReps(Pw, tag, trueLabels, trueCounts);

    figure; 
    bar_counts=bar(counts_for_bar); set(gca, 'FontSize', 14);
    xlabel('Activity Segment Number'); ylabel('Rep. Counts' + tag);
    legend(bar_counts, 'Predicted Counts', 'True Counts', ...
        'FontSize', 14, 'Orientation','horizontal');

    figure; hold on;
    pp = plot(counts_for_bar(:,1), 'color','red', 'LineWidth', 1);
    set(gca, 'FontSize', 14); xlabel('CSI Segments'); ylabel('Action Counts');
    pt = plot(counts_for_bar(:,2), 'color','blue', 'LineWidth', 1);
    legend([pp;pt], 'Predicted Counts', 'True Counts', ...
        'FontSize', 14, 'Orientation','horizontal');
end

% Obj. Function Signature: [objective,coupledconstraints,userdata] = fun(x)
function [objective] = objective_fun_overall_avg(x, H, wipt, ...
    true_labels, avgpca_acc_file)
    w = x.windowSize50Multiple * 50;
    y = x.smoothingFactor;
    [Pw, ~] = wipt.getAveragePCASeries(H, x.startPC, x.numPC, x.numRepAvg, w);
    Xp = wipt.SegmentByOverallAvg(Pw, w, y);
    f1score = wipt.logAccuracy(Xp, true_labels, ...
        [w x.startPC x.numPC x.numRepAvg y], avgpca_acc_file);
    objective = 1 - f1score;
end

function [objective] = objective_fun_local_avg(x, H, wipt, ...
    true_labels, localavg_acc_file)
    w = x.windowSize50Multiple * 50;
    wm = x.windowMultiplier;
    y = x.smoothingFactor;
    [Pw, ~] = wipt.getAveragePCASeries(H, x.startPC, x.numPC, x.numRepAvg, w);
    Xv = wipt.SegmentByLocalAvg(Pw, w, wm, y);
    f1score = wipt.logAccuracy(Xv, true_labels, ...
        [w wm x.startPC x.numPC x.numRepAvg y], localavg_acc_file);
    objective = 1 - f1score;
end

function [objective] = objective_fun_lookahead(x, H, wipt, ...
    true_labels, look_ahead_acc_file)
    w = x.windowSize50Multiple*50;
    wm = x.windowMultiplier;
    y = x.smoothingFactor;
    [Pw, ~] = wipt.getAveragePCASeries(H, x.startPC, x.numPC, x.numRepAvg, w);
    Xl = wipt.SegmentByLookAheadDelta(Pw, w, wm, y, x.delta);
    f1score = wipt.logAccuracy(Xl, true_labels, ...
        [w wm x.startPC x.numPC x.numRepAvg y x.delta], look_ahead_acc_file);
    objective = 1 - f1score;
end

function [objective] = objective_fun_counting(x, H, wipt, tag, ...
    true_action_labels, true_counts, counting_acc_file)
    w = x.windowSize50Multiple*50;
    [Pw, ~] = wipt.getAveragePCASeries(H, x.startPC, x.numPC, x.numRepAvg, w);
    [counting_accuracy, ~] = wipt.CountReps(Pw, tag, true_action_labels, true_counts);
    wipt.appendMetrics([w, tag, x.startPC, x.numPC, x.numRepAvg, ...
        counting_accuracy], counting_acc_file);
    objective = 100 - counting_accuracy;
end



