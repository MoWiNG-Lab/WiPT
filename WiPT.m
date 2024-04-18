classdef WiPT
    methods
function [Pw, explained] = getAveragePCASeries(~, H, pc_start, p, n, w)
    l = height(H);
    [coeff, ~, ~, ~, explained] = pca(H);
    Pw_temp = H * coeff(:, pc_start:pc_start+p-1); % retrieve PCs from H & resp. co-efficients
    if p>1
        x = Pw_temp';
        Pw_temp = (mean(x))'; %#ok<UDIM> 
    end
    Pw = zeros(height(Pw_temp), 1);
    for etao = 1:1:n
        for i=1:1:(l-w/2)
            a = max(1, i-w/2);
            b = min(l, i+w/2);
            Pw(i) = mean(Pw_temp(a:b));
            % Pw(min(l, i+w/2)) = mean(Pw_temp(i:min(l, i+w)));
        end
        Pw_temp = Pw;
    end
    Pw(isnan(Pw)) = mean(Pw(~isnan(Pw)));
end


function [X] = smoothByMergingSmallSegments(~, w, X, gamma)
global TAG_ACT TAG_NONACT; %#ok<*GVMIS> 
    l = height(X);
    for B=w:w:gamma*w+1
        for s_c = B+1 : 1 : l-2*B+1
                                e_c = s_c + B - 1;
            s_p = s_c - B;      e_p = s_p + B - 1;
            s_n = e_c + 1;      e_n = s_n + B - 1;
            
            cur_label = mean(X(s_c:e_c));
            prev_label = mean(X(s_p:e_p));
            next_label = mean(X(s_n:e_n));

%             if s_c>=3674 && e_c<3976 && B==300
%                 disp("Mergeable Window Found");
%             end

            %{
            %If any of the 3 comparable labels are mixed with activity &
            % non-activity, then we'll ignore the label. Otherwise, we'll
            % compare these 3 windows to merge the current label for 
            % 2 cases : (_-_) or (-_-).
            % Question: If the current-label is mixed, then are we
            % trying to merge a too long-window, so that we should stop 
            % with this smoothing-factor? --> NO. Anti-case: __-- can be
            % windowed as _-, so that the avg. of this current-window gets
            % mixed.
            %}
            if (cur_label==TAG_ACT || cur_label==TAG_NONACT) && ...
                (prev_label==TAG_ACT || prev_label==TAG_NONACT) && ...
                (next_label==TAG_ACT || next_label==TAG_NONACT) && ...
                (cur_label~=prev_label && cur_label~=next_label)
                    X(s_c:e_c) = prev_label;
            end 
        end
    end
end


function [Xp] = SegmentByOverallAvg(o, Pw, w, gamma)
global TAG_ACT TAG_NONACT;
    l = height(Pw);
    Xp = zeros(l ,1);
    tau = mean(Pw);
    for i=1:w:l
        s = max(1, i-w/2); e=min(l, i+w/2);
        if mean(Pw(s:e))>=tau
            Xp(s:e) = TAG_NONACT;
        else
            Xp(s:e) = TAG_ACT;
        end
    end
    Xp = smoothByMergingSmallSegments(o, w, Xp, gamma);
end


function [Xv] = SegmentByLocalAvg(o, Pw, w, psi, gamma)
global TAG_ACT TAG_NONACT;
    l = height(Pw);
    wm = psi * w; % larger windows
    Xv = zeros(l ,1);
    taus = zeros(l ,1);
    for i=1:wm:l
        s = max(1, i-wm/2); e=min(l, i+wm/2);
        taus(s:e) = mean(Pw(s:e));
    end
    for i=1:w:l
        s = max(1, i-w/2); e=min(l, i+w/2);
        if mean(Pw(s:e))>=taus(int32(floor((s+e)/2)))
            Xv(s:e) = TAG_NONACT;
        else
            Xv(s:e) = TAG_ACT;
        end
    end
    Xv = smoothByMergingSmallSegments(o, w, Xv, gamma);
end


function [Xl] = SegmentByLookAheadDelta(o, Pw, w, psi, gamma, delta)
global TAG_ACT TAG_NONACT;
    l = height(Pw);
    wm = psi * w; % larger windows
    Xl = zeros(l ,1);
    
    label = TAG_NONACT;
    for i=1:w:l
        s = max(1, i-w/2);  e = min(l, i+w/2);
        sl= max(1, i-wm/2); el = min(l, i+wm/2);%#ok<NASGU> %TODO: explore by changing these to sl=s & el=sl+wm;
        mid = int32((s+e)/2);

        d = mean(Pw(mid:el)) - mean(Pw(mid:e));
        if abs(d)>=delta
            if(d>0)
                label = TAG_NONACT;
            else
                label = TAG_ACT;
            end
        end
        Xl(mid:el) = label;
    end
    Xl = smoothByMergingSmallSegments(o, w, Xl, gamma);
end


function [counting_accuracy, count_bars] = CountReps(~, Pw, tag, labels, true_counts)
    l = height(labels);
    n = height(true_counts);
    pred_counts = zeros(n, 1);

    ct = 1; tct = 0;
    i=1;
    trimmed_true_counts = zeros(n, 1);
    while i < l
        if labels(i)==7 && labels(i)~=labels(i+1)
            % going from STATIC to any non-STATIC, so found the next true-segment
            tct = tct + 1;
        end
        if labels(i)==tag
            last = i;
            for j=i:1:l
                if labels(j)~=tag
                    last = j;
                    break;
                end
            end
            curr_seg = Pw(i:last);
            pred_counts(ct) = sum(islocalmin(curr_seg)) + sum(islocalmax(curr_seg));
            trimmed_true_counts(ct) = true_counts(max(1, tct));
            ct = ct +1;
            i=last+1;
        else 
            i = i + 1;
        end
    end
    pred_counts = pred_counts(1:ct-1);
    trimmed_true_counts = trimmed_true_counts(1:ct-1);


    count_accs = zeros(ct-1, 1);
    count_devs = zeros(ct-1, 1);
    count_bars = zeros(ct-1, 2);
    for i=1:ct-1
        count_devs(i) = abs(trimmed_true_counts(i) - pred_counts(i));
        count_accs(i) = max(0, 100 - (100.0 * count_devs(i))/trimmed_true_counts(i));
        count_bars(i,1)=pred_counts(i);
        count_bars(i,2)=trimmed_true_counts(i);
    end
    % avg_devs = mean(count_devs);
    % loss based on deviation?? --> think CLEARLY
    counting_accuracy = mean(count_accs);
end


function [counting_acc_file, avgpca_acc_file, local_avg_acc_file, look_ahead_acc_file] = initFiles(o, dateStr, singleFile)
constants;
    res_folder = fullfile(pwd, "results", dateStr);
    if ~isfolder(res_folder)
        mkdir(res_folder);
    end
    counting_acc_file = fullfile(res_folder, "counts.csv");
    avgpca_acc_file = fullfile(res_folder, "seg_avgpca.csv");
    local_avg_acc_file = fullfile(res_folder, "seg_localavg.csv");
    look_ahead_acc_file = fullfile(res_folder, "seg_lookahead.csv");
    
    headers_allavg = ["Window Size", "#StartPC", "#PC", "#RepAvg", ...
        "Smoothing Factor", "Forgiveness Offset", "Accuracy(%)", "F1-Score"];
    headers_localavg = ["Window Size", "Window Multiplier", "#StartPC", ...
        "#PC", "#RepAvg", "Smoothing Factor", "Forgiveness Offset", ...
        "Accuracy(%)", "F1-Score"];
    headers_lookahead = ["Window Size", "Window Multiplier", "#StartPC", ...
        "#PC", "#RepAvg", "Smoothing Factor", "Delta", ...
        "Forgiveness Offset", "Accuracy(%)", "F1-Score"];
    headers_counting = ["Window Size", "Action", "#StartPC", "#PC", ...
        "#RepAvg", "Counting Accuracy"];

    switch singleFile
        case METHOD_OVERALL_AVG_THRESHOLD
            o.overwriteMetrics(headers_allavg, avgpca_acc_file);
        case METHOD_LOCAL_AVG_THRESHOLD
            o.overwriteMetrics(headers_localavg, local_avg_acc_file);
        case METHOD_LOOK_AHEAD_DELTA
            o.overwriteMetrics(headers_lookahead, look_ahead_acc_file);
        case METHOD_BASE_COUNTING
            o.overwriteMetrics(headers_counting, counting_acc_file);
        otherwise
            o.overwriteMetrics(headers_allavg, avgpca_acc_file);
            o.overwriteMetrics(headers_localavg, local_avg_acc_file);
            o.overwriteMetrics(headers_lookahead, look_ahead_acc_file);
            o.overwriteMetrics(headers_counting, counting_acc_file);
    end
end


function [accuracy, f1score] = calcAccuracyF1(~, predLabels, trueLabels, ...
        fgv, warmup_offset, cooldown_offset)
global TAG_NONACT;
    tp = 0;
    tn = 0;
    fp = 0;
    fn = 0;
    l = height(trueLabels)-cooldown_offset;
    nframes = 0;
    c = fgv + 1 + warmup_offset;
    while c < l
        % Consider forgiveness range first, then count true/false pos./neg.
        e = min(l, c+2*fgv);
        if trueLabels(c)~=mean(trueLabels(c : e))
            % so, an edge is detected after `fgv` number of frames
            % forward c by 2*fgv (skipping the accuracy & F1 calculation)
            c = c + 2 * fgv;
        else
            lbt = trueLabels(c);
            lbp = predLabels(c);
            if lbt==lbp
                if lbt==TAG_NONACT; tn = tn + 1; else; tp = tp + 1; end
            else
                if lbt==TAG_NONACT; fn = fn + 1; else; fp = fp + 1; end
            end
            c = c + 1;
            nframes = nframes + 1;
        end
    end
    f1score = tp / (tp + (fp+fn)/2);
    accuracy = (100 * (tp+tn))/nframes;
end


function [bestF1] = logAccuracy(o, pred_labels, true_labels, params, file)
constants;
    bestF1 = 0;
    for fgv=0:forgiveness_inc:forgiveness_max
        [acc, f1] = o.calcAccuracyF1(pred_labels, true_labels, fgv, 0, 0);
        if f1>bestF1; bestF1 = f1; end
        o.appendMetrics([params, fgv, acc, f1], file);
    end
end


function appendMetrics(~, params, file)
    writematrix(params, file, 'Delimiter',',', 'WriteMode','append');
end


function overwriteMetrics(~, params, file)
    writematrix(params, file, 'Delimiter',',', 'WriteMode','overwrite');
end


    end
end