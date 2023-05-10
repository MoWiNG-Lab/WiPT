warning('off','all');
% fullcsv = readmatrix("wrist_finger_17Nov/csi_amplitudes.csv");  csirange=600:27500; fullcsv = fullcsv(csirange, :);
% fullcsv = readmatrix("wrist_finger_27Nov/csi_amps_horizontal.csv"); csirange = 250:49300; fullcsv = fullcsv(csirange, :);
fullcsv = readmatrix("wrist_finger_27Nov/csi_amps_vertical.csv"); csirange = 300:48600; fullcsv = fullcsv(csirange, :);
% fullcsv = readmatrix("wrist_finger_27Nov/csi_amps_train20x_vertical.csv"); csirange=3000:216000; fullcsv = fullcsv(:, :);
csicsv = fullcsv(:,3:66);

true_counts_17NovTestSeq = [5 5 5 5 5 5 5 5 5 6 5 5 5 5 5 5 5 5 5 5 5 5 5 6]';
true_counts_27NovTestRand = [5 6 6 13 10 10 13 14 16 5 5 5 9 8 8 10 12 15 6 11 9 7 10 7 5 5 7 6 6 8]';
true_counts_27NovTrain20x = []';


[coeff,score,latent,~,explained] = pca(csicsv);
% AvgPCA   --> p = avg_repi(avg_w(pca(H))) < tau
%          --> 70.34% (w=50, p=3, n=5, y=0, fg=200) --[wm=X, d=X]
%          --> 89.73% (w=150,p=1, n=3, y=0, fg=150~200)
% LocalAvg --> p = avg_repi(avg_w(pca(H))) < tau_series (=avg_mw(p)) 
%          --> 79.02% (w=150, wm=10,p=3, n=3, y=1, fg=170) --[d=X]
%          --> 90.48% (w=100, wm=1, p=1, n=5, y=0, fg=100~200)
% LookAhead--> p = avg_repi(avg_w(pca(H))); avg(larger_w)-avg(w) < delta 
%          --> 70.81% (w=50, wm=10, p=10, n=1, y=3, d=0.20, fg=200)
%          --> 63.01% (w=200, wm=3, p=10, n=1, y=2, d=0.15, fg=200)

nframes = height(csicsv);  % total num of CSI frames

%best-counting(80.9%) --> wsize=100, wmultiplier =5, num_pc=15, rep_avg = 3
%best-seq-count(72.92%)-> w=50, wm=1, p=10, n=5, \\y=0,d=0.15,fg=0

for wsize = 50 % [50 100 150 200 250]
    for wmultiplier = 10 % [1 3 5 10]
        for num_pc = 8 % [1 3 10]  % max num of PCs to be considered
            
            pcs = csicsv * coeff(:, 1:num_pc);
            if num_pc>1
                x = pcs';
                pcs_avg = (mean(x))';
            else
                pcs_avg = pcs;
            end
            
            for num_rep_avg = 13 % [1 3 5 7]
                for max_smoothing_factor = 0 % 0:1:3
                   for delta =0.15%:0.02:0.26
                        [pcs_avg_windowed, pred_counts] = main( ...
                            nframes, pcs_avg, fullcsv(:,2), csirange, ...
                            wsize, wmultiplier, num_pc, num_rep_avg, ...
                            max_smoothing_factor, delta, ... %true_counts_17NovTestSeq);
                            true_counts_27NovTestRand);
                            % true_counts_27NovTrain20x);
                    end
                end
            end
        end
    end
end



function [pcs_avg_windowed, pred_counts] = main(nframes, pcs_avg, true_action_labels, ...
    ~, window_size, window_multiplier, ...
    num_pc, n_rep_avg, max_smoothing_factor, delta, true_counts)

    n_pcs_frames = height(pcs_avg);
    pcs_avg_w_temp = pcs_avg;
    pcs_avg_windowed = zeros(n_pcs_frames, 1);
    
    % repetitive-average
    for n=1:1:n_rep_avg
        for i=1:1:n_pcs_frames
            e = min(n_pcs_frames, i + window_size);
            idx = min(n_pcs_frames, i+window_size/2);
            pcs_avg_windowed(idx) = mean(pcs_avg_w_temp(i : e));
        end
        pcs_avg_w_temp = pcs_avg_windowed;
    end
    
    % tau/threashold calculation
    tau = mean(pcs_avg_windowed);
    
    % calculate average-pca & look-ahead-ratio based segment predictions
    pred_labels_avgpca = zeros(n_pcs_frames, 1);
    pred_labels_lookahead = zeros(n_pcs_frames, 1);
    activity_val=tau+3; nonactivity_val=tau-3;
    last_label = nonactivity_val;
    for i=1:window_size:n_pcs_frames
        sw = max(1, i-window_size/2);    ew = min(n_pcs_frames, i+window_size/2);
        sw2 = max(1, i-window_multiplier*window_size/2); ew2 = min(nframes, i+window_multiplier*window_size/2); %#ok<NASGU> 
    
    
        % Activity-label (square-wave) calculation for mean of num_pc number of PCs
        m = mean(pcs_avg_windowed(sw:ew));
        if m >= tau; pred_labels_avgpca(sw:ew) = nonactivity_val;
        else; pred_labels_avgpca(sw:ew) = activity_val;
        end
        
        rm = int32((sw+ew)/2);
        % disp(ew2);
        a2 = abs(mean(pcs_avg_windowed(rm:ew2))-tau);
        a1 = abs(mean(pcs_avg_windowed(rm:ew))-tau);
        look_ahead_delta = a2 - a1;
        % look_back_ratio = mean(pcs_avg_windowed(sw2:rm)) / mean(pcs_avg_windowed(sw:rm));
        % s = sprintf('Look-ahead: %f ____a2=%f,  a1=%f', look_ahead_delta, a2, a1);
        % disp(s);

        if abs(look_ahead_delta) < delta % && look_ahead_delta> (1/delta)
            pred_labels_lookahead(rm:ew2) = last_label;
        elseif look_ahead_delta > 0
            pred_labels_lookahead(rm:ew2) = nonactivity_val;
            last_label = nonactivity_val;
        else
            pred_labels_lookahead(rm:ew2) = activity_val;
            last_label = activity_val;
        end
    end
    

    % calculate segment predictions by tau_series / local-average as the threshold
    pred_labels_localavg = zeros(n_pcs_frames, 1);
    tau_series = zeros(n_pcs_frames, 1);
    wm = window_multiplier*window_size;
    for i=1:wm:n_pcs_frames
        sw = max(1, i-wm/2);            ew = min(n_pcs_frames, i+wm/2);
        tau_series(sw:ew) = mean(pcs_avg_windowed(sw:ew));
    end
    for i=1:window_size:n_pcs_frames
        sw = max(1, i-window_size/2);   ew = min(n_pcs_frames, i+window_size/2);
        if mean(pcs_avg_windowed(sw:ew))>=tau_series(int32(floor((sw+ew)/2)))
            pred_labels_localavg(sw:ew)=nonactivity_val;
        else
            pred_labels_localavg(sw:ew)=activity_val;
        end
    end
    % Smooth-out `pred_labels_*3` => Merging single-windowed islands
    for s=1:1:max_smoothing_factor
        w = s * window_size;
        for i=3*w/2+1:1:n_pcs_frames-w/2+1
            scw = max(1, i-w/2);    ecw = min(n_pcs_frames, i+w/2);
            spw = max(1, i-3*w/2);  epw = max(1, scw-1);
            snw = min(n_pcs_frames,ecw+1);    enw = min(n_pcs_frames, i+3*w/2);
    
            cw_val = max(pred_labels_avgpca(scw:ecw)); % current-window's label
            pw_val = max(pred_labels_avgpca(spw:epw)); % previous-window's label
            nw_val = max(pred_labels_avgpca(snw:enw)); % next-window's label
            if cw_val~=pw_val && cw_val~=nw_val
                pred_labels_avgpca(scw:ecw) = nw_val; % smooth AvgPCA predictions
            end
    
            cw_val = max(pred_labels_localavg(scw:ecw)); % current-window's label
            pw_val = max(pred_labels_localavg(spw:epw)); % previous-window's label
            nw_val = max(pred_labels_localavg(snw:enw)); % next-window's label
            if cw_val~=pw_val && cw_val~=nw_val
                pred_labels_localavg(scw:ecw) = nw_val; % smooth LocalAvg predictions
            end
    
            cw_val = max(pred_labels_lookahead(scw:ecw)); % current-window's label
            pw_val = max(pred_labels_lookahead(spw:epw)); % previous-window's label
            nw_val = max(pred_labels_lookahead(snw:enw)); % next-window's label
            if cw_val~=pw_val && cw_val~=nw_val
                pred_labels_lookahead(scw:ecw) = nw_val; % smooth Look-Ahead predictions
            end
        end
    end

    % construct activity vs. non-activity true-label series by each CSI-frame
    labels = zeros(nframes, 1);
    for i=1:1:nframes
        l = true_action_labels(i);
        if l==0; l=5; elseif l==7; l=10; else; l = 15; end
        labels(i) = l;
    end
    pred_counts = zeros(height(true_counts), 1);

    ct = 1;
    i=1;
    while i<=nframes
        if labels(i)==15
            last = i;
            for j=i:1:nframes
                if labels(j)~=15
                    last = j;
                    break;
                end
            end
            curr_seg = pcs_avg_windowed(i:last);
            pred_counts(ct) = sum(islocalmin(curr_seg)) + sum(islocalmax(curr_seg));
            ct = ct +1;
            i=last+1;
        else 
            i = i + 1;
        end
    end
    count_accs = zeros(height(true_counts), 1);
    counts_for_bar = zeros(height(true_counts), 2);
    for i=1:height(true_counts)
        count_accs(i) = 100 - (100.0*abs(true_counts(i)-pred_counts(i)))/true_counts(i);
        counts_for_bar(i,1)=pred_counts(i);
        counts_for_bar(i,2)=true_counts(i);
    end
    counting_accuracy = mean(count_accs);

    %_{
    figure; 
    bar_counts=bar(counts_for_bar); set(gca, 'FontSize', 14);
    hold on;
    xlabel('Activity Segment Number'); ylabel('Rep. Counts');
    legend(bar_counts, 'Predicted Counts', 'True Counts',  'FontSize', 14, 'Orientation','horizontal');
    %}
    

    % Construct 3 methods' segmented-graphs for input to counting as
    % counts = sum(islocalmax+islocalmin)/2
    seg_graph_avgpca = zeros(n_pcs_frames, 1);
    seg_graph_localavg = zeros(n_pcs_frames, 1);
    seg_graph_lookahead = zeros(n_pcs_frames, 1);
    for i=1:1:n_pcs_frames
        if pred_labels_avgpca(i)<tau; seg_graph_avgpca(i) = tau;
        else; seg_graph_avgpca(i) = pcs_avg_windowed(i); end

        if pred_labels_localavg(i)<tau_series(i); seg_graph_localavg(i) = tau;
        else; seg_graph_localavg(i) = pcs_avg_windowed(i); end

        if pred_labels_lookahead(i)<tau; seg_graph_lookahead(i) = tau;
        else; seg_graph_lookahead(i) = pcs_avg_windowed(i); end
    end

    tmin_p = islocalmin(seg_graph_avgpca); tmax_p = islocalmax(seg_graph_avgpca);
    counts_avgpca = (sum(tmin_p)+sum(tmax_p))/2;

    tmin_l = islocalmin(seg_graph_localavg); tmax_l = islocalmax(seg_graph_localavg);
    counts_localavg = (sum(tmin_l) + sum(tmax_l))/2;

    tmin_la=islocalmin(seg_graph_lookahead); tmax_la=islocalmax(seg_graph_lookahead);
    counts_lookahead= (sum(tmin_la) + sum(tmax_la))/2;

    
    % Calculate the accuracy & F-1 score
    for fgv=170 % 0:10:201
        [accuracy, f1, pgp] = calcAccuracyAndF1Score( ...
            n_pcs_frames, tau, pred_labels_avgpca, true_action_labels, fgv); %#ok<ASGLU> 
        [accuracy_l, f1_l, pgl] = calcAccuracyAndF1Score( ...
            n_pcs_frames, tau, pred_labels_lookahead, true_action_labels, fgv); %#ok<ASGLU> 
        [accuracy_v, f1_v, pgv] = calcAccuracyAndF1Score( ...
            n_pcs_frames, tau_series, pred_labels_localavg, true_action_labels, fgv); %#ok<ASGLU> 
        s = sprintf(['w=%d,wm=%d,p=%d,n=%d,y=%d,d=%f, ' ...
            'al=%f,fl=%f, av=%f,fv=%f, ap=%f,fp=%f, ' ...
            'cl=%d, cv=%d, cp=%d, fgv=%d, c=%f'], ...
            window_size, window_multiplier, num_pc, n_rep_avg, ...
            max_smoothing_factor, delta, ...
            accuracy_l, f1_l, accuracy_v, f1_v, accuracy, f1, ...
            counts_lookahead, counts_localavg, counts_avgpca, fgv, counting_accuracy);
        disp(s); %#ok<DSPS> 
        % save("results.csv", s, '-append');
        writematrix([window_size, window_multiplier, num_pc, n_rep_avg, max_smoothing_factor, delta, ...
            accuracy_l, f1_l, accuracy_v, f1_v, accuracy, f1, ...
            counts_lookahead, counts_localavg, counts_avgpca, fgv, counting_accuracy], ...
            "seg_count_27NovTestVert.csv", 'Delimiter',',', 'WriteMode','append');
    end


    %_{
    figure;
    hold on;
    pw = plot(pcs_avg_windowed-tau, 'Color','blue', 'LineWidth', 1); set(gca, 'FontSize', 14);
    pt = plot(labels-mean(labels), 'color','black', 'LineWidth', 1);
    xlabel('CSI Frame Number'); ylabel('CSI Amplitude');
    % pv = plot(pred_labels_avgpca-mean(pred_labels_avgpca), 'Color','red', 'LineWidth', 1); % plot activity-labels (square-wave) for mean of num_pc number of PCs
    % pv = plot(pred_labels_lookahead-mean(pred_labels_lookahead), 'Color','red', 'LineWidth', 1); % plot activity-labels (square-wave) for mean of num_pc number of PCs
    pv = plot(pred_labels_localavg-mean(pred_labels_localavg), 'Color','red', 'LineWidth', 1);
    % plot(pgv, 'Color','green', 'LineWidth', 1);
    % legend(p2, 'Segments', 'Location',  'east');
    % csirange = csirange-csirange(1);
    % plot(csirange, seg_graph_localavg, csirange(tmin_p), seg_graph_localavg(tmin_p), 'r*');
    % plot(csirange, seg_graph_localavg, csirange(tmax_p), seg_graph_localavg(tmax_p), 'r*');
    legend([pw;pt;pv], 'Averaged PCA', 'True Segments', 'Predicted Segments', 'FontSize', 14,'Orientation','horizontal');%,'Location','south'
    %}
    % figure; plot(pred_counts);hold on; plot(true_counts);
end

function [accuracy, f1score, pred_graph] = calcAccuracyAndF1Score(nFrames, tau, predLabels, trueLabels, fgv)
    tp = 0;
    tn = 0;
    fp = 0;
    fn = 0;
    ff = zeros(nFrames, 1);
    ht = height(tau);
    pred_graph = zeros(nFrames, 1);
    total_forgiven=0;
    for c=1:1:nFrames
        true_label = trueLabels(c);
        if true_label==0 || true_label==7; true_label = 0; else; true_label = 1; end
    
        pred_label = predLabels(c);
        t=tau; if ht>1; t=tau(c);end
        if pred_label<t; pred_label = 0;else; pred_label = 1; end
        pred_graph(c) = pred_label;

        %{
        if c>2771 && c<2867
            a=1;
        end
        %}
    
        if true_label==pred_label
            if true_label==0; tn = tn + 1;ff(c)=60; else; tp = tp + 1;ff(c)=90; end
        else
            % Apply forgiveness-range first
            is_forgiven = 0;
            for cf=max(1, c-fgv) : min(nFrames, c+fgv+1) %+1-frame is always being forgiven
                if trueLabels(cf)==pred_label
                    is_forgiven = 1;
                    pred_graph(c:cf) = trueLabels(c:cf);
                    break;
                end
            end
            if is_forgiven==0
                if true_label==0; fn = fn + 1;ff(c)=70; else; fp = fp + 1;ff(c)=80; end
            else
                if true_label==0; tn = tn + 1;ff(c)=60; else; tp = tp + 1;ff(c)=90; end
            end
            total_forgiven = total_forgiven + is_forgiven;
        end
    end
    % total_forgiven
    f1score = tp / (tp + (fp+fn)/2);
    accuracy = (100 * (tp+tn))/nFrames;

    pred_graph = pred_graph * 5 + 37.2;
end



