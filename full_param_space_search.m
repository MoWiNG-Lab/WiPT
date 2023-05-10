[counting_acc_file, avgpca_acc_file, local_avg_acc_file, look_ahead_acc_file] = ...
    wipt.initFiles(string(datetime('now', 'Format', 'dMMMyy_HHmmss')), 0);

% Run all possible combinations of the parameters & log results
for w = window_sizes
    for ps = pc_start
        for p = max_pcas
            for n = rep_avgs
                [Pw, pca_explained_perc] = wipt.getAveragePCASeries(H, ps, p, n, w);
                [Pw_vld, ~] = wipt.getAveragePCASeries(H_vld, ps, p, n, w);
    
                % Counting Results
                % [counting_accuracy, count_bars] = wipt.CountReps(Pw, true_labels, true_counts);
                % wipt.appendMatrix([w, p_start, p, n, counting_accuracy], counting_acc_file);
                
                % Segmentations
                for y=0:2:max_smoothing_factor_gamma
                    % Xp = wipt.SegmentByOverallAvg(Pw, w, y);
                    % wipt.logAccuracy(Xp, true_labels, [w ps p n y], avgpca_acc_file);
    
                    for psi = window_multipliers
                        Xv = wipt.SegmentByLocalAvg(Pw, w, psi, y);
                        Xv_vld = wipt.SegmentByLocalAvg(Pw_vld, w, psi, y);

                        bestF1 = wipt.logAccuracy(Xv, true_labels, ...
                            [w psi ps p n y], local_avg_acc_file);
                        vldF1 = wipt.logAccuracy(Xv_vld, tls, ...
                            [w psi ps p n y], local_avg_acc_file);
                        fprintf("F1=%.4f (%.4f) <-- w=%d, wm=%d, " + ...
                            "ps=%d, p=%d, n=%d, y=%d\n\n", bestF1, vldF1, ...
                            w, psi, ps, p, n, y);
                    end
                end
            end
        end
    end
end