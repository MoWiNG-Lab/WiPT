filename = "results/26Dec22_223356_localavg_best/seg_localavg_best_26decTest.csv";

mat = readtable(filename); %mat=mat(:,1:9);


window_plot_table = mat(:,[1 2 8]);
pc_plot_table = mat(:,[3 4 8]);
rest_plot_table = mat(:,[5 6 7 8]);

figure; p = parallelplot(window_plot_table, 'GroupVariable','WindowSize');    set(gca, 'FontSize', 14); legend(p, 'Location', 'north');
figure; parallelplot(pc_plot_table, 'GroupVariable','x_StartPC');        set(gca, 'FontSize', 14);
figure; parallelplot(rest_plot_table, 'GroupVariable','x_RepAvg');      set(gca, 'FontSize', 14);