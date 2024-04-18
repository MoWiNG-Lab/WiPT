predictions = readmatrix("misc_data/predictions.99.csv");  % is_training, y_true, y_pred

y_test_true = predictions(predictions(:,1)==0, 2);
y_test_pred = predictions(predictions(:,1)==0, 3);
total_test_points = size(y_test_true);

percs = zeros(15, 1);
seg_count = 0; 
prev_c = y_test_true(1);
ttrue = 0;
tpred = 0;
for i=2:1:total_test_points
    ctrue = y_test_true(i);
    ttrue = ttrue + 1;

    if ctrue ~= prev_c
        seg_count = seg_count + 1;
        percs(seg_count) = 100 * tpred / ttrue;
        tpred = 0;
        ttrue = 0;
    end

    cpred = y_test_pred(i);
    if ctrue == cpred
        tpred = tpred + 1;
    end
    prev_c = ctrue;
end
seg_count = seg_count + 1;
percs(seg_count) = 100 * tpred / ttrue;
disp(percs);



x_axis = 1:total_test_points;
figure; hold on; set(gca, 'FontSize', 14);
plot(x_axis, y_test_true-1, 'Marker', '.');
scatter(x_axis, y_test_pred-1, 'Marker', 'x');
