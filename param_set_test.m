constants; wipt = WiPT;
DATA = DATA_5JanTest4xRandVert; % DATA_5JanTest4xRandVert;
fullcsv = readmatrix(csvfiles(DATA));
csirange = csiranges(2*DATA-1):csiranges(2*DATA); fullcsv = fullcsv(csirange, :);
H = fullcsv(:, 3:66); l = height(H);
global TAG_ACT TAG_NONACT; %#ok<GVMIS> 
true_action_labels = fullcsv(:,2); true_labels = zeros(l, 1);
for i=1:1:l
    l = true_action_labels(i);
    if l==0; l=5; elseif l==7; l=TAG_NONACT; else; l = TAG_ACT; end
    true_labels(i) = l;
end


csv = readmatrix('results/23Jan23_204811/seg_localavg.csv');
% 1. Window Size,
% 2. Window Multiplier,
% 3. #StartPC,
% 4. #PC,
% 5. #RepAvg,
% 6. Smoothing Factor,
% 7. Forgiveness Offset,
% 8. Accuracy(%),
% 9. F1-Score
bestSet = [0 0.0 0.0 0 0 0 0 0 0 0 0.0 0.0];
n = height(csv); i=1;
while i<=n
    row = csv(i, :);
    if row(8)>70
        [Pw, ~] = wipt.getAveragePCASeries(H, row(3), row(4), row(5), row(1));
        X = wipt.SegmentByLocalAvg(Pw, row(1), row(2), row(6));
        [acc, f1] = wipt.calcAccuracyF1(X, true_labels, row(7), 200, 200);
        
        fprintf('%d/%d ==> Accuracy: %f    F1-Score: %f\n', i, n, acc, f1);

        if acc>bestSet(11)
            bestSet = [i acc f1 row];
        end
        i = i + (150-row(7))/10 + 1;
    else
        i = i + 1;
    end
end




