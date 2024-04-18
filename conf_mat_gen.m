
Z3 = [94.24 3.05 2.71; 2.79 95.52 1.69; 3.36 5.12 91.52];
labels3 = ["Chest", "Hip", "Leg"];

Z6 = [92.68 3.88 0.26 0.00 1.64 1.55; ...
    3.12 90.71 3.46 0.00 2.28 0.42; ...
    0.34 3.88 93.45 1.81 0.26 0.26;
    0.51 0.00 4.26 92.50 2.47 0.26;
    0.33 0.08 0.33 4.11 94.16 0.99;
    6.14 0.17 0.26 0.09 7.59 85.76];
labels6 = ["Left-Chest", "Left-Hip", "Left-Leg", ...
    "Right-Chest", "Right-Hip", "Right-Leg"];

show_heatmap(Z3, labels3);
show_heatmap(Z6, labels6);


function show_heatmap(Z, labels)
    h = heatmap(Z, 'Colormap', (bone), 'CellLabelFormat', ...
        "%0.2f%%", 'XDisplayLabels', labels, 'YDisplayLabels', labels, ...
        'FontSize', 16);
    s = struct(h);
    s.YAxis.TickLabelRotation = 90;
    xlabel('Predicted Class')
    ylabel('True Class')
    caxis([0, 100])
end