clc; clear; constants; wipt = WiPT;  % close all; 
% Select the dataset: Data_RandTest4x, Data_RandTrain20x
DATA = Data_RandTrain20x;


% Construct the CSI series from the selected experiment's CSV file
fullcsv = readmatrix(csvfiles(DATA));
% csirange = csiranges(2*DATA-1):csiranges(2*DATA); fullcsv = fullcsv(csirange, :);
fullcsv = {fullcsv(89187:213560, :), fullcsv(268684:389840, :)}; % 14-april 1st 10x & last 10x for training
fullcsv = cat(1, fullcsv{:});
H = fullcsv(:, 3:66);

% Prepare the true-labels
h = height(H);
true_action_labels = fullcsv(:,2); true_labels = zeros(h, 1);
for i=1:1:h
    l = true_action_labels(i);
    if l==0; l=5; elseif l==7; l=TAG_NONACT; else; l = TAG_ACT; end
    true_labels(i) = l;
end


% 2. butterworth low-pass filter
b=butter(3, 0.1); % assuming the order n=3 and cut off freq. Wn=100Hz
filter(b, 1, H);


% 3. PCA
[coeff, ~, ~, ~, explained] = pca(H); % coeff = S
P = H * coeff(:, 1:3); % space of k eigenvectors & k=3


% 4. dynamic adaptive sliding window (T1/T1_new & T2/T2_new algo. impl.)

alpha = 0.5; %% TODO opimize it within 0.1~0.9
%(1) compute series of T1s as per T1=alpha * T1cur + (1-alpha) * T1new,
%length(T1) = 2 minutes = 120 seconds = 12000 CSI frames {T1cur=init? , T1new=var(W1)?}
%    as well as, T2s, where length(T2)=length(T1)/6
T1 = calculateThreshold(P, alpha, 1, 12000);
T2 = calculateThreshold(P, alpha, 1, 12000/6);

%(2) implement Algorithm-1
L = 300;
l = L/6;
startPoints = zeros(h, 1); sp = 1; 
endPoints = zeros(h, 1); ep = 1;
for i=L/2+1:L:h-L/2+1
    if(var(P(i-L/2:i))<mean(T1(i-L/2:i)) && ...
            var(P(i:i+L/2))>mean(T1(i:i+L/2)))
        if(var(P(i-l/2:i))<mean(T2(i-l/2:i)) && ...
                var(P(i:i+l/2))>mean(T2(i:i+l/2)))
            startPoints(sp) = i;
            sp = sp + 1;
        end
    end
    if(var(P(i-L/2:i))>mean(T1(i-L/2:i)) && ...
            var(P(i:i+L/2))<mean(T1(i:i+L/2)))
        if(var(P(i-l/2:i))>mean(T2(i-l/2:i)) && ...
                var(P(i:i+l/2))<mean(T2(i:i+l/2)))
            endPoints(ep) = i;
            ep = ep + 1;
        end
    end
end
true_segments = getSegments(true_labels);



function T = calculateThreshold(P, alpha, Tcur, t_length)
    T = zeros(height(P), 1);
    for t=1:t_length:height(T)
        t_end = t+t_length-1;
        % T1=alpha * T1cur + (1-alpha) * T1new,
        Tnew = mean(P(t:t_end)); % var(P(t:t_end));
        T(t:t_end) = alpha * Tcur + (1-alpha) * Tnew;
        Tcur = T(t:t_end);
    end
end


%%% Returns 2-D array consisting the starting & ending indices of the
%%% segments as specified by per frame in the input `tags` 1-D array
function segments = getSegments(tags)
global TAG_NONACT; %#ok<GVMIS> 
    l = height(tags);
    segments = zeros(1,2);

    ct = 1; i = 1;
    while i<=l
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
    segments(ct-1, 2) = l;
end

