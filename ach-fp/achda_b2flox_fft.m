s = struct;
tmp = {};
for x = 1:length(modAChDA)
    tmp{x,1} = modAChDA(x).rec;
    tmp{x,2} = modAChDA(x).FP; end
tmp = cell2struct(tmp, {'rec','fp'}, 2);
s(1).s = tmp; s(1).lbl = 'wt';

tmp = {};
for x = 1:length(b2flox)
    tmp{x,1} = b2flox(x).rec;
    tmp{x,2} = b2flox(x).FP; end
tmp = cell2struct(tmp, {'rec','fp'}, 2);
s(2).s = tmp; s(2).lbl = 'b2flox';

%%
p1_cell = cell(2,length(s));
for z = 1:length(s)
for y = 1:2 % Iterate over ACh, DA
    p1_mat = [];
for x = 1:length(s(z).s)
    vec = s(z).s(x).fp{y};
    % vec = vec - nanmean(vec);
    
    Fs = 50;                % Sampling frequency
    T = 1/Fs;               % Sampling period
    L = length(vec);        % Length of signal
    vec(isnan(vec)) = [];
    fftACh = fft(vec);      % Discrete Fourier Transform of photometry signal
    P2 = abs(fftACh/L);     % Two-sided spectrum P2
    P1 = P2(1:L/2+1);       % Single-sided spectrum P1 based on P2 and even-valued signal length L
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;     % Frequency domain vector
    P1 = medfilt1(P1);      % Median filter initial FFT
    p1_mat(:,x) = P1;
end
p1_cell{y,z} = p1_mat;
end
end

%% adjust by ANIMAL
p1_cell_an = cell(2,1);
for z = 1:length(s)
tmp = {}; for x = 1:length(s(z).s); tmp{x} = strtok(s(z).s(x).rec,'-'); end
uni = unique(tmp);
for y = 1:2; p1_cell_an{y,z} = [];
    for x = 1:length(uni)
        ii = find(strcmp(tmp,uni{x}));
        p1_cell_an{y,z}(:,x) = nanmean(p1_cell{y,z}(:,ii),2); % average all recordings for each animal
    end
end
end

%% PLOT
ds = 100; sm = 500;
y = 2; p1_mat = p1_cell{y};
vec = movmean(nanmean(p1_mat,2),sm); 
figure; loglog(f, vec);
title('Average FFT');
xlim([0.1 10]); xlabel('Frequency (Hz)'); ylabel('Power');

%% ADJUST
f_range = [0.01 10]; % frequency range to fit to
ds = 100; 
sm = 500;

vec = movmean(nanmean(p1_cell{1},2),sm);        % ACh average FFT
vec_log = log10(vec); vec_log = vec_log(:); % log10 of freq vector
f_log = log10(f); f_log = f_log(:);         % log10 of FFT vector
x1 = log10(f_range(1)); y1 = vec_log(f_log == x1);  % find FFT at freq = 0.1Hz
x2 = log10(f_range(2)); y2 = vec_log(f_log == x2);  % find FFT at freq = 10Hz
range = [find(f_log == x1):find(f_log == x2)];      % range of samples
coeff = polyfit([x1, x2], [y1, y2], 1);     % a*x+b coefficients for line connecting FFT power at 0.1Hz to power at 10Hz
f_fit = f_log(range); f_fit = f_fit(:);     % freq vector from [0.1 10] Hz
vec_sub = vec_log(range); vec_sub = vec_sub(:); % extract section corresponding to freq [0.1 10] Hz
vec_fit = (coeff(1)*f_fit) + coeff(2);      % fit line using freq vector and coefficients

% vec_flat = vec_sub - vec_fit;               % SUBTRACT away FFT power from aCSF average FFT fit above
% figure; loglog(10.^f_fit, 10.^vec_flat)

%% FLAT by RECORIDNG
vec_flat = [];
for y = 1:length(p1_cell_an)
    vec = movmean(nanmean([p1_cell_an{y}],2),sm);      % average FFT for this infusion condition
    vec_log = log10(vec); vec_log = vec_log(:);     % log10 of FFT vector
    vec_sub = vec_log(range); vec_sub = vec_sub(:); % extract section corresponding to freq [0.1 10] Hz
    vec_flat(:,y) = vec_sub - vec_fit;              % SUBTRACT away FFT power from aCSF average FFT fit above
end

figure;
loglog(10.^f_fit, 10.^vec_flat(:,1), 10.^f_fit, 10.^vec_flat(:,2));
% legend({s.lbl});
xlim([0.1 8.5]); xlabel('Frequency (Hz)'); ylabel('Power');
axis('square');
title('Average FFT - to ACh');

%% FLAT by ANIMAL
vec_flat_an = cell(2,length(s));
for y = 1:2
    f_range = [0.01 10]; % frequency range to fit to
    ds = 100; 
    sm = 500;

    vec = movmean(nanmean(p1_cell_an{y,1},2),sm);        % ACh average FFT
    vec_log = log10(vec); vec_log = vec_log(:); % log10 of freq vector
    f_log = log10(f); f_log = f_log(:);         % log10 of FFT vector
    x1 = log10(f_range(1)); y1 = vec_log(f_log == x1);  % find FFT at freq = 0.1Hz
    x2 = log10(f_range(2)); y2 = vec_log(f_log == x2);  % find FFT at freq = 10Hz
    range = [find(f_log == x1):find(f_log == x2)];      % range of samples
    coeff = polyfit([x1, x2], [y1, y2], 1);     % a*x+b coefficients for line connecting FFT power at 0.1Hz to power at 10Hz
    f_fit = f_log(range); f_fit = f_fit(:);     % freq vector from [0.1 10] Hz
    vec_sub = vec_log(range); vec_sub = vec_sub(:); % extract section corresponding to freq [0.1 10] Hz
    vec_fit = (coeff(1)*f_fit) + coeff(2);      % fit line using freq vector and coefficients

    for z = 1:length(s)
        for an = 1:size(p1_cell_an{y,z},2)
            vec = movmean([p1_cell_an{y,z}(:,an)],sm);      % average FFT for this infusion condition
            vec_log = log10(vec); vec_log = vec_log(:);     % log10 of FFT vector
            vec_sub = vec_log(range); vec_sub = vec_sub(:); % extract section corresponding to freq [0.1 10] Hz
            vec_flat_an{y,z}(:,an) = vec_sub - vec_fit; 
        end
    end
end
%%
figure;
loglog(10.^f_fit, 10.^nanmean(vec_flat_an{1,1},2)); hold on
for z = 1:length(s); clr = {'g','b';'m','r'};
    for y = 1:2
    shadederrbar(10.^f_fit, nanmean(10.^vec_flat_an{y,z},2), SEM(10.^vec_flat_an{y,z},2), clr{y,z});
    end
end
xlim([0.1 10]); xlabel('Frequency (Hz)'); ylabel('Power');
axis('square');
title('Average FFT - to ACh');
ylim([0 20]);

