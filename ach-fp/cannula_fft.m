z = 2; % ACh or DA

p1_cell = cell(3,1);
for y = 1:3
for x = 1:4 
    % vec = s(y).s(x).FP{z};
    vec = raw(y).raw(x).fp{z};
    
    Fs = 50;                % Sampling frequency
    T = 1/Fs;               % Sampling period
    L = length(vec);          % Length of signal
    vec(isnan(vec)) = [];
    fftACh = fft(vec);        % Discrete Fourier Transform of photometry signal
    P2 = abs(fftACh/L);     % Two-sided spectrum P2
    P1 = P2(1:L/2+1);       % Single-sided spectrum P1 based on P2 and even-valued signal length L
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;     % Frequency domain vector
    P1 = medfilt1(P1);      % Median filter initial FFT

    p1_cell{y} = [p1_cell{y}, P1];
end
end

%%
ds = 100; sm = 100;

figure;
vec = movmean(nanmean(p1_cell{1},2),sm);
loglog(f(1:ds:length(vec)-1), vec(1:ds:length(vec)-1),'--k');

clr = [0.05 0.75 0.45; 0.5 0.2 0.55; 0.85 0.35 0.1]; % teal, purple, orange
for y = 1:3
    vec = movmean(nanmean([p1_cell{y}],2),sm); z = movmean(SEM([p1_cell{y}],2),sm);
    plot(f(1:ds:length(vec)-1), vec(1:ds:length(vec)-1), 'Color', clr(y,:)); hold on
    % shadederrbar(f(1:ds:length(y)-1), y(1:ds:length(y)-1), z(1:ds:length(y)-1), clr(ii,:)); hold on
end
xlabel('Frequency'); ylabel('Power');

title('Average FFT');
xlim([0.1 10]); xlabel('Frequency (Hz)'); ylabel('Power');
legend({raw.inf});

%% PLOT
f_range = [0.01 24]; % frequency range to fit to
ds = 100; 
sm = 500;

vec = movmean(nanmean(p1_cell{1},2),sm);    % aCSF average FFT
vec_log = log10(vec); vec_log = vec_log(:); % log10 of freq vector
f_log = log10(f); f_log = f_log(:);         % log10 of FFT vector
x1 = log10(f_range(1)); y1 = vec_log(f_log == x1);  % find FFT at freq = 0.1Hz
x2 = log10(f_range(2)); y2 = vec_log(f_log == x2);  % find FFT at freq = 10Hz
range = [find(f_log == x1):find(f_log == x2)];      % range of samples
coeff = polyfit([x1, x2], [y1, y2], 1);     % a*x+b coefficients for line connecting FFT power at 0.1Hz to power at 10Hz
f_fit = f_log(range); f_fit = f_fit(:);     % freq vector from [0.1 10] Hz
vec_fit = (coeff(1)*f_fit) + coeff(2); % fit line using freq vector and coefficients
% vec_subtract = vec_log(range)' - vec_fit;

vec_flat = cell(3,1);
for y = 1:3
    for x = 1:4
        vec = movmean(p1_cell{y}(:,x),sm);      % average FFT for this infusion condition
        vec_log = log10(vec); vec_log = vec_log(:);     % log10 of FFT vector
        vec_sub = vec_log(range); vec_sub = vec_sub(:); % extract section corresponding to freq [0.1 10] Hz
        vec_flat{y}(:,x) = vec_sub - vec_fit;              % SUBTRACT away FFT power from aCSF average FFT fit above
    end
end

%%
figure;
loglog(10.^f_fit, 10.^vec_flat(:,1), 10.^f_fit, 10.^vec_flat(:,2), 10.^f_fit, 10.^vec_flat(:,3)); 
legend({raw.inf});
xlim([0.1 15]); xlabel('Frequency (Hz)'); ylabel('Power');
title('Average FFT - to aCSF');

%% SHADED BAR + DOWNSAMPLE
fig = figure; fig.Position([3 4]) = [560 365];
loglog(10.^f_fit(1:ds:end), 10.^nanmean(vec_flat{1}(1:ds:end,:),2)); hold on
for z = 1:length(vec_flat); clr = {'g','m','r'};
    shadederrbar(10.^f_fit(1:ds:end), nanmean(10.^vec_flat{z}(1:ds:end,:),2), SEM(10.^vec_flat{z}(1:ds:end,:),2), clr{z});
end
xlim([0.1 15]); xlabel('Frequency (Hz)'); 
ylim([0 5]); yticks([0:5]); ylabel('Power');
axis('square');
title('Average FFT - to ACh');

movegui(gcf,'center');

%% PLOT: ACh, DA
fig = figure; fig.Position(3) = 1375;
sp(1) = subplot(1,3,1); 
vec_flat = vec_flat_ach; % ACH FFT
loglog(10.^f_fit, 10.^vec_flat(:,1), 10.^f_fit, 10.^vec_flat(:,2), 10.^f_fit, 10.^vec_flat(:,3)); 
legend({raw.inf});
xlim([0.1 15]); xlabel('Frequency (Hz)'); ylabel('Power (a.u.)'); ylim([0 5])
axis('square');
title(sprintf('ACh - FFT to aCSF (r [%1.1f %d], ds %d, sm %d)',f_range(1),f_range(2),ds,sm));

sp(1) = subplot(1,3,2);
vec_flat = vec_flat_da; % ACH FFT
loglog(10.^f_fit, 10.^vec_flat(:,1), 10.^f_fit, 10.^vec_flat(:,2), 10.^f_fit, 10.^vec_flat(:,3)); 
legend({raw.inf});
xlim([0.1 15]); xlabel('Frequency (Hz)'); ylabel('Power (a.u.)'); ylim([0 15])
 axis('square');
title(sprintf('DA - FFT to aCSF (r [%1.1f %d], ds %d, sm %d)',f_range(1),f_range(2),ds,sm));
movegui(gcf,'center');

sp(3) = subplot(1,3,3);
loglog(10.^f_fit, 10.^vec_flat_ach(:,1), 10.^f_fit, 10.^vec_flat_da(:,1)); 
legend({raw.inf});
xlim([0.1 15]); xlabel('Frequency (Hz)'); ylabel('Power (a.u.)'); ylim([0 6])
axis('square');
title(sprintf('ACh + DA - FFT aCSF'));
movegui(gcf,'center');

linkaxes(sp,'y'); linkaxes(sp,'x');

%% SHADED BAR + DOWNSAMPLE -- ACH + DA
fig = figure; fig.Position([3 4]) = [560 365];
y = 1; % aCSF
loglog(10.^f_fit(1:ds:end), 10.^nanmean(vec_flat_ach{y}(1:ds:end,:),2)); hold on
shadederrbar(10.^f_fit(1:ds:end), nanmean(10.^vec_flat_ach{y}(1:ds:end,:),2), SEM(10.^vec_flat_ach{y}(1:ds:end,:),2), [0.05 0.75 0.45]);
ylim([0 10]); yticks([0:5]); ylabel('Power');
yyaxis right
shadederrbar(10.^f_fit(1:ds:end), nanmean(10.^vec_flat_da{y}(1:ds:end,:),2), SEM(10.^vec_flat_da{y}(1:ds:end,:),2), [1 0 1]);
xlim([0.1 15]); xlabel('Frequency (Hz)'); 
axis('square');
title('Average FFT - to ACh');
movegui(gcf,'center');