z = 1; % ACh or DA

p1_cell = cell(3,1);
for y = 1:3
for x = 1:4 
    vec = s(y).s(x).FP{z};
    
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
z = 1; % ACh or DA
p1_cell = cell(3,1);
for y = 1:3 
for x = 1:4
    fp = s(y).s(x).FP{z}; Fs = 50; % Extract photometry signal
    imm = s(y).s(x).onRest; imm(:,2) = s(y).s(x).offRest; % Extract rest samples
    winSize = 2; winOv = 1;
    [seg_fp, seg_samp] = segmentSignal(fp,winSize,winOv,Fs,imm);

    P1_out = [];
    for r = 1:size(seg_fp,2)
        vec = seg_fp(:,r);  % Photometry signal segment
        if any(isnan(vec)); continue; end % Skip any segments with nans
        Fs = 50;                % Sampling frequency
        T = 1/Fs;               % Sampling period
        L = length(vec);    % Length of signal
        vec(isnan(vec)) = [];
        fftACh = fft(vec);  % Discrete Fourier Transform of photometry signal
        P2 = abs(fftACh/L);     % Two-sided spectrum P2
        P1 = P2(1:L/2+1);       % Single-sided spectrum P1 based on P2 and even-valued signal length L
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;     % Frequency domain vector
        P1 = medfilt1(P1);      % Median filter initial FFT
        P1_out(:,r) = P1;
    end
    
    p1_cell{y} = [p1_cell{y}, nanmean(P1_out,2)];
end
end