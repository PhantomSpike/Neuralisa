%Here we will make the two types of filters, notch and inverse notch
clear all; close all;
%Load the data
filename = 'sample_data.wav';
[data.input,fs] = audioread(filename);
t_start_s = 5;
dur_s = 10;
t_end_s = t_start_s+dur_s;
data.input = data.input(t_start_s*fs:t_end_s*fs, 1);

%Here we design the notch filter
%>> Input params
f0_Hz = 12e3; %Center freq of the filter
% fs = 44100; %Sampling rate
bw_oct = 1/2; %The width of the filter in Octaves
f1_Hz = f0_Hz*2^-bw_oct;
f2_Hz = f0_Hz*2^bw_oct;

%2nd order notch filter
Q = f0_Hz/(f2_Hz-f1_Hz);
wo = f0_Hz/(fs/2);  
bw = wo/Q;
[b,a] = iirnotch(wo,bw);

%Make white-noise sound 
% dur_s = 5; %Duration in s

% data.input = randn(dur_s*fs,1);

%Filter the input
data.output = filtfilt(b,a,data.input);

fnames = fieldnames(data);


%Loop over the resulrs
figure('units','normalized','outerposition',[0 0 1 1]);
lwidth = 2;
dt_spec_ms = 20; %The bin of the cochleagram in ms
f_min = 500;
f_max = 20000;
n_f = 64;
n_tlab = 10;



for j = 1:length(fnames)
    [P1,f] = fft_run(data.(fnames{j}),fs);
    %Plot the results
    subplot(2,length(fnames),j);
    plot(f,P1,'LineWidth',lwidth);
    set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
    xlabel('Freqeuncy [kHz]','FontSize',16,'FontWeight','bold');
    ylabel('|P1(f)|','FontSize',16,'FontWeight','bold');
    title(['Amplitude Spectrum of ',fnames{j}],'FontSize',16,'FontWeight','bold');
    
    [X_ft, t, params] = cochleagram_spec_log(data.(fnames{j}), fs, dt_spec_ms, 'log', f_min, f_max, n_f);
    skip_t = round(length(t)/n_tlab);
    freqs = fliplr(params.freqs);
    skip_f = 3;
    freqs = freqs(1:skip_f:end);
    
    for f = 1:numel(freqs)
        y_labels{f} = num2str(freqs(f)./1000,'%.1f');
    end
    
    for tm = 1:n_tlab
        x_labels{tm} = num2str(t((tm-1)*skip_t +1),'%.0f');
    end
    
    subplot(2,length(fnames),j+2);
    imagesc(flipud(X_ft));
    yticks([1:skip_f:n_f]);
    yticklabels(y_labels);
    xticks([1:skip_t:length(t)]);
    xticklabels(x_labels);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Power Spectrogram of ',fnames{j}],'FontSize',16,'FontWeight','bold');
end

%% Play the sounds
sound(data.input,fs); pause(dur_s); sound(data.output,fs);