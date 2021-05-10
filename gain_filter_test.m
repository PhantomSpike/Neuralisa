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
f0_Hz = 12420; %Center freq of the filter
% fs = 44100; %Sampling rate
bw_oct = 1/10; %The width of the filter in Octaves
f1_Hz = f0_Hz*2^-bw_oct;
f2_Hz = f0_Hz*2^bw_oct;

%Bandpass filter
order = 4;
nyquist = fs/2;
[b,a] = butter(order,[f1_Hz f2_Hz]/nyquist);


%Make white-noise sound 
x_white = randn(length(data.input),1);

%Filter the input
x_filt = filtfilt(b,a,x_white);
db_input = get_db(data.input);
db_noise = db_input - 5;

adj_coeff = db_adjust(x_filt,db_noise);
x_filt = x_filt.*adj_coeff;
data.output = data.input + x_filt;

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