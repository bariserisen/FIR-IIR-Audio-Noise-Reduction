%% =====================================================================
%  EE0857 DSP PROJECT — Removal of Noise from a Music Signal Using
%                        IIR and FIR Digital Filters
%  ---------------------------------------------------------------------
%  This script requires NO toolbox and NO external file.
%  Only core MATLAB functions are used:
%     filter, conv, fft, roots, poly, residue, unwrap, sound, plot, stem
%  (Signal Processing Toolbox functions such as butter, fir1, freqz, zplane,
%   hamming, impz and filtfilt are re-implemented at the end of this file.)
%
%  Section map (course schedule):
%     Section 1  -> Week 2  : Discrete-time signals
%     Section 2  -> Week 3  : Discrete-time systems (linearity, convolution)
%     Section 3  -> Week 4  : Discrete Fourier Transform (own DFT vs fft)
%     Section 4  -> Week 5  : Sampling and aliasing (Nyquist)
%     Section 5  -> Week 6  : z-transform, poles/zeros, stability
%     Section 6  -> Week 7  : Inverse z-transform (partial fractions)
%     Section 7  -> Week 10 : System and design tools (own freqz/zplane)
%     Section 8  -> Week 11 : IIR filter design (Butterworth + bilinear, notch)
%     Section 9  -> Week 12 : Linear-phase FIR (symmetry, constant group delay)
%     Section 10 -> Week 13 : FIR design by the window method
%     Section 11 -> Week 9/14: Application — noise removal, SNR, listening
%     Section 12 -> Report outputs: specifications, coefficients,
%                   impulse responses, saving figures as PNG
%
%  Run: open this file in MATLAB and press F5 (R2016b or newer).
% ======================================================================

clear; close all; clc;
rng(0);                 % reproducible noise
PLAY_SOUND = true;      % set to false if you do not want to listen

%% ===================== SECTION 1: DISCRETE-TIME SIGNALS ================
fprintf('\n=== SECTION 1: Discrete-time signals ===\n');

n = -10:20;
delta = double(n == 0);                 % unit impulse  delta[n]
u     = double(n >= 0);                 % unit step     u[n]
expo  = (0.8.^n) .* u;                  % exponential   (0.8)^n u[n]
sinus = cos(2*pi*0.1*n);                % discrete cosine, w0 = 0.2*pi

figure('Name','Section 1 - Basic sequences');
subplot(2,2,1); stem(n,delta,'filled'); title('\delta[n]'); xlabel('n'); grid on;
subplot(2,2,2); stem(n,u,'filled');     title('u[n]'); xlabel('n'); grid on;
subplot(2,2,3); stem(n,expo,'filled');  title('(0.8)^n u[n]'); xlabel('n'); grid on;
subplot(2,2,4); stem(n,sinus,'filled'); title('cos(0.2\pi n)'); xlabel('n'); grid on;

% --- Project signal: a melody generated in code instead of an audio file ---
Fs    = 8000;                                % sampling frequency [Hz]
notes = [262 294 330 349 392 440 494 523];   % C4 D4 E4 F4 G4 A4 B4 C5
dur   = 0.35;                                % duration of each note [s]
clean = [];
for f0 = notes
    t    = (0:round(dur*Fs)-1)/Fs;
    env  = exp(-3*t) .* min(1, t/0.01);      % piano-like envelope (attack + decay)
    note = env .* ( sin(2*pi*f0*t) + 0.5*sin(2*pi*2*f0*t) + 0.25*sin(2*pi*3*f0*t) );
    clean = [clean note];                    %#ok<AGROW>
end
clean = clean / max(abs(clean));
N  = length(clean);
t  = (0:N-1)/Fs;

% --- Noise: 50 Hz mains hum + 3 kHz whistle + white noise ---
hum     = 0.30*sin(2*pi*50*t);
whistle = 0.25*sin(2*pi*3000*t);
white   = 0.10*randn(1,N);
noisy   = clean + hum + whistle + white;

figure('Name','Section 1 - Project signal');
subplot(2,1,1); plot(t,clean); title('Clean melody x[n]'); xlabel('t [s]'); grid on;
subplot(2,1,2); plot(t,noisy); title('Noisy melody (50 Hz + 3 kHz + white noise)');
xlabel('t [s]'); grid on;

%% ===================== SECTION 2: DISCRETE-TIME SYSTEMS ================
fprintf('\n=== SECTION 2: Discrete-time systems ===\n');

% Example system: 5-point moving average  y[n] = (1/5) sum_{k=0}^{4} x[n-k]
h_ma = ones(1,5)/5;

% Convolution: own implementation vs MATLAB conv
x_test = [1 2 3 4 5 4 3 2 1];
y_own  = my_conv(x_test, h_ma);
y_mat  = conv(x_test, h_ma);
fprintf('Convolution error (my_conv - conv): %.2e\n', max(abs(y_own - y_mat)));

% Linearity test:  T{a x1 + b x2} == a T{x1} + b T{x2} ?
x1 = randn(1,50); x2 = randn(1,50); a1 = 2; b1 = -3;
lhs = filter(h_ma,1,a1*x1 + b1*x2);
rhs = a1*filter(h_ma,1,x1) + b1*filter(h_ma,1,x2);
fprintf('Linearity test error: %.2e  -> system is linear\n', max(abs(lhs-rhs)));

% Time-invariance test: delay the input by 7 samples, output must be delayed too
k0 = 7;
y1 = filter(h_ma,1,x1);
y2 = filter(h_ma,1,[zeros(1,k0) x1(1:end-k0)]);
fprintf('Time-invariance error: %.2e -> system is LTI\n', ...
        max(abs(y2(k0+1:end) - y1(1:end-k0))));

figure('Name','Section 2 - Convolution');
subplot(3,1,1); stem(0:length(x_test)-1, x_test,'filled'); title('x[n]'); grid on;
subplot(3,1,2); stem(0:4, h_ma,'filled'); title('h[n] (5-point moving average)'); grid on;
subplot(3,1,3); stem(0:length(y_own)-1, y_own,'filled'); title('y[n] = x[n] * h[n]');
xlabel('n'); grid on;

%% ===================== SECTION 3: DISCRETE FOURIER TRANSFORM ===========
fprintf('\n=== SECTION 3: DFT ===\n');

Nd = 64;
nd = 0:Nd-1;
xd = cos(2*pi*8*nd/Nd) + 0.5*cos(2*pi*20*nd/Nd);   % peaks expected at k=8 and k=20
X_own = my_dft(xd);
X_fft = fft(xd);
fprintf('DFT error (my_dft - fft): %.2e\n', max(abs(X_own - X_fft)));

figure('Name','Section 3 - DFT');
subplot(2,1,1); stem(nd, xd,'filled'); title('x[n] (N = 64)'); xlabel('n'); grid on;
subplot(2,1,2); stem(nd, abs(X_own),'filled'); hold on;
plot(nd, abs(X_fft),'r.','MarkerSize',12);
title('|X[k]|  (blue: my\_dft, red: fft)'); xlabel('k'); grid on;

figure('Name','Section 3 - Spectrum of the project signal');
amplitude_spectrum(clean, Fs, 'b'); hold on;
amplitude_spectrum(noisy, Fs, 'r');
legend('Clean','Noisy'); title('Amplitude spectrum: 50 Hz and 3 kHz peaks are noise');

%% ===================== SECTION 4: SAMPLING =============================
fprintf('\n=== SECTION 4: Sampling and aliasing ===\n');

f_analog = 1000;                         % 1 kHz "analog" sinusoid
tc = 0:1e-6:5e-3;                        % continuous-time approximation
xc = cos(2*pi*f_analog*tc);

Fs_good = 8000;                          % Nyquist satisfied (8000 > 2*1000)
Fs_bad  = 1200;                          % Nyquist violated  (1200 < 2000)
f_alias = abs(f_analog - round(f_analog/Fs_bad)*Fs_bad);
fprintf('With Fs = %d Hz a 1000 Hz sinusoid appears at %d Hz (aliasing)\n', Fs_bad, f_alias);

figure('Name','Section 4 - Sampling');
for i = 1:2
    if i == 1, Fsx = Fs_good; else, Fsx = Fs_bad; end
    ts = 0:1/Fsx:5e-3;
    subplot(2,1,i); plot(tc*1e3, xc,'Color',[.7 .7 .7]); hold on;
    stem(ts*1e3, cos(2*pi*f_analog*ts),'filled');
    if i == 2
        plot(tc*1e3, cos(2*pi*f_alias*tc),'r--');
        legend('1 kHz analog','Samples',sprintf('%d Hz alias',f_alias));
    end
    title(sprintf('Fs = %d Hz', Fsx)); xlabel('t [ms]'); grid on;
end

%% ===================== SECTION 5: Z-TRANSFORM ==========================
fprintf('\n=== SECTION 5: z-transform ===\n');

% Example system: H(z) = 1 / ((1 - 0.9 z^-1)(1 - 0.5 z^-1))
b_z = 1;
a_z = conv([1 -0.9],[1 -0.5]);          % = [1 -1.4 0.45]
poles_z = roots(a_z);
fprintf('Poles: %s  -> all |z|<1, so the system is STABLE\n', mat2str(poles_z.',3));

% Unstable example: pole at z = 1.1
a_unstable = [1 -1.1];
h_unstable = filter(1, a_unstable, [1 zeros(1,39)]);

figure('Name','Section 5 - z-transform');
subplot(2,2,1); my_zplane(b_z, a_z); title('Stable: poles at 0.9 and 0.5');
subplot(2,2,2); stem(0:39, filter(b_z,a_z,[1 zeros(1,39)]),'filled');
title('h[n] decays'); xlabel('n'); grid on;
subplot(2,2,3); my_zplane(1, a_unstable); title('Unstable: pole at 1.1');
subplot(2,2,4); stem(0:39, h_unstable,'filled'); title('h[n] grows without bound');
xlabel('n'); grid on;

%% ===================== SECTION 6: INVERSE Z-TRANSFORM ==================
fprintf('\n=== SECTION 6: Inverse z-transform (partial fractions) ===\n');

% With H(z) given in powers of z^-1, multiply numerator and denominator by z^M
% and expand H(z)/z in partial fractions (residue is a core MATLAB function):
%   H(z)/z = sum r_k / (z - p_k)  =>  h[n] = sum r_k (p_k)^n u[n]
[r, p] = inv_z_partial(b_z, a_z);
for k = 1:length(r)
    fprintf('  r = %+.4f , p = %+.4f\n', real(r(k)), real(p(k)));
end
nn = 0:29;
h_closed = inv_z_eval(r, p, nn);
h_recur  = filter(b_z, a_z, [1 zeros(1,29)]);      % recursion of the difference equation
fprintf('Closed-form h[n] vs filter() difference: %.2e\n', max(abs(h_closed - h_recur)));
fprintf('Closed form: h[n] = (%.3f)(0.9)^n + (%.3f)(0.5)^n,  n>=0\n', ...
        real(r(abs(p-0.9)<1e-6)), real(r(abs(p-0.5)<1e-6)));

figure('Name','Section 6 - Inverse z-transform');
stem(nn, h_recur,'filled'); hold on; plot(nn, h_closed,'ro');
legend('Recursion (filter)','Partial fractions'); title('Inverse z-transform check');
xlabel('n'); grid on;

%% ===================== SECTION 7: SYSTEM AND DESIGN TOOLS ==============
% my_freqz and my_zplane (end of file) replace the toolbox functions
% freqz and zplane. Here the moving-average filter is analysed.
fprintf('\n=== SECTION 7: Frequency-response tool ===\n');
[Hm, w] = my_freqz(h_ma, 1, 1024);
figure('Name','Section 7 - Moving-average frequency response');
subplot(2,1,1); plot(w/pi, 20*log10(abs(Hm)+eps)); ylim([-60 5]);
title('|H(e^{j\omega})| [dB]'); xlabel('\omega/\pi'); grid on;
subplot(2,1,2); plot(w/pi, unwrap(angle(Hm))); title('Phase [rad]'); xlabel('\omega/\pi'); grid on;

%% ===================== SECTION 8: IIR FILTER DESIGN ====================
fprintf('\n=== SECTION 8: IIR design ===\n');

% (a) Butterworth low-pass — analog prototype + bilinear transformation
fc_lp = 2000;   N_iir = 6;
[b_lp, a_lp] = my_butter(N_iir, fc_lp, Fs);
fprintf('Butterworth (N=%d, fc=%d Hz) pole magnitudes: %s\n', N_iir, fc_lp, ...
        mat2str(abs(roots(a_lp)).',3));

% (b) 50 Hz notch filter — zeros on the unit circle, poles just inside
[b_n, a_n] = notch_filter(50, Fs, 0.995);

% Cascade: H(z) = H_lp(z) * H_notch(z)  (polynomial product = conv)
b_iir = conv(b_lp, b_n);
a_iir = conv(a_lp, a_n);

[H_lp, w] = my_freqz(b_lp, a_lp, 4096);
[H_n,  ~] = my_freqz(b_n,  a_n,  4096);
f_axis    = w/pi*Fs/2;

figure('Name','Section 8 - IIR filters');
subplot(2,2,1); plot(f_axis, 20*log10(abs(H_lp)+eps)); ylim([-100 5]);
my_xline(fc_lp); title(sprintf('Butterworth LPF, N = %d',N_iir));
xlabel('f [Hz]'); ylabel('dB'); grid on;
subplot(2,2,2); my_zplane(b_lp, a_lp); title('Butterworth pole-zero plot');
subplot(2,2,3); plot(f_axis, 20*log10(abs(H_n)+eps)); xlim([0 300]); ylim([-60 5]);
title('50 Hz notch filter'); xlabel('f [Hz]'); ylabel('dB'); grid on;
subplot(2,2,4); my_zplane(b_n, a_n); title('Notch pole-zero plot (zoom in!)');

% Check the Butterworth impulse response with the inverse z-transform (Section 6)
[r_lp, p_lp] = inv_z_partial(b_lp, a_lp);
h_lp_closed  = inv_z_eval(r_lp, p_lp, 0:99);
h_lp_recur   = filter(b_lp, a_lp, [1 zeros(1,99)]);
fprintf('Butterworth h[n]: partial fractions vs filter() difference = %.2e\n', ...
        max(abs(h_lp_closed - h_lp_recur)));

%% ===================== SECTION 9: LINEAR-PHASE FIR =====================
fprintf('\n=== SECTION 9: Linear-phase FIR ===\n');

M = 101;                                  % filter length (odd -> Type I)
h_fir = windowed_fir(M, 2*pi*fc_lp/Fs, 'hamming');
sym_err = max(abs(h_fir - fliplr(h_fir)));
fprintf('Symmetry error h[n] = h[M-1-n]: %.2e -> linear phase\n', sym_err);
fprintf('Expected group delay: (M-1)/2 = %d samples = %.2f ms\n', (M-1)/2, (M-1)/2/Fs*1e3);

[H_fir, w] = my_freqz(h_fir, 1, 4096);
gd_fir = group_delay(H_fir, w);
gd_iir = group_delay(H_lp,  w);

figure('Name','Section 9 - Linear phase');
subplot(3,1,1); stem(0:M-1, h_fir, '.'); title('FIR h[n] — symmetric about n = 50');
xlabel('n'); grid on;
subplot(3,1,2); plot(f_axis, unwrap(angle(H_fir)), f_axis, unwrap(angle(H_lp)));
xlim([0 fc_lp]); legend('FIR (linear)','IIR (nonlinear)');
title('Phase in the passband'); xlabel('f [Hz]'); ylabel('rad'); grid on;
subplot(3,1,3); plot(f_axis(2:end), gd_fir, f_axis(2:end), gd_iir);
xlim([0 fc_lp]); ylim([0 80]); legend('FIR','IIR');
title('Group delay [samples]: FIR constant, IIR frequency-dependent'); xlabel('f [Hz]'); grid on;

%% ===================== SECTION 10: WINDOW METHOD =======================
fprintf('\n=== SECTION 10: FIR design by the window method ===\n');

windows = {'rectangular','hann','hamming','blackman'};
colors  = {'k','b','r','m'};
figure('Name','Section 10 - Window comparison');
for i = 1:numel(windows)
    h_i = windowed_fir(M, 2*pi*fc_lp/Fs, windows{i});
    [H_i, w] = my_freqz(h_i, 1, 4096);
    HdB = 20*log10(abs(H_i)+eps);
    idx = find(f_axis > fc_lp);                      % find the first null after the cutoff,
    k1  = find(diff(HdB(idx)) > 0, 1);               % then the highest sidelobe after it
    stop_peak = max(HdB(idx(k1):end));
    fprintf('  %-11s : highest stopband sidelobe = %6.1f dB\n', windows{i}, stop_peak);
    subplot(2,1,1); plot(f_axis, HdB, colors{i}); hold on;
    subplot(2,1,2); plot(0:M-1, make_window(M, windows{i}), colors{i}); hold on;
end
subplot(2,1,1); ylim([-120 5]); my_xline(fc_lp); legend(windows);
title(sprintf('Windowed FIR magnitude responses (M = %d)',M)); xlabel('f [Hz]'); ylabel('dB'); grid on;
subplot(2,1,2); legend(windows); title('Window functions w[n]'); xlabel('n'); grid on;

% Effect of the filter length
figure('Name','Section 10 - Effect of filter length');
for Mi = [21 51 101 201]
    [H_i, w] = my_freqz(windowed_fir(Mi, 2*pi*fc_lp/Fs,'hamming'), 1, 4096);
    plot(f_axis, 20*log10(abs(H_i)+eps)); hold on;
end
ylim([-100 5]); legend('M=21','M=51','M=101','M=201'); grid on;
title('Hamming: the transition band narrows as M increases'); xlabel('f [Hz]'); ylabel('dB');

%% ===================== SECTION 11: APPLICATION — NOISE REMOVAL =========
fprintf('\n=== SECTION 11: Noise removal ===\n');

% FIR path: a high-pass FIR for 50 Hz (spectral inversion of a low-pass),
% followed by the 2 kHz low-pass FIR. The convolution of two FIRs is one FIR.
M_hp = 801;                                          % 50 Hz is very close -> long filter
h_lp_short = windowed_fir(M_hp, 2*pi*120/Fs, 'blackman');
h_hp = -h_lp_short;  h_hp((M_hp+1)/2) = h_hp((M_hp+1)/2) + 1;   % delta[n-c] - h_lp[n]
h_fir_total = conv(h_hp, h_fir);
delay = (length(h_fir_total)-1)/2;                   % linear phase -> constant delay

y_fir = filter(h_fir_total, 1, [noisy zeros(1,delay)]);
y_fir = y_fir(delay+1:end);                          % compensate the delay

% IIR path: notch + Butterworth
y_iir      = filter(b_iir, a_iir, noisy);                    % causal (phase distortion)
y_iir_zero = zero_phase_filter(b_iir, a_iir, noisy);         % forward-backward (zero phase)

snr = @(ref, y) 10*log10(sum(ref.^2) / sum((ref - y).^2));
fprintf('  SNR noisy input              : %6.2f dB\n', snr(clean, noisy));
fprintf('  SNR FIR (delay compensated)  : %6.2f dB\n', snr(clean, y_fir));
fprintf('  SNR IIR (causal)             : %6.2f dB  <- phase distortion lowers the SNR\n', snr(clean, y_iir));
fprintf('  SNR IIR (forward-backward)   : %6.2f dB\n', snr(clean, y_iir_zero));
fprintf('  FIR coefficients: %d   |   IIR coefficients: %d\n', ...
        length(h_fir_total), length(b_iir)+length(a_iir));

figure('Name','Section 11 - Spectra');
amplitude_spectrum(noisy, Fs, 'r'); hold on;
amplitude_spectrum(y_fir, Fs, 'b');
amplitude_spectrum(y_iir_zero, Fs, 'g');
legend('Noisy','FIR output','IIR output'); title('Spectrum before and after filtering');

figure('Name','Section 11 - Time domain');
rng_idx = round(1.0*Fs):round(1.05*Fs);              % 50 ms close-up
subplot(3,1,1); plot(t(rng_idx), noisy(rng_idx)); title('Noisy input'); grid on;
subplot(3,1,2); plot(t(rng_idx), clean(rng_idx),'k', t(rng_idx), y_fir(rng_idx),'b--');
legend('Clean','FIR'); title('FIR output'); grid on;
subplot(3,1,3); plot(t(rng_idx), clean(rng_idx),'k', t(rng_idx), y_iir(rng_idx),'r--', ...
                     t(rng_idx), y_iir_zero(rng_idx),'g:');
legend('Clean','IIR causal','IIR forward-backward'); title('IIR output'); xlabel('t [s]'); grid on;

%% ===================== SECTION 12: REPORT OUTPUTS ======================
% Specification table, filter coefficients and impulse responses required
% by the grading criteria; all figures are saved as PNG for the report.
fprintf('\n=== SECTION 12: Report outputs ===\n');

% --- 12.1 Filter specification table ---
fprintf('\nFilter specifications (Fs = %d Hz)\n', Fs);
fprintf('%-24s %-14s %-12s %-22s %-24s\n','Filter','Passband','Stopband','Requirement','Design');
fprintf('%-24s %-14s %-12s %-22s %-24s\n','IIR notch','not 50+-6 Hz','50 Hz','null at 50 Hz','order 2, r = 0.995');
fprintf('%-24s %-14s %-12s %-22s %-24s\n','IIR Butterworth LPF','0-1569 Hz','>= 3000 Hz','Ap <= 0.1, As >= 40 dB','N = 6, fc = 2000 Hz');
fprintf('%-24s %-14s %-12s %-22s %-24s\n','FIR high-pass','>= 190 Hz','<= 50 Hz','Ap <= 0.1, As >= 60 dB','Blackman, order 800');
fprintf('%-24s %-14s %-12s %-22s %-24s\n','FIR low-pass','<= 1800 Hz','>= 2200 Hz','Ap <= 0.1, As >= 50 dB','Hamming, order 100');

% Minimum orders from the lecture design formulas
N_min    = log10(1e4 - 1) / (2*log10(tan(0.375*pi)/tan(pi/4)));    % Butterworth, As = 40 dB
M_lp_min = 6.64*pi  / (0.55*pi   - 0.45*pi);                         % Hamming,  dw = 0.1*pi
M_hp_min = 11.12*pi / (0.0475*pi - 0.0125*pi);                       % Blackman, dw = 0.035*pi
fprintf('Minimum orders: Butterworth N >= %.2f, Hamming M >= %.1f, Blackman M >= %.1f\n', N_min, M_lp_min, M_hp_min);

% Verify the specifications by measurement
[H_iirT, wT] = my_freqz(b_iir, a_iir, 8192);
[H_firT, ~ ] = my_freqz(h_fir_total, 1, 8192);
fT = wT/pi*Fs/2;
meas = @(H,f0) 20*log10(abs(H(find(fT >= f0, 1))) + eps);
fprintf('\nMeasured magnitude responses [dB]:\n');
fprintf('%-10s %8s %8s %8s %8s %8s\n','Freq.','50 Hz','262 Hz','1000 Hz','1569 Hz','3000 Hz');
fprintf('%-10s %8.1f %8.1f %8.1f %8.1f %8.1f\n','IIR', meas(H_iirT,50), meas(H_iirT,262), meas(H_iirT,1000), meas(H_iirT,1569), meas(H_iirT,3000));
fprintf('%-10s %8.1f %8.1f %8.1f %8.1f %8.1f\n','FIR', meas(H_firT,50), meas(H_firT,262), meas(H_firT,1000), meas(H_firT,1569), meas(H_firT,3000));

% --- 12.2 Filter coefficients (screen + text files) ---
fprintf('\nButterworth b = %s\n', mat2str(b_lp, 6));
fprintf('Butterworth a = %s\n', mat2str(a_lp, 6));
fprintf('Notch       b = %s\n', mat2str(b_n, 6));
fprintf('Notch       a = %s\n', mat2str(a_n, 6));

% Output folder next to this script (full path, so MATLAB does not pick up
% another 'outputs' folder on the search path). If that folder is not
% writable, fall back to the user's temporary folder.
out_dir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
if isempty(fileparts(mfilename('fullpath'))), out_dir = fullfile(pwd, 'outputs'); end
if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end
fid = fopen(fullfile(out_dir, 'write_test.txt'), 'w');
if fid < 0
    out_dir = fullfile(tempdir, 'dsp_project_outputs');
    if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end
    warning('Script folder is not writable; outputs are saved to %s', out_dir);
else
    fclose(fid); delete(fullfile(out_dir, 'write_test.txt'));
end
fprintf('Output folder: %s\n', out_dir);
fid = fopen(fullfile(out_dir,'iir_coefficients.txt'), 'w');
fprintf(fid, 'k\tb_butter\ta_butter\tb_notch\ta_notch\n');
for k = 1:length(b_lp)
    if k <= length(b_n), bn_k = sprintf('%.8f', b_n(k)); an_k = sprintf('%.8f', a_n(k));
    else, bn_k = '-'; an_k = '-'; end
    fprintf(fid, '%d\t%.8f\t%.8f\t%s\t%s\n', k-1, b_lp(k), a_lp(k), bn_k, an_k);
end
fclose(fid);

fid = fopen(fullfile(out_dir,'fir_lowpass_coefficients.txt'), 'w');   % 101 coefficients
fprintf(fid, 'n\th[n]\n');
fprintf(fid, '%d\t%.8f\n', [0:M-1; h_fir]);
fclose(fid);

fid = fopen(fullfile(out_dir,'fir_highpass_coefficients.txt'), 'w');  % 801 coefficients
fprintf(fid, 'n\th[n]\n');
fprintf(fid, '%d\t%.8f\n', [0:M_hp-1; h_hp]);
fclose(fid);
fprintf('Coefficients written to %s\n', out_dir);

% --- 12.3 Impulse responses ---
Lb = 60;  Ln = 1500;
h_butter = filter(b_lp, a_lp, [1 zeros(1, Lb-1)]);
h_notch  = filter(b_n,  a_n,  [1 zeros(1, Ln-1)]);
h_iirT   = filter(b_iir, a_iir, [1 zeros(1, Ln-1)]);

figure('Name','Section 12 - IIR impulse responses');
subplot(3,1,1); stem(0:Lb-1, h_butter, 'filled', 'MarkerSize', 3); grid on;
title('Butterworth LPF impulse response h[n]'); xlabel('n');
subplot(3,1,2); plot(0:Ln-1, h_notch); grid on;
title('50 Hz notch impulse response h[n] (h[0] \approx 1, then slow decay)'); xlabel('n');
subplot(3,1,3); plot(0:Ln-1, h_iirT); grid on;
title('IIR cascade (notch + Butterworth) impulse response'); xlabel('n');

figure('Name','Section 12 - FIR impulse responses');
subplot(3,1,1); stem(0:M-1, h_fir, 'filled', 'MarkerSize', 3); grid on;
title(sprintf('FIR low-pass (Hamming, M = %d) h[n]', M)); xlabel('n');
subplot(3,1,2); plot(0:M_hp-1, h_hp); grid on;
title(sprintf('FIR high-pass (Blackman, M = %d) h[n]', M_hp)); xlabel('n');
subplot(3,1,3); plot(0:length(h_fir_total)-1, h_fir_total); grid on;
title(sprintf('FIR cascade h_{HP} * h_{LP} (L = %d, delay = %d samples)', ...
      length(h_fir_total), delay)); xlabel('n');

% --- 12.4 Frequency response of the final filters ---
figure('Name','Section 12 - Final filter frequency responses');
subplot(2,1,1);
semilogx(fT(2:end), 20*log10(abs(H_iirT(2:end))+eps), fT(2:end), 20*log10(abs(H_firT(2:end))+eps));
xlim([10 Fs/2]); ylim([-120 5]); grid on; legend('IIR cascade','FIR cascade','Location','southwest');
title('Magnitude response'); xlabel('f [Hz]'); ylabel('dB');
subplot(2,1,2);
plot(fT, unwrap(angle(H_iirT)), fT, unwrap(angle(H_firT))); grid on; xlim([0 Fs/2]);
legend('IIR cascade','FIR cascade','Location','southwest');
title('Phase response'); xlabel('f [Hz]'); ylabel('rad');

% --- 12.5 Save every figure as PNG for the report ---
figs = findobj('Type', 'figure');
fig_dir = fullfile(out_dir, 'figures');
if exist(fig_dir, 'dir') ~= 7, mkdir(fig_dir); end
for i = 1:numel(figs)
    fh = figs(i);
    name = regexprep(get(fh,'Name'), '[^A-Za-z0-9]+', '_');
    print(fh, fullfile(fig_dir, sprintf('fig_%02d_%s', get(fh,'Number'), name)), '-dpng', '-r200');
end
fprintf('%d figures saved to %s\n', numel(figs), fig_dir);

% --- Listening ---
if PLAY_SOUND
    normalize = @(x) x / max(abs(x));
    fprintf('\nPlaying: 1) clean  2) noisy  3) FIR  4) IIR\n');
    sound(normalize(clean),      Fs); pause(N/Fs + 0.5);
    sound(normalize(noisy),      Fs); pause(N/Fs + 0.5);
    sound(normalize(y_fir),      Fs); pause(N/Fs + 0.5);
    sound(normalize(y_iir_zero), Fs);
end
fprintf('\nDone.\n');


%% =====================================================================
%                 HELPER FUNCTIONS (replacing toolbox functions)
% ======================================================================

function y = my_conv(x, h)
% Direct convolution sum: y[n] = sum_k x[k] h[n-k]
    Lx = length(x); Lh = length(h);
    y = zeros(1, Lx+Lh-1);
    for nI = 1:Lx+Lh-1
        for k = max(1, nI-Lh+1):min(nI, Lx)
            y(nI) = y(nI) + x(k) * h(nI-k+1);
        end
    end
end

function X = my_dft(x)
% X[k] = sum_n x[n] e^{-j 2 pi k n / N}   (matrix product, O(N^2))
    x  = x(:).';
    Nx = length(x);
    k  = (0:Nx-1).';
    W  = exp(-1j*2*pi*k*(0:Nx-1)/Nx);
    X  = (W * x.').';
end

function amplitude_spectrum(x, Fs, col)
% Magnitude spectrum of x computed with the DFT (fft), in dB
    L = length(x);
    X = fft(x);
    X = abs(X(1:floor(L/2)+1)) / L;
    X(2:end-1) = 2*X(2:end-1);
    f = (0:floor(L/2)) * Fs / L;
    plot(f, 20*log10(X + 1e-8), col);
    xlabel('f [Hz]'); ylabel('Amplitude [dB]'); grid on; ylim([-100 0]);
end

function [H, w] = my_freqz(b, a, Nf)
% H(e^{jw}) = B(e^{jw}) / A(e^{jw}),  w in [0, pi]
    w = linspace(0, pi, Nf).';
    B = exp(-1j*w*(0:length(b)-1)) * b(:);
    A = exp(-1j*w*(0:length(a)-1)) * a(:);
    H = B ./ A;
end

function my_zplane(b, a)
% Pole-zero plot: convert z^-1 polynomials to z polynomials and find the roots
    L = max(length(b), length(a));
    b = [b(:).' zeros(1, L-length(b))];
    a = [a(:).' zeros(1, L-length(a))];
    z = roots(b); p = roots(a);
    th = linspace(0, 2*pi, 400);
    plot(cos(th), sin(th), 'k:'); hold on;
    plot(real(z), imag(z), 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(real(p), imag(p), 'rx', 'MarkerSize', 10, 'LineWidth', 1.5);
    axis equal; grid on; xlabel('Re\{z\}'); ylabel('Im\{z\}');
    lim = max([1.2; abs(z)*1.1; abs(p)*1.1]);
    axis([-lim lim -lim lim]); hold off;
end

function [b, a] = my_butter(Nb, fc, Fs)
% Butterworth LPF: analog prototype + prewarping + bilinear transformation
%   s = 2Fs (1 - z^-1)/(1 + z^-1)
    Wc = 2*Fs*tan(pi*fc/Fs);                         % prewarped analog cutoff [rad/s]
    k  = 1:Nb;
    s_poles = Wc * exp(1j*pi*(2*k + Nb - 1)/(2*Nb)); % left-half-plane poles
    z_poles = (2*Fs + s_poles) ./ (2*Fs - s_poles);  % bilinear mapping
    a = real(poly(z_poles));
    b = real(poly(-ones(1, Nb)));                    % N zeros at z = -1
    b = b * sum(a) / sum(b);                         % unity DC gain
end

function [b, a] = notch_filter(f0, Fs, r)
% Second-order notch: zeros at e^{+-jw0}, poles at r e^{+-jw0}
    w0 = 2*pi*f0/Fs;
    b = [1, -2*cos(w0), 1];
    a = [1, -2*r*cos(w0), r^2];
    b = b * sum(a) / sum(b);                         % unity DC gain
end

function w = make_window(M, type)
    nI = 0:M-1;
    switch lower(type)
        case 'rectangular', w = ones(1, M);
        case 'hann',        w = 0.5  - 0.5 *cos(2*pi*nI/(M-1));
        case 'hamming',     w = 0.54 - 0.46*cos(2*pi*nI/(M-1));
        case 'blackman',    w = 0.42 - 0.5 *cos(2*pi*nI/(M-1)) + 0.08*cos(4*pi*nI/(M-1));
        otherwise, error('Unknown window: %s', type);
    end
end

function h = windowed_fir(M, wc, type)
% Ideal LPF impulse response h_d[n] = sin(wc m)/(pi m), m = n-(M-1)/2, then window
    m  = (0:M-1) - (M-1)/2;
    hd = sin(wc*m) ./ (pi*m);
    hd(m == 0) = wc/pi;
    h  = hd .* make_window(M, type);
    h  = h / sum(h);                                 % unity DC gain
end

function gd = group_delay(H, w)
% tau(w) = -d(phase)/dw  (numerical derivative)
    gd = -diff(unwrap(angle(H))) ./ diff(w);
end

function [r, p] = inv_z_partial(b, a)
% H(z) = B(z^-1)/A(z^-1). After padding to equal length and multiplying by z^L,
% b and a are polynomials in descending powers of z. One extra z is added to the
% denominator to expand H(z)/z.
    L = max(length(b), length(a));
    b = [b(:).' zeros(1, L-length(b))];
    a = [a(:).' zeros(1, L-length(a))];
    [r, p] = residue(b, [a 0]);
end

function h = inv_z_eval(r, p, nn)
% h[n] = sum r_k p_k^n   (the p = 0 term contributes only at n = 0; 0^0 = 1)
    h = zeros(size(nn));
    for k = 1:length(r)
        h = h + r(k) * p(k).^nn;
    end
    h = real(h);
end

function y = zero_phase_filter(b, a, x)
% Filter forward, reverse, filter again, reverse -> zero phase, magnitude |H|^2
    y = filter(b, a, x);
    y = fliplr(filter(b, a, fliplr(y)));
end

function my_xline(f)
% Vertical line compatible with older MATLAB versions
    yl = ylim; hold on; plot([f f], yl, 'k--'); hold off;
end
