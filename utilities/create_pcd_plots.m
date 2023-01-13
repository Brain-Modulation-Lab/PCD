function XCORR = create_pcd_plots(z, Xpcd, vlen, sf, out_figures, title_fig, fo_, idxplot)
XCORR = [];
bands= [50, 250];
[bf, af] = butter(5, bands./(sf/2));

bands_h= 1;
[bh, ah] = butter(5, bands_h./(sf/2), 'high');
npoints = 0.2*sf;

if nargin > 7
    ki=idxplot;
    ntoplot = idxplot;
else
    ki=1;
    ntoplot = size(Xpcd,1); 
end
    
for k =ki:ntoplot
xpcd = Xpcd(k,:);
pcd_filt = filtfilt(bf, af, xpcd);
pcd_filt = (pcd_filt-mean(pcd_filt(:)))./std(pcd_filt(:));

xpcd_h = filtfilt(bh, ah,xpcd);
z_h = filtfilt(bh, ah,z);

% power spectrum of the reconstructed signal (denoised)
[pxx_d,f] = pwelch(detrend(xpcd_h,0),hanning(npoints),[],2^nextpow2(npoints)*5,sf);
% power spectrum of the ground truth signal
[pxx_gt,f] = pwelch(detrend(z_h,0),hanning(npoints),[],2^nextpow2(npoints)*5,sf);

%%

ax = figure('units','normalized','outerposition',[0 0 1 1], "renderer","painters");
title_fig_f = [title_fig ' PCD component ', int2str(k)];
sgtitle(title_fig_f)


subplot(2,4,1)
plot(z, 'r')
title('raw Audio')
xlabel(" Time [s] ")

subplot(2,4,5)
plot(pcd_filt, 'k')
title('PCD component')
xlabel(" Time [s] ")



subplot(2,4,6)
semilogy(f,pxx_d,'linestyle','-','color','k','linewidth',1.4)
hold on
semilogy([fo_ fo_],[1e-6 max(pxx_d)]);

grid on
grid minor
xlabel(" Frequency [Hz] ")
ylabel( " Power [dB] ")
title('Power PCD');
xlim([0 200])

subplot(242)
semilogy(f,pxx_gt,'linestyle','-','color','k','linewidth',1.4)
hold on
semilogy([fo_ fo_],[1e-6 max(pxx_gt)]);

grid on
grid minor
xlim([0 200])
xlabel(" Frequency [Hz] ")
ylabel( " Power [dB] ")
title("Power Audio");


%% performing time frequency analysis

nwin = 100; nvlp = 80; 
fint = 2:2:250;  
[Bz,f,T] = spectrogram(z_h, hamming(nwin), nvlp,fint,sf, 'power');
Bz = 20*log10(abs(Bz));
B_idx = Bz < max(max(Bz))-60;
Bz(B_idx) = 0;
subplot(2,4,3)
imagesc(T, f, Bz);
axis xy;
colorbar
xlabel('Time (s)');ylabel('Frequency (Hz)'); 
title("Spectogram Audio");


%% spectrogram
[Bx,f,T] = spectrogram(xpcd_h, hamming(nwin), nvlp,fint,sf, 'power');
Bx = 20*log10(abs(Bx));
B_idx = Bx < max(max(Bx))-60;
Bx(B_idx) = 0;
subplot(2,4,7)

imagesc(T, f, Bx);
axis xy;
colorbar
xlabel('Time (s)');ylabel('Frequency (Hz)'); 
title('Spectogram PCD');


%% phase
ph = angle(hilbert(z));
ph_c = angle(hilbert(xpcd)); 
dif_ph = ph_c - ph;
subplot(2,4,4)
% circ_plot(dif_ph','pretty','ro',true,'linewidth',2,'color','r');
polarhistogram(dif_ph)
title('Phase Difference. Vlen: ', num2str(vlen(k)));
%% xspect
subplot(2,4,8)
% circ_plot(dif_ph','pretty','ro',true,'linewidth',2,'color','r');
xpec = get_xsimilarities(Bz, Bx, 'pixel-xcorr');
xcorr = get_xsimilarities(Bz, Bx, 'corr2');
XCORR = [XCORR; xcorr];

imagesc(T, f, xpec);
axis xy;
colorbar
xlabel('Time (s)');ylabel('Frequency (Hz)'); 
title(['Xspectrogram. xcorr: ', num2str(xcorr)]);
if ~isempty(out_figures)
    fig_name = strcat(out_figures, '_PCD ', int2str(k));
    saveas(ax, fig_name, 'png')
end
end
end