function [score, idx_] = select_ica_components(z, Xica, sf, fo_, naive, doplot, out_figures)
% selects ICA components based on the PLV wrt to the audio signal z

if nargin < 7
  out_figures = [];
end
if nargin < 6
   doplot = 0;
end

XCORR = [];
PLV=[];
bands= [50, 250];
[bf, af] = butter(5, bands./(sf/2));

bands_h= 1;
[bh, ah] = butter(5, bands_h./(sf/2), 'high');
npoints = 0.2*sf;


ntoplot = size(Xica,1); 
    
for k =1:ntoplot
xpcd = Xica(k,:);
ica_filt = filtfilt(bf, af, xpcd);
ica_filt = (ica_filt-mean(ica_filt(:)))./std(ica_filt(:));

z_n = (z-mean(z(:)))./std(z(:));


%% phase
ph = angle(hilbert(z_n));
ph_c = angle(hilbert(ica_filt)); 
dif_ph = ph_c - ph;
plv_art = circ_r(dif_ph');

PLV = [PLV, plv_art];
%% pearson corr
xcorr = corrcoef(ica_filt, z_n);
XCORR = [XCORR, abs(xcorr(1,2))];

if naive
    %sort components
    [score, idx_] = sort(XCORR, 'descend');
else
    [score, idx_] = sort(PLV, 'descend');
end

end

%% plot
if doplot

% sort sources to plot
Xica_sort=Xica(idx_,:);
for k =1:ntoplot
    xpcd = Xica_sort(k,:);

    xpcd_h = filtfilt(bh, ah,xpcd);
    z_h = filtfilt(bh, ah,z);

    % power spectrum of the reconstructed signal (denoised)
    [pxx_d,f] = pwelch(detrend(xpcd_h,0),hanning(npoints),[],2^nextpow2(npoints)*5,sf);
    % power spectrum of the ground truth signal
    [pxx_gt,f] = pwelch(detrend(z_h,0),hanning(npoints),[],2^nextpow2(npoints)*5,sf);


    %% performing time frequency analysis

    nwin = 100; nvlp = 80; 
    fint = 2:2:250;  
    [Bz,f_sp,T] = spectrogram(z_h, hamming(nwin), nvlp,fint,sf, 'power');
    Bz = 20*log10(abs(Bz));
    B_idx = Bz < max(max(Bz))-60;
    Bz(B_idx) = 0;



    %% spectrogram
    [Bx,f_sp,T] = spectrogram(xpcd_h, hamming(nwin), nvlp,fint,sf, 'power');
    Bx = 20*log10(abs(Bx));
    B_idx = Bx < max(max(Bx))-60;
    Bx(B_idx) = 0;

    %% phase
    ph = angle(hilbert(z_n));
    ph_c = angle(hilbert(xpcd)); 
    dif_ph = ph_c - ph;

    %% pearson corr
    xpec = get_xsimilarities(Bz, Bx, 'pixel-xcorr');
    ica_filt = filtfilt(bf, af, xpcd);
    ica_filt = (ica_filt-mean(ica_filt(:)))./std(ica_filt(:));


    ax = figure('units','normalized','outerposition',[0 0 1 1], "renderer","painters");
    title_fig_f = [' ICA component ', int2str(k)];
    sgtitle(title_fig_f)


    subplot(2,4,1)
    plot(z, 'r')
    title('raw Audio')
    xlabel(" Time [s] ")

    subplot(2,4,5)
    plot(ica_filt, 'k')
    title('ICA component')
    xlabel(" Time [s] ")

    subplot(2,4,6)
    semilogy(f,pxx_d,'linestyle','-','color','k','linewidth',1.4)
    hold on
    semilogy([fo_ fo_],[1e-6 max(pxx_d)]);

    grid on
    grid minor
    xlabel(" Frequency [Hz] ")
    ylabel( " Power [dB] ")
    title('Power ICA');
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
    
    subplot(2,4,3)
    imagesc(T, f, Bz);
    axis xy;
    colorbar
    xlabel('Time (s)');ylabel('Frequency (Hz)'); 
    title("Spectogram Audio");
    
    subplot(2,4,7)

    imagesc(T, f, Bx);
    axis xy;
    colorbar
    xlabel('Time (s)');ylabel('Frequency (Hz)'); 
    title('Spectogram ICA');

    subplot(2,4,4)
    % circ_plot(dif_ph','pretty','ro',true,'linewidth',2,'color','r');
    polarhistogram(dif_ph)
    title('Phase Difference');

    subplot(2,4,8)
    imagesc(T, f, xpec);
    axis xy;
    colorbar
    xlabel('Time (s)');ylabel('Frequency (Hz)'); 
    title('Xspectrogram ');
    if ~isempty(out_figures)
        fig_name = strcat(out_figures, '_ICA ', int2str(k));
        saveas(ax, fig_name, 'png')
    end
end
end
end