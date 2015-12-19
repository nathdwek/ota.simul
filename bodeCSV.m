function [] = bodeCSV(magFile,phFile, sys, titleStr)
    mag = importdata(magFile);
    ph = importdata(phFile);
    w = 2*pi*mag.data(:,1);
    [magth, phth] = bode(sys,w);
    figure;subplot(2,1,1);
    semilogx(mag.data(:,1),mag.data(:,2), mag.data(:,1), 20*log10(squeeze(magth(1,1,:))));grid on;
    title(sprintf('Bodeplot: %s', titleStr));
    legend('Spectre', 'Matlab')
    xlabel('Frequency [Hz]');ylabel('Magnitude [dB]');
    subplot(2,1,2);
    semilogx(ph.data(:,1),ph.data(:,2), mag.data(:,1), squeeze(phth(1,1,:)));grid on;
    xlabel('Frequency [Hz]');ylabel('Phase [degs]');
    legend('Spectre','Matlab')
end
