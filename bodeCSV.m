function [] = bodeCSV(magFile,phFile, sys, legendStr)
    mag = importdata(magFile);
    ph = importdata(phFile);
    w = 2*pi*mag.data(:,1);
    [magth, phth] = bode(sys,w);
    subplot(2,1,1);
    semilogx(mag.data(:,1),mag.data(:,2), 'DisplayName',sprintf('%s: Spectre', legendStr));
    hold all;
    semilogx(mag.data(:,1), 20*log10(squeeze(magth(1,1,:))), 'DisplayName', sprintf('%s: Matlab', legendStr));
    grid on;
    xlabel('Frequency [Hz]');ylabel('Magnitude [dB]');
    subplot(2,1,2);
    semilogx(ph.data(:,1),ph.data(:,2), 'DisplayName',sprintf('%s: Spectre', legendStr));
    hold all;
    semilogx(mag.data(:,1), squeeze(unwrap(wrapTo180(phth(1,1,:)))), 'DisplayName', sprintf('%s: Matlab', legendStr));
    grid on;
    xlabel('Frequency [Hz]');ylabel('Phase [degs]');
    %legend('Spectre','Matlab')
    set(findall(gcf,'-property','FontSize'),'FontSize',16) 
end
