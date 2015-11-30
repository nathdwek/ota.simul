close all;
%% Loading MOS tables

addpath(genpath('circuitDesign'));
clear all;
close all;
load ('UMC65_RVT.mat');
%% Initialize everything
designkitName = 'umc65';
circuitTitle= 'Analog Design - Session 2';

%Declaration of the circuit components
elementList.nmos = {'Mn3','Mn4','Mn6'};
elementList.pmos = {'Mp1', 'Mp2', 'Mp5', 'Mp7', 'Mp8'};
 
spec=0;
choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 1e-6;
simulator='spectre';
simulFile=0;
simulSkelFile=0;
analog = cirInit('analog', circuitTitle, 'top', elementList, spec , choice, ...
    designkitName, NRVT, PRVT, simulator, simulFile,simulSkelFile);
analog = cirCheckInChoice(analog, choice);

%% Specs
fprintf('\n--- Project: 2 stage diff amp ---\n');
fGBW.spec       = 70e6;      % [Hz] GBW frequency
Cl.spec         = 50e-12;    % [F] load capacitance
VDD.spec        = 1.1;       % [V] Power supply voltage
DCGain.spec     = 10^(48/20);% []  DC voltage gain
mPhi.spec       = pi/3;      % [rad] Phase margin
swing.spec      = 0.8;       % [V] Output swing

%% size mn6 for gm
Mn6.vds = VDD.spec/2; %First, maximize output swing.
Mn6.vsb = 0;
Mn6.gm = 2*pi*2.5*fGBW.spec*Cl.spec;
Mn6.vov = 0.05;% To ensure gm/Ids = 10 knowing gm = Ids/2Vov
%Also: low inversion is good for low vdsat and thus high output swing.
Mn6.lg = 150e-9;%low lg because we need current to have gain
%and we don't want a super wide transistor. Lower limit is set by the fact
%that we need gds to not be that high or else stage 2 does not give any
%gain whatsoever.
Mn6.nFingers = 2;

Mn6.vth = tableValueWref('vth',NRVT,Mn6.lg,0,Mn6.vds,Mn6.vsb);
Mn6.vgs = Mn6.vov + Mn6.vth;

Mn6.w = mosWidth('gm', Mn6.gm, Mn6);
Mn6 = mosOpValues(Mn6);
stage2Current = Mn6.ids

if mosCheckSaturation(Mn6)
	fprintf('\nMn6 is in sat\n')
end

%% size Mp1,2 for gm
Cm = 10e-12;%First try value. Too big: load stage1 too much.
%Too small: load Mn6 too much
Mp2.nFingers = 2;
Mp2.gm = 2*pi*fGBW.spec*Cm*1.1;
%We need a bit some more gain than predicted here
Mp2.lg = 200e-9;%pif
Mp2.vov = -0.05;% To ensure gm/Ids = 10 knowing gm = Ids/2Vov
Vdb2 = Mn6.vgs - VDD.spec
Vgb2 = Vdb2 + 0.1
Mp2.vsb = mosVsbBody(Mp2, Vgb2, Vdb2, Mp2.vov, -0.2);
Mp2.vgs = Vgb2 - Mp2.vsb;
Mp2.vds = Vdb2 - Mp2.vsb;

Mp2.w = mosWidth('gm',Mp2.gm, Mp2);
Mp2 = mosOpValues(Mp2);
stage1Current = 2*Mp2.ids

Mp1 = cirElementCopy(Mp2,Mp1);

if mosCheckSaturation(Mp2)
 	fprintf('\nMp1,2 in sat\n')
end

%% Size Mn3,4 for ids
Mn4.vds = Mn6.vgs;
Mn4.vgs = Mn4.vds;
Mn4.vsb = 0;
Mn4.lg = 300e-9;%try to not have gds too much higher than Mp2.gds
Mn4.ids = Mp2.ids;
Mn4.nFingers = 2;

Mn4.w = mosWidth('ids', Mn4.ids, Mn4);
Mn4 = mosOpValues(Mn4);

Mn3 = cirElementCopy(Mn4,Mn3);

if mosCheckSaturation(Mn4)
 	fprintf('\nMn3,4 in sat\n')
end

%% size Mp7 for ids, high output impedance
Mp7.lg = 900e-9;%Really long for high output impedance
Mp7.nFingers = 2;
Mp7.vsb = 0;
Mp7.vds = Mn6.vgs - Mp2.vds - VDD.spec;

Mp7.vov = -0.3;%increase |vov| to achieve smaller w
%and thus smaller parasitic capacitances
%But this sets Mp5.vgs and we need output swing.
Mp7.vth = tableValueWref('vth', PRVT, Mp7.lg, 0, Mp7.vds, Mp7.vsb);
Mp7.vgs = Mp7.vov + Mp7.vth;

Mp7.ids = stage1Current;
Mp7.w = mosWidth('ids', Mp7.ids, Mp7);
fprintf('\nMp7.w = %f\n',Mp7.w);

Mp7 = mosOpValues(Mp7);

if mosCheckSaturation(Mp7)
 	fprintf('\nMp7 in sat\n');
end

%% size Mp5 for ids
Mp5.lg = 120e-9;%Smallest lg that matches Mn6 in terms of gds.
Mp5.nFingers = 2;
Mp5.vsb = 0;
Mp5.vds = Mn6.vds - VDD.spec;
Mp5.vgs = Mp7.vgs;

Mp5.ids = stage2Current;
Mp5.w = mosWidth('ids', Mp5.ids, Mp5);
fprintf('\nMp5.w = %f\n', Mp5.w);

Mp5 = mosOpValues(Mp5);

if mosCheckSaturation(Mp5)
 	fprintf('\nMp5 in sat\n')
end

%% Size Mp8
%No real constraints on small signal parameter of Mp8?
Mp8.lg = 200e-9;
Mp8.vsb = 0;
Mp8.vgs = Mp7.vgs;
Mp8.vds = Mp8.vgs;
Mp8.nFingers = 2;

Mp8.ids = max(stage1Current,stage2Current)/4;
%B of a current mirror never higher than 5
Mp8.w = mosWidth('ids', Mp8.ids, Mp8);

Mp8 = mosOpValues(Mp8);
fprintf('\nMp8.w = %f\n', Mp8.w);

if mosCheckSaturation(Mp8)
 	fprintf('\nMp8 in sat\n')
end

%% calculate gain, dominant pole, GBW = gain* dominant pole (verify)
gain = Mp2.gm/(Mp2.gds + Mn4.gds) * Mn6.gm/(Mn6.gds + Mp5.gds);

p1 = -(Mp2.gds + Mn4.gds)/((Cm + Mn6.cgd)*Mn6.gm/(Mn6.gds + Mp5.gds) + Mn6.cgs + Mn6.cgb + Mn4.cdb + Mp2.cdb + Mp2.cgd);
p2 = -(Mp5.gds + Mn6.gds + Mn6.gm * (Cm/(Cm + Mp5.cdb + Mn6.cdb + Mp5.cgd))) / (Cl.spec + Cm + Mp5.cdb + Mn6.cdb + Mp5.cgd);
Rm = -1/(Cm*p2)*(1 - p2*Cm/Mn6.gm);
fprintf('\n Miller resistor = %f ohms', Rm);
z1 = 1/(Cm*(1/Mn6.gm - Rm));
p4 = -1/(Rm*(Mn6.cgs + Mn6.cgb + Mn4.cdb + Mp2.cdb + Mp2.cgd));

sys = tf([-gain/z1 gain],[1/(p1*p2) -(p1+p2)/(p1*p2) 1]);
sys = series(sys, tf(1, [-1/p4 1]));
bodeplot(sys);grid on;

domPole.real = p1/(2*pi);
fGBW.real = -gain*p1/(2*pi);
p2OverGBW.real = p2/(gain*p1);
totCurrent.real = stage1Current + stage2Current;
swing.real = VDD.spec + Mp5.vdsat - Mn6.vdsat;
[~,mPhi.real,~,~] = margin(sys);

fprintf('\n\n--<OTA values>--\n\n')
fprintf('\ngain = %f dB\n', 20*log10(gain))
fprintf('\nDominant pole = %f kHz\n', -domPole.real/1e3)
fprintf('\nfGBW = %f MHz\n', fGBW.real/1e6);
fprintf('\np2/gbw = %f\n', p2OverGBW.real)
fprintf('\nCurrent consumption = %f mA\n', totCurrent.real/1e-3);
fprintf('\nfom = %d MHz*pF/mA\n', (fGBW.real/1e6)*(Cl.spec/1e-12)/(totCurrent.real/1e3));
fprintf('\nOutput swing = %f V\n', swing.real);
fprintf('\nPhase margin = %f°\n', mPhi.real);
