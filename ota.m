%% Loading MOS tables

addpath(genpath('circuitDesign'));
clear all;
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
mPhi.spec       = 60;        % [deg] Phase margin
swing.spec      = 0.8;       % [V] Output swing

%% size mn6 for gm
Mn6.vds = VDD.spec/2; %First, maximize output swing.
Mn6.vds = Mn6.vds + 0.1;%But it is Mn6 who produces gain, so give him a tad more vds
Mn6.vsb = 0;
Mn6.gm = 2*pi*2.5*fGBW.spec*Cl.spec;
Mn6.gm = Mn6.gm/1.7;%Reduce the spec on Mn6.gm thanks to nulling resistor;
Mn6.vov = 0;%Reduce Mn6.ids as much as possible. Also good for swing
%Also: low inversion is good for low vdsat and thus high output swing.
Mn6.lg = 132e-9;%Lg rises: W rises: upper limit W<500e-6
%Lg decreases: gds increases. lower limit set by gain.
%Lg rises: better intrinsic gain
Mn6.nFingers = 2;
%For some reason, we need nFingers to be set before calling transistor
%functions. So we set it to 2, and call those
%But nFingers need to be unset to be able to call mosNFingers. So
%afterwards we delete the field and call mosNfingers.
%This is done for every transistor.
%The choice of nFingers before calling the transistor functions doesn't
%seem to affect the results.

Mn6.vth = tableValueWref('vth',NRVT,Mn6.lg,0,Mn6.vds,Mn6.vsb);
Mn6.vgs = Mn6.vov + Mn6.vth;

Mn6.w = mosWidth('gm', Mn6.gm, Mn6);
fprintf('\nMn6.w = %f um\n', Mn6.w*1e6);

Mn6 = rmfield(Mn6, 'nFingers');
Mn6 = mosNfingers(Mn6);
Mn6 = mosOpValues(Mn6);
stage2Current = Mn6.ids;
fprintf('\nStage 2: %f mA\n', stage2Current*1e3);

mosCheckSaturation(Mn6);

%% size Mp1,2 for gm
Cm = Cl.spec/5;%First try value. Too big: load stage1 too much.
%Too small: load Mn6 too much
Mp2.nFingers = 2;
Mp2.gm = 2*pi*fGBW.spec*Cm*1.2;
%We need a bit more gain than predicted here
Mp2.lg = 212e-9;%Minimize gds while staying in width bounds
Mp2.vov = 0.05;%weak inversion
Vdb2 = Mn6.vgs - VDD.spec;
Vgb2 = Vdb2 + 0.25;
%Saturation: vds<vov (pmos). <=> vgb > vdb + vth where vth is found
%iteratively. (here ~ -0.3V).
%if we go deeper in saturation: vgb > vdb - 0.3V: more intrinsic gain and
%smaller w. (or a longer transistor if w is desired just under 500um)
Mp2.vsb = mosVsbBody(Mp2, Vgb2, Vdb2, Mp2.vov, -0.2);
Mp2.vgs = Vgb2 - Mp2.vsb;
Mp2.vds = Vdb2 - Mp2.vsb;

Mp2.w = mosWidth('gm',Mp2.gm, Mp2);
fprintf('\nMp2.w = %f um\n', Mp2.w*1e6);

Mp2 = rmfield(Mp2, 'nFingers');
Mp2 = mosNfingers(Mp2);
Mp2 = mosOpValues(Mp2);
stage1Current = 2*Mp2.ids;
fprintf('\nStage 1: %f mA\n', stage1Current*1e3);

Mp1 = cirElementCopy(Mp2,Mp1);

mosCheckSaturation(Mp2);


%% Size Mn3,4 for ids
Mn4.vds = Mn6.vgs;
Mn4.vgs = Mn4.vds;
Mn4.vsb = 0;
Mn4.lg = 400e-9;%minimize gds
Mn4.ids = Mp2.ids;
Mn4.nFingers = 2;

Mn4.w = mosWidth('ids', Mn4.ids, Mn4);
fprintf('\nMn4.w = %f um\n', Mn4.w*1e6);

Mn4 = rmfield(Mn4, 'nFingers');
Mn4 = mosNfingers(Mn4);
Mn4 = mosOpValues(Mn4);

Mn3 = cirElementCopy(Mn4,Mn3);

mosCheckSaturation(Mn4);

%% size Mp5 for ids
Mp7.lg = 500e-9;%Mp7 should be an ideal current source
%Maximize lg of Mp5, Mp7 until W of any <= 500 um
%And better linearity if all transistor in a given mirror have the same L
%But actually increasing lg beyond 400 doesn't seem to lower gds 
Mp5.lg = Mp7.lg;
Mp5.nFingers = 2;
Mp5.vsb = 0;
Mp5.vds = Mn6.vds - VDD.spec;

Mp5.vov = -0.2;%increase |vov| to achieve smaller w
%and thus smaller parasitic capacitances
%But this sets Mp5.vgs and we need output swing.
%Also decrease |vov| to decrease gds
%Overall this doesn't change the final vales of the ota much except the
%output swing
Mp5.vth = tableValueWref('vth', PRVT, Mp5.lg, 0, Mp5.vds, Mp5.vsb);
Mp5.vgs = Mp5.vov + Mp5.vth;

Mp5.ids = stage2Current;
Mp5.w = mosWidth('ids', Mp5.ids, Mp5);
fprintf('\nMp5.w = %f um\n', Mp5.w*1e6);

Mp5 = rmfield(Mp5, 'nFingers');
Mp5 = mosNfingers(Mp5);
Mp5 = mosOpValues(Mp5);

mosCheckSaturation(Mp5);


%% size Mp7 for ids, high output impedance
%We chose Mp7.lg previously since it's the same for Mp5, Mp7, Mp8
Mp7.nFingers = 2;
Mp7.vsb = 0;
Mp7.vds = Mn6.vgs - Mp2.vds - VDD.spec;
Mp7.vgs = Mp5.vgs;

Mp7.ids = stage1Current;
Mp7.w = mosWidth('ids', Mp7.ids, Mp7);
fprintf('\nMp7.w = %f um\n',Mp7.w*1e6);

Mp7 = rmfield(Mp7, 'nFingers');
Mp7 = mosNfingers(Mp7);
Mp7 = mosOpValues(Mp7);

mosCheckSaturation(Mp7);


%% Size Mp8
Mp8.lg = Mp5.lg;%Better linearity if all transistor in a given mirror have the same L
Mp8.vsb = 0;
Mp8.vgs = Mp7.vgs;
Mp8.vds = Mp8.vgs;
Mp8.nFingers = 2;

Mp8.ids = min(stage1Current,stage2Current)/100;
%Mp8 should not contribute significantly to current consumption!
Mp8.w = mosWidth('ids', Mp8.ids, Mp8);

Mp8 = rmfield(Mp8, 'nFingers');
Mp8 = mosNfingers(Mp8);
Mp8 = mosOpValues(Mp8);

mosCheckSaturation(Mp8);



%% calculate gain, dominant pole, GBW = gain* dominant pole (verify)
% If exercise
%Cm = 0;
%Cl.spec = 1.5 * Cl.spec;
%Cm = 1.2 * Cm;
%
gain = Mp2.gm/(Mp2.gds + Mn4.gds) * Mn6.gm/(Mn6.gds + Mp5.gds);
p1 = -(Mp2.gds + Mn4.gds)/((Cm + Mn6.cgd)*Mn6.gm/(Mn6.gds + Mp5.gds) + Mn6.cgs + Mn6.cgb + Mn4.cdb + Mp2.cdb + Mp2.cgd);
p2 = -(Mp5.gds + Mn6.gds + Mn6.gm * (Cm/(Cm + Mp5.cdb + Mn6.cdb + Mp5.cgd))) / (Cl.spec + Cm + Mp5.cdb + Mn6.cdb + Mp5.cgd);
Rm = -1/(Cm*p2)*(1 - p2*Cm/Mn6.gm);
%If exercise
%Rm = 0;
%Rm = 213.6;
%
fprintf('\n Miller resistor = %f ohms', Rm);
z1 = 1/(Cm*(1/Mn6.gm - Rm));
p3 = -Mn3.gm/(Mp1.cdb+Mp1.cgd+Mn3.cdb+Mn3.cgs+Mn3.cgb);
z3 = 2*p3;
p4 = -1/(Rm*(Mn6.cgs + Mn6.cgb + Mn4.cdb + Mp2.cdb + Mp2.cgd));

%If exercise
%sys = tf(gain,[1/(p1*p2) -(p1+p2)/(p1*p2) 1]);%This is with no nulling
%resistor
%
sys = tf([-gain/z1 gain],[1/(p1*p2) -(p1+p2)/(p1*p2) 1]);
sys = series(sys, tf([-1/z3 1],[-1/p3 1]));
sys = series(sys, tf(1, [-1/p4 1]));
bodeplot(sys);grid on;

domPole.real = p1/(2*pi);
fGBW.real = -gain*p1/(2*pi);
p2OverGBW.real = p2/(gain*p1);
totCurrent.real = stage1Current + stage2Current + Mp8.ids;
swing.real = VDD.spec + Mp5.vdsat - Mn6.vdsat;
[~,mPhi.real,~,~] = margin(sys);
inputRange = VDD.spec + Mp7.vov + Mp1.vov + Mp1.vth - Mp1.vth - Mn3.vth - Mn3.vov;


fprintf('\n\n--<OTA values>--\n\n')
fprintf('\ngain = %f dB\n', 20*log10(gain))
fprintf('\nDominant pole = %f kHz\n', -domPole.real/1e3)
fprintf('\nfGBW = %f MHz\n', fGBW.real/1e6);
fprintf('\np2/gbw = %f\n', p2OverGBW.real)
fprintf('\nCurrent consumption = %f mA\n', totCurrent.real/1e-3);
fprintf('\nfom = %f MHz*pF/mA\n', (fGBW.real/1e6)*(Cl.spec/1e-12)/(totCurrent.real/1e-3));
fprintf('\nOutput swing = %f V\n', swing.real);
fprintf('\nPhase margin = %f°\n', mPhi.real);
fprintf('\nIbias = %fmA\n',Mp8.ids*1e3);
fprintf('\nInput Common Mode Range: %.0fmV\n', inputRange*1e3)

%% Noise Calculations

oneTwoFN = 2*Mp1.di2_fn*(log(100e9) - log(1))/Mp1.gm^2;
threeFourFN = 2*Mn3.di2_fn*(log(100e9) - log(1))/Mp1.gm^2;
threeFourID = 2*Mn3.di2_id*(100e9 - 1)/Mp1.gm^2;
oneTwoID = 2*Mp1.di2_id*(100e9 - 1)/Mp1.gm^2;
fprintf('\nApproximated input referred voltage noise power:%fVrms\n',sqrt(oneTwoFN+threeFourFN+threeFourID+oneTwoID) );
