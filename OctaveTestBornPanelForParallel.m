% The parallelized Octave version of testBornPanelForParallel
Octave_driver_timer = tic();
addpath('../pointbem');
addpath('../panelbem');
loadConstants
pkg load parallel

epsIn  =  4;
epsOut = 80;
conv_factor = 332.112;
q=1;
asymParams = struct('alpha',0.0, 'beta', -60.0,'EfieldOffset',-0.5);
kappa_list = [0.01];% 0.500];

pqrData = struct('xyz',[0 0 0],'q',q,'R',0);



srfNoSternFile = 'born/born_1A_4.srf'; 
srfData = loadSrfIntoPanels(srfNoSternFile);

srfSternFile   = 'born/stern_1A_2.srf';
srfSternData = loadSternSrfIntoPanels(srfSternFile);

% part 1: no-salt reference calculation: PCM/NLBC formulation
% part 1 currently uses Octave versions of files
disp("\n");
disp("start part 1 in driver");
OctaveBemPcm = OctaveMakePanelBemEcfQualMatrices(srfData, pqrData, epsIn, epsOut);
disp("from driver, after OctaveMakePanelBemEcfQualMatrices");
asymBemPcm = makePanelAsymEcfCollocMatrices(srfData, OctaveBemPcm, pqrData);
[phiReac, sigma] = solvePanelConsistentAsymmetric(srfData, OctaveBemPcm, ...
						        epsIn, epsOut, ...
						        conv_factor, ...
						        pqrData, ...
						        asymParams,asymBemPcm);
dG_nosalt_pcm = 0.5 * pqrData.q' * phiReac;

disp("part 1 in driver is done \n");


% part 2: no-salt reference calculation: YL/NLBC formulation
% part 2 currently using original files meant for matlab parfor,
% not using all Octave versions of files
disp("start part 2 in driver");
OctaveBemYoonDiel = OctaveMakePanelBemYoonDielMatrices(srfData,pqrData,epsIn,epsOut);
asymBemYL = asymBemPcm; % we're reusing asymBemPcm
[phiReacYL, phiBndy, dphiDnBndy] = ...
    solvePanelConsistentYoonNoSternAsym(srfData, OctaveBemYoonDiel, ...
					epsIn, epsOut, ...
					conv_factor, pqrData, ...
					asymParams, asymBemYL);
dG_nosalt_yl = 0.5 * pqrData.q' * phiReacYL;
disp("part 2 in driver is done \n");

Octave_driver_timerVal = toc(Octave_driver_timer);
disp("Octave_driver_timerVal:"), disp(Octave_driver_timerVal);
return
%{
for l=1:length(kappa_list)
        kappa = kappa_list(l)
	% part 3: salt without Stern, YL/NLBC formulation
	bemYoonNoStern = makePanelBemYoonLPBMatrices(srfData, pqrData, ...
						     epsIn, epsOut, ...
						     kappa);
	% can re-use asymBemYL from no-salt YL/NLBC case because it
        % only involves the solute cavity
	[phiReacNoStern, phiBndy, dphiDnBndy] = ...
	    solvePanelConsistentYoonNoSternAsym(srfData, bemYoonNoStern, ...
						epsIn, epsOut, ...
						conv_factor, pqrData, ...
						asymParams, ...
						asymBemYL);
	dG_nostern = 0.5 * pqrData.q' * phiReacNoStern;

	% part 4: salt with Stern, YL/NLBC formulation
	bemYoonStern = makePanelBemSternMatrices(srfSternData, ...
						 pqrData, epsIn, ...
						 epsOut, kappa);
	% again, can re-use asymBemYL
	[phiReacStern, phiBndy, dphiDnBndy] = ...
	    solvePanelConsistentSternAsym(srfSternData.dielBndy(1), ...
					  srfSternData.sternBndy(1), ...
					  pqrData, bemYoonStern, ...
					  epsIn, epsOut, kappa, ...
					  conv_factor, asymParams, ...
					  asymBemYL);
	dG_stern = 0.5 * pqrData.q' * phiReacStern;
	
end
%}
