driver_timer = tic();
addpath('../pointbem');
addpath('../panelbem');
loadConstants

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
disp("\n");
disp("start part 1 in driver");
bemPcm = makePanelBemEcfQualMatrices(srfData, pqrData, epsIn, epsOut);
disp("from driver, after makePanelBemEcfQualMatrices");
asymBemPcm = makePanelAsymEcfCollocMatrices(srfData, bemPcm, pqrData);
[phiReac, sigma] = solvePanelConsistentAsymmetric(srfData, bemPcm, ...
						        epsIn, epsOut, ...
						        conv_factor, ...
						        pqrData, ...
						        asymParams,asymBemPcm);
dG_nosalt_pcm = 0.5 * pqrData.q' * phiReac;

disp("part 1 in driver is done \n");

% part 2: no-salt reference calculation: YL/NLBC formulation

disp("start part 2 in driver");
bemYoonDiel = makePanelBemYoonDielMatrices(srfData,pqrData,epsIn,epsOut);
asymBemYL   = asymBemPcm; % we're reusing asymBemPcm
[phiReacYL, phiBndy, dphiDnBndy] = ...
    solvePanelConsistentYoonNoSternAsym(srfData, bemYoonDiel, ...
					epsIn, epsOut, ...
					conv_factor, pqrData, ...
					asymParams, asymBemYL);
dG_nosalt_yl = 0.5 * pqrData.q' * phiReacYL;
disp("part 2 in driver is done \n");

driver_timerVal = toc(driver_timer);
disp("driver_timerVal:"), disp(driver_timerVal);
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
