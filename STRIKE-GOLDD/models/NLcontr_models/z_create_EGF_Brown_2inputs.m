% File that creates the model from Brown et al. 2004

clear all;

%--------------------------------------------------------------------------
% Model features common to all experiments:
%--------------------------------------------------------------------------



%% 48 PARAMETERS
syms global_par_krbEGF global_par_kruEGF global_par_krbNGF ...
    global_par_kruNGF global_par_kEGF global_par_KmEGF global_par_kNGF...
    global_par_KmNGF global_par_kdSos global_par_KmdSos global_par_kSos...
	global_par_KmSos global_par_kRasGap global_par_KmRasGap...
    global_par_kRasToRaf1 global_par_KmRasToRaf1 global_par_kpRaf1...
	global_par_KmpRaf1 global_par_kpBRaf global_par_KmpBRaf...
	global_par_kdMek global_par_KmdMek global_par_kpMekCytoplasmic...
	global_par_KmpMekCytoplasmic global_par_kdErk global_par_KmdErk...
	global_par_kpP90Rsk global_par_KmpP90Rsk global_par_kPI3K...
	global_par_KmPI3K global_par_kPI3KRas global_par_KmPI3KRas...
    global_par_kAkt global_par_KmAkt global_par_kdRaf1ByAkt...
	global_par_KmRaf1ByAkt global_par_kC3GNGF global_par_KmC3GNGF...
    global_par_kC3G global_par_KmC3G global_par_kRapGap...
    global_par_KmRapGap global_par_kRap1ToBRaf global_par_KmRap1ToBRaf...
	global_par_kdRaf1 global_par_KmdRaf1 global_par_kdBRaf...
	global_par_KmdBRaf
    
p = [global_par_krbEGF 
    global_par_kruEGF 
    global_par_krbNGF 
    global_par_kruNGF
    global_par_kEGF
    global_par_KmEGF
    global_par_kNGF
    global_par_KmNGF
    global_par_kdSos
    global_par_KmdSos
    global_par_kSos
    global_par_KmSos
    global_par_kRasGap
    global_par_KmRasGap
    global_par_kRasToRaf1
    global_par_KmRasToRaf1
    global_par_kpRaf1
    global_par_KmpRaf1
    global_par_kpBRaf
    global_par_KmpBRaf
    global_par_kdMek
    global_par_KmdMek
    global_par_kpMekCytoplasmic
    global_par_KmpMekCytoplasmic
    global_par_kdErk
    global_par_KmdErk
    global_par_kpP90Rsk
    global_par_KmpP90Rsk
    global_par_kPI3K
    global_par_KmPI3K
    global_par_kPI3KRas
    global_par_KmPI3KRas
    global_par_kAkt
    global_par_KmAkt
    global_par_kdRaf1ByAkt
    global_par_KmRaf1ByAkt
    global_par_kC3GNGF
    global_par_KmC3GNGF
    global_par_kC3G
    global_par_KmC3G
    global_par_kRapGap
    global_par_KmRapGap
    global_par_kRap1ToBRaf
    global_par_KmRap1ToBRaf
    global_par_kdRaf1
    global_par_KmdRaf1
    global_par_kdBRaf
    global_par_KmdBRaf
];


%% INPUT AND CONSTANTS
syms compartment_cell

compartment_cell=1.0;
% Species:   id = RasGapActive, name = RasGapActive, constant	
const_species_RasGapActive=120000.0;
% Species:   id = RapGapActive, name = RapGapActive, constant	
const_species_RapGapActive=120000.0;
% Species:   id = PP2AActive, name = PP2AActive, constant	
const_species_PP2AActive=120000.0;
% Species:   id = Raf1PPtase, name = Raf1PPtase, constant	
const_species_Raf1PPtase=120000.0;

k = [
    compartment_cell;
    const_species_RasGapActive;
    const_species_RapGapActive;
    const_species_PP2AActive;
    const_species_Raf1PPtase
    ];

%% 26 STATES, 2 INPUTS
syms x
x = sym(char(x),[28 1]);    

% Reaction: id = EGFBindingReaction, name = EGF binding
	reaction_EGFBindingReaction=compartment_cell*global_par_krbEGF*x(1)*x(3);
% Reaction: id = EGFUnbindingReaction, name = EFG unbinding
	reaction_EGFUnbindingReaction=compartment_cell*global_par_kruEGF*x(4);
% Reaction: id = NGFBindingReaction, name = NGF binding
	reaction_NGFBindingReaction=global_par_krbNGF*x(2)*x(5)*compartment_cell;
% Reaction: id = NGFUnbindingReaction, name = NGF unbinding
	reaction_NGFUnbindingReaction=global_par_kruNGF*x(6)*compartment_cell;
% Reaction: id = SosActivationByEGFReaction, name = SOS activation by EGF
	reaction_SosActivationByEGFReaction=compartment_cell*global_par_kEGF*x(4)*x(7)/(x(7)+global_par_KmEGF);
% Reaction: id = SosActivationByNGFReaction, name = SOS activation by NGF
	reaction_SosActivationByNGFReaction=compartment_cell*global_par_kNGF*x(6)*x(7)/(x(7)+global_par_KmNGF);
% Reaction: id = SosDeactivationReaction, name = SOS deactivation
	reaction_SosDeactivationReaction=compartment_cell*global_par_kdSos*x(10)*x(8)/(x(8)+global_par_KmdSos);
% Reaction: id = RasActivationReaction, name = Ras activation
	reaction_RasActivationReaction=compartment_cell*global_par_kSos*x(8)*x(11)/(x(11)+global_par_KmSos);
% Reaction: id = RasDeactivationReaction, name = Ras deactivation
	reaction_RasDeactivationReaction=compartment_cell*global_par_kRasGap*const_species_RasGapActive*x(12)/(x(12)+global_par_KmRasGap);
% Reaction: id = Raf1ByRasActivationReaction, name = Raf1 activation by Ras
	reaction_Raf1ByRasActivationReaction=compartment_cell*global_par_kRasToRaf1*x(12)*x(13)/(x(13)+global_par_KmRasToRaf1);
% Reaction: id = MekbyRaf1ActivationReaction, name = Mek activation by Raf1
	reaction_MekbyRaf1ActivationReaction=compartment_cell*global_par_kpRaf1*x(14)*x(17)/(x(17)+global_par_KmpRaf1);
% Reaction: id = MekbyBRafActivationReaction, name = Mek activation by B-Raf
	reaction_MekbyBRafActivationReaction=compartment_cell*global_par_kpBRaf*x(16)*x(17)/(x(17)+global_par_KmpBRaf);
% Reaction: id = ErkActivationReaction, name = Erk activation
	reaction_ErkActivationReaction=compartment_cell*global_par_kpMekCytoplasmic*x(18)*x(19)/(x(19)+global_par_KmpMekCytoplasmic);
% Reaction: id = MekDeactivationReaction, name = Mek deactivation
	reaction_MekDeactivationReaction=compartment_cell*global_par_kdMek*const_species_PP2AActive*x(18)/(x(18)+global_par_KmdMek);
% Reaction: id = ErkDeactivationReaction, name = Erk deactivation
	reaction_ErkDeactivationReaction=compartment_cell*global_par_kdErk*const_species_PP2AActive*x(20)/(x(20)+global_par_KmdErk);
% Reaction: id = Raf1byPPtaseDeactivationReaction, name = Raf1 deactivation by PPase
	reaction_Raf1byPPtaseDeactivationReaction=compartment_cell*global_par_kdRaf1*const_species_Raf1PPtase*x(14)/(x(14)+global_par_KmdRaf1);
% Reaction: id = P90RskActivationReaction, name = P90Rsk activation
	reaction_P90RskActivationReaction=compartment_cell*global_par_kpP90Rsk*x(20)*x(9)/(x(9)+global_par_KmpP90Rsk);
% Reaction: id = PI3KbyEGFRActivationReaction, name = PI3K activation by EGFR
	reaction_PI3KbyEGFRActivationReaction=compartment_cell*global_par_kPI3K*x(4)*x(21)/(x(21)+global_par_KmPI3K);
% Reaction: id = PI3KbyRasActivationReaction, name = PI3K activation by Ras
	reaction_PI3KbyRasActivationReaction=compartment_cell*global_par_kPI3KRas*x(12)*x(21)/(x(21)+global_par_KmPI3KRas);
% Reaction: id = AktActivationReaction, name = Akt activation
	reaction_AktActivationReaction=compartment_cell*global_par_kAkt*x(22)*x(23)/(x(23)+global_par_KmAkt);
% Reaction: id = Raf1ByAktDeactivationReaction, name = Raf1 deactivation by Akt
	reaction_Raf1ByAktDeactivationReaction=compartment_cell*global_par_kdRaf1ByAkt*x(24)*x(14)/(x(14)+global_par_KmRaf1ByAkt);
% Reaction: id = C3GActivationReaction, name = C3G activation
	reaction_C3GActivationReaction=compartment_cell*global_par_kC3GNGF*x(6)*x(25)/(x(25)+global_par_KmC3GNGF);
% Reaction: id = Rap1ActivationReaction, name = Rap1 activation
	reaction_Rap1ActivationReaction=compartment_cell*global_par_kC3G*x(26)*x(27)/(x(27)+global_par_KmC3G);
% Reaction: id = Rap1DeactivationReaction, name = Rap1 deactivation
	reaction_Rap1DeactivationReaction=compartment_cell*global_par_kRapGap*const_species_RapGapActive*x(28)/(x(28)+global_par_KmRapGap);
% Reaction: id = BRafByRap1ActivationReaction, name = BRaf activation by Rap1
	reaction_BRafByRap1ActivationReaction=compartment_cell*global_par_kRap1ToBRaf*x(28)*x(15)/(x(15)+global_par_KmRap1ToBRaf);
% Reaction: id = BRafbyPPtaseDeactivationReaction, name = BRaf deactivation by PPase
	reaction_BRafbyPPtaseDeactivationReaction=compartment_cell*global_par_kdBRaf*const_species_Raf1PPtase*x(16)/(x(16)+global_par_KmdBRaf);

xdot = sym(zeros(numel(x),1));
% Species:   id = EGF, name = EGF, affected by kineticLaw
xdot(1) = (1/(compartment_cell))*((-1.0 * reaction_EGFBindingReaction) + ( 1.0 * reaction_EGFUnbindingReaction));	
% Species:   id = NGF, name = NGF, affected by kineticLaw
xdot(2) = (1/(compartment_cell))*((-1.0 * reaction_NGFBindingReaction) + ( 1.0 * reaction_NGFUnbindingReaction));	
% Species:   id = freeEGFReceptor, name = freeEGFReceptor, affected by kineticLaw
xdot(3) = (1/(compartment_cell))*((-1.0 * reaction_EGFBindingReaction) + ( 1.0 * reaction_EGFUnbindingReaction));	
% Species:   id = boundEGFReceptor, name = boundEGFReceptor, affected by kineticLaw
xdot(4) = (1/(compartment_cell))*(( 1.0 * reaction_EGFBindingReaction) + (-1.0 * reaction_EGFUnbindingReaction));	
% Species:   id = freeNGFReceptor, name = freeNGFReceptor, affected by kineticLaw
xdot(5) = (1/(compartment_cell))*((-1.0 * reaction_NGFBindingReaction) + ( 1.0 * reaction_NGFUnbindingReaction));	
% Species:   id = boundNGFReceptor, name = boundNGFReceptor, affected by kineticLaw
xdot(6) = (1/(compartment_cell))*(( 1.0 * reaction_NGFBindingReaction) + (-1.0 * reaction_NGFUnbindingReaction));	
% Species:   id = SosInactive, name = SosInactive, affected by kineticLaw
xdot(7) = (1/(compartment_cell))*((-1.0 * reaction_SosActivationByEGFReaction) + (-1.0 * reaction_SosActivationByNGFReaction) + ( 1.0 * reaction_SosDeactivationReaction));	
% Species:   id = SosActive, name = SosActive, affected by kineticLaw
xdot(8) = (1/(compartment_cell))*(( 1.0 * reaction_SosActivationByEGFReaction) + ( 1.0 * reaction_SosActivationByNGFReaction) + (-1.0 * reaction_SosDeactivationReaction));	
% Species:   id = P90RskInactive, name = P90RskInactive, affected by kineticLaw
xdot(9) = (1/(compartment_cell))*((-1.0 * reaction_P90RskActivationReaction));	
% Species:   id = P90RskActive, name = P90RskActive, affected by kineticLaw
xdot(10) = (1/(compartment_cell))*(( 1.0 * reaction_P90RskActivationReaction));
% Species:   id = RasInactive, name = RasInactive, affected by kineticLaw
xdot(11) = (1/(compartment_cell))*((-1.0 * reaction_RasActivationReaction) + ( 1.0 * reaction_RasDeactivationReaction));
% Species:   id = RasActive, name = RasActive, affected by kineticLaw
xdot(12) = (1/(compartment_cell))*(( 1.0 * reaction_RasActivationReaction) + (-1.0 * reaction_RasDeactivationReaction));	
% Species:   id = Raf1Inactive, name = Raf1Inactive, affected by kineticLaw
xdot(13) = (1/(compartment_cell))*((-1.0 * reaction_Raf1ByRasActivationReaction) + ( 1.0 * reaction_Raf1byPPtaseDeactivationReaction) + ( 1.0 * reaction_Raf1ByAktDeactivationReaction));	
% Species:   id = Raf1Active, name = Raf1Active, affected by kineticLaw
xdot(14) = (1/(compartment_cell))*(( 1.0 * reaction_Raf1ByRasActivationReaction) + (-1.0 * reaction_Raf1byPPtaseDeactivationReaction) + (-1.0 * reaction_Raf1ByAktDeactivationReaction));	
% Species:   id = BRafInactive, name = BRafInactive, affected by kineticLaw
xdot(15) = (1/(compartment_cell))*((-1.0 * reaction_BRafByRap1ActivationReaction) + ( 1.0 * reaction_BRafbyPPtaseDeactivationReaction));	
% Species:   id = BRafActive, name = BRafActive, affected by kineticLaw
xdot(16) = (1/(compartment_cell))*(( 1.0 * reaction_BRafByRap1ActivationReaction) + (-1.0 * reaction_BRafbyPPtaseDeactivationReaction));	
% Species:   id = MekInactive, name = MekInactive, affected by kineticLaw
xdot(17) = (1/(compartment_cell))*((-1.0 * reaction_MekbyRaf1ActivationReaction) + (-1.0 * reaction_MekbyBRafActivationReaction) + ( 1.0 * reaction_MekDeactivationReaction));	
% Species:   id = MekActive, name = MekActive, affected by kineticLaw
xdot(18) = (1/(compartment_cell))*(( 1.0 * reaction_MekbyRaf1ActivationReaction) + ( 1.0 * reaction_MekbyBRafActivationReaction) + (-1.0 * reaction_MekDeactivationReaction));	
% Species:   id = ErkInactive, name = ErkInactive, affected by kineticLaw
xdot(19) = (1/(compartment_cell))*((-1.0 * reaction_ErkActivationReaction) + ( 1.0 * reaction_ErkDeactivationReaction));	
% Species:   id = ErkActive, name = ErkActive, affected by kineticLaw
xdot(20) = (1/(compartment_cell))*(( 1.0 * reaction_ErkActivationReaction) + (-1.0 * reaction_ErkDeactivationReaction));	
% Species:   id = PI3KInactive, name = PI3KInactive, affected by kineticLaw
xdot(21) = (1/(compartment_cell))*((-1.0 * reaction_PI3KbyEGFRActivationReaction) + (-1.0 * reaction_PI3KbyRasActivationReaction));	
% Species:   id = PI3KActive, name = PI3KActive, affected by kineticLaw
xdot(22) = (1/(compartment_cell))*(( 1.0 * reaction_PI3KbyEGFRActivationReaction) + ( 1.0 * reaction_PI3KbyRasActivationReaction));	
% Species:   id = AktInactive, name = AktInactive, affected by kineticLaw
xdot(23) = (1/(compartment_cell))*((-1.0 * reaction_AktActivationReaction));	
% Species:   id = AktActive, name = AktActive, affected by kineticLaw
xdot(24) = (1/(compartment_cell))*(( 1.0 * reaction_AktActivationReaction));	
% Species:   id = C3GInactive, name = C3GInactive, affected by kineticLaw
xdot(25) = (1/(compartment_cell))*((-1.0 * reaction_C3GActivationReaction));	
% Species:   id = C3GActive, name = C3GActive, affected by kineticLaw
xdot(26) = (1/(compartment_cell))*(( 1.0 * reaction_C3GActivationReaction));	
% Species:   id = Rap1Inactive, name = Rap1Inactive, affected by kineticLaw
xdot(27) = (1/(compartment_cell))*((-1.0 * reaction_Rap1ActivationReaction) + ( 1.0 * reaction_Rap1DeactivationReaction));
% Species:   id = Rap1Active, name = Rap1Active, affected by kineticLaw
xdot(28) = (1/(compartment_cell))*(( 1.0 * reaction_Rap1ActivationReaction) + (-1.0 * reaction_Rap1DeactivationReaction));
%%
x0=zeros(28,1);
x0(1) = 1.0002E7; % EGF stimulation
x0(2) = 456000.0; % NGF stimulation
x0(3) = 80000.0;
x0(4) = 0.0;
x0(5) = 10000.0;
x0(6) = 0.0;
x0(7) = 120000.0;
x0(8) = 0.0;
x0(9) = 120000.0;
x0(10) = 0.0;
x0(11) = 120000.0;
x0(12) = 0.0;
x0(13) = 120000.0;
x0(14) = 0.0;
x0(15) = 120000.0;
x0(16) = 0.0;
x0(17) = 600000.0;
x0(18) = 0.0;
x0(19) = 600000.0;
x0(20) = 0.0;
x0(21) = 120000.0;
x0(22) = 0.0;
x0(23) = 120000.0;
x0(24) = 0.0;
x0(25) = 120000.0;
x0(26) = 0.0;
x0(27) = 120000.0;
x0(28) = 0.0;    


%%
h   = [x(12),x(14),x(16),x(18),x(20),x(28)].';
f   = [];
ics = [];
for i=3:28
    f=[f;xdot(i)];
    ics=[ics,x0(i)];
end
u = [x(1);x(2)];	
x = [x(3);x(4);x(5);x(6);x(7);x(8);x(9);x(10);x(11);x(12);x(13);x(14);x(15);x(16);x(17);x(18);x(19);x(20);x(21);x(22);x(23);x(24);x(25);x(26);x(27);x(28)];
known_ics=zeros(26,1);

save('EGF_Brown_2inputs','x','h','u','p','f','ics','known_ics'); 
