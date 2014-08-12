
Enl_solv = Enl_twosweep-Eref_twosweep;
El_solv = El_twosweep-Eref_twosweep;

Enl_burial = squeeze(Enl_solv(1,:,:)-Enl_solv(2,:,:))
El_burial = squeeze(El_solv(1,:,:)-El_solv(2,:,:))
