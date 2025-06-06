%% update the nomenclature of the metabolites
load('SPD39_modified.mat')
ids=findMetIDs(SPD39,{'glc-D[c]';'gln-L[c]';'trp-L[c]';'tyr-L[c]';'rib-D[c]';'cys-L[c]';'val-L[c]';'tyr-L[c]';'trp-L[c]';'ser-L[c]';'rib-D[c]'});
SPD39.mets(ids)={'glc_D[c]';'gln_L[c]';'trp_L[c]';'tyr_L[c]';'rib_D[c]';'cys_L[c]';'val_L[c]';'tyr_L[c]';'trp_L[c]';'ser_L[c]';'rib_D[c]'};
