mem_tr={'M02050[e] <=> M02050[c]'
		 'prgstrn[e] <=> prgstrn[c]'
		 'Lcystin[e] <=> Lcystin[c]'
		 'avite1[e] <=> avite1[c]'
		 'ca2[e] <=> ca2[c]'
		 'chol[e] <=> chol[c]'
		 'cortsn[e] <=> cortsn[c]'
		 'creat[e] <=> creat[c]'
		 'fe3[e] <=> fe3[c]'
		 'ile_L[e] <=> ile_L[c]'
		 'ncam[e] <=> ncam[c]'
		 'prostgf2[e] <=> prostgf2[c]'
		 'pydxn[e] <=> pydxn[c]'
		 'thm[e] <=> thm[c]'
         'retinol[e] <=> retinol[c]'
         'tststerone[e] <=> tststerone[c]'
		 'tyr_L[e] <=> tyr_L[c]'
		 'val_L[e] <=> val_L[c]'
		 'prostge2[e] <=> prostge2[c]'
		 'pyr[e] <=> pyr[c]'
		 'bilirub[e] <=> bilirub[c]'};
rxnNamesnew={'U1';'U2';'U3';'U4';'U5';'U6';'U7';'U8';'U9';'U10';'U11';'U12';'U13';'U14';'U15';'U16';'U17';'U18';'U19';'U20';'U21';'U22';'U23';'U24';'U25';'U26';'U27';'U28';'U29';'U30';'U31';'U32';'U33';'U34';'U35';'U36';'U37';'U38';'U39';'U40';'U41';'U42';'U43';'U44';'U45';'U46';'U47';'U48';'U49';'U50'};
memsubsytem={'mem_transports'};
%% add the transport reactions to the mem nutreints
for i = 1:21
model3 = addReaction(model3,rxnNamesnew{i,1},'reactionFormula',mem_tr{i,1},'subSystem',memsubsytem{1,1});
end
optimizeCbModel(model3)
