# Host-pathogen-interaction-models
Codes used to build host-pathogen interaction models based on Recon 3D and Dual RNA-seq of pathogen and host. Two cases were considered. First case was pneumonia and the second case was gastric cancer. The RNA seq data of the infected host tissues belonging to the two cases were taken from the following sources 
1. Aprianto, Rieza, et al. "Time-resolved dual RNA-seq reveals extensive rewiring of lung epithelial and pneumococcal transcriptomes during early infection." _Genome biology_ 17 (2016): 1-16.
2. Sharafutdinov, Irshad, et al. "Early and late genome-wide gastric epithelial transcriptome response during infection with the human carcinogen Helicobacter pylori." _Cell Insight_ 1.3 (2022): 100032.

## Description of the codes in the following folders
1. Pathogen GEMs
   Two pathogen genome-scale metabolic models were used in this project. They are _Streptococcus pneumoniae_ D39 and _Helicobacter pylori_. The models were taken from the following sources.
   1. Pedram, Narges, Hamid Rashedi, and Ehsan Motamedian. "A systematic strategy using a reconstructed genome-scale             metabolic network for pathogen Streptococcus pneumoniae D39 to find novel potential drug targets." Pathogens and           Disease 78.6 (2020): ftaa051.
   2. Thiele, Ines, et al. "Expanded metabolic reconstruction of helicobacter pylori (i it341 gsm/gpr): an in silico             genome-scale characterization of single-and double-deletion mutants." Journal of bacteriology 187.16 (2005): 5818-         5830.

   These models were further modified by adding a few additional reactions and removing thermodynamically infeasible cycles. Further the nomenclature of the metabolites were changed according to the nomenclature of the Rexon 3D derived host GEMs. These models were constrained using the lunglike media for _Streptococcus pneumoniae_ D39 and a minimal media for _Helicobacter pylori_ to support their growth respectively. The nutrient media were taken from the sources of their respective GEMs.

3. Infected host model
     Following steps were taken to build the infected host models from the CCLE cancer cell line GEMs as both the dual RNA seq experiments were done on the cell lines. The detailed methodology common to both the scenarious is shown in the figure below.
   ![image-url](https://github.com/Subasree514/Host-pathogen-interaction-models/blob/main/hp.png)
4. HP model building
      Codes used to reconstruct the host-pathogen interaction models for pneumonia and gastric cancer.
