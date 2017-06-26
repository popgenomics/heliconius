
## Locus specific model comparison  
### statistics
**output_loci_pIso_pMig_v1.txt** contains:  
**geneName** = locus name (_scaffoldID_ _ _#gene_)  
**contig** = scaffold ID  
**start** = position of the first nucleotide on the scaffold  
**end** = position of the last nucleotide on the scaffold  
**bialsites_avg** = number of synonymous polymorphic positions in the interspecific alignement  
**sf_avg** = number of synonymous positions with a alleles differentially fixed between the two studied species  
**sxA_avg** = number of synonymous positions with a polymorphism exclusive to species A  
**sxB_avg** = number of synonymous positions with a polymorphism exclusive to species A  
**ss_avg** = number of synonymous positions with a shared polymorphism between A and B  
**Wald_avg** = spatial arrangement of fixed differences, and shared and exclusive polymorphisms (DOI: 10.1186/1471-2148-14-89)  
**piA_avg** = Tajima's estimator of \theta for species A  
**piB_avg** = Tajima's estimator of \theta for species B  
**thetaA_avg** = Watterson's estimator of \theta for species A  
**thetaB_avg** = Watterson's estimator of \theta for species B  
**DtajA_avg** = Tajima's _D_ for species A  
**DtajB_avg** = Tajima's _D_ for species B  
**divAB_avg** = raw divergence between A and B (_d<sub>xy</sub>_)  
**netdivAB_avg** = net divergence between A and B (_d<sub>a</sub>_)  
**minDivAB_avg** = minimum divergence between A and B among all pairs of individuals  
**maxDivAB_avg** = maximum divergence between A and B among all pairs of individuals  
**Gmin_avg** = minDivAB / divAB  
**Gmax_avg** = maxDivAB / divAB  
**FST_avg** = Fst computed by 1 - (**\pi A** + **\pi B**)/(2 * **\pi Tot**)  

### caracteristics of loci  
**nSamA** = number of gene copies for species A  
**nSamB** = number of gene copies for species B  
**Lsyno** = synonymous length  
**P_iso** = relative posterior probability of model with no migration  
**P_mig** = relative posterior probability of model with migration  

### compared models  
**M0** = model of strict-isolation between surveyed species.
**M1** = model of uni-lateral gene flow, only from species 2 (_rosina_) to species 1 (_chioneus_)
**M2** = model of uni-lateral gene flow, only from species 1 (_chioneus_) to species 1 (_rosina_)
**M3** = model of bi-lateral gene flow, at rates that can be different in both directions.


