(**
#Gibbs-Sampler
###Motif search de novo or with prior knowledge


#####In this documentation the basic principles of the Gibbs-Sampling will be explained with a focus on its use in computational biology.

First the basic idea of the Gibbs-Sampling-Algorithm will be explained, together with the biological-background and goal of implementing this algorithm. 
After that, the working principle of the Gibbs-Algorithm will be explained with examples and then how to use the code. 
In the end is a short summarizing with the advanteages and disadvanteages of the algorithm.

Gibbs-Sampling is an algorithm based on the Markov-Chain-Monte-Carlo method, which enables the calculation of a joint distribution
solely based on conditional distributions without prior knowledge. In computational biology it is commonly used to find 
conserved segments, so called motifs, that are part of different sequences. Motifs can be, for example, transcription-factor-binding-sites in DNA or
ligand-binding-sites in proteins. These binding-sites enable the binding of those elements, in order to regulate the activity of the target. 
Transcription-factor-binding-sites are often located in the promoter and 5'UTR of genes, which enable the regulation of the gene activity downstreams. 

My research-group researches the heath-shock-response of the alga Chlamydomonas reinhardtii by evaluating the proteome compositions of different mutants and 
targeting the corresponding genes. The goal of implementing this algorithm is to detect the transcription-factor-binding-sites of Heat-Shock-Factors (HSFs), proteins that
regulate the gene-activity during heat-stress. These binding-sites are called Heat-Shock-Elements (HSEs) and it is suspected that they are located in many genes. 
Sadly, they are poorly conserved compared to the HSEs of other organims and we hope that we can detect more HSE regulated genes by implementing a Gibbs-Sampling-Algorithm!


*)
