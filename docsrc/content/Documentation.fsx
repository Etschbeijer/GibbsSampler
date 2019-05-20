(**
#Gibbs-Sampler
###motive search de novo or with prior knowledge


#####In this documentation the basic principles of the Gibbs-Sampling will be explained with a focus on its use in computational biology.

First the basic idea of the Gibbs-Sampling-Algorithm will be explained, together with the biological-background and goal of implementing this algorithm. 
After that, the working principle of the Gibbs-Algorithm will be explained with examples and then how to use the code. 
In the end is a short summarizing with the advanteages and disadvanteages of the algorithm.

Gibbs-Sampling is an algorithm based on the Markov-Chain-Monte-Carlo method, which enables the calculation of a joint distribution
solely based on conditional distributions without prior knowledge. In computational biology it is commonly used to find 
conserved segments, so called motives, that are part of different sequences. motives can be, for example, transcription-factor-binding-sites in DNA or
ligand-binding-sites in proteins. These binding-sites enable the binding of those elements, in order to regulate the activity of the target. 
Transcription-factor-binding-sites are often located in the promoter and 5'UTR of genes, which enable the regulation of the gene activity downstreams. 

My research-group researches the heath-shock-response of the alga Chlamydomonas reinhardtii by evaluating the proteome compositions of different mutants and 
targeting the corresponding genes. The goal of implementing this algorithm is to detect the transcription-factor-binding-sites of Heat-Shock-Factors (HSFs), proteins that
regulate the gene-activity during heat-stress. These binding-sites are called Heat-Shock-Elements (HSEs) and it is suspected that they are located in many genes. 
Sadly, they are poorly conserved compared to the HSEs of other organims and we hope that we can detect more HSE regulated genes by implementing a Gibbs-Sampling-Algorithm!


There are two distinct types of Gibbs-Sampling-Algorithms, the Site-Sampler and motive-Sampler. Both of them share the basic idea but while the Site-Sampler requires
a motive in every sequence to work properly, the motive-Sampler does not. The motive-Sampler can also find multiple motives per sequence, while the Site-Sampler can only find one
but works faster and is more sensitive for less conserved structures then the motive-Sampler. In order for both to wok proberly, the length of the sequence must be known.
The core idea of the Site-Sampler will now be explained with an example:

![Figure 1](img/sequences1.png)

Figure 1.                       

You have four sequences (see Figure 1) and you pick one at random and put it aside. 

![Figure 2](img/sequences2.png) 

Figure 2.

After that you choose one segment with the length of your motive in every sequence you
did not choose (see Figure 2). 

![Figure 3](img/pfm.png)

Figure 3.

![Figure 4](img/ppm.png)

Figure 4.

Align these segments with eacht other and count the amount of elements at each position in order to create a so called Position-Frequency-Matrix (PFM) (see Figure 3).
In the next step ypu create a so called position probability matrix by adding a pseudocount to each element at each position and normalize by dividing each element at 
each position through the summ of the column (see Figure 4). You add a pseudocount in order to avoid dividing through zero because often times the segments are rather short and do not 
contain every element and it is still possible to encounter one in the rest of the sequence.

![Figure 5](img/cfv.png)

Figure r.

Then you create a so called Frequency-Composite-Vector (FCV) by counting the amount of elements in the unchosen sequences without the elements that are part of the segments 
that are used to create the PFM.

*)
