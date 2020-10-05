Analysis Script is designed to run Single CLIP sample at a time from the CLIP_Pipeline.
The script will generate two Report files per sample:
	SampleName: Structure Analysis.html
	SampleName: LINE/SINE Distance.html

1. Run CLIP_Pipeline
2. Run Random CLIP data : RandomCLIP_V2.sh
    - Arg1 - Number of permutation
    - Arg2 - window size
    - Arg3 - minimum average chuck basepair probability
    - Arg4 - minimum Max chuck basepair probability
    - Arg5 - number of CPUs for individual tasks
    - Arg6 - batch number
    - Arg7 - CLIPtype
    - Arg8 - RunMode
    - Arg9 - TrasncLoc
    -
3. Place folder “CLIPoutput/CLIP_Structure_Analsysis” into Pipeline output Folder
4. Edit CLIP_Analysis.sh script for settings related to output
    - Note: Generating MotifAnalysis Figures and number of permutations for LINE/SINE distance will greatly increase run time
    - LINE/SINE Distance analysis should only have to be run once per sample as it generates all permutations in a single output




## Set Analysis mode  
**CLIPtype:**   
    Introns - Select peaks that are classified as "Protein Coding : Exons"   
    Exons - Select peaks that are classified as "Protein Coding : Introns"   
    Introns_LINESINETLR - Select peaks that are classified as Protein Coding : Introns and Intron peaks that are classified as "Repeat Elements" but overelap with LINE SINE or TLR   
       
**RunMode:**   
    Genomic - Include Intronic and Exonic Sequences when finding 100nt region Upsteram and Downstream of CLIP peak   
    Transcriptomic - Only include Exonic Sequences when finding 100nt region Upsteram and Downstream of CLIP peak. If Upstream/Downstream region crosses Intronic region, the intronic sequence is skpped and the next exon in included.   
       
**TranscLoc:** select subset of Protein Coding CLIP 
 
**window** Select the number of nts upstream and downstream of the 5’ CLIP peak in which you would like to interrogate.


## Nucleotide Content 
	Set the sliding window size and number of clusters for the Nucleotide content Heatmap

## Local Structure Settings
	**meanBPP** select the minimum, average base paring probability of a 10nt chuck in the Upstream region of the CLIP peak to cluster on the Heatmap
	**LScluster** Select number of clusters to Run