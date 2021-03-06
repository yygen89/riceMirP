The our training and testing datasets, and some datasets from other approaches are included in this folder.


1. train_positive.fa & test_positive.fa :  
Pre-miRNA sequences of rice were downloaded from the miRBase database (Release 22.1). 
After removing sequences containing non-AGCU characters, 604 known pre-miRNAs of rice were collected as positive dataset.
Then, 422 known pre-miRNAs are randomly selected as positive training dataset, and the other remaining 182 known pre-miRNAs 
are collected as independent positive testing dataset.


2. train_negative.fa & test_negative.fa :  
Pseudo pre-miRNAs we selected from cDNAs.  The CDS sequences of rice were obtained from 
PlantGDB database (http://www.plantgdb.org/download/Download/xGDB/OsGDB/Osativa_193_cds.fa.bz2), 
and then fragmented into non-overlapped segments under a constraint condition 
that the length distribution of extracted segments was identical with that of known plant pre-miRNAs. 
Further, they should satisfy some criteria. The criteria are determined by observing real plant pre-miRNAs. 
The criteria for selecting the pseud pre-miRNA are: minimum of 14 base pairings in the hairpins 
and maximum of −9.7 kcal/mol free energy of secondary structures (including GU wobble pairs).
Finally, 502 and 216 pseudo pre-miRNAs were randomly selected as negative training dataset 
and testing dataset, respectively.


3. datasets_from_miPlantPreMat:
The training and testing datasets of miPlantPreMat, which included in the software of miPlantPreMat, were downloaded from website: https://github.com/kobe-liudong/miPlantPreMat.


4. datasets_from_PlantMiRNAPred:
The training and testing datasets of PlantMiRNAPred were downloaded from web site of PlantMiRNAPred: http://nclab.hit.edu.cn/PlantMiRNAPred/.
