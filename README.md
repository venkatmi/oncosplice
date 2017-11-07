**Oncosplice**

Unsupervised Detection of Splicing Subtypes from RNA-Seq

![workflow](images/workflow.png)

Oncosplice is an automated pipeline to identify sample subtypes in an unsupervised manner for splicing data. It is mainly written in python and makes some calls to R scripts. It usually takes input files generated from the AltAnalyze splicing pipeline. 

***How to run the workflow***

```javascript
python CompleteWorkflow.py -—Inputfile psifile.txt —-EventAnnotation Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt
```

where psifile and eventannot are files generated from AltAnalyze, see below for description of these files.

How to generate the input files

AltAnalyze is a multi-functional and easy-to-use software package for automated single-cell and bulk gene and splicing analyses. Easy-to-use precompiled graphical user-interface versions available from our website (https://github.com/AltAnalyze/altanalyze). Using the GUI, users can enter the path where their RNA-Seq BAMs are stored and run the tool with default parameters. The tool will generate all the splicing results inside yourexperimentfolder/AltResults/AlternateOutput. The eventannot file is the Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt in the folder. The simplest way to generate the psifile is filter the eventannot file for the UID column and sample columns. The psifile contains samples as columns and events as rows where the first column is the eventid name provided in the UID column of the event annotation file.