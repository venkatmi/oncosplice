# Oncosplice # 

Unsupervised Detection of Splicing Subtypes from RNA-Seq

![workflow](https://github.com/venkatmi/oncosplice/wiki/images/workflow.png)

Oncosplice is an automated pipeline to identify sample subtypes in an unsupervised manner from splicing quantification data. Oncosplice is currently dependent on input files from AltAnalyze using the recently develop MultiPath-PSI algorithm. 

 # Dependencies # 

  * Python 2.7
  * numpy
  * scipy
  * matplotlib
  * R 3.0 or greater

 # How to run the workflow # 

```javascript
python MetaSpliceCompleteWorkflow.py —-EventAnnotation Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt
```

 # Generating the input files #

The input to this workflow can be obtained through AltAnalyze (http://www.altanalyze.org). The option to run the Oncosplice pipeline directly from BAM files will be added in the release version of this pipeline. Using the GUI, users can enter the path where their RNA-Seq BAMs are stored and run the tool with default parameters. The input PSI results file can be found in the folder *YourExpDirectory*/AltResults/AlternateOutput with the name “Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt” in the folder.

 # Additional Details Coming Soon # 