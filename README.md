# Oncosplice # 

Unsupervised Detection of Splicing Subtypes from RNA-Seq

![workflow](https://github.com/venkatmi/oncosplice/wiki/images/workflow.png)

Oncosplice is an automated pipeline to identify sample subtypes in an unsupervised manner from splicing quantification data. Oncosplice is currently dependent on input files from AltAnalyze using the recently developed MultiPath-PSI algorithm. 

 # Dependencies # 

  * Python 2.7
  * numpy
  * scipy
  * matplotlib
  * R 3.0 or greater with Hopach library included
  * nimfa

 # How to run the workflow # 

```javascript
python MetaSpliceWorkflow.py --EventAnnotation Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt"
```

Before running this workflow, you will need to install the appropriate species RNA-Seq database using the associated command-line:

```javascript
python AltAnalyze.py --species Hs --platform RNASeq --update Official --version EnsMart72
```

 # Generating the input files #

The input to this workflow can be obtained through AltAnalyze version 2.1.0 or greater (http://www.altanalyze.org) using the default analysis workflow described [here](http://altanalyze.readthedocs.io/en/latest/Algorithms/#multipath-psi-splicing-algorithm). The option to run the Oncosplice pipeline directly from RNA-Seq aligned BAM files will be added in the release version of this pipeline. If using the AltAnalyze GUI, users can enter the path where their RNA-Seq BAMs are stored and run the tool with default parameters. If the using the command-line version AltAnalyze, see the BAM file analysis instructions [here](https://github.com/nsalomonis/altanalyze/wiki/CommandLineMode). For AltAnalyze, at least two groups need to specified (can be random - won’t impact OncoSplice results). The input PSI results file can be found in the folder *YourExpDirectory*/AltResults/AlternateOutput with the name “Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt” in the folder.

 # Additional Information # 

See our [wiki](https://github.com/venkatmi/oncosplice/wiki) pages for additional details on the pipeline algorithms and available functions.