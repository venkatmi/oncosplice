# Oncosplice # 

Unsupervised Detection of Splicing Subtypes from RNA-Seq

![workflow](https://github.com/venkatmi/oncosplice/wiki/images/workflow.png)

Oncosplice is an automated pipeline to identify sample subtypes in an unsupervised manner from splicing quantification data. Oncosplice is currently runs on input files generated using MultiPath-PSI algorithm developed recently in AltAnalyze.

 # Dependencies # 

You will need to install the below packages (majorly python packages) before running OncoSplice.
  * Python 2.7
  * numpy
  * scipy
  * matplotlib
  * R 3.0 or greater with Hopach library included
  * nimfa
  * scikit-learn
  
Before running the OncoSplice workflow, you will need to install the appropriate species RNA-Seq database using the associated command-line:

```javascript
python AltAnalyze.py --species Hs --platform RNASeq --update Official --version EnsMart72
```

Replace the existing AltDatabase folder in oncosplice/scripts with the downloaded AltDatabase folder. This will update all the reference files required for running Oncosplice pipeline.

 # Steps involved # 
 
1. Generating the input files  

The input to this workflow can be obtained through AltAnalyze version 2.1.0 or greater (http://www.altanalyze.org) using the default analysis workflow described [here](http://altanalyze.readthedocs.io/en/latest/Algorithms/#multipath-psi-splicing-algorithm). If using the AltAnalyze GUI, users can enter the path where their RNA-Seq BAMs are stored and run the tool with default parameters.If using the command-line version AltAnalyze, see the installation instructions [here](https://github.com/nsalomonis/altanalyze/wiki/CommandLineMode).  
Generating bed files from BAM files (generated using STAR/Tophat)

- Junction Bed
```javascript
python BAMtoJunctionBED.py --i BAM_dir --species Hs --r software/AltAnalyze/AltDatabase/EnsMart72/ensembl/Hs/Hs_Ensembl_exon.txt
```

- Intron Bed
```javascript
python BAMtoExonBED.py --i BAM_dir --r exonrefdir/Hs.bed --s Hs
```
*the exon reference bed is generated in the Bam directory

Running AltAnalyze on Bed files (AltAnalyze generated/TCGA generated)
```javascript
python AltAnalyze.py --species Hs --platform RNASeq --bedDir bed_file_dir --output output_dir --groupdir /output_dir/ExpressionInput/groups_file.txt --compdir /output_dir/ExpressionInput/comps_file.txt --expname Exp_Name --runGOElite no
```
*For AltAnalyze, at least two groups need to specified (can be random - won’t impact OncoSplice results). Please read AltAnalyze manual to understand the naming and format of the groups and comps files [here](https://github.com/nsalomonis/altanalyze/wiki/ManualGroupsCompsCreation).


2. Running Oncosplice on AltAnalyze generated input PSI files  

The input PSI results file can be found in the folder *YourExpDirectory*/AltResults/AlternateOutput with the name “Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt” in the folder.

```javascript
python Oncosplice.py --EventAnnotation Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt"
```

*For the Event Annotation file provide the full path.

Please click [here](https://github.com/venkatmi/oncosplice/wiki/Demo-Example-with-scripts) for detailed description on how to run OncoSplice on the demo files. Please read the additional information section for detailed description of the pipeline.

 # Additional Information # 

1. Please click [here](https://github.com/venkatmi/oncosplice/wiki/Feature-Selection-and-Sample-Subtype-Classification) for more information on the computational pipeline and algorithms used in Oncosplice. 

2. Please click [here](https://github.com/venkatmi/oncosplice/wiki/Demo-Example-with-scripts) for detailed description on how to run the entire OncoSplice on the demo files.

3. Please click [here](https://github.com/venkatmi/oncosplice/wiki/Input-Output-Folder-Structure) for detailed information on the output files and folders generated from OncoSplice.

4. Please click [here](https://github.com/venkatmi/oncosplice/wiki/OncoSplice-Software-Functions-Description) for information on additional input parameters, python scripts and modules used by OncoSplice. 

See our [wiki](https://github.com/venkatmi/oncosplice/wiki) pages for additional details on the pipeline algorithms and available functions.


