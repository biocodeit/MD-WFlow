# MD-WFlow

## Description
 This repo consist of scripts that will aid in MD (Molecular Dynamic ) Simulation preperation for Gromacs. A pre-requisite to run a simulations are that there are various files that need to be prepared. These scripts will help to prepare the systems with reduced effort, provided that input files and correct prompts are provided. 

 ## Features
 - Command line scripts
 - Flags systems e.g: python3 MDSysPrep.py --c complexFile.pdb
 - Resumable i.e will keep track of where the script last stopped and resume from there
 - Resume from any step by editing log file.

## Scripts

- **1. MDSysPrep.py**
  \<ins> </ins>
       <ins>Why:</ins>\
           Helps inpreparing the MD of Protein Ligand System in water solvent. The script is interactive and guide you through the steps.\
           For Amber itp:prm can be made from prmtop\
       <ins>Usage:</ins>\
           1. python3 MDSysPrep.py --c complex.pdb\
           2. Follow the prompts given.\
       <ins>Requirements:</ins>\
           1. Script:mdp files:Complex PDB:Ligand forcefield:Ligand itp:\
           2. The ligand files should begin from same name. eg. ATP.itp,ATP.ff,ATP.prm,ATP.pdb etc.\

- **2. MutantSystemPrep.py**
  \<ins> </ins>
       <ins>Why:</ins>\
           Creates sytems for mutants of already prepared system, .\
           Large number of Mutants can be creates simulatneously\
           Useful for Single, Double Mutant System Preperations.\
       <ins>Usage:</ins>\
           1. python3 MutatnSystemPrep.py --in MutantInput.in\
           2. Paths of all the things mentioned must be provided in Input file.\
       <ins>Requirements:</ins>\
           1. Script:Input.in:Folder with already Preped sytem
           2. Paths of files mentioned in input.in.\

  
