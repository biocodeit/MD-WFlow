import sys
sys.path.append("/home/spurge/scripts")

import logging, os, pathlib, shutil, argparse

import MDAnalysis as mda
import nglview as nv
import numpy as np
import parmed as pmd
import gromacs as gmx

gmx.config.get_configuration()

#loggging
requests_logger = logging.getLogger('INFO')
requests_logger.setLevel(logging.INFO)


#argument parsing
parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')

parser= argparse.ArgumentParser(prog='systemPreper', description='This will prepare system for you', )


def check_log(fileName, word):
    steps=[]
    with open(fileName, 'r') as readFile:
        for line in readFile.readlines():
            if line.startswith(word):
                step=line.split(':')[-1].strip('\n')
                steps.append(step)
            
    if len(steps)!=0:
        last_step=steps[-1]
    else:
        last_step=None
    return last_step

def ligFromSys(sys, ligand_resName):
    a=ligand_resName
    ats=sys.select_atoms('resname '+a)
    if len(ats.residues.resids) > 1:
        print('More than one residue exists for resname '+a)
        print(f'creating {len(ats.residues.resids)} pdbs')
        for resid in ats.residues.resids:
            ats.select_atoms('resid '+ str(resid)).write(a+str(resid)+'.pdb')
            print(f'made {a} into PDB')
    else:
        ligand=sys.select_atoms('resname '+a)
        if len(ligand)==0:
            raise ValueError('Atomgroup has 0 atoms')
        else:    
            ligand.write(a+'.pdb')
        print(f'made {a} into PDB')

def ligFromPDB(pdb, ligand_resName):
    sys=mda.Universe(pdb)
    a=ligand_resName
    ats=sys.select_atoms('resname '+a)
    if len(ats.residues.resids) > 1:
        print('More than one residue exists for resname '+a)
        print(f'creating {len(ats.residues.resids)} pdbs')
        for resid in ats.residues.resids:
            ats.select_atoms('resid '+ str(resid)).write(a+str(resid)+'.pdb')
            print(f'made {a} into PDB')
    else:
        ligand=sys.select_atoms('resname '+a)
        if len(ligand)==0:
            raise ValueError('Atomgroup has 0 atoms')
        else:    
            ligand.write(a+'.pdb')
        print(f'made {a} into PDB')

def choseLigPDB(**kwargs):
    try:
        kwargs['sys']
    except:
        print('sys not given checking for pdb')
    else:
        sys=kwargs['sys']
    try:
        kwargs['pdb']
    except:
        print('give parameters sys=, or pdb=')
    else:
        sys=mda.Universe(kwargs['pdb'])
    hetRes=np.unique(sys.select_atoms('not protein').residues.resnames)
    water_ions=['SOL', 'WAT', 'NA', 'SOD', 'CLA', 'HOH', 'CL']
    
    print('Following will be not be considered, Water_Ions:')
    print(water_ions)
    non_aa=input('any non-standard a.a involved:\n this will not be made into PDB ')
    water_ions.append(non_aa)
    for a in hetRes:
        if a not in water_ions:
            ligFromSys(sys=sys, ligand_resName=a)

def get_file_type(extension, exclude='####'):
    '''returns list of files of given extension'''
    files=[]
    for a in os.listdir():
        if a.endswith(extension) and not exclude in a:
            files.append(a)

    return files
    
def complexToPDB(sys, **kwargs):
    ##---- writing protein pdb
    try:
        kwargs['hetRes']
    except:
        residue=input('Any non-standard residue to be added to protein = ')
        kwargs['hetRes']=[a.strip(' ') for a in residue.split(',')]
    else:
        if not type(kwargs['hetRes'])==list:
             kwargs['hetRes']=kwargs['hetRes'].split(',')
        
    protein=sys.select_atoms('protein')        
    
    if residue !='n':
        for a in kwargs['hetRes']:
            ligand=sys.select_atoms('resname '+a)
            if len(ligand)==0:
                raise ValueError('Atom group Empty')
            protein=protein+ligand
    protein=protein.sort(key='ids')
    protein.write('protein.pdb')

def prmtopToGmxTop(fileName):
    x=pmd.load_file(f"{fileName}")
    x.save(f"{fileName[:-7]}.top",overwrite=True)
    print(f'{fileName} has been processed')

def topToPrm(fileName):
    with open(f'{fileName}', 'r') as readFile:
        content=readFile.read()
        sections=content.split('\n\n')
        with open(f'{fileName[:-4]}.prm', 'w') as writeFile:
            writeFile.write('[ defaults ]\n')
            writeFile.write(sections[2])
        with open(f'{fileName[:-4]}.itp', 'w') as writeFile:
            writeFile.write('\n\n'.join(sections[3:-2]))
        print(f'{fileName} is processed')
        
def ligAmberToGmx():
    is_parmeterized=input('Is the hetatms already parameterized (y/n)')
    if is_parmeterized == 'n':
        raise ValueError('Parameterize the Ligands')
    
    prmtops=get_file_type('.prmtop')
    
    top=get_file_type('.top') 

    itps=get_file_type('.itp')  
    
    if len(itps)>0:
        skip_lig_param=input('itps have been found , do we skip(y/n)')

    try:
        skip_lig_param
    except :
        skip_lig_param=None

    if skip_lig_param==None or skip_lig_param.lower()[0]=='n':
        if len(prmtops)+len(top)==0:
            raise ValueError('prmtop or top files not found!')
            
        # give option if top already present
        if len(top)> 1:
            proceed=input('prmtop files found. But .top files also present, do you want to proceed(y/n)')

        elif len(top)== 0:
            proceed='y'
            
        if proceed=='y' or len(top)== 0:                                    
            if len(prmtops) > 0:                                               ## if prmtop files are present
                print('Found following prmtops')
                print(prmtops)
                fileName=input('*'*5+'do we proceed(y/n)'+'*'*5)
                choices=['y','n']
                if fileName not in choices:
                    raise ValueError('Please put y/n')
                    
                while fileName != 'n':
                    fileName=input('which one to be used: (\'n\' to stop)')
                    if fileName != 'n':
                        prmtopToGmxTop(fileName)
        
                                                                            ## if top files are present
        top=get_file_type('.top')
                
        if len(top)> 1:
            proceed=input('The .top files are present, do you want to proceed(y/n)')
            if proceed=='y':   
                print('Found top Files. Are these for cof, lig:')
                print(top)
                fileName= input('*'*5+'do we proceed (y/n)'+'*'*5)
                choices=['y','n']
                if fileName not in choices:
                    raise ValueError('Please put y/n')
                while fileName != 'n':
                    fileName=input('which one to be used : (\'n\' for stop)' )
                    if fileName !='n':
                        topToPrm(fileName)

def add_to_topol(fileName, ff_extension, resName, mols, readFN='topol.top', 
                 writeFN='topol_complex.top'):
    fileName=fileName[:-4]
    print('Checking if itp and prm present')
    if fileName+'.itp' in os.listdir() and fileName+'.'+ff_extension in os.listdir():
        print(fileName+'.itp', fileName+ff_extension, 'present')
        with open(readFN, 'r') as readFile, open(writeFN, 'w') as writeFile:
            content=readFile.read()
            sections=content.split('\n\n')
            header=sections[0]
            ffBlock=sections[1]
            itpBlock=sections[2]
            ionWaterBlock='\n\n'.join(sections[3:6])
            sysName=sections[6]
            noOFmolecules=sections[7]
        
            ffBlock=ffBlock+f'\n#include "{fileName}.{ff_extension}"'
            itpBlock=itpBlock+f'\n#include "{fileName}.itp"'
            noOFmolecules=noOFmolecules+f'{resName}         {mols}\n'
        
            content=[header,ffBlock,itpBlock,ionWaterBlock,sysName,noOFmolecules]
            writeFile.write('\n\n'.join(content))
        print(fileName, 'processed')
        shutil.copy(writeFN, 'prev'+writeFN)
    elif fileName+'.itp' not in os.listdir() and fileName+ff_extension not in os.listdir():
        print(fileName+'.itp', fileName+ff_extension, 'not present')
        exit()

def update_topol():
    print('Will be writing the FolLowing into the topol')
    anyLigand=input('Are there any ligands to be included(y/n)')
    if anyLigand=='y':
        shutil.copy('topol.top', 'prevtopol_complex.top')      #backup the topol
        ff_extension=input('\n what is ff extension [prm, ff]')
        itps=get_file_type('itp', exclude='Protein')
        print(itps)
        
        fileName='y'
        while fileName != 'n':
            fileName=input('which of the follwoing to be included(\'n\' to stop)')
            
            if fileName != 'n':
                resName=input('Resname to be given')
                mols=input('No. of Molecules')
                
                add_to_topol(fileName, ff_extension, resName, mols, readFN='prevtopol_complex.top')

def make_complex_gmx_PDB(complexFile, proteinFile, ):
    with open(complexFile, 'w') as writeFile:                         # make atom only file
        with open(proteinFile, 'r') as readFile:
            for line in readFile.readlines():
                if line.startswith(('ATOM','TER')):
                    writeFile.write(line)
    
    pdbs=get_file_type('.pdb')
    print('The following will be added to protein, into complex')
    print('Please note the order should be same as appearing in the topol')
    
    fileName='y'
    print(pdbs)
    while fileName != 'n':
        fileName=input('Choose a pdb file(\'n\' to stop)')
        if fileName !='n':
            with open(complexFile, 'a') as writeFile:
                with open(f'{fileName}', 'r') as readFile:
                    for line in readFile.readlines():
                        if line.startswith(('ATOM','HETATM')):
                            writeFile.write(line)

def make_gmx_box(**kwargs):
    print('Complex will be boxed')

    try:
        kwargs['boxType']
    except:
        print('What box to be used. [triclinic, cubic, dodecahedron, octahedron]')
        kwargs['boxType']=input('Use default(dodecahedron) with \'y\', or specifiy')
        if kwargs['boxType']=='y':
            kwargs['boxType']='dodecahedron'
    try:
        kwargs['edgeDist']
    except:
        kwargs['edgeDist']=input('Mention the edge distance to be used. \'y\' for 1.0, or specify')
        if kwargs['edgeDist']=='y':
            kwargs['edgeDist']=1.0
        
    gmx.editconf(f='complex.pdb',
                o='newbox.gro',
                bt=kwargs['boxType'],
                c=True,
                d=kwargs['edgeDist'])

def gmxSolvate(**kwargs):
    options={'cp':'newbox.gro',
                 'p':'topol_complex.top',
                 'o':'solv.gro'}

    options.update(kwargs)
    
    gmx.solvate(cp=options['cp'],
           p=options['p'],
           o=options['o'],
            cs='spc216.gro',)

    
def gmxIonTpr(**kwargs):
    options={'f':'ions.mdp','c':'solv.gro','r':'solv.gro','p':'topol_complex.top','o':'ion.tpr'}

    options.update(kwargs)
    
    gmx.grompp(f=options['f'],
           c=options['c'],
           r=options['c'],
           p=options['p'],
           o=options['o'])

    
def gmxAdd_Ions(**kwargs):
    options={'s':'ion.tpr', 'p':'topol_complex.top','o':'solv_ions.gro', 'pname':'NA', 'nname':'CL',
          'neutral':True, 'input':['SOL'],'rmin':0.6}

    options.update(kwargs)

    gmx.genion(s=options['s'],
           p=options['p'],
           o=options['o'],
          pname=options['pname'],
          nname=options['nname'],
          neutral=options['neutral'],
          input=options['input'],
          rmin=options['rmin'])

def gmxIndexGroup(**kwargs):
    options={'f':'solv_ions.gro', 'o':'index.ndx'}

    options.update(kwargs)
    print('Coupling groups can be made')
    
    gmx.make_ndx(f= options['f'],
            o = options['o'],
            input=['q'])
    try:
        kwargs['group']
    except:
        index=gmx.fileformats.ndx.uniqueNDX(options['o'])
        
        groupName=input('Input group name: ')
        group_input=input('write input, (\'n\' to stop), (\'g\' another group)')
        index[groupName]={}
        
        while group_input !='n':
            if group_input=='g':
                groupName=input('Input group name: ')
                index[groupName]={}
                group_input=input('write input')
            index[groupName]=index[groupName]+index[group_input]
            group_input=input('write input')
        index.write(options['o'])
        gmx.check(n=options['o'])
    else:
        index=gmx.fileformats.ndx.uniqueNDX(options['o'])
        for key, value in kwargs['group'].items():
            index[key]={}
            for item in value:
                index[key]=index[key]+index[item]

        index.write(options['o'])
        gmx.check(n=options['o'])

def gmx_insert_mols(**kwargs):
    options={'f':'newbox.gro','o':'newbox_inserted.gro'}
    options.update(kwargs)

    gmx.insert_molecules(f=options['f'],
    o=options['o'],
    ci=options['ci'],
    nmol=options['nmol']
    )

def insert_mols_interactive(choice):
    if choice.lower()[0]!='n':
        ci=input('input the gro of mol to be added')
        nmol=input('how many molecules to be added')
        gmx_insert_mols(ci=ci, nmol=nmol)

        

if __name__=="__main__":
    ##---- READ THE LOG FILE TO GET THE LAST STEP
    if os.path.exists('MDSysPrep.log'):
        last_step= check_log('MDSysPrep.log', 'WARNING:MD')
        if last_step==None:
            last_step='S1.ProteinRead'
    
    if not os.path.exists('MDSysPrep.log'):
        last_step='S1.ProteinRead'
        
    logger = logging.getLogger('MD')
    logging.basicConfig(filename='MDSysPrep.log', encoding='utf-8', level=logging.WARNING, )
    
    #DIFINE the ARGUMENT complex
    parser.add_argument('-c', '--complex', type=str, help='pdb file complex', )
    args = parser.parse_args()
    
    print(f'\n {args.complex} will be used as complex')
    
    #LOADING OF THE SYSTEM
    print('\n loading the system...')
    sys=mda.Universe(args.complex)
    print(len(sys.atoms))
    
    
    if last_step=='S1.ProteinRead':
        ##---- identifying non proteins
        print('non standard residues =>', np.unique(sys.select_atoms('not protein').residues.resnames))
        
        ##---- writing protein pdb
        complexToPDB(sys)
        
        ##---- making hetatom pdbs
        choseLigPDB(sys=sys)
    
        ##any particular pdb to be made
        yesPDB=input('Do you have any resName pdb to be make(resname/(\'n\' to stop))')
        if not yesPDB=='n':
            ligFromSys(sys,yesPDB)
            
        logger.warning('S2.Protein_processed')
        last_step='S2.Protein_processed'
    
        
    ##---- PARAMETERIZATION OF SMALL MOLECULES
    if last_step=='S2.Protein_processed':
        ligAmberToGmx()
        
        logger.warning('S3.Lig_paramerterized')
        last_step='S3.Lig_paramerterized'
        
    ##---- PROTEIN PARAMETERIZATION
    if last_step=='S3.Lig_paramerterized':
        os.system('gmx pdb2gmx -f protein.pdb -o protein.gro -ignh -ter')              #amber forcefield is used
        
        gmx.editconf(f='protein.gro',
                o='protein_gmx.pdb')
        
        last_step='S4.Protein_parameterized'
        logger.warning('S4.Protein_parameterized')
    
    ##---- TOPOL.TOP ADDITION OF LIGANDS
    if last_step=='S4.Protein_parameterized':
        update_topol()
                
        last_step='S5.topol_modified'
        logger.warning('S5.topol_modified')
        
    
    ##---- MAKING COMPLEX PDB
    
    if last_step=='S5.topol_modified':
        make_complex_gmx_PDB('complex.pdb', 'protein_gmx.pdb')
    
        last_step='S6.complex_made'
        logger.warning('S6.complex_made')
    
    ##---- MAKING BOX
    
    if last_step=='S6.complex_made':
        make_gmx_box()
        
        last_step='S7.Box_made'
        logger.warning('S7.Box_made')
    
    if last_step=='S7.Box_made':
        choice=input('Is there molecules to be inserted')
        insert_mols_interactive(choice)
    
        last_step='S7.1.Inserted_mols'
        logger.warning('S7.1.Inserted_mols')
        if choice=='n':
            last_step=='S7.Box_made'
            logger.warning('S7.Box_made')


    ##---- SOLVATING
    if last_step=='Box_made':
        gmxSolvate()
    
        last_step='S8.Solvated'
        logger.warning('S8.Solvated')

    if last_step=='S7.1.Inserted_mols':
        gmxSolvate(cp='newbox_inserted.gro')
        last_step='S8.Solvated'
        logger.warning('S8.Solvated')

    
    ##---- ADDING IONS
    if last_step=='S8.Solvated':
        ionError=input('Are you getting any eroor in ion adding, use rmin?(y/n)')
        if ionError.lower()[0]=='y':
            rmin=input('rmin of your choice')
            gmxSolvate(cp='newbox_inserted.gro')
            gmxIonTpr()
            gmxAdd_Ions(rmin=rmin)
        if ionError.lower()[0]=='n':
            gmxIonTpr()
            gmxAdd_Ions()
        
        last_step='S9.Ions_added'
        logger.warning('S9.Ions_added')
    
    
    ##---- MAKE COUPLING GROUPS
    if last_step=='S9.Ions_added': 
        gmxIndexGroup()
    
        last_step='S10.Sys_prepared'
        logger.warning('S10.Sys_prepared')
    ##---- 
    
    if last_step=='S10.Sys_prepared':
        print('System is already Prepared')
    
