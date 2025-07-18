import sys
sys.path.append("/home/spurge/scripts/Gwork")
import logging, argparse
import sys, os, shutil, time

import gromacs as gmx
from modeller import *
from modeller.optimizers import MolecularDynamics, ConjugateGradients
from modeller.automodel import autosched
from pymsaviz import MsaViz, get_msa_testdata
import MDAnalysis as mda

from MD_WFlow import MDSysPrep

'''
#Input file for Mutant Prep

PATHS:
path_parent_folder== /home/spurge/work/sdproj/ongoing/MD/system_preping/MutantPrep/
template==/home/spurge/work/sdproj/ongoing/MD/system_preping/MutantPrep/reordered.pdb

REQUIREMENTS:
path_ligand_pdb== LIG.pdb
path_cofactor_pdb== COF.pdb
path_ligand_itp== LIG.itp
path_cofactor_itp== COF.itp
path_ligand_prm== LIG.prm
path_cofactor_prm== COF.prm
path_ligand_posre== posre_LIG.itp
path_cofactor_posre== posre_COF.itp
#path_forcefield_file==

MUTANTS:
#maintain the pattern of mutations
1A31==ASP 234 LYS, GLU 545 ARG
1A32==GLY 182 ARG, LEU 132 LYS, LEU 172 ARG, LEU 236 ASP


TOTOPOL:
#add the itps to be added in order
#mol==ffFile:itpFile:resName:mols
mol1==COF.prm:COF.itp:posre_COF.itp:COF:2
mol2==LIG.prm:LIG.itp:posre_LIG.itp:LIG:1

PDBTOCOMPLEX:
#put in order as above as topol
#mol1:mol2
COF.pdb:LIG.pdb

INDEXGROUPS:
#name the indexgroups
Complex==Protein:COF:LIG
Solvent==SOL:Ion
'''

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

def complexToPDB(sys, **kwargs):
    ##---- writing protein pdb
    if len(kwargs['hetRes'])==0:
        residue='n'
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

class make_mutant:
    def __init__(self, mutant, mutations, template):
        '''
        mutant='mutant name':str
        mutations='ASP 122 LYS, GLU 583 ARG' : str
        non-standard residues will not work
        '''
        #create list of mutations & remove spaces
        #seperate mutations and restyp
        self.mutant=mutant
        sep_res_num=[mutations.split(',')[i].strip().split(' ') for i in range(len(mutations.split(',')))]
        #get orginal residues present
        self.orgRes=[sep_res_num[i][0] for i in range(len(sep_res_num))]
        #get the renumber of mutations
        self.mutResNum=[sep_res_num[i][1] for i in range(len(sep_res_num))] 
        # get the res type of mutattion
        self.mutResTyp=[sep_res_num[i][2] for i in range(len(sep_res_num))] 
        self.template=template
        


    def run(self):
        #first argument
        #template, respos, restyp, chain, = 'WT_chimera','93','ARG','A'
        def optimize(atmsel, sched):
            #conjugate gradient
            for step in sched:
                step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
            #md
            refine(atmsel)
            cg = ConjugateGradients()
            cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
        
        
        #molecular dynamics
        def refine(atmsel):
            # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
            md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0,
                                   md_return='FINAL')
            init_vel = True
            for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                        (200, 600,
                                         (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
                for temp in temps:
                    md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                                 max_iterations=its, equilibrate=equil)
                    init_vel = False
        
        
        #use homologs and dihedral library for dihedral angle restraints
        def make_restraints(mdl1, aln):
           rsr = mdl1.restraints
           rsr.clear()
           s = Selection(mdl1)
           for typ in ('stereo', 'phi-psi_binormal'):
               rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
           for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
               rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                        spline_dx=0.3, spline_min_points = 5, aln=aln,
                        spline_on_site=True)
        
        
        log.verbose()
        
        # Set a different value for rand_seed to get a different final model
        env = Environ(rand_seed=-49837)
        
        env.io.hetatm = True
        #soft sphere potential
        env.edat.dynamic_sphere=False
        #lennard-jones potential (more accurate)
        env.edat.dynamic_lennard=True
        env.edat.contact_shell = 4.0
        env.edat.update_dynamic = 0.39
        
        # Read customized topology file with phosphoserines (or standard one)
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        
        # Read customized CHARMM parameter library with phosphoserines (or standard one)
        env.libs.parameters.read(file='$(LIB)/par.lib')
        
        
        # Read the original PDB file and copy its sequence to the alignment array:
        mdl1 = Model(env, file=self.template)
        ali = Alignment(env)
        ali.append_model(mdl1, atom_files=self.template, align_codes=self.template)

        #check for correct residue positions
        for res, s in zip(self.orgRes, self.mutResNum):
            for chn in mdl1.chains:
                if s in [residue.num for residue in mdl1.chains[chn.name].residues]:
                    if res != mdl1.chains[chn.name].residues[s].name:
                        raise ValueError(f"the orginal residues {chn.name}{s}{mdl1.chains[chn.name].residues[s].name} dont match given {res}, please jalview")
                # elif s in [residue.num for residue in mdl1.chains['B'].residues]:
                #     if res != mdl1.chains['B'].residues[s].name:
                #         raise ValueError(f"the orginal residues {s} dont match, please jalview")     
            
        #set up the mutate residue selection segment
        
        selections={}
        
        for s in self.mutResNum:
            selections[f's{s}']=Selection()
            for chn in mdl1.chains:
                if s in [residue.num for residue in mdl1.chains[chn.name].residues]:
                    selections[f's{s}'].add(mdl1.chains[chn.name].residues[s])
    
        
        #perform the mutate residue operation
        for mutRes,s in zip(self.mutResTyp,selections):
                    selections[s].mutate(residue_type=mutRes)
        #get two copies of the sequence.  A modeller trick to get things set up
        ali.append_model(mdl1, align_codes=self.template)
        
        # Generate molecular topology for mutant
        mdl1.clear_topology()
        mdl1.generate_topology(ali[-1])
        
        
        # Transfer all the coordinates you can from the template native structure
        # to the mutant (this works even if the order of atoms in the native PDB
        # file is not standard):
        #here we are generating the model by reading the template coordinates
        mdl1.transfer_xyz(ali)
        
        # Build the remaining unknown coordinates
        mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
        
        #yes model2 is the same file as model1.  It's a modeller trick.
        mdl2 = Model(env, file=self.template)
        
        #required to do a transfer_res_numb
        #ali.append_model(mdl2, atom_files=template, align_codes=template)
        #transfers from "model 2" to "model 1"
        mdl1.res_num_from(mdl2,ali)
        
        #It is usually necessary to write the mutated sequence out and read it in
        #before proceeding, because not all sequence related information about MODEL
        #is changed by this command (e.g., internal coordinates, charges, and atom
        #types and radii are not updated).
        
        mdl1.write(file=self.mutant+'.tmp')
        mdl1.read(file=self.mutant+'.tmp')
        
        #set up restraints before computing energy
        #we do this a second time because the model has been written out and read in,
        #clearing the previously set restraints
        make_restraints(mdl1, ali)
        
        #a non-bonded pair has to have at least as many selected atoms
        mdl1.env.edat.nonbonded_sel_atoms=1
        
        sched = autosched.loop.make_for_model(mdl1)
        
        #only optimize the selected residue (in first pass, just atoms in selected
        #residue, in second pass, include nonbonded neighboring atoms)
        #set up the mutate residue selection segment
        s=Selection()
        
        for res in self.mutResNum:
            for chn in mdl1.chains:
                if res in [residue.num for residue in mdl1.chains[chn.name].residues]:
                    s.add(mdl1.chains[chn.name].residues[res])
        
        
        mdl1.restraints.unpick_all()
        mdl1.restraints.pick(s)
        
        s.energy()
        
        s.randomize_xyz(deviation=4.0)
        
        mdl1.env.edat.nonbonded_sel_atoms=2
        optimize(s, sched)
        
        #feels environment (energy computed on pairs that have at least one member
        #in the selected)
        mdl1.env.edat.nonbonded_sel_atoms=1
        optimize(s, sched)
        
        s.energy()
        
        #give a proper name
        mutFileName=self.mutant+'_'.join(self.mutResTyp)+'-'+'_'.join(self.mutResNum)+'.pdb'
        mdl1.write(file=mutFileName)
        shutil.copy(mutFileName, 'modelled.pdb')
        #delete the temporary file
        os.remove(self.mutant+'.tmp') 

def setup_mutant_folders(dictionary, parent_directory, template):
    """creates folder of named with specified mutation supplied as dictionary, 
    puts a mutant pdb using modeller,
    copies necessaty files to this folder for system prep.
    requirements:
    1. dictionary => Eg. mutants['variant1']='ASP 122 LYS, GLU 583 ARG'
    2. parent directory ==> path of already system prepared folder
    3. necessary_files ==> a list of paths of all the files needed to be in newly created,
        variant folder, eg. ['./topol.top','COF.itp','charmm.ff','./toppar/posre.itp']
    """
    for mutant, mutation in dictionary.items():
        os.chdir(parent_directory)
        # parent_directory=os.getcwd()
        mutant_folder='../'+mutant
        destin=os.path.abspath(mutant_folder)
        os.makedirs(mutant_folder, exist_ok=True)
        os.chdir(mutant_folder)
        mutations = make_mutant(mutant, mutation, template)
        mutations.run()
        os.chdir(parent_directory)
    
        # for file in necessary_files:
        #     copy_selectFiles(file, parent_directory, mutant_folder)

def copy_selectFiles(file, source, dest):
    '''source us parent folder
    '''
    sourceFile=os.path.join(source, file)
    destFile=os.path.join(dest, file)

    if os.path.isdir(sourceFile):
        shutil.copytree(sourceFile, destFile)
    else:
        shutil.copy(sourceFile,destFile)

def reorder_atoms(self):
    self.rename_segments(segment_ids=['A', 'B'])

def add_to_topol(ff_File, itpFile, posreFile, resName, mols, readFN='topol.top', 
                 writeFN='topol_complex.top'):
    
    print('Checking if itp and ff present')
    if itpFile in os.listdir() and ff_File in os.listdir():
        print(itpFile, ff_File, 'present!!!')
        with open(readFN, 'r') as readFile:
            content=readFile.read()
        with open(writeFN, 'w') as writeFile:
            sections=content.split('\n\n')
            header=sections[0]
            ffBlock=sections[1]
            itpBlock=sections[2]
            ionWaterBlock='\n\n'.join(sections[3:6])
            sysName=sections[6]
            noOFmolecules=sections[7]
        
            ffBlock=ffBlock+f'\n#include "{ff_File}"'
            itpBlock=itpBlock+f'\n#include "{itpFile}"'
            itpBlock=itpBlock+f'\n#ifdef POSRES\n#include "{posreFile}"\n#endif'
            noOFmolecules=noOFmolecules+f'{resName}         {mols}\n'
        
            content=[header,ffBlock,itpBlock,ionWaterBlock,sysName,noOFmolecules]
            writeFile.write('\n\n'.join(content))
        print(itpFile, 'processed')
        shutil.copy(writeFN, 'prev'+writeFN)
    else:
        raise FileNotFoundError(ff_File, itpFile, 'not found')

def make_complex_gmxPDB(complexFile, proteinFile, hetRes):
    '''give list of pdbs as hetRes in order to added to protein'''
    with open(complexFile, 'w') as writeFile:
        with open(proteinFile, 'r') as readFile:
            for line in readFile.readlines():
                if line.startswith('ATOM'):
                    writeFile.write(line)
        for ligand in hetRes:
            with open(ligand, 'r') as readFile:
                for line in readFile.readlines():
                     if line.startswith(('ATOM','HETATM')):
                        writeFile.write(line)

def copy_mdps(source, dest):
    mdps=['ions.mdp','steep1.mdp', 'steep2.mdp', 'nvt.mdp', 'npt.mdp', 'md.mdp']
    for mdp in mdps:
        file=os.path.join(source,mdp)
        if not os.path.exists(file):
            raise ValueError(mdp,'does not exist')
        shutil.copy(file, dest)
        

#argument parsing
parser= argparse.ArgumentParser(prog='MutantSystemPreper',
                                description='This will prepare system for you')
parser.add_argument('--inp',type=str, help='put input.in file')
args=parser.parse_args()
args.inp

if __name__=="__main__":
    if os.path.exists('MutPrep.log'):
        last_step=check_log('MutPrep.log', 'WARNING:MtD')
        if last_step==None:
                last_step='S1.beginning'
    
    if not os.path.exists('MutPrep.log'):
        last_step='S1.beginning'
    
    logger = logging.getLogger('MtD')
    logging.basicConfig(filename='MutPrep.log', encoding='utf-8', level=logging.WARNING, )
    
    with open(args.inp , 'r') as readFile:
        paths, mutants, extra_files, ff_files, indexGroups, sysVal, inserts={}, {}, {}, {}, {}, {}, {}
        pdbsToTopol=[]
        trigger='noPrint'
        for line in readFile.readlines():
            if not line.startswith('#'):
                if line.startswith('REQUIREMENTS:'):
                    trigger='requirements'
                if line.startswith('MUTANTS:'):
                    trigger='mutants'
                if line.startswith('PATHS:'):
                    trigger='path'
                if line.startswith('TOTOPOL:'):
                    trigger='topol'
                if line.startswith('PDBTOCOMPLEX:'):
                    trigger='pdbtocomplex'
                if line.startswith('INDEXGROUPS:'):
                    trigger='indexgroups'
                if line.startswith('INSERTMOLECULE:'):
                    trigger='insertMols'
                if line.startswith('SYSTEMSETUP:'):
                    trigger='systemSetup'
                if line.startswith('\n'):
                    trigger='noPrint'
                if not trigger=='noPrint':
                    line=line.split('==')
                    if trigger=='path' and len(line)>1:
                        paths[line[0]]=line[1].strip('\n').strip(' ')
                    if trigger=='mutants' and len(line)>1:
                        mutants[line[0]]=line[1].strip('\n').strip(' ')
                    if trigger=='requirements' and len(line)>1:
                        extra_files[line[0]]=line[1].strip('\n').strip(' ')
                    if trigger=='topol' and len(line)>1:
                        ff_adds=[a.strip() for a in line[1].split(':')]
                        ff_files[line[0]]=ff_adds
                    if trigger=='pdbtocomplex' and len(line)>1:
                        for a in line[1].split(':'):
                            pdbsToTopol.append(a.strip())
                    if trigger=='indexgroups' and len(line)>1:
                        indexGroups[line[0]]=[a.strip() for a in line[1].split(':')] 
                    if trigger=='insertMols' and len(line)>1:
                        inserts[line[0]]=[a.strip() for a in line[1].split(':')]
                    if trigger=='systemSetup' and len(line)>1:
                        sysVal[line[0]]=line[1].strip()              
    need_files=[]
    for value in ff_files.values():
        for a in value[:-2]:
            need_files.append(a)
    for pdb in pdbsToTopol:
        need_files.append(pdb)
    for value in inserts.values():
        need_files.append(value[0])
        print(pdbsToTopol)  
    for value in extra_files.values():
        need_files.append(value)    
    if last_step=='S1.beginning':
        setup_mutant_folders(mutants, paths['path_parent_folder'], paths['template'] )
    
        last_step='S2.MutantsMade'
        logger.warning('S2.MutantsMade')
        
        print(100*'*')
        print(' '*40+'Mutants Made Succesfully')
        print(100*'*')
    
    if last_step=='S2.MutantsMade':
        necesarryFiles=[]
        for value in need_files:
            if (not value.isspace()) and value:
                necesarryFiles.append(value)
            else:
                raise ValueError(f'Not found {key}')
        
        for mutant in mutants.keys():
            dest=os.path.join(paths['path_parent_folder'],'../', mutant)
            for file in necesarryFiles:
                print(file)
                copy_selectFiles(file, paths['path_parent_folder'], dest)
        
        print(100*'*')
        print(' '*40+'Copied Files Succesfully')
        print(100*'*')
        last_step='S3.FilesCopied'
        logger.warning('S3.FilesCopied')
    
    if last_step=='S3.FilesCopied':
        for mutant in mutants.keys():
            destination=os.path.join(paths['path_parent_folder'],'../', mutant)
            os.chdir(destination)
            print('Now in', destination)
            sys=mda.Universe('modelled.pdb')
            complexToPDB(sys, hetRes=[])
        
            
            gmx.pdb2gmx(f='protein.pdb',
                  o='protein.gro',
                  ter=True,
                    ignh=True,
                   input=['1','1','0','0','0','0'])              #amber forcefield is used
            
            gmx.editconf(f='protein.gro',
                    o='protein_gmx.pdb')
            
            last_step='S4.Protein_parameterized'
            logger.warning('S4.Protein_parameterized')
    
    ##---- TOPOL.TOP ADDITION OF LIGANDS
    if last_step=='S4.Protein_parameterized':
        for mutant in mutants.keys():
            destination=os.path.join(paths['path_parent_folder'],'../', mutant)
            os.chdir(destination)
            print('Now in', destination)
            shutil.copy('topol.top', 'topol_complex.top')
            print('ok')
            time.sleep(3)

            for block in ff_files.values():
                add_to_topol(block[0], 
                        block[1],
                        block[2],
                        block[3],
                        block[4],
                            readFN='topol_complex.top', writeFN='topol_complex.top')
                    
        last_step='S5.topol_modified'
        logger.warning('S5.topol_modified')
        print(100*'*')
        print(' '*40+'Topol Modified Succesfully')
        print(100*'*')
        
    
    ##---- MAKING COMPLEX PDB
    
    if last_step=='S5.topol_modified':
        for mutant in mutants.keys():
            destination=os.path.join(paths['path_parent_folder'],'../', mutant)
            os.chdir(destination)
            print('Now in', destination)
            make_complex_gmxPDB('complex.pdb', 'protein_gmx.pdb', pdbsToTopol)
    
        last_step='S6.complex_made'
        logger.warning('S6.complex_made')
        print(100*'*')
        print(' '*40+'Complex Made Succesfully')
        print(100*'*')
    
    ##---- MAKING BOX
    
    if last_step=='S6.complex_made':
        for mutant in mutants.keys():
            destination=os.path.join(paths['path_parent_folder'],'../', mutant)
            os.chdir(destination)
            print('Now in', destination)
            make_gmx_box(boxType=sysVal['box'], edgeDist=float(sysVal['edgeDist']))
        
        last_step='S7.Box_made'
        logger.warning('S7.Box_made')
        print(100*'*')
        print(' '*40+'Complex Made Succesfully')
        print(100*'*')

    ##---- INSERT MOLECULES
    if last_step=='S7.Box_made':
        ##---- SOLVATING
        if len(inserts)==0:
            for mutant in mutants.keys():
                destination=os.path.join(paths['path_parent_folder'],'../', mutant)
                os.chdir(destination)
                print('Now in', destination)
                gmxSolvate()
        
            last_step='S8.Solvated'
            logger.warning('S8.Solvated')
            print(100*'*')
            print(' '*40+'Solvated Succesfully')
            print(100*'*')
        if len(inserts)>0:
            for mutant in mutants.keys():
                destination=os.path.join(paths['path_parent_folder'],'../', mutant)
                os.chdir(destination)
                print('Now in', destination)
                for key, value in inserts.items():
                    gmx_insert_mols(ci=value[0], nmol=value[1])
                    gmxSolvate(cp='newbox_inserted.gro')
        
            last_step='S8.Solvated'
            logger.warning('S8.Solvated')
            print(100*'*')
            print(' '*40+'Solvated Succesfully')
            print(100*'*')
            
    ##---- ADDING IONS
    if last_step=='S8.Solvated':
        # write_mdps()
        for mutant in mutants.keys():
            destination=os.path.join(paths['path_parent_folder'],'../', mutant)
            os.chdir(destination)
            print('Now in', destination)
            copy_mdps(paths['path_parent_folder'], destination)
            gmxIonTpr()
            gmxAdd_Ions()
        
        last_step='S9.Ions_added'
        logger.warning('S9.Ions_added')
        print(100*'*')
        print(' '*40+ 'Ions Added SuccessFully')
        print(100*'*')
    
    ##---- MAKE COUPLING GROUPS
    if last_step=='S9.Ions_added': 
        for mutant in mutants.keys():
            destination=os.path.join(paths['path_parent_folder'],'../', mutant)
            os.chdir(destination)
            print('Now in', destination)
            MDSysPrep.gmxIndexGroup(group=indexGroups)
    
        last_step='S10.Sys_prepared'
        logger.warning('S10.Sys_prepared')
        print(100*'*')
        print(' '*40+ 'System Prepared SuccessFully')
        print(100*'*')
    ##---- 
    
    if last_step=='S10.Sys_prepared':
        print('System is already Prepared')
    
