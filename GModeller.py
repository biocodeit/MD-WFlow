def reorder_atoms(self):
    self.rename_segments(segment_ids=['A', 'B'])


query='WT'
template='WT_docked.pdb'
## can be only 'PDB' or 'SEQ'
queryType='SEQ'
Sequence='''AHHHHHHMVNMETNETFALLDATRAFLAKPKQMLIGAEWSDAASGRQLDVVNPADGTVIARVPEADERDVQQAVAAARRAFDAGPWRTAKTTDRERLMLVLADLIEANARELAEIESLDNGKPVMVAQGLDVAMAAQCFRYMAGWATKIEGSVIDAGMPYLPDSEIFAYTRKEPVGVVGAIIPWNFPLLMAAWKIAPALATGCTVVLKPAEDTPLSALRLGELIQAAGFPDGVVNIVTGYGHTAGAALSRDPRIDKIAFTGSTQTGKTIGHAALDNMTRMSLELGG'''


from modeller import *
# Get the sequence of tahe 1qg8 PDB file, and write to an alignment file
code = query
log.verbose()

env = Environ()
env.io.hetatm=True
aln = Alignment(env)

#if query is pdb then,
if queryType=='PDB':
    m = Model(env, file=code)
    aln.append_model(m, align_codes=code)
    
#if sequence is given then
elif queryType=='SEQ':
    aln.append_sequence(Sequence)
    aln[0].code=query
aln.write(file=code+'.ali')

from modeller import *

env = Environ()
env.io.hetatm=True
aln = Alignment(env)
# if 'LAST' not included you will get error
mdl = Model(env, file= template, model_segment=('FIRST:X','LAST:X'))
aln.append_model(mdl, align_codes=template, atom_files=f'{template}.pdb')
aln.append(file=f'{query}.ali', align_codes=query)
aln.align2d(max_gap_length=10)
aln.write(file=f'{query}AND{template}.ali', alignment_format='PIR')
aln.write(file=f'{query}AND{template}.pap', alignment_format='PAP')

from modeller import *
from modeller.automodel import *

env = Environ()
env.io.hetatm = True
a = AutoModel(env, alnfile=f'{query}AND{template}.ali',
        knowns=(template), sequence=query)
a.starting_model = 1
a.ending_model = 1
a.make()


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
            
            for residue in mdl1.residues:
                print(residue.num)    
                ?if res != mdl1.chains[chn.name].residues[s].name:
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
    
