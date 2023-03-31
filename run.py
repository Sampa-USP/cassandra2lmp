from utils.pdb2lmp import parse_mol_info
from utils.utils import write_last
from utils.rename_pdb import rename_oxygens_bonded_to_hydrogens

from ase.io import read

import glob
import os
import shutil 

c_files = glob.glob("./*.*")
paths = glob.glob("../cristobalite/*/gcmc_run*.xyz")
#paths = glob.glob("../../GCMC/cristobalite/*/gcmc_run*.xyz")


sio2_pdb_out = "sio2.pdb"
h2o_pdb_out = "h2o.pdb"
sio2_pdb_out_renamed = "sio2_tmp.pdb"

charges_file = './files/charges.txt'
pdb_file = "r_pore.pdb"

charges_sio2 = './files/charges_sio2.txt'
charges_h2o = './files/charges_h2o.txt'

axis = "z"
buffer = 2.0
buffer_length_orthogonal = 0.0
pbc_bonds = True
printdih = ""
ignore_bonds_solute = False
ignore_impropers = True
box_size = ""
supress_coeffs = ""
atom_labels_pdb = True
do_not_split_molecules = True
extra_impropers = ""

for gcmc_path in paths :

    print(gcmc_path)

    all_frames = read(gcmc_path, index=':')
    last_frame = all_frames[-1]

    lista = ['','','','','','']
    count = 0

    for indice,atom in enumerate(last_frame):
        _symbol = atom.symbol
        lista.append(_symbol)

        lista = lista[-4:]
        
        if lista == ['O','H','H','O']:
            break
    sio2 = last_frame[:indice - 3]
    h2o = last_frame[indice - 3:]

    h2o.write(h2o_pdb_out)
    sio2.write(sio2_pdb_out)

    rename_oxygens_bonded_to_hydrogens(sio2_pdb_out,sio2_pdb_out_renamed,7.5,False,2,False)

    lmpstring_sio2 = parse_mol_info(sio2_pdb_out_renamed, charges_sio2,axis,buffer,buffer_length_orthogonal,
                            pbc_bonds,printdih,ignore_bonds_solute,ignore_impropers,box_size,
                            supress_coeffs,atom_labels_pdb,do_not_split_molecules,extra_impropers)

    lmpstring_h2o = parse_mol_info(h2o_pdb_out, charges_h2o,axis,buffer,buffer_length_orthogonal,
                            pbc_bonds,printdih,ignore_bonds_solute,ignore_impropers,box_size,
                            supress_coeffs,atom_labels_pdb,do_not_split_molecules,extra_impropers)

    with open("sio2.lmp",'w') as file:
        file.write(lmpstring_sio2)

    with open("h2o.lmp",'w') as file:
        file.write(lmpstring_h2o)


    end_path = '/'.join(gcmc_path.split("/")[1:-1]).replace("GCMC","MD") + '/'
    os.makedirs(os.path.dirname(end_path), exist_ok=True)
    
    shutil.copy("sio2.lmp",end_path)
    shutil.copy("h2o.lmp",end_path)

    break



"""
for file in glob.glob("./*.*"):
    if file in c_files:
        continue
    else:
        os.remove(file)
"""