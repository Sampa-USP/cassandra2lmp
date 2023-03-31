"""
Este código em Python tem como objetivo renomear os átomos de Oxigênio ligados a Hidrogênio em uma molécula. 
Ele faz uso da biblioteca Open Babel para ler o arquivo PDB contendo a estrutura da molécula e identificar os átomos que precisam ser renomeados. 
O usuário pode escolher entre renomear o rótulo do átomo ou o nome completo do átomo. O código também pode detectar átomos de Oxigênio ionizados na superfície do poro e renomeá-los de acordo.
A função retorna a carga final da molécula após a renomeação dos átomos de Oxigênio. O código é útil para a análise de porosidade em materiais e para a simulação de sistemas químicos.
"""

import os
import numpy as np
from openbabel import openbabel

def rename_oxygens_bonded_to_hydrogens(pdbfile,ofile, radius, detect_ionized_o=False, buffer_radius=2.0, rename_label=False):
    """
    Renames the Oxygens that are bonded to hydrogens.

    Parameters:
        pdbfile (str): Path of the PDB file containing the structure.
        radius (float): Pore radius in Angstroms to consider the surface.
        detect_ionized_o (bool): Flag to detect ionized oxygens on the surface of the pore. Default is False.
        buffer_radius (float): Buffer added to pore radius to find surface atoms. Default is 2.0.
        rename_label (bool): Flag to rename the label instead of atom name. Default is False.

    Returns:
        final_charge (float): The final charge of the molecule after renaming the oxygens.

    Author: Henrique Musseli Cezar
    Date: OCT/2021
    """

    charge_SiB = 1.1
    charge_OB = -0.55
    charge_SiI = 0.725
    charge_OI = -0.9
    charge_H = 0.4
    charge_OH = -0.675

    rbufsq = (radius + buffer_radius) * (radius + buffer_radius)

    # set openbabel file format
    base, ext = os.path.splitext(pdbfile)
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(ext[1:], "pdb")
    # trick to disable ring perception and make the ReadFile waaaay faster
    # Source: https://sourceforge.net/p/openbabel/mailman/openbabel-discuss/thread/56e1812d-396a-db7c-096d-d378a077853f%40ipcms.unistra.fr/#msg36225392
    obConversion.AddOption("s", openbabel.OBConversion.INOPTIONS)
    obConversion.AddOption("c", openbabel.OBConversion.INOPTIONS)  # ignore the CONECTs

    # read molecule to OBMol object
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, pdbfile)
    mol.SetFlag(openbabel.OB_PERIODIC_MOL)

    #mol.ConnectTheDots()  # necessary because of the 'b' INOPTION

    # for each bond, check if it is an hydrogen connected to an oxygen
    oxygensh = []

    counter = 0

    atypes = []

    sil = 0
    hyd = 0
    oxy = 0

    for atom in openbabel.OBMolAtomIter(mol):
        if atom.GetAtomicNum() == 14:
            sil += 1
        elif atom.GetAtomicNum() == 8:
            oxy += 1
        else:
            hyd += 1

    for bond in openbabel.OBMolBondIter(mol):
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        b1 = a1.GetAtomicNum()
        b2 = a2.GetAtomicNum()

        atypes.append(b1)
        atypes.append(b2)

        if ((not (b1 in [8, 14])) and (b2 == 8)):
            oxygensh.append(a2.GetId() + 1)
        elif (b1 == 8) and (not(b2 in [8,14])):
            oxygensh.append(a1.GetId()+1)

    #print(np.unique(atypes, return_counts=True))

    # for each oxygen, check if it's an ionized oxygen of a silanol at the pore surface
    ionoxygen = []
    ionsilicon = []
    if detect_ionized_o:
        for atom in openbabel.OBMolAtomIter(mol):
            if (atom.GetAtomicNum() == 8) and (atom.GetExplicitDegree() == 1):
                z = atom.GetZ()
                if ((z * z) <= rbufsq):
                    ionoxygen.append(atom.GetId() + 1)
                    for neigh in openbabel.OBAtomAtomIter(atom):
                        ionsilicon.append(neigh.GetId() + 1)

    #print("Found {} OH, {} OI and {} SiI atoms".format(len(oxygensh), len(ionoxygen), len(ionsilicon)))

    #print(f"Total number of atoms : {oxy + sil + hyd}")

    final_charge = len(oxygensh) * charge_OH + len(ionoxygen) * charge_OI + len(ionsilicon) * charge_SiI + (sil - len(
        ionsilicon)) * charge_SiB + (oxy - len(oxygensh) - len(ionoxygen)) * charge_OB + hyd * charge_H

    print(f"Final charge : {final_charge}")
    # now rename based on the list
    anum = 0
    fout = open(ofile, "w")
    with open(pdbfile, "r") as f:
        for line in f:
            if "HETATM" in line or "ATOM" in line:
                anum += 1
                if anum in oxygensh:
                    llist = list(line.rstrip())
                    if rename_label:
                        llist[12] = "O"
                        llist[13] = "H"
                        line = "".join(llist)
                    else:
                        line = "".join(llist) + " OH"
                elif anum in ionoxygen:
                    llist = list(line.rstrip())
                    if rename_label:
                        llist[12] = "O"
                        llist[13] = "I"
                        line = "".join(llist)
                    else:
                        line = "".join(llist) + " OI"
                elif anum in ionsilicon:
                    llist = list(line.rstrip())
                    if rename_label:
                        llist[12] = "S"
                        llist[13] = "i"
                        llist[14] = "I"
                        line = "".join(llist)
                    else:
                        line = "".join(llist) + " SiI"
                else:
                    if rename_label:
                        line = line.rstrip()
                    else:
                        line = line.rstrip() + " %s" % (line.split()[-1])
                fout.write(line + "\n")
            elif "CONECT" in line:
                pass
            else:
                fout.write(line)

    fout.close()
    return final_charge < 1e-8




