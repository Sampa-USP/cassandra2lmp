"""
Receives a .pdb and returns the lammps data file for
a system with only atoms, bonds and angles.

Author: Henrique Musseli Cezar
Date: MAY/2020
"""

import sys
import argparse
import os

from openbabel import openbabel
openbabel.obErrorLog.StopLogging()

from openbabel import pybel
ob3 = True
import numpy as np
import math
import itertools

# from https://stackoverflow.com/a/11541495
def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def parse_mol_info(fname, fcharges, axis, buffa, buffo, pbcbonds, printdih, ignorebonds, ignoreimproper, boxl, suppcoeff, labelpdb, notsplit, extraimp):
  iaxis = {"x": 0, "y": 1, "z": 2}
  if axis in iaxis:
    repaxis = iaxis[axis]
  else:
    print("Error: invalid axis")
    sys.exit(0)

  if fcharges:
    chargesLabel = {}
    with open(fcharges, "r") as f:
      for line in f:
        chargesLabel[line.split()[0]] = float(line.split()[1])

# set openbabel file format
  base, ext = os.path.splitext(fname)
  obConversion = openbabel.OBConversion()
  obConversion.SetInAndOutFormats(ext[1:], "pdb")
  # trick to disable ring perception and make the ReadFile waaaay faster
  # Source: https://sourceforge.net/p/openbabel/mailman/openbabel-discuss/thread/56e1812d-396a-db7c-096d-d378a077853f%40ipcms.unistra.fr/#msg36225392
  obConversion.AddOption("s", openbabel.OBConversion.INOPTIONS)
  obConversion.AddOption("c", openbabel.OBConversion.INOPTIONS)  # ignore the CONECTs

  # read molecule to OBMol object
  mol = openbabel.OBMol()
  obConversion.ReadFile(mol, fname)
  mol.SetFlag(openbabel.OB_PERIODIC_MOL)


  # split the molecules
  if notsplit:
    molecules = [mol]
  else:
    molecules = mol.Separate()

  # detect the molecules types
  mTypes = {}
  mapmTypes = {}
  atomIdToMol = {}
  nty = 0
  for i, submol in enumerate(molecules, start=1):
    atomiter = openbabel.OBMolAtomIter(submol)
    atlist = []
    for at in atomiter:
      atlist.append(at.GetAtomicNum())
      atomIdToMol[at.GetId()] = i
    foundType = None

    for ty in mTypes:
      # check if there's already a molecule of this type
      if atlist == mTypes[ty]:
        foundType = ty

    # if not, create a new type
    if not foundType:
      nty += 1
      foundType = nty
      mTypes[nty] = atlist

    mapmTypes[i] = foundType

  # get atomic labels from pdb
  idToAtomicLabel = {}
  if labelpdb:
    for res in openbabel.OBResidueIter(mol):
      for atom in openbabel.OBResidueAtomIter(res):
        idToAtomicLabel[atom.GetId()] = res.GetAtomID(atom).strip() 
  elif ext[1:] == "pdb":
    for res in openbabel.OBResidueIter(mol):
      for atom in openbabel.OBResidueAtomIter(res):
        if (atomIdToMol[atom.GetId()] > 1) and (len(mTypes) > 1):
          idToAtomicLabel[atom.GetId()] = res.GetAtomID(atom).strip()+str(mapmTypes[atomIdToMol[atom.GetId()]])
        else:
          idToAtomicLabel[atom.GetId()] = res.GetAtomID(atom).strip()
  else:
    if not ob3:
      etab = openbabel.OBElementTable()
    for atom in openbabel.OBMolAtomIter(mol):
      if (atomIdToMol[atom.GetId()] > 1) and (len(mTypes) > 1):
        if ob3:
          idToAtomicLabel[atom.GetId()] = openbabel.GetSymbol(atom.GetAtomicNum())+str(mapmTypes[atomIdToMol[atom.GetId()]])
        else:
          idToAtomicLabel[atom.GetId()] = etab.GetSymbol(atom.GetAtomicNum())+str(mapmTypes[atomIdToMol[atom.GetId()]])
      else:
        if ob3:
          idToAtomicLabel[atom.GetId()] = openbabel.GetSymbol(atom.GetAtomicNum())
        else:
          idToAtomicLabel[atom.GetId()] = etab.GetSymbol(atom.GetAtomicNum())

  # print(idToAtomicLabel)

  # identify atom types and get masses
  outMasses = "Masses\n\n"

  massTypes = {}
  mapTypes = {}
  nmassTypes = 0
  atomIterator = openbabel.OBMolAtomIter(mol)
  for atom in atomIterator:
    i = atom.GetId()
    if idToAtomicLabel[i] not in massTypes:
      nmassTypes += 1
      mapTypes[nmassTypes] = idToAtomicLabel[i]
      massTypes[idToAtomicLabel[i]] = nmassTypes
      outMasses += "\t%d\t%.3f\t# %s\n" % (nmassTypes, atom.GetAtomicMass(), idToAtomicLabel[i])

  # create atoms list
  outAtoms = "Atoms # full\n\n"

  xmin = float("inf")
  xmax = float("-inf")
  ymin = float("inf")
  ymax = float("-inf")
  zmin = float("inf")
  zmax = float("-inf")
  natoms = 0
  for mnum, imol in enumerate(molecules, start=1):
    atomIterator = openbabel.OBMolAtomIter(imol)
    for atom in sorted(atomIterator, key=lambda x: x.GetId()):
      natoms += 1
      i = atom.GetId()
      apos = (atom.GetX(), atom.GetY(), atom.GetZ())

      # look for the maximum and minimum x for the box (improve later with numpy and all coordinates)
      if apos[0] > xmax:
        xmax = apos[0]
      if apos[0] < xmin:
        xmin = apos[0]
      if apos[1] > ymax:
        ymax = apos[1]
      if apos[1] < ymin:
        ymin = apos[1]
      if apos[2] > zmax:
        zmax = apos[2]
      if apos[2] < zmin:
        zmin = apos[2]

      if fcharges:
        outAtoms += "\t%d\t%d\t%d\t%.6f\t%.4f\t%.4f\t%.4f\t# %s\n" % (i+1, mnum, massTypes[idToAtomicLabel[i]], chargesLabel[idToAtomicLabel[i].upper()], atom.GetX(), atom.GetY(), atom.GetZ(), idToAtomicLabel[i])
      else:
        outAtoms += "\t%d\t%d\t%d\t0.00000\t%.4f\t%.4f\t%.4f\t# %s\n" % (i+1, mnum, massTypes[idToAtomicLabel[i]], atom.GetX(), atom.GetY(), atom.GetZ(), idToAtomicLabel[i])

  # define box shape and size
  if boxl:
    fromBounds = False
    v1 = [boxl[0], 0., 0.]
    v2 = [0., boxl[1], 0.]
    v3 = [0., 0., boxl[2]]
    orthogonal = True
  else:
    try:
      fromBounds = False
      rcell = mol.GetData(12)
      cell = openbabel.toUnitCell(rcell)
      v1 = [cell.GetCellVectors()[0].GetX(), cell.GetCellVectors()[0].GetY(), cell.GetCellVectors()[0].GetZ()]
      v2 = [cell.GetCellVectors()[1].GetX(), cell.GetCellVectors()[1].GetY(), cell.GetCellVectors()[1].GetZ()]
      v3 = [cell.GetCellVectors()[2].GetX(), cell.GetCellVectors()[2].GetY(), cell.GetCellVectors()[2].GetZ()]
      boxinfo = [v1,v2,v3]
      orthogonal = True
      for i, array in enumerate(boxinfo):
        for j in range(3):
          if i == j:
            continue
          if not math.isclose(0., array[j], abs_tol=1e-6):
            orthogonal = False
    except:
      fromBounds = True
      v1 = [xmax - xmin, 0., 0.]
      v2 = [0., ymax - ymin, 0.]
      v3 = [0., 0., zmax - zmin]
      orthogonal = True

  # add buffer
  if orthogonal:
    buf = []
    boxinfo = [v1,v2,v3]
    for i, val in enumerate(boxinfo[repaxis]):
      if i == repaxis:
        buf.append(val+buffa)
      else:
        buf.append(val)
    boxinfo[repaxis] = buf
    for i in range(3):
      if i == repaxis:
        continue
      buf = []
      for j, val in enumerate(boxinfo[i]):
        if j == i:
          buf.append(val+buffo)
        else:
          buf.append(val)
      boxinfo[i] = buf
    
    # set mol unitcell
    cell = openbabel.OBUnitCell()
    cell.SetData(boxinfo[0][0], boxinfo[1][1], boxinfo[2][2], 90.0, 90.0, 90.0)
    # remove old info
    if mol.HasData(12):
      mol.DeleteData(12)
    for bond in openbabel.OBMolBondIter(mol):
      mol.DeleteBond(bond)
    # set new cell and create the bonds
    mol.CloneData(cell)
    mol.ConnectTheDots()

  # print(boxinfo)

  # identify bond types and create bond list
  outBonds = "Bonds # harmonic\n\n"

  bondTypes = {}
  mapbTypes = {}
  nbondTypes = 0
  nbonds = 0
  bondIterators = []
  if ignorebonds:
    sepmols = mol.Separate()
    for smol in sepmols[1:]:
      bondIterators.append(openbabel.OBMolBondIter(smol))
  else:
    bondIterators.append(openbabel.OBMolBondIter(mol))

  lastidx = 1
  for iterator in bondIterators:
    for i, bond in enumerate(iterator, lastidx):
      b1 = bond.GetBeginAtom().GetId()    
      b2 = bond.GetEndAtom().GetId()
      if bond.GetBeginAtom().GetAtomicNum() == bond.GetEndAtom().GetAtomicNum():
          continue
      # identify bond type
      btype1 = "%s - %s" % (idToAtomicLabel[b1],idToAtomicLabel[b2])
      btype2 = "%s - %s" % (idToAtomicLabel[b2],idToAtomicLabel[b1])

      if btype1 in bondTypes:
        bondid = bondTypes[btype1]
        bstring = btype1
      elif btype2 in bondTypes:
        bondid = bondTypes[btype2]
        bstring = btype2
      else:
        nbondTypes += 1
        mapbTypes[nbondTypes] = btype1
        bondid = nbondTypes
        bondTypes[btype1] = nbondTypes
        bstring = btype1

      nbonds += 1
      outBonds += "\t%d\t%d\t%d\t%d\t# %s\n" % (nbonds, bondid, b1+1, b2+1, bstring)

    lastidx = i

  # identify angle types and create angle list
  angleTypes = {}
  mapaTypes = {}
  nangleTypes = 0
  nangles = 0
  angleIterators = []

  if ignorebonds:
    sepmols = mol.Separate()
    for smol in sepmols[1:]:
      smol.FindAngles()
      angleIterators.append(openbabel.OBMolAngleIter(smol))
    prevnumatoms = sepmols[0].NumAtoms()
  else:
    mol.FindAngles()
    angleIterators.append(openbabel.OBMolAngleIter(mol))
  
  outAngles = "Angles # harmonic\n\n"

  lastidx = 1
  for j, iterator in enumerate(angleIterators, 1):
    for i, angle in enumerate(iterator, lastidx):
      if ignorebonds:
        a1 = angle[1] + prevnumatoms
        a2 = angle[0] + prevnumatoms
        a3 = angle[2] + prevnumatoms
      else:
        a1 = angle[1]
        a2 = angle[0]
        a3 = angle[2]
      

      if (idToAtomicLabel[a1]==idToAtomicLabel[a2]) or (idToAtomicLabel[a2]==idToAtomicLabel[a3]):
        continue
    
      atype1 = "%s - %s - %s" % (idToAtomicLabel[a1],idToAtomicLabel[a2],idToAtomicLabel[a3])
      atype2 = "%s - %s - %s" % (idToAtomicLabel[a3],idToAtomicLabel[a2],idToAtomicLabel[a1])

      if atype1 in angleTypes:
        angleid = angleTypes[atype1]
        astring = atype1
      elif atype2 in angleTypes:
        angleid = angleTypes[atype2]
        astring = atype2
      else:
        nangleTypes += 1
        mapaTypes[nangleTypes] = atype1
        angleid = nangleTypes
        angleTypes[atype1] = nangleTypes
        astring = atype1

      nangles += 1
      outAngles += "\t%d\t%d\t%d\t%d\t%d\t# %s\n" % (nangles, angleid, a1+1, a2+1, a3+1, astring)

    lastidx = i
    if ignorebonds:
      prevnumatoms += sepmols[j].NumAtoms()

  # identify dihedral types and create dihedral list
  if printdih:
    dihedralTypes = {}
    mapdTypes = {}
    ndihedralTypes = 0
    ndihedrals = 0
    dihedralIterators = []

    if ignorebonds:
      sepmols = mol.Separate()
      for smol in sepmols[1:]:
        smol.FindTorsions()
        dihedralIterators.append(openbabel.OBMolTorsionIter(smol))
    else:
      mol.FindTorsions()
      dihedralIterators.append(openbabel.OBMolTorsionIter(mol))

    outDihedrals = "Dihedrals # charmmfsw\n\n"

    lastidx = 1
    for iterator in dihedralIterators:
      for i, dihedral in enumerate(iterator, lastidx):
        a1 = dihedral[0]
        a2 = dihedral[1]
        a3 = dihedral[2]
        a4 = dihedral[3]

        dtype1 = "%s - %s - %s - %s" % (idToAtomicLabel[a1],idToAtomicLabel[a2],idToAtomicLabel[a3],idToAtomicLabel[a4])
        dtype2 = "%s - %s - %s - %s" % (idToAtomicLabel[a4],idToAtomicLabel[a3],idToAtomicLabel[a2],idToAtomicLabel[a1])

        if dtype1 in dihedralTypes:
          dihedralid = dihedralTypes[dtype1]
          dstring = dtype1
        elif dtype2 in dihedralTypes:
          dihedralid = dihedralTypes[dtype2]
          dstring = dtype2
        else:
          ndihedralTypes += 1
          mapdTypes[ndihedralTypes] = dtype1
          dihedralid = ndihedralTypes
          dihedralTypes[dtype1] = ndihedralTypes
          dstring = dtype1

        ndihedrals += 1
        outDihedrals += "\t%d\t%d\t%d\t%d\t%d\t%d\t# %s\n" % (ndihedrals, dihedralid, a1+1, a2+1, a3+1, a4+1, dstring)

      lastidx = i

    if not ignoreimproper:
      # look for the improper dihedrals
      improperDihedralTypes = {}
      mapiDTypes = {}
      niDihedralTypes = 0
      niDihedrals = 0
      mollist = []

      if ignorebonds:
        sepmols = mol.Separate()
        for smol in sepmols[1:]:
          smol.PerceiveBondOrders()
          mollist.append(smol)
      else:
        mol.PerceiveBondOrders()
        mollist.append(mol)

      outImpropers = "Impropers # harmonic\n\n"

      for imol in mollist:
        atomIterator = openbabel.OBMolAtomIter(imol)
        for atom in atomIterator:
          try:
            # print(atom.GetHyb(), atom.GetAtomicNum(), atom.GetValence())
            expDegree = atom.GetValence()
          except:
            # print(atom.GetId()+1, atom.GetHyb(), atom.GetAtomicNum(), atom.GetExplicitDegree(), atom.GetExplicitValence(), atom.GetHvyDegree(), atom.GetHeteroDegree())
            expDegree = atom.GetExplicitDegree()

          # returns impropers for atoms with connected to other 3 atoms and SP2 hybridization
          if atom.GetHyb() == 2 and expDegree == 3:
            impAtom = True
          elif extraimp and (atom.GetAtomicNum() == 6 and expDegree == 3):
            impAtom = True
          else:
            impAtom = False

          if impAtom:
            connectedAtoms = []
            for atom2, depth in openbabel.OBMolAtomBFSIter(imol, atom.GetId()+1):
              if depth == 2:
                connectedAtoms.append(atom2)

            torsional = [atom.GetId()+1, connectedAtoms[0].GetId()+1, connectedAtoms[1].GetId()+1, connectedAtoms[2].GetId()+1]

            a1 = torsional[0]-1
            a2 = torsional[1]-1
            a3 = torsional[2]-1
            a4 = torsional[3]-1

            alldtypes = []
            for perm in itertools.permutations([idToAtomicLabel[a2], idToAtomicLabel[a3], idToAtomicLabel[a4]]):
              alldtypes.append("%s - %s - %s - %s" % (idToAtomicLabel[a1], perm[0], perm[1], perm[2]))

            founddtype = False
            for dtype in alldtypes:
              if dtype in improperDihedralTypes:
                idihedralid = improperDihedralTypes[dtype]
                dstring = dtype
                founddtype = True
                break

            if not founddtype:
              niDihedralTypes += 1
              mapiDTypes[niDihedralTypes] = alldtypes[0]
              idihedralid = niDihedralTypes
              improperDihedralTypes[alldtypes[0]] = niDihedralTypes
              dstring = alldtypes[0]

            niDihedrals += 1
            outImpropers += "\t%d\t%d\t%d\t%d\t%d\t%d\t# %s\n" % (niDihedrals, idihedralid, a1+1, a2+1, a3+1, a4+1, dstring)

  # print header
  header = "LAMMPS topology created from %s using pdb2lmp.py - By Henrique Musseli Cezar, 2021\n\n" % fname
  header += "\t%d atoms\n" % natoms
  if nbonds > 0:
    header += "\t%d bonds\n" % nbonds
  if nangles > 0:
    header += "\t%d angles\n" % nangles
  if printdih and (ndihedrals > 0):
    if ignoreimproper or (niDihedrals == 0):
      header += "\t%d dihedrals\n" % ndihedrals
    else:
      header += "\t%d dihedrals\n\t%d impropers\n" % (ndihedrals, niDihedrals)

  header += "\n\t%d atom types\n" % nmassTypes

  if nbondTypes > 0:
    header += "\t%d bond types\n" % nbondTypes
  if nangleTypes > 0:
    header += "\t%d angle types\n" % nangleTypes
  if printdih and (ndihedralTypes > 0):
    if ignoreimproper or (niDihedralTypes == 0):
      header += "\t%d dihedral types\n" % ndihedralTypes
    else:
      header += "\t%d dihedral types\n\t%d improper types\n" % (ndihedralTypes, niDihedralTypes)

  header += "\n"

  # add box info
  if fromBounds:
    boxsize = [(xmin,xmax),(ymin,ymax),(zmin,zmax)]
    boxsize[repaxis] = (boxsize[repaxis][0]-buffa/2., boxsize[repaxis][1]+buffa/2.)
    for i in range(3):
      if i == repaxis:
        continue
      boxsize[i] = (boxsize[i][0]-buffo/2., boxsize[i][1]+buffo/2.)
    header += "\t%.8f\t%.8f\t xlo xhi\n\t%.8f\t%.8f\t ylo yhi\n\t%.8f\t%.8f\t zlo zhi\n" % (boxsize[0][0], boxsize[0][1], boxsize[1][0], boxsize[1][1], boxsize[2][0], boxsize[2][1])
  else:
    if orthogonal:
      header += "\t%.8f\t%.8f\t xlo xhi\n\t%.8f\t%.8f\t ylo yhi\n\t%.8f\t%.8f\t zlo zhi\n" % (0., boxinfo[0][0], 0., boxinfo[1][1], 0., boxinfo[2][2])
    else:
      header += "\t%.8f\t%.8f\t xlo xhi\n\t%.8f\t%.8f\t ylo yhi\n\t%.8f\t%.8f\t zlo zhi\n\t%.8f\t%.8f\t%.8f\t xy xz yz\n" % (0., boxinfo[0][0], 0., boxinfo[1][1], 0., boxinfo[2][2], boxinfo[1][0], boxinfo[2][0], boxinfo[2][1])

  # print Coeffs
  outCoeffs = "Pair Coeffs\n\n"

  for i in range(1,nmassTypes+1):
    outCoeffs += "\t%d\teps\tsig\t# %s\n" % (i, mapTypes[i])

  if nbonds > 0:
    outCoeffs += "\nBond Coeffs\n\n"

    for i in range(1,nbondTypes+1):
      outCoeffs += "\t%d\tK\tr_0\t# %s\n" % (i, mapbTypes[i])

  if nangles > 0:
    outCoeffs += "\nAngle Coeffs\n\n"

    for i in range(1,nangleTypes+1):
      outCoeffs += "\t%d\tK\ttetha_0 (deg)\t# %s\n" %(i, mapaTypes[i])

  if printdih and (ndihedrals > 0):
    outCoeffs += "\nDihedral Coeffs\n\n"

    for i in range(1,ndihedralTypes+1):
      outCoeffs += "\t%d\tK\tn\tphi_0 (deg)\tw\t# %s\n" % (i, mapdTypes[i])

    if not ignoreimproper and (niDihedralTypes > 0):
      outCoeffs += "\nImproper Coeffs\n\n"

      for i in range(1,niDihedralTypes+1):
        outCoeffs += "\t%d\tK\txi_0 (deg)\t# %s\n" % (i, mapiDTypes[i])

  lmpstring = header+"\n"+outMasses+"\n"
  if not suppcoeff:
    lmpstring += outCoeffs+"\n"+outAtoms
  else:
    lmpstring += outAtoms

  if nbonds > 0:
    lmpstring += "\n"+outBonds
  if nangles > 0:
    lmpstring += "\n"+outAngles
  if printdih and (ndihedrals > 0):
    if ignoreimproper or (niDihedralTypes == 0):
      lmpstring += "\n"+outDihedrals
    else:
      lmpstring += "\n"+outDihedrals+"\n"+outImpropers


  with open("file.lmp",'w') as file:
    file.write(lmpstring)

  return lmpstring
