import pandas as pd
from pathlib import Path
import hjson
import csv
from rdkit import Chem
import numpy as np
import progressbar
import re
from openbabel import openbabel
from rdkit.Chem import AllChem

converter = openbabel.OBConversion()
converter.SetInAndOutFormats('smi', 'can')


def canonicalize(can_smiles: str):
    molecule = openbabel.OBMol()
    converter.ReadString(molecule,((can_smiles.replace(" ","")).replace("'","")).strip())
    return (converter.WriteString(molecule)).strip()

this_Path = Path.cwd()

metal_halide = '[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79]~[ClX1,BrX1,IX1,FX1,O]'
no_C_N = '[c,C,n,N]~[ClX1,BrX1,IX1,FX1,OX1X2]' #exclude Carbon Nitrogen
Metal_Covalent = '[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79]~[C,c,N,n,P]'
metal = '[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79]' #für M-Ligand-Schleife
#unbound_metal_SMARTS = Chem.MolFromSmarts(unbound_metal)
metal_halide_SMARTS = Chem.MolFromSmarts(metal_halide)
no_C_N_SMARTS = Chem.MolFromSmarts(no_C_N)
M_Cov_SMARTS = Chem.MolFromSmarts(Metal_Covalent)
metal_SMARTS = Chem.MolFromSmarts(metal)

with open(this_Path / 'final_output'/'BHR_sorted_reduced.csv') as f:
    dataframe = pd.read_csv(f)
dataframe = dataframe.assign(Catalyst=(''), Ligand=(''))
print(dataframe['Additive'])
a0 = dataframe
print(a0)

#Palladium[Pd-]-Reiniger: output ist kanonisch
for i in progressbar.progressbar(range(0, len(dataframe.index))):
    smiles_import = dataframe.at[i, 'Additive'].split('.')
    smiles = []
    [smiles.append(x) for x in smiles_import if x not in smiles]

    for j in range(0, len(smiles)):
        smiles_mol = Chem.MolFromSmiles(smiles[j])
        for z in range(1, 5):
            #try:
            negative_metal_SMARTS = Chem.MolFromSmarts(f'[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;-{z}]')
            if len(smiles_mol.GetSubstructMatches(negative_metal_SMARTS)) >= 1:
                rxn = AllChem.ReactionFromSmarts(f'[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;-{z}:1]>>[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;+0:1]')
                product = [Chem.MolToSmiles(y) for y in rxn.RunReactants((Chem.MolFromSmiles(smiles[j]),))[0]] #Änderung 4
                #for c in range(0, len(product)):
                print(i)
                print(product[0])
                for l in range(0, len(product)):
                    product[l] = canonicalize(product[l])
                print(product[0])
                dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'].replace(smiles[j], '')
                dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'] + '.' + product[0] #Platz getauscht mit replace
print('Test1')
a1 = dataframe
print(a1)
M_Lig_double_SMARTS = Chem.MolFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79]=[#6]')
M_Lig_single_SMARTS = Chem.MolFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79]-[C,c,N,n,P,Cl,Br,I,F]')
unbound_metal = '[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X0]'
unbound_metal_SMARTS = Chem.MolFromSmarts(unbound_metal)
#M=Lig breaker ->Cat./Add. output ist kanonisch
for i in progressbar.progressbar(range(0, len(dataframe.index))):
    smiles_import = dataframe.at[i, 'Additive'].split('.')
    smiles = []
    [smiles.append(x) for x in smiles_import if x not in smiles]

    for j in range(0, len(smiles)):
        product = [smiles[j]]
        smiles_mol = Chem.MolFromSmiles(smiles[j])
        while len(Chem.MolFromSmiles(product[0]).GetSubstructMatches(M_Lig_double_SMARTS)) >= 1:
            try:
                rxn = AllChem.ReactionFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79:1]=[#6:2]>>[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79:1].[#6:2]')
                product = [Chem.MolToSmiles(y) for y in rxn.RunReactants((Chem.MolFromSmiles(product[0]),))[0]] #Änderung2: 1 rausnehmen
                print(i)
                product[0] = canonicalize(product[0])
                product[1] = canonicalize(product[1])
                #if len(Chem.MolFromSmiles(product[1]).GetSubstructMatches(metal_SMARTS)) == 0:
                    #dataframe.at[i, 'Ligand'] = dataframe.at[i, 'Ligand'] + '.' + product[1] #Kohlenstoff in rxn-SMARTS anpassen!
                if len(Chem.MolFromSmiles(product[0]).GetSubstructMatches(unbound_metal_SMARTS)) >= 1:
                    dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'] + '.' + product[0]
                    dataframe.at[i, 'Ligand'] = dataframe.at[i, 'Ligand'] + '.' + product[1]
                if len(Chem.MolFromSmiles(product[1]).GetSubstructMatches(metal_SMARTS)) >= 1: #für Ferrocen????
                    dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'] + '.' + product[1]
                elif len(Chem.MolFromSmiles(product[0]).GetSubstructMatches(M_Lig_single_SMARTS)) >= 1:
                    dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'] + '.' + product[0] #neue Zeile
                    dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'] + '.' + product[1]
            except ValueError or AttributeError:
                print(i)
                continue
        if len(smiles_mol.GetSubstructMatches(M_Lig_double_SMARTS)) >= 1:
            dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'].replace(smiles[j], '') #PROBLEMZEILE

print('Hello')
a2 = dataframe
print('Test2')


M_Lig_single = '[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79]-[C,c,N,n,P,Cl,Br,I,F]'
M_Lig_single_SMARTS = Chem.MolFromSmarts(M_Lig_single)


print('Test3')
print(dataframe['Additive'])
print(dataframe['Catalyst'])
print(dataframe['Ligand'])
a3 = dataframe

Ligi_SMARTS = Chem.MolFromSmarts('[c,C,n,N,P]')
Addi_SMARTS = Chem.MolFromSmarts('[F,Cl,Br,I;X0]')

Cyclopentadienyl_SMARTS = Chem.MolFromSmarts('[#6;-1]~1~[#6]~[#6]~[#6]~[#6]~1')
#Cp-Fix: output ist kanonisch
for i in progressbar.progressbar(range(0, len(dataframe.index))):
    smiles_import = dataframe.at[i, 'Additive'].split('.')
    smiles = []
    [smiles.append(x) for x in smiles_import if x not in smiles]
    # Cp canonicalisieren:
    try:
        for j in range(0, len(smiles)):
            if len(Chem.MolFromSmiles(smiles[j]).GetSubstructMatches(Cyclopentadienyl_SMARTS)) >= 1:
                #print('Cp:')
                #print(i)
                dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'] + '.' + canonicalize(smiles[j])
                dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'].replace(smiles[j], '')
    except AttributeError:
        continue
#M-Lig -> Catalyst/Ligand/Additive
for i in progressbar.progressbar(range(0, len(dataframe.index))):
    smiles_import = dataframe.at[i, 'Additive'].split('.')
    smiles = []
    [smiles.append(x) for x in smiles_import if x not in smiles]
    # Cp canonicalsieren:
    for j in range(0, len(smiles)):
        product = [smiles[j]]
        smiles_mol = Chem.MolFromSmiles(smiles[j])
        if Chem.MolFromSmiles(product[0]) is None:
            #print(f'Smiles is not convertable {product[0]}')
            #print(i)
            continue
        #if len(Chem.MolFromSmiles(product[0]).GetSubstructMatches(Cyclopentadienyl_SMARTS)) >= 1:
            #product[0] = canonicalize(product[0])
        while len(Chem.MolFromSmiles(product[0]).GetSubstructMatches(M_Lig_single_SMARTS)) >= 1:
            try:
                #a = Chem.rdmolops.GetFormalCharge(smiles_mol)
                a = Chem.rdmolops.GetFormalCharge(Chem.MolFromSmiles(product[0]))
                print(i)
                if '-]' in product[0]: #ändern durch abfrage nach negativer Ladung im SMARTS!
                    a = a+1
                rxn = AllChem.ReactionFromSmarts(f'[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;+{a},X1:1]-[c,C,n,N,P,Cl,Br,I,F:2]>>[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;+{a+1}:1].[c,C,n,N,P,Cl,Br,I,F;-1:2]')
                product = [Chem.MolToSmiles(y) for y in rxn.RunReactants((Chem.MolFromSmiles(product[0]),))[0]] #(y, 1) ist böse!
                #product = [Chem.MolToSmiles(y, 1, kekuleSmiles=True) for y inrxn.RunReactants((Chem.MolFromSmiles(product[0]),))[0]]
                #if i == 9429:
                    #print(product[0])
                    #print(product[1])
                    #product[1] = canonicalize(product[1])
                if len(Chem.MolFromSmiles(canonicalize(product[1])).GetSubstructMatches(Cyclopentadienyl_SMARTS)) >= 1:
                    product[1] = canonicalize(product[1])
                    print('neu')
                    #print(product[1])
                if len(Chem.MolFromSmiles(product[0]).GetSubstructMatches(unbound_metal_SMARTS)) >= 1:
                    dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'] + '.' + canonicalize(product[0])
                if len(Chem.MolFromSmiles(product[1]).GetSubstructMatches(Ligi_SMARTS)) >= 1: #Error bei 9429 (+2)
                    dataframe.at[i, 'Ligand'] = dataframe.at[i, 'Ligand'] + '.' + canonicalize(product[1])
                elif len(Chem.MolFromSmiles(product[1]).GetSubstructMatches(Addi_SMARTS)) >= 1:
                    dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'] + '.' + canonicalize(product[1])
            except ValueError:
                continue
        if len(smiles_mol.GetSubstructMatches(M_Lig_single_SMARTS)) >= 1:
            dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'].replace(smiles[j], '')

print('Test4')
print(dataframe['Additive'])
print(dataframe['Catalyst'])
print(dataframe['Ligand'])
a4 = dataframe
print('a4')
#unbound metal
for i in progressbar.progressbar(range(0, len(dataframe.index))):
    smiles_import = dataframe.at[i, 'Additive'].split('.')
    smiles = []
    [smiles.append(x) for x in smiles_import if x not in smiles]

    for j in range(0, len(smiles)):
        drop_unboundsmiles = False
        #drop_boundmetal = False #unnöööööööööötig du Lauch
        smiles_mol = Chem.MolFromSmiles(smiles[j])

        if smiles_mol is None:
            #print(f'Smiles is not convertable {smiles[j]}')
            continue
        if len(smiles_mol.GetSubstructMatches(unbound_metal_SMARTS)) >= 1: #len(smiles_mol.GetSubstructMatches(no_C_N_SMARTS)) == 0 and #condition für einzelnes Metall(Salz) und Lig+Metall schreiben
            drop_unboundsmiles = True
        #elif len(smiles_mol.GetSubstructMatches(M_Cov_SMARTS)) >= 1 or len(smiles_mol.GetSubstructMatches(metal_halide_SMARTS)) >= 1:
        #    drop_boundmetal = True
        if drop_unboundsmiles == True:
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'] + '.' + smiles[j]
            dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'].replace(smiles[j], '') #Änderung 5, war auskommentiert

dataframe.reset_index(drop=True, inplace=True)

dataframe.replace(to_replace='^((\.)*)|\.((\.)*)$',value='', inplace=True, regex=True)
dataframe.replace(to_replace='\.\.((\.)*)',value='.', inplace=True, regex=True)

print('Test5')
print(dataframe['Additive'])
print(dataframe['Catalyst'])
print(dataframe['Ligand'])
a5 = dataframe
#dataframe.to_csv(this_Path / 'output_additive' / 'Cat#3.csv', index=False)

#with open(this_Path / 'output_additive' / 'Cat#3.csv') as f:
#dataframe = pd.read_csv(f, keep_default_na=False)


SMARTS_ligands = {
    'phosphane_PR3':['P([C,c])([C,c])[C,c]'],
    'NHCs':['nc(=*)n', '[C,c]1[N,n][C,c][C,c][N,n]1[c,C]', 'C1N(CCN1[C,c])[C,c]', 'C1-,=N(CCN1[C,c])[C,c]', '[C,c]1[N,n]([C,c]~[C,c][N,n]1[C,c])[c,C]'], #3 ist sehr allgemein (aliph. NHC core)
    'dba':['C(C=Cc1ccccc1)(C=Cc2ccccc2)=O'],
    'BINOL+SH':['c1cccc2c1c(c(cc2)[O,S])c3c(ccc4c3cccc4)[O,S]'],
    'bipyridines':['c1cc(ncc1)c2ccccn2'],
    '1,3-dicarbonyls+acac+coumarin':['C(C=O)C=O','CC(=CC(C)=O)[O;-1]','[C,c]([C,c]=O)[C,c]=O'],#noch prüfen
    '1,4-dicarbonyls':['C(=O)CCC=O'],#noch prüfen
    'ethylendiamines':['[#6]([#6]~[#7X3;H2,H1])~[#7X3;H2,H1]', '[#6]([#6][#7]=[#6])[#7]=[#6]', '[#6]([#6]=[#7][#6,O])=[#7][#6,O]'],#noch prüfen#für Kupfer
    'sugars':['[#6]([OX2H1])[#6][OX2H1]'],#noch prüfen
    'Proline_4-OH':['[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]','[*$([NX3H1,NX4H2+1]),*$([NX3](C)(C)C)]1[CH2]C([CH2][CX4H1]1[CX3](=[OX1])[OX2H1,OX1-1,N])O[H]'],#noch prüfen#alle statt nur Prolin?
    'COD':['C1CC=CCCC=C1'],#noch prüfen
    'NCCO': ['[#7]~[#6]~[#6]~[#8]'],  # noch prüfen, !sehr allgemein!
    'Allyl_Pd': ['C(C=C)[#46]'],  # noch prüfen
    'P=O-OR-OR-R': ['P(O[#6])(O[#6])(=O)[#6]'],  # noch prüfen
    'P=S-OR-OR-SR': ['P(O[#6])(O[#6])(=S)S[#6]'],  # noch prüfen
    #'biphenol': ['c1cc(c(cc1)[OX2H1])c2c(cccc2)[OX2H1]'],  # noch prüfen
    'all_phenols': ['c[OX2H1]'],
    'Pyridine-N-Oxid':['[O-][n+]1ccccc1'],
    'Cyclopentadienyl':['*[#6-1]1[#6]=[#6][#6]=[#6]1'],
    'Methylendiamin':['[#7H2;X3][#6][#7H2;X3]'],
    'PX3': ['[P;X3]([#6,#7,#5,#8])([#6,#7,#5,#8])[#6,#7,#5,#8]'],
    'c(nn)3': ['[#6]([#7][#7])([#7][#7])[#7][#7]'],
    #'Pd-Ar': ['[#46]c'], #aufpassen bei Zuordnung später! #aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaah
    'benzodiox_anisol_lig': ['c1(c2c(ccc1[#6])OCO2)[#6]'],
    'TH-quinoline':['c1ccnc2c1CCCC2~[#6,#7]'],
    'PN-(Bn)2':['c1cc(ccc1)C[NX3]Cc2ccccc2'],
    'phenanthroline++': ['c2:1c3c([c,n][c,n]c1[c,n][c,n][c,n]:[n,c]2)[c,n][c,n][c,n]:[n,c]3'],
    'thiophene-2C': ['c1c(scc1)C(=O)[O;-1]'],
    'nccn': ['C(C[NX3])[NX3]'],
    '1,2-dicarbonyl': ['C(=O)C=O'],
    'glyme': ['[CH3]OCCO[CH3]'],
    'tetralin=NOH':['c1c2c(ccc1)C(CCC2)=N[OX2H1]'],
    'beta-alanine': ['[NH2]CC[CX3](=O)[OX2H1]'], #check
    'xanthene': ['c31ccccc1Oc2ccccc2C3'], #check
    'SteglichCat+': ['n1ccc(cc1)N(C)C'], #check
    'letzte Ausnahmen': ['C1CN(C(N(C1)C)=O)C', 'C(C)(NC(C)=O)=O', 'C1N(C=C[N;+1]1)CCC', 'c1c(nccc1)[N;+1]2(CCCN(C2)C)', 'C(c1c(nnc1C)C)(c2c(nnc2C)C)c3c(nnc3C)C'],

}



#Ligandenzuordnung
for key in SMARTS_ligands:
    smartslist = SMARTS_ligands[key]
    for e in range(0, len(smartslist)):
        smarts_mol = Chem.MolFromSmarts(smartslist[e])
        #print('\n' + smartslist[e])
        for i in progressbar.progressbar(range(0, len(dataframe.index))):
            smiles_import = dataframe.at[i, 'Additive'].split('.')
            smiles = []
            [smiles.append(x) for x in smiles_import if x not in smiles]
            move_ligand = False
            for j in range(0, len(smiles)):
                smiles_mol = Chem.MolFromSmiles(smiles[j])

                if smiles_mol is None:
                    #print(i)
                    #print(f'Smiles is not convertable {smiles[j]}')
                    continue
                if len(smiles_mol.GetSubstructMatches(smarts_mol)) >= 1:
                    dataframe.at[i, 'Ligand'] = dataframe.at[i, 'Ligand'] + '.' + smiles[j]
                    dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'].replace(smiles[j], '') #Änderung1, war vorher auskommentiert

dataframe.reset_index(drop=True, inplace=True)

dataframe.replace(to_replace='^((\.)*)|\.((\.)*)$',value='', inplace=True, regex=True)
dataframe.replace(to_replace='\.\.((\.)*)',value='.', inplace=True, regex=True)

print('Test6')
print(dataframe['Additive'])
print(dataframe['Catalyst'])
print(dataframe['Ligand'])
a6 = dataframe

metal_Add_SMARTS = Chem.MolFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79]')

#Metals from Additive to Catalyst:
for i in progressbar.progressbar(range(0, len(dataframe.index))):
    smiles_import = dataframe.at[i, 'Additive'].split('.')
    smiles = []
    [smiles.append(x) for x in smiles_import]  # if x not in smiles rausgenommen

    for j in range(0, len(smiles)):
        smiles_mol = Chem.MolFromSmiles(smiles[j])
        try:
            if len(smiles_mol.GetSubstructMatches(metal_Add_SMARTS)) >= 1:
                dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst']+'.'+smiles[j]
                dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'].replace(smiles[j], '')
        except AttributeError:
            print('Cleverhans')
            print(i)
            continue

#Cat Refinement
Trihal_SMARTS = Chem.MolFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X3]([Cl,Br,I,F,O])([Cl,Br,I,F,O])[Cl,Br,I,F,O]')
Dihal_SMARTS = Chem.MolFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X2]([Cl,Br,I,F,O])[Cl,Br,I,F,O]')
Monohal_SMARTS = Chem.MolFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X1]-[Cl,Br,I,F,O]')
I_Oxides_SMARTS = Chem.MolFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X1]-[OX2]-[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X1]')
II_Oxides_SMARTS = Chem.MolFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79]=[O]')
Pd_corr_SMARTS = Chem.MolFromSmarts('[#46X0;+1]')
Pd_corr_SMARTS_3 = Chem.MolFromSmarts('[#46X0;+3]')
hal_hydr_oxide_SMARTS = Chem.MolFromSmarts('[O,Cl,Br,I,F]')

#Catalyst refinement
for i in progressbar.progressbar(range(0, len(dataframe.index))):
    smiles_import = dataframe.at[i, 'Catalyst'].split('.')
    smiles = []
    [smiles.append(x) for x in smiles_import] #if x not in smiles rausgenommen

    for j in range(0, len(smiles)):
        smiles_mol = Chem.MolFromSmiles(smiles[j])
        if smiles_mol is None:
            #print(f'Smiles is not convertable {smiles[j]}')
            continue
        if len(smiles_mol.GetSubstructMatches(Trihal_SMARTS)) >= 1:
            rxn = AllChem.ReactionFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X3:1]([Cl,Br,I,F,O:2])([Cl,Br,I,F,O:3])[Cl,Br,I,F,O:4]>>[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X3;+3:1].[Cl,Br,I,F,O;-1:2]')
            trihalogenides = [Chem.MolToSmiles(x) for x in rxn.RunReactants((Chem.MolFromSmiles(smiles[j]),))[0]] #Änderung5
            for c in range(0,len(trihalogenides)):
                dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst']+'.'+trihalogenides[c]
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'].replace(smiles[j], '')
        if len(smiles_mol.GetSubstructMatches(Dihal_SMARTS)) >= 1:
            rxn = AllChem.ReactionFromSmarts('[Cl,Br,I,F,O:2]-[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X2:1]-[Cl,Br,I,F,O:3]>>[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;+2:1].[Cl,Br,I,F,O;-1:2]')
            dihalogenides = [Chem.MolToSmiles(x) for x in rxn.RunReactants((Chem.MolFromSmiles(smiles[j]),))[0]]
            for c in range(0, len(dihalogenides)):
                dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst']+'.'+dihalogenides[c]
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'].replace(smiles[j], '')
        if len(smiles_mol.GetSubstructMatches(Monohal_SMARTS)) >= 1:
            rxn = AllChem.ReactionFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X1:1]-[Cl,Br,I,F,O:2]>>[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;+1:1].[Cl,Br,I,F,O;-1:2]')
            monohalogenides = [Chem.MolToSmiles(x) for x in rxn.RunReactants((Chem.MolFromSmiles(smiles[j]),))[0]]
            for c in range(0,len(monohalogenides)):
                dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst']+'.'+monohalogenides[c]
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'].replace(smiles[j], '')
        if len(smiles_mol.GetSubstructMatches(I_Oxides_SMARTS)) >= 1:
            rxn = AllChem.ReactionFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X1:1]-[OX2:2]-[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;X1:3]>>[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;+1:1].[O;-2:2]')
            I_Oxides = [Chem.MolToSmiles(x) for x in rxn.RunReactants((Chem.MolFromSmiles(smiles[j]),))[0]]
            for c in range(0, len(I_Oxides)):
                dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst']+'.'+I_Oxides[c]
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'].replace(smiles[j], '')
        if len(smiles_mol.GetSubstructMatches(II_Oxides_SMARTS)) >= 1:
            rxn = AllChem.ReactionFromSmarts('[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79:1]=[O:2]>>[#26,#27,#28,#29,#44,#45,#46,#47,#76,#77,#78,#79;+2:1].[O;-2:2]')
            II_Oxides = [Chem.MolToSmiles(x) for x in rxn.RunReactants((Chem.MolFromSmiles(smiles[j]),))[0]]
            for c in range(0, len(II_Oxides)):
                dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst']+'.'+II_Oxides[c]
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'].replace(smiles[j], '')
print('Test7')
print(dataframe['Additive'])
print(dataframe['Catalyst'])
print(dataframe['Ligand'])
a7 = dataframe
#Pd+ - Cleanup
for i in progressbar.progressbar(range(0, len(dataframe.index))):
    smiles_import = dataframe.at[i, 'Catalyst'].split('.')
    smiles = []
    [smiles.append(x) for x in smiles_import]  # if x not in smiles rausgenommen

    for j in range(0, len(smiles)):
        smiles_mol = Chem.MolFromSmiles(smiles[j])
        if len(smiles_mol.GetSubstructMatches(Pd_corr_SMARTS)) >= 1:
            rxn = AllChem.ReactionFromSmarts('[#46X0;+1]>>[#46X0;+2]')
            Pd_corr = [Chem.MolToSmiles(x) for x in rxn.RunReactants((Chem.MolFromSmiles(smiles[j]),))[0]] #Änderung6
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst']+'.'+Pd_corr[0]
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'].replace(smiles[j], '')
        if len(smiles_mol.GetSubstructMatches(hal_hydr_oxide_SMARTS)) >= 1:
            #dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'] + '.' + smiles[j]
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'].replace(smiles[j], '')
#pd3+
for i in progressbar.progressbar(range(0, len(dataframe.index))):
    smiles_import = dataframe.at[i, 'Catalyst'].split('.')
    smiles = []
    [smiles.append(x) for x in smiles_import]  # if x not in smiles rausgenommen

    for j in range(0, len(smiles)):
        smiles_mol = Chem.MolFromSmiles(smiles[j])
        if len(smiles_mol.GetSubstructMatches(Pd_corr_SMARTS_3)) >= 1:
            rxn = AllChem.ReactionFromSmarts('[#46X0;+3]>>[#46X0;+2]')
            Pd_corr = [Chem.MolToSmiles(x) for x in rxn.RunReactants((Chem.MolFromSmiles(smiles[j]),))[0]] #Änderung6
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst']+'.'+Pd_corr[0]
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'].replace(smiles[j], '')
        if len(smiles_mol.GetSubstructMatches(hal_hydr_oxide_SMARTS)) >= 1:
            #dataframe.at[i, 'Additive'] = dataframe.at[i, 'Additive'] + '.' + smiles[j]
            dataframe.at[i, 'Catalyst'] = dataframe.at[i, 'Catalyst'].replace(smiles[j], '')

dataframe.reset_index(drop=True, inplace=True)


print('Cleverhans!')
#Unify NHCs/Ligands
NHC_aromatic_SMARTS = Chem.MolFromSmarts('[c]1[n](ccn1[C,c])[c,C]')
NHC_unaromatic_db_SMARTS = Chem.MolFromSmarts('[C]1[N]([C]=[C][N]1[c,C])[C,c]')
NHC_unaromatic_nodb_SMARTS = Chem.MolFromSmarts('[C]1[N]([C]-[C][N]1[c,C])[C,c]')

for i in progressbar.progressbar(range(0, len(dataframe.index))):
    smiles_import = dataframe.at[i, 'Ligand'].split('.')
    smiles = []
    [smiles.append(x) for x in smiles_import if x not in smiles]
    for j in range(0, len(smiles)):
        if smiles[j] != 'ValueNotFound':
            smiles_mol = Chem.MolFromSmiles(smiles[j])
            if len(smiles_mol.GetSubstructMatches(NHC_aromatic_SMARTS)) >= 1:
                rxn = AllChem.ReactionFromSmarts("[#6:1]1[#7:2][#6:3][#6:4][#7:5]1>>[CH0:1]1[NH0;+0:2][C:3]=[C:4][NH0;+0:5]1")
                product = [Chem.MolToSmiles(x) for x in rxn.RunReactants((smiles_mol,))[0]]
                p = canonicalize(product[0])
                p = re.sub('N1C=CN\(C1\)', 'N1C=CN([C]1)', p)
                dataframe.at[i, 'Ligand'] = dataframe.at[i, 'Ligand'].replace(smiles[j], '')
                dataframe.at[i, 'Ligand'] = dataframe.at[i, 'Ligand'] + '.' + p
            if len(smiles_mol.GetSubstructMatches(NHC_unaromatic_db_SMARTS)) >= 1:
                    rxn = AllChem.ReactionFromSmarts("[#6:1]1[#7:2][#6:3]-,=[#6:4][#7:5]1>>[CH0:1]1[NH0;+0:2][C:3]=[C:4][NH0;+0:5]1")
                    product = [Chem.MolToSmiles(x) for x in rxn.RunReactants((smiles_mol,))[0]]
                    p = canonicalize(product[0])
                    p = re.sub('N1C=CN\(C1\)', 'N1C=CN([C]1)', p)
                    dataframe.at[i, 'Ligand'] = dataframe.at[i, 'Ligand'].replace(smiles[j], '')
                    dataframe.at[i, 'Ligand'] = dataframe.at[i, 'Ligand'] + '.' + p
            if len(smiles_mol.GetSubstructMatches(NHC_unaromatic_nodb_SMARTS)) >= 1:
                rxn = AllChem.ReactionFromSmarts("[#6:1]1[#7:2][#6:3][#6:4][#7:5]1>>[CH0:1]1[NH0;+0:2][C:3][C:4][NH0;+0:5]1")
                product = [Chem.MolToSmiles(x) for x in rxn.RunReactants((smiles_mol,))[0]]
                p = canonicalize(product[0])
                p = re.sub('N1CCN\(C1\)', 'N1CCN([C]1)', p)
                dataframe.at[i, 'Ligand'] = dataframe.at[i, 'Ligand'].replace(smiles[j], '')
                dataframe.at[i, 'Ligand'] = dataframe.at[i, 'Ligand'] + '.' + p

#alles Kanonikalisieren!
for i in progressbar.progressbar(range(0, len(dataframe.index))):
    for k in 0, 2, 4, 5, 6, 7, 8:
        smiles_import = dataframe.iat[i, k].split('.')
        smiles = []
        [smiles.append(x) for x in smiles_import]
        for j in range(0, len(smiles)):
            dataframe.iat[i, k] = dataframe.iat[i, k].replace(smiles[j], canonicalize(smiles[j]))

dataframe.replace(to_replace='^((\.)*)|\.((\.)*)$',value='', inplace=True, regex=True)
dataframe.replace(to_replace='\.\.((\.)*)',value='.', inplace=True, regex=True)
dataframe.replace('', 'ValueNotFound', inplace=True)
#Reaction SMILES entfernen
dataframe.drop(columns='Reaction SMILES', inplace=True)

dataframe.to_csv(this_Path / 'full_program' / 'FRS_004.csv', index=False)