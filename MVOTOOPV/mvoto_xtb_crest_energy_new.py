import pandas.util  # Assuming 'util' is an alias for pandas.util
import pandas as pd

df=pd.read_pickle('filtered_molecules.pkl')



from pyscf.data import nist
import time

#conversion en eV
au2ev = nist.HARTREE2EV

def find_homo_lumo(myhf, au2ev):
    """Function that returns the HOMO and LUMO index and the HOMO energy in eV

    Args:
        mf_pyscf (pyscf object): pyscf meam-field object of the molecule to be evaluated.
    """
    # Index of HOMO and LUMO
    lumo_idx = myhf.mo_occ.tolist().index(0.)
    homo_idx = lumo_idx - 1

    # Calculate the HOMO Homo-LUMO
    E_HOMO = myhf.mo_energy[homo_idx]*au2ev
    E_LUMO = myhf.mo_energy[lumo_idx]*au2ev
    E_g = abs(E_HOMO - E_LUMO)

    return E_HOMO, E_LUMO, E_g


from pyscf import dft,scf,gto
def pyscf_calculation(mol_xtb_xyz, directory):

    mol_pyscf = gto.Mole(
        atom=mol_xtb_xyz,
        charge=0,
        spin=0,
        basis="6-31G(2df,p)",
        symmetry=True,
        unit='Angstrom'
    )
    mol_pyscf.build()

   # DFT calculation with B3LYP functional
    mdf = dft.RKS(mol_pyscf, xc="B3LYP").run()


    Eks1_homo, Eks1_lumo, Eks1_g = find_homo_lumo(mdf, au2ev)


    return Eks1_homo, Eks1_lumo, Eks1_g


import os
from rdkit import Chem

from rdkit.Chem import AllChem
import rdkit
import subprocess as sp

import re

def log_1error(file_path, process, error_message):
    with open("error_log.txt", "a") as error_file:
        error_file.write(f"Error processing {file_path} with {process}: {error_message}\n")

def log_error(smiles_key, method_gfn, stage, error_message, error_file="error_log.txt"):
    with open(error_file, "a") as f:
        f.write(f"SMILES: {smiles_key}, GFN: {method_gfn}, Stage: {stage}, Error: {error_message}\n")

"""###  <a id='toc1_1_'></a>[Fonction pour le Calcul XTB : Optimisation et Pré-optimisation](#toc0_)

Calcul d'optimisation et pré-optimisation avec `XTB` en utilisant les modèles `GFN2, GFN1, et GFN0`.
"""

def clean_xtb_files1():
    #---------------------------------------------------------------------
        # Clean up output files from CREST processes
        #FIXME To call after each xtb function
        #---------------------------------------------------------------------
        sp.run(['rm', 'bondlengths', 'charges', 'coord', 'coord.original', 'cregen_0.tmp',
                'cregen_1.tmp', 'cre_members', 'crest_best.xyz', 'crest_conformers.xyz',
                'crest.energies', 'crest_rotamers.xyz', 'gfnff_charges', 'gfnff_topo',
                '.history.xyz', 'wbo', 'crest_property.xyz', 'gfnff_adjacency', '.UHF',
                'ensemble_energies.log', 'charges3', 'charges', 'molden.input', 'crest_0.mdrestart',
                'crest_dynamics.trj', 'crestopt.log', 'crest.restart','crest.engrad','crestopt.xyz', 'crest_input_copy.xyz','g98.out'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        # For folder
        sp.run(['rm', '-r', 'calculation.level.1'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)

def clean_xtb_files():
    #---------------------------------------------------------------------
    # Clean up output files from xtb
    #FIXME To call after running crest function and leave crest_bestxyz file
    #---------------------------------------------------------------------
    sp.run(['rm', 'bondlengths', 'charges', 'coord', 'coord.original' , 'vibspectrum', 'hessian', 'gfnff_charges', 'gfnff_topo', 'wfn.xtb', 'xtbhess.xyz',
            '.history.xyz', 'struc.xyz', 'wbo', 'xtbopt.xyz', 'xtbopt.log', '.xtboptok','g98.out',
            'xtbrestart', 'xtbtopo.mol', 'xtblast.xyz', 'gfnff_adjacency', '.UHF',
            'ensemble_energies.log', 'charges3', 'charges', 'molden.input','test.smi','pat.xyz','sample.sdf'],
           stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    # For folder
    sp.run(['rm', '-r', 'calculation.level.1', 'PROP'], stdout=sp.DEVNULL, stderr=sp.DEVNULL)

def generate_3d_conformation(smiles_key, smiles, working_dir,path_xyz):
    """
    Generate a 3D conformation using Open Babel from SMILES and store the coordinates in a file.
    """
    # Define path for the output XYZ file
    path_xyz = working_dir / f'{smiles_key}.xyz'

    # Check if the file already exists
    if not path_xyz.exists():
        try:
            # Execute the Open Babel command
                # optimisation avec openbabel

            with open('test.smi', 'w') as f:
                f.writelines([smiles])

            OB_process= sp.run('obabel  test.smi -O  pat.xyz --gen3d -h --best --errorlevel 2 --minimize --ff MMFF94s  --steps 15000 --sd --crit 1e-9 --log', shell=True, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True, text=True)
            sp.run(['cp', 'pat.xyz', str(path_xyz)], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
            return path_xyz
        except sp.CalledProcessError as e:
            print(f"An error occurred while processing molecule {smiles_key}: {e}")
        return None

def get_xtb_energy_1(path_xtb1_opt_log, smiles_key, working_dir):
        path_xtb1_opt_log = working_dir / f'{smiles_key}_xtb1_opt.log'

        if path_xtb1_opt_log.exists():
            with open(path_xtb1_opt_log, 'r') as f:
                text_content = f.readlines()

            # Read the output (implementation details omitted)

            output_index = [i for i in range(len(text_content)) if 'Property Printout' in text_content[i]]
            text_content = text_content[output_index[0]:]
            homo_data = [x for x in text_content if '(HOMO)' in x]
            lumo_data = [x for x in text_content if '(LUMO)' in x]
            homo_lumo_gap = [x for x in text_content if 'HOMO-LUMO GAP' in x]
            total = [x for x in text_content if 'TOTAL ENERGY' in x]
            total0= float(total[0].split(' ')[-5])
            total1=total0*au2ev
            total2=total0*627.509474
            lumo_val = float(lumo_data[0].split(' ')[-2])
            homo_val = float(homo_data[0].split(' ')[-2])
            homo_lumo_val = float(homo_lumo_gap[0].split(' ')[-5])

             # Write the properties to a single file (modify as needed)
            with open(os.path.join(working_dir, f'{smiles_key}_properties.txt'), 'a') as f:
                 f.write(f'LUMO: {lumo_val}\n')
                 f.write(f'HOMO: {homo_val}\n')
                 f.write(f'HOMO-LUMO GAP: {homo_lumo_val}\n')

        return  homo_lumo_val, homo_val, lumo_val


import time
from pathlib import Path


def run_xtb_process(path_xyz, output_xyz, log_path, path_output_sdf, smiles_key, method_gfn,charge):
    try:
        start_xtb = time.time()
        xtb_process = sp.run(
            ["xtb", str(path_xyz) ,"--opt", "vtight", "--gfn",str(method_gfn), "--chrg",str(charge),"--parallel", "4", "--alpb", "toluene"],
            stdout=sp.PIPE, stderr=sp.PIPE, text=True, check=True
        )
        xtb_time = time.time() - start_xtb

        if os.path.exists("xtbopt.xyz"):
            sp.run(['cp', 'xtbopt.xyz', str(output_xyz)], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        if os.path.exists("xtbopt.xyz"):
            sp.run(['obabel', 'xtbopt.xyz', '-O', str(path_output_sdf)], stdout=sp.DEVNULL, stderr=sp.DEVNULL)

        with open(log_path, "w") as fl:
            fl.write(xtb_process.stdout)
        return xtb_time, output_xyz, log_path
    except sp.CalledProcessError as e:
        log_error(smiles_key, method_gfn, "xtb", str(e))
        return None, None, None

def run_crest_process(path_xyz, output_xyz, log_path, smiles_key, method_gfn,charge):
    try:
        start_crest = time.time()
        crest_process = sp.run(
            ["crest", str(path_xyz), "--gfn", str(method_gfn), "--mquick", "--noreftopo", "--opt", "vtight","--chrg",str(charge),"--parallel", "4", "--alpb", "toluene"],
            stdout=sp.PIPE, stderr=sp.PIPE, text=True, check=True
        )
        crest_time = time.time() - start_crest

        if os.path.exists("crestopt.xyz"):
            sp.run(['cp', 'crestopt.xyz', str(output_xyz)], stdout=sp.DEVNULL, stderr=sp.DEVNULL)

        with open(log_path, "w") as fl1:
            fl1.write(crest_process.stdout)
        return crest_time, output_xyz
    except sp.CalledProcessError as e:
        log_error(smiles_key, method_gfn, "crest", str(e))
        return None, None

import pandas as pd
from pathlib import Path
import gc
import os

def calculate_properties_xtb_crest(df, working_dir):
    """
    Évalue l'énergie avec xTB sur les coordonnées xyz générées pour chaque entrée du DataFrame,
    en traitant les données par lots de 100 molécules pour éviter les problèmes de mémoire.

    Paramètres:
    df : Pandas DataFrame contenant les colonnes 'SMILES' et 'smiles_key'.
    working_dir : Répertoire où tous les fichiers seront sauvegardés et traités.

    Retourne:
    df_ENERGY : DataFrame contenant les propriétés calculées.
    """
    # Définition des colonnes pour le résultat final
    columns = ["smiles_key", "SMILES", "HOMO(eV)", "HOMO_xtb(eV)", "HOMO_DFT(eV)", 
               "LUMO(eV)", "LUMO_xtb(eV)", "LUMO_DFT(eV)", "GAP(eV)", "GAP_xtb(eV)", "GAP_DFT(eV)"]
    results = []
    df = df.copy()

    # Créer le répertoire de travail s'il n'existe pas
    working_dir = Path(working_dir) if isinstance(working_dir, str) else working_dir
    working_dir.mkdir(parents=True, exist_ok=True)

    # Traitement par lots (100 molécules par lot)
    batch_size = 100
    num_batches = (len(df) + batch_size - 1) // batch_size  # Calcul du nombre total de lots

    for batch_index in range(num_batches):
        # Sélection du lot courant
        batch_df = df.iloc[batch_index * batch_size:(batch_index + 1) * batch_size]

        for _, row in batch_df.iterrows():
            smiles_key = row["smiles_key"]
            smiles = row["smiles"]
            HOMO, LUMO, GAP = row["HOMO(eV)"], row["LUMO(eV)"], row["Gap(eV)"]
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            if mol == None: 
                print("INVALID")
        
            charge = Chem.rdmolops.GetFormalCharge(mol)
            atom_number = mol.GetNumAtoms()
            # Définition des chemins pour les fichiers
            paths = {key: working_dir / f"{smiles_key}_{key}.xyz" for key in ["xyz", "xtb1_pre_opt", "xtb1_opt", "crest1_opt"]}
            logs = {key: working_dir / f"{smiles_key}_{key}.log" for key in ["xtb1_pre_opt", "crest1", "xtb1_opt"]}
            sdf_paths = {key: working_dir / f"{smiles_key}_{key}.sdf" for key in ["xtb_pre_opt", "xtb_opt"]}

            # Génération de la structure 3D
            result = generate_3d_conformation(smiles_key, smiles, working_dir, paths["xyz"])
            if result is None:
                log_error(smiles_key, "N/A", "3D generation", "Failed to generate 3D structure")
                continue

            # Essayer plusieurs méthodes pour xTB et CREST
            retry_methods = ['2', '1', '0']
            for method in retry_methods:
                xtb_result = run_xtb_process(result, paths["xtb1_pre_opt"], logs["xtb1_pre_opt"], sdf_paths["xtb_pre_opt"], smiles_key, method,charge)
                clean_xtb_files()

                if xtb_result[0] is not None:
                    crest_result = run_crest_process(paths["xtb1_pre_opt"], paths["crest1_opt"], logs["crest1"], smiles_key, method,charge)
                    clean_xtb_files1()
                    if crest_result[0] is not None:
                        break
                    log_error(smiles_key, method, "conformer search", "Topological changes detected or CREST failed")
                else:
                    log_error(smiles_key, method, "preoptimization", "Topological changes detected or xTB failed")
            else:
                continue

            # Optimisation finale avec xTB
            final_xtb_result = run_xtb_process(paths["crest1_opt"], paths["xtb1_opt"], logs["xtb1_opt"], sdf_paths["xtb_opt"], smiles_key, method,charge)
            if final_xtb_result[0] is None or not os.path.exists(paths["xtb1_opt"]):
                log_error(smiles_key, method, "optimization", "Topological changes detected or xTB failed")
                continue
            clean_xtb_files()

            # Lire les résultats optimisés
            with open(paths["xtb1_opt"], 'r') as f:
                lines = f.readlines()

            if len(lines) < 3:
                log_error(smiles_key, method, "file_read", "xtb1_opt.xyz file is empty or corrupted")
                continue

            mol_xtb1_xyz = '\n'.join(lines[2:])
            Eks1_homo, Eks1_lumo, Eks1_g = pyscf_calculation(mol_xtb1_xyz, working_dir)
            homo_lumo_val, homo_val, lumo_val = get_xtb_energy_1(logs["xtb1_opt"], smiles_key, working_dir)

            # Ajouter les résultats au tableau
            results.append([
                smiles_key, smiles, HOMO, homo_val, Eks1_homo, LUMO, lumo_val, Eks1_lumo, GAP, homo_lumo_val, Eks1_g
            ])

        # Forcer la libération de la mémoire après chaque lot
        gc.collect()

    # Retourner le DataFrame final contenant les résultats
    return pd.DataFrame(results, columns=columns)

MY_crest1_GDB9 = Path(os.getcwd()) / 'MY_crest1_GDB9'
MY_crest1_GDB9.mkdir(exist_ok=True)
df_ENERGY1 = calculate_properties_xtb_crest(df, MY_crest1_GDB9)
df_ENERGY1.to_pickle('dataset_for_ml.pkl')
