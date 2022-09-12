"""
Program aimed to calculate the Solvent-Accessible Surface Area (SASA) of a protein 
by using the functions of the module tools.
"""


#Required modules:
import tools as tl
import argparse
import math
import time
from tqdm import tqdm
from tabulate import tabulate


def parse_args():  
    
    parser = argparse.ArgumentParser(description = "Program aimed to calculate the Solvent-Accessible Surface Area (SASA) of a protein")  #Create an ArgumentParser object.
    parser.add_argument("-pdb", type=str, dest="pdbFile", help="PDB file directory", required="TRUE")  #Add positional argument 1: -pdb.
    parser.add_argument("-n", type=int, dest="pointsNum", help="Number of points evenly distributed on a sphere", required="TRUE")  #Add positional argument 2: -n.

    return parser.parse_args()


def SASA():  #The main.
    
    start_time = time.time() #Time when the program execution starts.
    
    args = parse_args() 
    pointsNum = args.pointsNum 
    df = tl.atomic_coordinates(args.pdbFile) 
    
    VDW_RADIUS = {'H': 1.2, 'C': 1.5, 'N': 1.5, 'O': 1.4, 'P': 1.8, 'S': 1.8} 
    SASA_MAX = {"ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0, "GLU": 223.0,	"GLN": 225.0, "GLY": 104.0, "HIS":224.0, "ILE":197.0, "LEU":201.0, "LYS":236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0, "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0}                       
 
    dico = {}
    
    s_exp_abs_res = 0
    s_tot_res = 0
    s_exp_prot = 0
    s_tot_prot = 0
    
    res_count = int(df.loc[1, ["res_num"]][-1])  #The first residue number is different from 1 in some PDB files.
    pos = 1
    
    for i, atom1 in tqdm(df.iterrows(), total=len(df)):  #Iterate over DataFrame rows.
        name_atom1, coord_atom1, name_res, num_res = atom1["atom_name"], [atom1["x"], atom1["y"], atom1["z"]], atom1["res_name"], atom1["res_num"]
        r_vdw_atom1 = VDW_RADIUS[name_atom1]
        points_coordinates = tl.sphere(pointsNum, coord_atom1, r_vdw_atom1)
        points_nonexp = []
        for j, atom2 in df.iterrows():
            if i == j : continue
            name_atom2, coord_atom2 = atom2["atom_name"], [atom2["x"], atom2["y"], atom2["z"]]
            r_vdw_atom2 = VDW_RADIUS[name_atom2]
            distance = math.dist(coord_atom1, coord_atom2)  #Distance between the 2 atoms.
            if distance < 10: 
                for k, p in enumerate(points_coordinates): 
                    if k not in points_nonexp:
                        coord_point_atom1 = [p[0], p[1], p[2]] 
                        distance_patom1_atom2 = math.dist(coord_point_atom1, coord_atom2)
                        if tl.check_distance(r_vdw_atom2, distance_patom1_atom2):  #Check if dist < threshold.
                            points_nonexp.append(k)  #Add the hidden/non exposed points.
                                              
        s_exp_atom, s_tot_atom = tl.surface_exp(r_vdw_atom1, pointsNum-len(points_nonexp), pointsNum)
        s_exp_prot += s_exp_atom
        s_tot_prot += s_tot_atom
        
        if res_count != int(num_res):  #If the residue changes, reinitialize the variables.      
            s_exp_abs_res = 0
            s_tot_res = 0
            res_count = int(num_res)
            pos += 1

        s_exp_abs_res += s_exp_atom 
        s_tot_res += s_tot_atom
        s_exp_rel_res = s_exp_abs_res/SASA_MAX[name_res]
        perc_res = (s_exp_abs_res/s_tot_res)*100
                    
        dico[pos] = (name_res, num_res, round(s_exp_abs_res, 2), round(s_exp_rel_res, 2), round(perc_res, 2))  #In homodimers, residue names and numbers are the same in both chains. For this reason, pos (position) is used as key.
    
    headers = ["POS", "RES", "NUM", "S-ABS", "S-REL", "PERC"]
    print(tabulate([(k,) + v for k, v in dico.items()], headers=headers))

    print(f"{round(s_exp_prot, 2)} ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)")
    print(f"{round((s_exp_prot/s_tot_prot)*100, 2)} PERCENTAGE OF ACCESSIBILITY")
    
    print(f"--- {(time.time() - start_time):.2f} seconds ---")  #Time taken for the program execution.

#NB: The main has no return                   
#Example of use : python3 SASA.py -pdb data/pdb/1bzv.pdb -n 92
                      
                                       
if __name__ == "__main__":  #If the program SASA.py is executed as a script in a shell, the result of the if test will be True and the corresponding instruction block will be executed.
    SASA()                  #Otherwise, if the program SASA.py is imported as a module, then the result of the test if will be False (and the corresponding instruction block will not be executed).