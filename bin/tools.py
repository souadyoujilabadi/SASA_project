"""
Module containing the functions used by the program SASA aimed to calculate 
the Solvent-Accessible Surface Area (SASA) of proteins: 

    - atomic_coordinates() 
    - sphere() 
    - check_distance() 
    - surface_exp() 
"""


#Required python modules: 
import numpy as np
import pandas as pd
#from matplotlib import pyplot as pp  #If we want to visualize the sphere.


#Functions used by the program SASA: 
def atomic_coordinates(pdb):
    """
    Function that reads the .PDB file of a protein to extract its atomic three-dimensional coordinates.
    It returns as a result a DataFrame containing the name and 3D coordinates associated to each atom, 
    as well as the residue name and number.
        Arguments: 
            pdb: .PDB file.
        Return: 
            df: DataFrame containing the atomic coordinates and the other mentioned associated data.
    """
    with open(pdb, "r") as pdb_file:  #Open the PDB file in read "r" mode.                                                                  
        atom_data = {"res_name" : [], "res_num" : [], "atom_name" : [], "x" : [], "y" : [], "z" : []}  #Declare a dictionary, {key : value}, in which values = empty lists.
        atoms_names, res_names, res_nums, x, y, z = [], [], [], [], [], []  #Create empty lists. 
        #count_atoms = 0 #If we want to count the number of atoms. 
        for line in pdb_file:
            if line.startswith("ATOM"):  #Take only the lines starting with ATOM. These lines contain the protein atoms. (It is necessary to exclude the other lines, e.g. the lines starting with HETATM (= correspond to atoms that belong to other molecules, e.g. ions and H2O).                                                               
                atoms_names.append(line[13].strip())  #Take only the first letter of an atom name (e.g. C instead of CA, CB, ...), remove any leading (spaces at the beginning) and trailing (spaces at the end) characters with the method .strip(), and add it to the list atoms_names using the method .append().                                               
                res_names.append(line[17:20].strip())  #Extracte the residue name.
                res_nums.append((line[22:26]))  #Extracte the residue number.
                x.append(float(line[30:38]))  #Extracte the x coordinate, and add it to the list called x.
                y.append(float(line[38:46]))
                z.append(float(line[46:54]))
                #count_atoms += 1    
        #"This protein contains {} atoms".format(count_atoms)
        atom_data["res_name"] = res_names  #Assign the list res_names to the key "res_name".
        atom_data["res_num"] = res_nums 
        atom_data["atom_name"] = atoms_names
        atom_data["x"], atom_data["y"], atom_data["z"] = x, y, z
        df = pd.DataFrame(atom_data)  #Create a DataFrame with the atom_data dictionary, (the keys become the column names).
    return df


def sphere(n, center, r_vdw):  
    """
    Function that evenly distributes points on a sphere using the Golden spiral algorithm.
    The sphere has as center the 3D coordinates of an atom and as size the vdw radius of this atom + the vdw radius of Oxygen.
        Arugments: 
            n: Number of points.
            center: Center of the sphere (3D coordinates of an atom).
            r_vdw: Van der Waals radius of an atom.
        Return: 
            points_coordinates: Array containing the points coordinates. 
    """
    ind = np.arange(n)  #Create an array (vector) object of n integers. 
    phi = np.pi * (3 - np.sqrt(5))  #Golden angle (= 2.39996 rad, 137,51 °)
    theta = phi * ind  #Golden angle increment/Theta angle 
    z = np.linspace(start=1 , stop=-1, num=n)  #Create the variable z with the function linspace() which returns an array of n evenly spaced numbers over the specified interval, [1, -1]. (=Z axis)  
    r = np.sqrt(1 - z**2)  #Circles radius 
    r_sphere = r_vdw + 1.4  #Sphere defined radius (= sphere size): Van der Waals radius of an atom + Van der Waals radius of Oxygen (= 1.4).
    #Points coordinates: 
    z = np.add(z * r_sphere, center[2])
    x = np.add(r * np.cos(theta) * r_sphere, center[0])
    y = np.add(r * np.sin(theta) * r_sphere, center[1]) 
    points_coordinates = np.zeros(shape=(n, 3))  #Create an array of given shape (= n sequences containing 3 zeros each).
    points_coordinates[:, 0] = x  #Replace the first zero of each sequence by the x coordinates.
    points_coordinates[:, 1] = y 
    points_coordinates[:, 2] = z
    #pp.figure().add_subplot(projection="3d").scatter(x, y, z)  #If we want to visualize the sphere.
    return points_coordinates#, pp.show()


def check_distance(r_vdw, distance):
    """
    Function that checks if the given distance is inferior to a fixed threshold.
        Arugments: 
            r_vdw: Van der Waals radius of an atom.
            distance: Distance between a point on sphere 1 to the center of sphere 2.
        Return: 
            distance < r_vdw + 1.4: True or False. 
    """
    return distance < (r_vdw + 1.4)  #True (= hidden points) or False (= solvent-exposed points).


def surface_exp(r_vdw, points_exp, n):
    """
    Function that converts the solvent-accessible surface area expressed in term of points to Å**2.
        Arugments: 
            r_vdw: Van der Waals radius of an atom.
            points_exp: Solvent-exposed points.
            n: Total number of points.
        Return: 
            s_exp: Sphere's exposed surface area.
            s_tot: Sphere's total surface area (exposed + hidden).
    """
    s_tot = 4 * np.pi * (r_vdw + 1.4)**2  #Surface area of a sphere.
    ratio_exp = points_exp / n  #Accessibility ratio. 
    s_exp = s_tot * ratio_exp  #Cross-multiplication.
    return s_exp, s_tot  