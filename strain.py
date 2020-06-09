import sys
#from importlib import reload
#reload(sys)
#sys.setdefaultencoding("ISO-8859-1")
import os
import numpy as np
from numpy import linalg as LA
from array import *

def array_list(array_num): 
    num_list = array_num.tolist() # list 
    #print(num_list)

def format(value):
    return "%.3f" % value

def extractArrayInSameLine(filename,beginReadingString,linesToRead,columnsToRead,startColumnIndex):
    extracted_arr = []
    with open(filename, "r") as f:
        beginReading = False
        read_counter = 0
        for line in f:
            line = line.strip()

            if read_counter == linesToRead:
                beginReading = False
                break

            if beginReading==False:
                line_no_spaces = "".join(line.split())

                if line_no_spaces.startswith(beginReadingString):
                    beginReading = True

            if beginReading == True:
                splits = line.split()
                for i in range(columnsToRead):
                    extracted_arr.append(float(splits[startColumnIndex + i]))
                read_counter += 1

    return extracted_arr;


def extractArray(filename,beginReadingString,linesToRead,columnsToRead,startColumnIndex,startRowIndex,includeAdditionalData):
    extracted_arr = [[]]
    for i in range(linesToRead-1):
        extracted_arr.append([])

    with open(filename, "r") as f:
        beginReading = False
        read_counter = 0
        for line in f:
            line = line.strip()

            if read_counter == linesToRead+startRowIndex:
                beginReading = False
                break

            if beginReading == True:
                if read_counter>=startRowIndex:

                    splits = line.split()
                    if includeAdditionalData==True:
                        extracted_arr[read_counter-startRowIndex].append(splits[1]);
                    for i in range(columnsToRead):
                        extracted_arr[read_counter-startRowIndex].append(float(splits[startColumnIndex + i]))
                read_counter += 1
            else:
                line_no_spaces = "".join(line.split())

                if line_no_spaces == beginReadingString:
                    beginReading = True
    return extracted_arr;





def get_all_coord(filename):
    """@returns all_coords
    """
    
    rcut=10.0 # cutoff radius of neighbors
    natoms0=extractArrayInSameLine(filename,"numberofatoms/cell",1,1,4) # how many atoms per simualtion 
    natoms=int(natoms0[0])
    #fractional cooridates
    native_frac_1=extractArray(filename, 
                    "siten.atompositions(cryst.coord.)",
                    natoms,3,6,0, True
                 );
    native_frac_2=np.asarray(native_frac_1)
    native_frac=np.transpose(native_frac_2)

    symbols=native_frac[0,:]
    native_frac_pos=native_frac[1:4,:]
    native_frac_pos=native_frac_pos.astype('float64');



    #box-vector 
    crystal_1=extractArray(filename, "".join("crystal axes: (cart. coord. in units of alat)".split()),3,3,3,0, False);
    crystal=np.asarray(crystal_1)
    
    lattice_param_EXTRACT=extractArrayInSameLine(filename,"latticeparameter(alat)=",1,1,4)

    lattice_param=lattice_param_EXTRACT[0]
    lattice_param=lattice_param*0.529177249
    
    
    boxVec=lattice_param*crystal

    avec=boxVec[0,:]
    bvec=boxVec[1,:]
    cvec=boxVec[2,:]


    #position=lattice_param*np.matmul(crystal,native_frac_pos)
    position=lattice_param*np.matmul(np.transpose(native_frac_pos),crystal)
    position = np.transpose(position)
    periods=3 # how many periodic images we want to replicate
    # if we need 5 periodic images than range will change (-2,3) : -2,-1,0,+1,+2
    atomPerCell=natoms
    cellPerDim=3
    totalCell=periods**3 

    totalAtom=totalCell*atomPerCell

    replicated_atom_names=[];

    superCell_position=[]
    replicated_atom_parent_index=[];
    id=0

    superCell_position=np.zeros((3,totalAtom))

    for na in range(-1,2):
        for nb in range(-1,2):
            for nc in range(-1,2):
                for i in range(atomPerCell):
                    #print 'Hello'
                    temp_pos=position[:,i]+na*avec+nb*bvec+nc*cvec
                    superCell_position[0,id]=temp_pos[0]
                    superCell_position[1,id]=temp_pos[1]
                    superCell_position[2,id]=temp_pos[2]

                    replicated_atom_names.append(symbols[i]);
                    replicated_atom_parent_index.append(i);
                   # print 'C', superCell_position[:,id] 
                    id=id+1
   

    unique_atoms=list(set(symbols))

#    force_1=extractArray(filename, "Forcesactingonatoms(cartesianaxes,Ry/au):", natoms,3,6,1, False);
#    force=np.asarray(force_1)
#    conversion=25.7110
#    force_ev_angstrom=conversion*force


    total_1=extractArrayInSameLine(filename,"!totalenergy=",1,1,4)
    total=np.asarray(total_1)
    e=(total[0]/ natoms)*13.6056080659


    native=lattice_param*np.matmul(crystal,native_frac_pos)
    
    
    id=0
    all_neigh=[]
    atom_combos=[];
    all_types=[];
    original_types=[];
    all_neigh_index=[];

    
    for i in range(atomPerCell):

        if symbols[i] == 'Al':
            original_types.append(13); # atomic number of central atom
        else:
            original_types.append(30); # atomic number of neighbor atom

        neigh = []
        indx=[];
        types = []

        for j in range(totalAtom):
            relDist = []
            for k in range(3):
                relDist.append(position[k, i] - superCell_position[k, j])
            inorm = LA.norm(relDist)
            # print inorm
            if (inorm < rcut and inorm > 0.00001):
                neigh.append(relDist)
                indx.append(replicated_atom_parent_index[j]);
                if replicated_atom_names[j] == 'Al':
                    types.append(13);
                else:
                    types.append(30);

        all_neigh.append(neigh)
        all_neigh_index.append(indx)
        all_types.append(types);

    for i in range(len(all_neigh)):
       for j in range(len(all_neigh[i])-1):
          for k in range(j+1,len(all_neigh[i])):
             if(LA.norm(all_neigh[i][j])>LA.norm(all_neigh[i][k])):
                temp=all_neigh[i][j]
                all_neigh[i][j]=all_neigh[i][k]
                all_neigh[i][k]=temp

                temp=all_types[i][j];
                all_types[i][j]=all_types[i][k];
                all_types[i][k]=temp;

                temp = all_neigh_index[i][j];
                all_neigh_index[i][j] = all_neigh_index[i][k];
                all_neigh_index[i][k] = temp;


    length=[]
    all_coord=[]
    all_name=[]


    #flatten list
    for i in range(len(all_neigh)):
        icoords=[]
        names=[];
        for j in range(len(all_neigh[i])):
        #for j in range(30):
            for k in range(len(all_neigh[i][j])):
                icoords.append(all_neigh[i][j][k])
        icoords.extend(all_neigh_index[i])
        #icoords.extend(all_types[i]); # for alloy, print types of neighbors
        all_coord.append(icoords)



#    force_ev_angstrom_list=force_ev_angstrom.tolist()

    for i in range(len(all_coord)):
        all_coord[i].extend([e,original_types[i],natoms])
#        all_coord[i].extend([force_ev_angstrom_list[i],e,original_types[i],natoms])

    file = open("strain_fccEV.dat", "a")
    
    for i in range(len(all_neigh)):
        file.write(str(all_coord[i][:]).replace('[','').replace(']','').replace("'", "").replace(',','')+'\n')
    return all_coord

if __name__=="__main__":
# directory name where the data files are located
    dir_name="/MN/"
    input_dir="{0}{1}".format(os.getcwd(),dir_name )
    input_data=[f for f in os.listdir(input_dir) if ((os.path.isfile(os.path.join(input_dir,f))) and (os.stat(os.path.join(input_dir,f)).st_size > 0))]
    output=[]
    for fname in input_data:
        #print fname
        data="{0}{1}".format(input_dir,fname)
        output.extend(get_all_coord(data))
