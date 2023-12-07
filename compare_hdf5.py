import os
import h5py

def are_hdf5_files_equal(file1, file2):
    with h5py.File(file1, 'r') as f1:
        with h5py.File(file2, 'r') as f2:
            if f1.keys() != f2.keys():
                return False
            for key in f1.keys():
                if isinstance(f1[key], h5py.Dataset):
                    if not isinstance(f2[key], h5py.Dataset):
                        return False
                    if not (f1[key].shape == f2[key].shape and
                            f1[key].dtype == f2[key].dtype and
                            (f1[key][()] == f2[key][()]).all()):
                        return False
                elif isinstance(f1[key], h5py.Group):
                    if not isinstance(f2[key], h5py.Group):
                        return False
                    if not are_hdf5_files_equal(f1[key], f2[key]):
                        return False
                else:
                    raise ValueError('Unknown HDF5 object type')
            return True

def compare_hdf5_files(dir1, dir2):
    list_of_files_1 = sorted(os.listdir(dir1))
    list_of_files_2 = sorted(os.listdir(dir2))
    if len(list_of_files_1) != len(list_of_files_2):
      raise Exception('Number of files in directories are not equal')

    for i in range(len(list_of_files_1)):
      if list_of_files_1[i].endswith('.h5') and list_of_files_2[i].endswith('.h5'):
        file1_path = os.path.join(dir1, list_of_files_1[i])
        file2_path = os.path.join(dir2, list_of_files_2[i])
        print('Comparing:')
        print(f'  {file1_path}')
        print(f'  {file2_path}')
        if are_hdf5_files_equal(file1_path, file2_path):
            print('Files are equal\n')
        else:
            raise Exception('Files are not equal')

# Example usage:
if __name__ == '__main__':
    # Compare 2_nrTrans_instS_PorMedia
    dir1 = '/home/kck540/OpenWQ-Projects/Chris_Files/Summa-openWQ/bin/Output_OpenWQ_2_comp/HDF5'
    dir2 = '/home/kck540/OpenWQ-Projects/OpenWQ_Reference_Setup/Summa-openWQ/bin/Output_OpenWQ_2_Ref/HDF5'
    compare_hdf5_files(dir1, dir2)
    
    # Compare 4_nrTrans_contS_PorMedia
    dir1 = '/home/kck540/OpenWQ-Projects/Chris_Files/Summa-openWQ/bin/Output_OpenWQ_4_comp/HDF5'
    dir2 = '/home/kck540/OpenWQ-Projects/OpenWQ_Reference_Setup/Summa-openWQ/bin/Output_OpenWQ_4_Ref/HDF5'
    compare_hdf5_files(dir1, dir2)

    # Compare 6_nrTrans_instS_PorMedia_linDecay
    dir1 = '/home/kck540/OpenWQ-Projects/Chris_Files/Summa-openWQ/bin/Output_OpenWQ_6_Comp/HDF5'
    dir2 = '/home/kck540/OpenWQ-Projects/OpenWQ_Reference_Setup/Summa-openWQ/bin/Output_OpenWQ_6_Ref/HDF5'
    compare_hdf5_files(dir1, dir2)
    
    # Compare 8_nrTrans_contS_PorMedia_linDecay
    dir1 = '/home/kck540/OpenWQ-Projects/Chris_Files/Summa-openWQ/bin/Output_OpenWQ_8_Comp/HDF5'
    dir2 = '/home/kck540/OpenWQ-Projects/OpenWQ_Reference_Setup/Summa-openWQ/bin/Output_OpenWQ_8_Ref/HDF5'
    compare_hdf5_files(dir1, dir2)

    # Compare 9_batch_singleSp_1storder
    dir1 = '/home/kck540/OpenWQ-Projects/Chris_Files/Summa-openWQ/bin/Output_OpenWQ_9_Comp/HDF5'
    dir2 = '/home/kck540/OpenWQ-Projects/OpenWQ_Reference_Setup/Summa-openWQ/bin/Output_OpenWQ_9_Ref/HDF5'
    compare_hdf5_files(dir1, dir2)

    # Compare 10_batch_singleSp_2ndorder
    dir1 = '/home/kck540/OpenWQ-Projects/Chris_Files/Summa-openWQ/bin/Output_OpenWQ_10_Comp/HDF5'
    dir2 = '/home/kck540/OpenWQ-Projects/OpenWQ_Reference_Setup/Summa-openWQ/bin/Output_OpenWQ_10_Ref/HDF5'
    compare_hdf5_files(dir1, dir2)
    
    # Compare 11_1_batch_3species
    dir1 = '/home/kck540/OpenWQ-Projects/Chris_Files/Summa-openWQ/bin/Output_OpenWQ_11_1_Comp/HDF5'
    dir2 = '/home/kck540/OpenWQ-Projects/OpenWQ_Reference_Setup/Summa-openWQ/bin/Output_OpenWQ_11_1_Ref/HDF5'
    compare_hdf5_files(dir1, dir2)

    # Compare 11_batch_2species
    dir1 = '/home/kck540/OpenWQ-Projects/Chris_Files/Summa-openWQ/bin/Output_OpenWQ_11_Comp/HDF5'
    dir2 = '/home/kck540/OpenWQ-Projects/OpenWQ_Reference_Setup/Summa-openWQ/bin/Output_OpenWQ_11_Ref/HDF5'
    compare_hdf5_files(dir1, dir2)

    # Compare 12_batch_nitrogencycle
    dir1 = '/home/kck540/OpenWQ-Projects/Chris_Files/Summa-openWQ/bin/Output_OpenWQ_12_Comp/HDF5'
    dir2 = '/home/kck540/OpenWQ-Projects/OpenWQ_Reference_Setup/Summa-openWQ/bin/Output_OpenWQ_12_Ref/HDF5'
    compare_hdf5_files(dir1, dir2)
    
    # Compare 13_batch_oxygenBODcycle
    dir1 = '/home/kck540/OpenWQ-Projects/Chris_Files/Summa-openWQ/bin/Output_OpenWQ_13_Comp/HDF5'
    dir2 = '/home/kck540/OpenWQ-Projects/OpenWQ_Reference_Setup/Summa-openWQ/bin/Output_OpenWQ_13_Ref/HDF5'
    compare_hdf5_files(dir1, dir2)