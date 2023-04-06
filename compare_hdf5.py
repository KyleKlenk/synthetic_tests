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
    dir1 = '/code/case_studies/Output_OpenWQ_13/HDF5'
    dir2 = '/code/case_studies/Output_OpenWQ_13_Reference/HDF5'
    compare_hdf5_files(dir1, dir2)