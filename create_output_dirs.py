import os

def create_output_dirs(top_level_path):
    # Define the subdirectories to be created
    subdirs = ['summa', 'openwq']
    subsubdirs = ['openwq_2', 'openwq_4', 'openwq_6', 'openwq_8', 'openwq_9',
                  'openwq_10', 'openwq_11', 'openwq_11_1', 'openwq_12',
                  'openwq_13']
    
    path = os.path.join(top_level_path, 'output')
    try:
        os.makedirs(path, exist_ok=True)
        print(f"Directory '{path}' created successfully.")
    except OSError as error:
        print(f"Error creating directory '{path}': {error}")
    
    top_level_path = path

    # Create each subdirectory
    for subdir in subdirs:
        path = os.path.join(top_level_path, subdir)
        try:
            os.makedirs(path, exist_ok=True)
            print(f"Directory '{path}' created successfully.")
        except OSError as error:
            print(f"Error creating directory '{path}': {error}")
        for subsubdir in subsubdirs:
            path = os.path.join(top_level_path, subdir, subsubdir)
            try:
                os.makedirs(path, exist_ok=True)
                print(f"Directory '{path}' created successfully.")
            except OSError as error:
                print(f"Error creating directory '{path}': {error}")

if __name__ == "__main__":
    # Take the top-level path as input
    top_level_path = input("Enter the top-level path: ")
    create_output_dirs(top_level_path)