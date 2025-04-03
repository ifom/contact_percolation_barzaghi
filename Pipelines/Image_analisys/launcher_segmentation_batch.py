import time
import os
import sys
import yaml

def process_image(filename, params):
    mixed_pop = params.get('mixed_pop')
    print('Processing:', filename)
    print('mixed_pop = ' + str(mixed_pop), flush=True)

    if mixed_pop:
        import mixedsegmentation
        mixedsegmentation.mixedsegmentation(filename, params)
    else:
        import segmentation
        segmentation.segmentation(filename, params)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_yaml>")
        sys.exit(1)

    yaml_file = sys.argv[1]

    # Read the YAML file
    try:
        with open(yaml_file, 'r') as f:
            params = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"File {yaml_file} not found.")
        sys.exit(1)
    except yaml.YAMLError as e:
        print("Error parsing the YAML file:", e)
        sys.exit(1)
        
    folder_path = params.get('folder_path')
    if not folder_path:
        print("Error: 'folder_path' not found in the YAML file.")
        sys.exit(1)
        
    # Process each image in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith(('.tif', '.tiff')):
            file_path = os.path.join(folder_path, filename)
            process_image(file_path, params)

