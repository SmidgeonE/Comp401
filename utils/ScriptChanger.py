import os
import sys

def replace_line_content(file_path, line_number, new_content):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    if len(lines) >= line_number:
        parts = lines[line_number - 1].split('=')
        if len(parts) > 1:
            parts[-1] = new_content + '\n'
            lines[line_number - 1] = '='.join(parts)
    
    with open(file_path, 'w') as file:
        file.writelines(lines)

def main(arg1, arg2):
    directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../bc4/')
    for filename in os.listdir(directory):
        if filename.endswith('.sh'):
            file_path = os.path.join(directory, filename)
            replace_line_content(file_path, 9, arg1)
            replace_line_content(file_path, 11, arg2)
            replace_line_content(file_path, 20, arg2)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: scriptChanger.py <numNodes> <numThreads> -> This util changes all of the bash scripts for each job quickly.")
        sys.exit(1)
    
    arg1 = sys.argv[1]
    arg2 = sys.argv[2]
    main(arg1, arg2)