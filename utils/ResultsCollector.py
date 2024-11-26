import os
import sys

def process_files():
    result = {}
    current_directory = os.getcwd()
    for filename in os.listdir(current_directory):
        if filename.endswith(".out"):
            with open(filename, 'r') as file:
                lines = file.readlines()
                if len(lines) < 3:
                    continue
                fourth_last_line = lines[-4]
                last_word = int(fourth_last_line.split()[-1])
                
                last_line = lines[-1]
                time_part = float(last_line.split()[-2])
                
                if last_word not in result:
                    result[last_word] = [time_part]
                else: 
                    result[last_word] += [time_part]
    return result

if __name__ == "__main__":
    if "-h" in sys.argv:
        print("This util collects the data from the output files and prints it in a readable format.")
        print("Note: It assumes the output files are in the parent directory.")
        sys.exit(0)

    collected_data = process_files()
    for key in sorted(collected_data.keys()):
        print(f"{key}: {collected_data[key]}\n")
