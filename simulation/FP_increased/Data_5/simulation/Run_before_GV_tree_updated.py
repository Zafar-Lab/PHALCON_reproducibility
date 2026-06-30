# THIS can handle all cases. two digit number, three digit number
# Fool-proof file
import re
benchmark = 'fp_inc'
with open('data_'+benchmark+'_2000_10_inferred_tree.nw','r') as f:
        string = f.read()

print(string)
string = re.sub('C','',string)
print(string)
with open('inferred_format.nw','w') as f:
        f.write(string)


def decrease_numbers(text):
    def replace(match):
        num = match.group(0)
        return str(int(num) - 1)

    pattern = r'\b\d+\b'
    return re.sub(pattern, replace, text)

# Read the text file
file_path = 'inferred_format.nw'  # Update this with the path to your text file
with open(file_path, 'r') as file:
    original_text = file.read()

# Decrease numbers in the text
modified_text = decrease_numbers(original_text)

# Write the modified text back to the file
with open(file_path, 'w') as file:
    file.write(modified_text)

with open('inferred_format.nw','r') as f:
        content = f.read()
        print(content)
print("Done")

        
    