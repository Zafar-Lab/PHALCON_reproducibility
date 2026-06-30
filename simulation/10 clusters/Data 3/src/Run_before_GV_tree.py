import re
with open('sc_2000_inferred_tree.nw','r') as f:
        string = f.read()

print(string)
string = string[3:-2]
string = re.sub('c','',string)
print(string)
with open('inferred_format.nw','w') as f:
        for i in string:
                if i.isdigit():
                        i = int(i)-1
                f.write(str(i))
        f.write(";")
        
    