import re
string = "((c3,(c4,c5)),(c6,(c1,(c7,c2))))"
print(string)
string = re.sub('c','',string)
print(string)
with open('inferred_format.nw','w') as f:
        for i in string:
                if i.isdigit():
                        i = int(i)-1
                f.write(str(i))
        f.write(";")
        
    