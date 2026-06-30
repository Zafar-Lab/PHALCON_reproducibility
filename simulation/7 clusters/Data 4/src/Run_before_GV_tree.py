import re
string = "((((c4,(c6,c3)),c7),(c2,c5)),c1)"
print(string)
string = re.sub('c','',string)
print(string)
with open('inferred_format.nw','w') as f:
        for i in string:
                if i.isdigit():
                        i = int(i)-1
                f.write(str(i))
        f.write(";")
        
    