import re
string = "((c5,(c1,c6)),(((c3,c4),c2),c7))"
print(string)
string = re.sub('c','',string)
print(string)
with open('inferred_format.nw','w') as f:
        for i in string:
                if i.isdigit():
                        i = int(i)-1
                f.write(str(i))
        f.write(";")
        
    