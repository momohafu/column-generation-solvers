p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
print('m L',m,L)
for frac in [2,3,4,5]:
    th=L/frac
    cnt=sum(d for l,d in items if l>th)
    print('>L/',frac,'th',th,'count',cnt,'min bins ceil(cnt/(frac-1))', (cnt+(frac-2))//(frac-1))
# pairs/triples
