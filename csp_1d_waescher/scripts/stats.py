p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
total=sum(l*d for l,d in items); dem=sum(d for _,d in items)
print('m',m,'L',L,'total_length',total,'ceil_LB',(total+L-1)//L,'dem',dem,'max_len',max(l for l,d in items))
