p='/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt'
lines=open(p).read().splitlines()
N,T=map(int,lines[0].split())
out=[0]*(N+1); inc=[0]*(N+1)
for line in lines[1+N:]:
    i,j,c=map(int,line.split()); out[i]+=1; inc[j]+=1
noinc=[i for i in range(1,N+1) if inc[i]==0]
noout=[i for i in range(1,N+1) if out[i]==0]
print('no incoming (must start path)',len(noinc),noinc)
print('no outgoing (must end path)',len(noout),noout)
print('max of those',max(len(noinc),len(noout)))
