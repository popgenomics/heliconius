# Python script to compute statistics from Hudson's ms output  
## mscalc.py  
**computes statistics from a written ms file, or a _fifo_ file.**  
mknod myfifo p  
mscalc.py < myfifo &  
priorgen_GeneralModel_beta.py bpfile=bpfile n1=0 n1=800 n1past=0 n1past=2 n2=0 n2=200 n2past=0 n2past=2 nA=0 nA=800 tau=0 tau=100 M1=0 M1=10 M2=0 M2=10 shape1=0 shape1=40 shape2=0 shape2=400 model={0} nreps=10000 Nvariation=hetero Mvariation={1} symMig=asym parameters=priorfile | msnsam tbs 50000000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs >myfifo &  
  
  
## mscalc_stdin.py  
**computes statistics from a pipe (_stdin_).**  
priorgen_GeneralModel_beta.py bpfile=bpfile n1=0 n1=800 n1past=0 n1past=2 n2=0 n2=200 n2past=0 n2past=2 nA=0 nA=800 tau=0 tau=100 M1=0 M1=10 M2=0 M2=10 shape1=0 shape1=40 shape2=0 shape2=400 model={0} nreps=100 Nvariation=hetero Mvariation={1} symMig=asym parameters=priorfile | msnsam tbs 50000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs | /home/cr600/scratch/ABC/mscalc_stdin.py  

