import sys
from os import path,makedirs,_exit,system
from collections import Counter
import discover
import numpy as np
import pandas as pd
import argparse
from itertools import combinations
from itertools import product
sys.setrecursionlimit(100000)

parser = argparse.ArgumentParser()
parser.add_argument('-i','--maf',type=str,required=True)
parser.add_argument('-f','--freq',type=float,required=False,default=0.05)
parser.add_argument('-n','--num',type=int,required=False,default=30)
parser.add_argument('-k','--geneSetNum',type=int,required=False,default=2)
parser.add_argument('-o','--outdir',type=str,required=True)
parser.add_argument('-p1','--cop',type=float,required=False,default=0.05)
parser.add_argument('-p2','--mup',type=float,required=False,default=0.01)
parser.add_argument('-m','--mod',type=str,required=False,default='vcfilter',choices=['vcfilter','vcall'])
parser.add_argument('-fdr','--fdr',type=float,required=False,default=1)
args = parser.parse_args(sys.argv[1:])

mafs=args.maf
freq=args.freq
geneShownumber=args.num
outdir=args.outdir
mod=args.mod
fdr=args.fdr
k=args.geneSetNum
vc=["Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Translation_Start_Site","Nonsense_Mutation","Nonstop_Mutation","In_Frame_Del","In_Frame_Ins","Missense_Mutation"]

if not path.exists(outdir):
        makedirs(outdir)
if not path.exists(outdir+'/datafile'):
	makedirs(outdir+'/datafile')
if not path.exists(outdir+'/discover'):
        makedirs(outdir+'/discover')
if not path.exists(outdir+'/wext'):
        makedirs(outdir+'/wext')
if not path.exists(outdir+'/maftools'):
        makedirs(outdir+'/maftools')

def clean(lst):
	outlst=lst[:]
	for i in lst:
		if len(i) > 2:
			for j in range(2,len(i)):
				for c in combinations(i,j):
					c=sorted(list(c))
					c=list(map(lambda x:x.strip(),c))
					if c in outlst:
						del outlst[outlst.index(c)]
	return outlst

def detect_comb(dic,comb):
	for g1,g2 in product(comb,comb):
		if g1 != g2:
			if g1 not in dic[g2]:
				return False
	return True

def cluster(lst):
	relationship={}
	for i in lst:
		for gene1,gene2 in product(i,i):
			if gene1 != gene2:
				gene1=gene1.strip();gene2=gene2.strip()
				relationship.setdefault(gene1,[]).append(gene2)
	combinationList=[]
	for gene in relationship.keys():
		netgene=relationship[gene]
		netgene.append(gene)
		for i in range(2,(len(netgene)+1)):
			for c in combinations(netgene,i):
				genes=sorted(list(c))
				genes=list(map(lambda x:x.strip(),genes))
				if genes not in combinationList:
					combinationList.append(genes)
	outlist=[]
	for comb in combinationList:
		if detect_comb(relationship,comb):
			outlist.append(comb)
	out=clean(outlist)
	return(out)

def myfilter(f,p_threshold,mod):
	group=[]
	out=[]
	flag=0
	for line in open(f,'r'):
		line=line.strip()
		if flag == 0:
			flag += 1
			out.append(line)
			continue
		if mod == 'discover':
			gene1,gene2,p=line.split('\t')[1:4]
			gene1=gene1.strip()
			gene2=gene2.strip()
			genes=gene1+','+gene2
		elif mod == 'wext':
			genes,p=line.split('\t')[0:2]
		if float(p) <= float(p_threshold):
			gene_out=[]
			for g in genes.split(','):
				g=g.strip()
				gene_out.append('\''+str(g)+'\'')
			group.append(gene_out)
			out.append(line)
	return group,out

def detect_overlap(g1,g2):
	out=[]
	flag=0
	for g in g1:
		if g in g2:
			flag=1
			break
	if flag == 1:
		out=g1[:]
		for g in g2:
			if g not in out:
				out.append(g)
	return out

def mygroup(group):
	out=[]
	if len(group) == 1:
		return group
	overlap=0
	flag=0
	for i1,g1 in enumerate(group):
		for g2 in group[i1+1:]:
			comba=detect_overlap(g1,g2)
			if len(comba) > 0:
				overlap += 1
				out.append(comba)
				for g in group:
					if (g is not g1) and (g is not g2):
						out.append(g)
						flag=1
			if flag == 1:
				break
		if flag == 1:
			break
	if overlap == 0:
		return group
	else:
		return mygroup(out)

allmaf=open(outdir+'/datafile/allSNV_INDEL.maf','w')
gene=[]
info={}
flag=0
for maf in open(mafs,'r').readlines():
	maf=maf.strip()
	tmp=[]
	for i in open(maf,'r',encoding='ISO-8859-1').readlines():
		if i.startswith('#'):
			continue
		if i.startswith('Hugo_Symbol') and (flag == 0):
			allmaf.writelines(i)
			flag = 1
			continue
		if i.startswith('Hugo_Symbol') and (flag == 1):
			continue
		if mod == 'vcfilter':
			if i.split('\t')[8] not in vc:
				continue
		allmaf.writelines(i)
		g=i.split('\t')[0]
		samplename=i.split('\t')[15]
		if info.get(samplename):
			if g not in info[samplename]:
				info[samplename].append(g)
		else:
			info.setdefault(samplename,[g])
		if g not in tmp:
			tmp.append(g)
	gene=gene+tmp
allmaf.close()

count_result=Counter(gene)
genelen=len(info)
gene=[]
if (1/float(genelen)) >= freq:
	print("frequence of one mutation is higher than the setting %.2f" % freq)
	print("this will lead to error, so we break this pipeline!")
	_exit(0)
for i,n in count_result.items():
	if (n/float(genelen)) >= freq:
		gene.append(i)

wext=open(outdir+'/datafile/wext.input','w')
for samp,genelst in info.items():
	lst=[samp]
	for g in gene:
		if g in genelst:
			lst.append(g)
	if len(lst) > 1:
		wext.writelines('\t'.join(lst)+'\n')
wext.close()

dis=open(outdir+'/datafile/discover.input','w')
samp=info.keys()
dis.writelines("\t"+'\t'.join(samp)+'\n')
for g in gene:
	lst=[g]
	for s in samp:
		if g in info[s]:
			lst.append('1')
		else:
			lst.append('0')
	dis.writelines('\t'.join(lst)+'\n')
dis.close()

#######DISCOVER RUNNING#################
arr=np.loadtxt(outdir+"/datafile/discover.input",str,delimiter='\t')
mut=pd.DataFrame(arr[1:,1:].astype(np.int32),index=arr[1:,0],columns=arr[0,1:])
events = discover.DiscoverMatrix(mut)
result_mutex = discover.pairwise_discover_test(events,alternative='less',correct=True)
result_co = discover.pairwise_discover_test(events,alternative='greater',correct=True)
result_mutex.significant_pairs().to_csv(outdir+'/discover/discover_mutex.txt',sep='\t')
result_co.significant_pairs().to_csv(outdir+'/discover/discover_co.txt',sep='\t')
########################################

#######wext RUNNING#####################
cmd1="/hwfssz4/BC_COM_FP/bc_tumor/Pipeline/Cancer_Genomics_Interpretation_Analysis/CGIA_v2.0/bin/python /zfssz5/BC_PS/chengyuanfang/software/wext-master/process_mutations.py -m %s/datafile/wext.input -ct NA -o %s/wext/data.json" % (outdir,outdir)
cmd2="/hwfssz4/BC_COM_FP/bc_tumor/Pipeline/Cancer_Genomics_Interpretation_Analysis/CGIA_v2.0/bin/python /zfssz5/BC_PS/chengyuanfang/software/wext-master/compute_mutation_probabilities.py -mf %s/wext/data.json -np 1000 -nc 4 -wf %s/wext/weights.npy -v 1" % (outdir,outdir)
cmd3="/hwfssz4/BC_COM_FP/bc_tumor/Pipeline/Cancer_Genomics_Interpretation_Analysis/CGIA_v2.0/bin/python /zfssz5/BC_PS/chengyuanfang/software/wext-master/find_sets.py -mf %s/wext/data.json -wf %s/wext/weights.npy -s exclusivity -fdr %f -k %d -c 4 -f 2 -o %s/wext/exclusivity_results -v 0" % (outdir,outdir,fdr,k,outdir)
cmd5="/hwfssz4/BC_COM_FP/bc_tumor/Pipeline/Cancer_Genomics_Interpretation_Analysis/CGIA_v2.0/bin/python /zfssz5/BC_PS/chengyuanfang/software/wext-master/find_sets.py -mf %s/wext/data.json -wf %s/wext/weights.npy -s all-co-occurrence -fdr %f -k %d -c 4 -f 2 -o %s/wext/all-co-occurrence_results -v 0" % (outdir,outdir,fdr,k,outdir)
system(cmd1)
system(cmd2)
system(cmd3)
system(cmd5)
##########################################


###############process outfile###########
group1,outline1=myfilter(outdir+'/wext/exclusivity_results-k'+str(k)+'.tsv',args.mup,'wext')
#result1=mygroup(group1)
result1=cluster(group1)
group2,outline2=myfilter(outdir+'/wext/all-co-occurrence_results-k'+str(k)+'.tsv',args.cop,'wext')
#result2=mygroup(group2)
result2=cluster(group2)
f1=open(outdir+'/exclusivity.filtered.txt','w')
f2=open(outdir+'/cooccurrence.filtered.txt','w')
f1.writelines('\n'.join(outline1))
#f1.writelines('\ncluster these gene:\n')
f2.writelines('\n'.join(outline2))
#f2.writelines('\ncluster these gene:\n')
str1=''
str2=''
for i,lst in enumerate(result1):
#	f1.writelines(','.join(lst)+'\n')
	str1 +='''
		pdf(paste(outdir,'/exclusivity_oncostrip%d.pdf',sep=''),width=10,height=10)
		oncostrip(maf = f, genes = c(%s), removeNonMutated = TRUE, colors=col)
		dev.off()
		''' % (i,','.join(lst))
for i,lst in enumerate(result2):
#	f2.writelines(','.join(lst)+'\n')
	str2 +='''
		pdf(paste(outdir,'/co_oncostrip%d.pdf',sep=''),width=10,height=10)
		oncostrip(maf = f, genes = c(%s), removeNonMutated = TRUE, colors=col)
		dev.off()
		''' % (i,','.join(lst))
f1.close()
f2.close()
##########################################

##############maftools RUNNING############
maftools=open(outdir+'/maftools/runMaftools.r','w')
outstr='''
library(maftools)
library(ulimit)
ulimit::memory_limit(4000)

args=commandArgs(T)
outdir=args[2]
'''
if mod == 'vcall':
	outstr +='''
f<-read.maf(args[1],removeDuplicatedVariants = FALSE, vc_nonSyn=c("De_novo_Start_InFrame","De_novo_Start_OutOfFrame","lincRNA","IGR","Intron","RNA","Silent","3'UTR","5'Flank","5'UTR","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Site","Start_Codon_Ins","Start_Codon_Del","Start_Codon_SNP","Stop_Codon_Del","Stop_Codon_Ins","Translation_Start_Site"))

col = c(RColorBrewer::brewer.pal(8,name = "Dark2"), RColorBrewer::brewer.pal(8,name = "Accent"),RColorBrewer::brewer.pal(8,name = "Set2"),'black', 'violet', 'royalblue')
col = grDevices::adjustcolor(col = col, alpha.f = 1)
names(col) = names = c("De_novo_Start_InFrame","De_novo_Start_OutOfFrame","lincRNA","IGR","Intron","RNA","Silent","3'UTR","5'Flank","5'UTR","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Site","Start_Codon_Ins","Start_Codon_Del","Start_Codon_SNP","Stop_Codon_Del","Stop_Codon_Ins","Translation_Start_Site","Multi_Hit",'Amp','Del')
'''
else:
	outstr +='''
f<-read.maf(args[1],removeDuplicatedVariants = FALSE)
col = c(RColorBrewer::brewer.pal(12,name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'violet', 'royalblue')
col = grDevices::adjustcolor(col = col, alpha.f = 1)
names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','IGR','Missense_Mutation','Silent','Nonsense_Mutation','RNA','Splice_Site','Intron','Frame_Shift_Ins','Nonstop_Mutation','In_Frame_Del','ITD','In_Frame_Ins','Translation_Start_Site',"Multi_Hit", 'Amp', 'Del')
col
'''

outstr +='''
#pdf(paste(outdir,'/somaticInteractions.pdf',sep=''),width=10,height=10)
#somaticInteractions(maf = f, top = %s, pvalue = c(0.05, 0.1))
#dev.off()

pdf(paste(outdir,'/oncoplot.pdf',sep=''),width=10,height=10)
oncoplot(maf = f, top = %s, colors=col)
dev.off()

%s
%s
''' % (geneShownumber,geneShownumber,str1,str2)
maftools.writelines(outstr)
maftools.close()
system("export R_LIBS=/zfssz5/BC_PS/chengyuanfang/lib/R:/hwfssz4/BC_COM_FP/bc_tumor/Pipeline/Cancer_Genomics_Interpretation_Analysis/CGIA_v2.0/lib/R_LIBS:/zfssz5/BC_PS/chengyuanfang/lib/R:$R_LIBS")
cmd6="/hwfssz4/BC_COM_FP/bc_tumor/Pipeline/Cancer_Genomics_Interpretation_Analysis/CGIA_v2.0/bin/Rscript %s/maftools/runMaftools.r %s/datafile/allSNV_INDEL.maf %s" % (outdir,outdir,outdir)
system(cmd6)
