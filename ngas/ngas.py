from pylab import *
import subprocess

def lercsv(arq):
	dados=loadtxt('cromatografias.csv',skiprows=1,usecols=range(4,18),delimiter='\t')
	X=[]
	rho=[]
	for i in range(len(dados)):
		X.append(dados[i,:]/sum(dados[i,:]))
	return X	

def writeinfile(p,T,X):
	sep='   '
	template=open('ngas.template','r')
	infile=open('ngas.dat','w')
	infile.write(template.read())
	infile.write('\n\n\n\n**=-=-=Composition\n**REM\n*COMPOSITION   *PRIMARY')
	infile.write('\n'+str(X[0])+sep+str(X[1])+sep+str(X[2])+sep+str(X[3])+sep+str(X[4]))	
	infile.write('\n'+str(X[5])+sep+str(X[6])+sep+str(X[7])+sep+str(X[8])+sep+str(X[9]))	
	infile.write('\n'+str(X[10])+sep+str(X[11])+sep+str(X[12])+sep+str(X[13]))	
	infile.write('\n\n')
	infile.write('**=-=-=Single-phase Calculation\n*SINGLE\n*LABEL'+sep+'\'\'\n'+'*FEED *MIXED 1.0\n*PRES '+str(p)+'\n*TEMP '+str(T)+'\n\n**=-=-=END')
	infile.close()
	

def density(p,T,X):
	writeinfile(p,T,X)
	
	winprop='C:\\Program Files (x86)\\CMG\\WINPROP\\2014.10\\Win_x64\\EXE\\pr201410.exe'
	infile='C:\\Users\\bfv8\\Documents\\github\\ngas\\ngas.dat'

#	subprocess.call([winprop,'-f',infile])
	p = subprocess.Popen([winprop,'-f',infile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()

	r=-1
	outfile=file('ngas.out')
	for line in outfile:
		if 'Density' in line:
			r=float(line.strip().split(' ')[-1])
	return r

####------------------------------------------------------------------------------------------------------

X=lercsv('cromatografias.csv')
p=10. # kgf/cm2
T=70. # oC

rho=[]
for i in range(len(X)):
	print str(i+1)+'/'+str(len(X))
	rho.append(density(p,T,X[i]))

print '\n\n\n\n Densidades, kg/m3 = ',rho
savetxt('rho.txt',array(rho))
