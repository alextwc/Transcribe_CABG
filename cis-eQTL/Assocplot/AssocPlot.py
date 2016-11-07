
import sys

#############################################################################################################
*************************************************************************************************************
In Jupyter Notebook Home, open a New Python3 shell ....
*************************************************************************************************************
#############################################################################################################

import os; import sys; import math; import types;
os.getcwd()
os.chdir("D:\\paths\\Meetings02\\2016-08-15")
ls
from assocplots.qqplot import * ;
from assocplots.manhattan import * ;
from assocplots.interactive import * ;
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
%matplotlib inline

Baseline_Log2_NoReg    = np.genfromtxt('MEQTL_NoReg.Baseline.Log2.imputed121.txt', dtype=None);
Baseline_Log2_RegOut   = np.genfromtxt('MEQTL_RegOut.Baseline.Log2.imputed121.txt', dtype=None);
Baseline_PEER10_NoReg  = np.genfromtxt('MEQTL_NoReg.Baseline.PEER10.imputed121.txt', dtype=None);
Baseline_PEER10_RegOut = np.genfromtxt('MEQTL_RegOut.Baseline.PEER10.imputed121.txt', dtype=None);
Ischemia_Log2_NoReg    = np.genfromtxt('MEQTL_NoReg.Ischemia.Log2.imputed121.txt', dtype=None);
Ischemia_Log2_RegOut   = np.genfromtxt('MEQTL_RegOut.Ischemia.Log2.imputed121.txt', dtype=None);
Ischemia_PEER10_NoReg  = np.genfromtxt('MEQTL_NoReg.Ischemia.PEER10.imputed121.txt', dtype=None);
Ischemia_PEER10_RegOut = np.genfromtxt('MEQTL_RegOut.Ischemia.PEER10.imputed121.txt', dtype=None);

lambda_Baseline_Log2_NoReg    = get_lambda(Baseline_Log2_NoReg['f3'], definition = 'median');
lambda_Baseline_Log2_RegOut   = get_lambda(Baseline_Log2_RegOut['f3'], definition = 'median');
lambda_Baseline_PEER10_NoReg  = get_lambda(Baseline_PEER10_NoReg['f3'], definition = 'median');
lambda_Baseline_PEER10_RegOut = get_lambda(Baseline_PEER10_RegOut['f3'], definition = 'median');
lambda_Ischemia_Log2_NoReg    = get_lambda(Ischemia_Log2_NoReg['f3'], definition = 'median');
lambda_Ischemia_Log2_RegOut   = get_lambda(Ischemia_Log2_RegOut['f3'], definition = 'median');
lambda_Ischemia_PEER10_NoReg  = get_lambda(Ischemia_PEER10_NoReg['f3'], definition = 'median');
lambda_Ischemia_PEER10_RegOut = get_lambda(Ischemia_PEER10_RegOut['f3'], definition = 'median');

LAMBDA=[lambda_Baseline_Log2_NoReg,lambda_Baseline_Log2_RegOut,lambda_Baseline_PEER10_NoReg,lambda_Baseline_PEER10_RegOut,lambda_Ischemia_Log2_NoReg,lambda_Ischemia_Log2_RegOut,lambda_Ischemia_PEER10_NoReg,lambda_Ischemia_PEER10_RegOut];
LAMBDA
[1.2755015584185816,
 1.3253828247438937,
 1.5075929863128841,
 1.4956062686095888,
 1.2630205643665071,
 1.3120597064714479,
 1.4892568965957846,
 1.4884128058714166]

mpl.rcParams['figure.dpi']=300
mpl.rcParams['savefig.dpi']=300
mpl.rcParams['figure.figsize']=5.375, 5.375
qqplot([Baseline_Log2_NoReg['f3'], Baseline_Log2_RegOut['f3']], 
       ['-RegOut(1.27)', '+RegOut(1.32)'], 
       color=['b','r'], 
       fill_dens=[0.2,0.2], 
       error_type='experimental', 
       distribution='beta',
       title='QQ Plot of Baseline_Log2 eQTLs')
plt.savefig('QQ2_Baseline_Log2.png', dpi=300)

mpl.rcParams['figure.dpi']=300
mpl.rcParams['savefig.dpi']=300
mpl.rcParams['figure.figsize']=5.375, 5.375
qqplot([Baseline_PEER10_NoReg['f3'], Baseline_PEER10_RegOut['f3']], 
       ['-RegOut(1.51)', '+RegOut(1.49)'], 
       color=['b','r'], 
       fill_dens=[0.2,0.2], 
       error_type='experimental', 
       distribution='beta',
       title='QQ Plot of Baseline_PEER10 eQTLs')
plt.savefig('QQ2_Baseline_PEER10.png', dpi=300)

mpl.rcParams['figure.dpi']=300
mpl.rcParams['savefig.dpi']=300
mpl.rcParams['figure.figsize']=5.375, 5.375
qqplot([Ischemia_Log2_NoReg['f3'], Ischemia_Log2_RegOut['f3']], 
       ['-RegOut(1.26)', '+RegOut(1.31)'], 
       color=['b','r'], 
       fill_dens=[0.2,0.2], 
       error_type='experimental', 
       distribution='beta',
       title='QQ Plot of Ischemia_Log2 eQTLs')
plt.savefig('QQ2_Ischemia_Log2.png', dpi=300)

mpl.rcParams['figure.dpi']=300
mpl.rcParams['savefig.dpi']=300
mpl.rcParams['figure.figsize']=5.375, 5.375
qqplot([Ischemia_PEER10_NoReg['f3'], Ischemia_PEER10_RegOut['f3']], 
       ['-RegOut(1.49)', '+RegOut(1.49)'], 
       color=['b','r'], 
       fill_dens=[0.2,0.2], 
       error_type='experimental', 
       distribution='beta',
       title='QQ Plot of Ischemia_PEER10 eQTLs')
plt.savefig('QQ2_Ischemia_PEER10.png', dpi=300)

mpl.rcParams['figure.dpi']=2000
mpl.rcParams['savefig.dpi']=2000
mpl.rcParams['figure.figsize']=17, 7
chrs = [str(i) for i in range(1,23)]
chrs_names = np.array([str(i) for i in range(1,23)])
chrs_names[1::2] = ''
cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in [0.0,0.33,0.67,0.90]]
manhattan(     Baseline_PEER10_RegOut['f3'], Baseline_PEER10_RegOut['f1'], Baseline_PEER10_RegOut['f0'].astype(str), 'Baseline_RegOut_PEER10',
               p2=Baseline_Log2_RegOut['f3'], pos2=Baseline_Log2_RegOut['f1'], chr2=Baseline_Log2_RegOut['f0'].astype(str), label2='Baseline_RegOut_Log2',
               type='inverted',
               chrs_plot=[str(i) for i in range(1,23)],
               chrs_names=chrs_names,
               cut = 0,
               title='Baseline eQTLs of TRANSCRiBE 121 patients',
               xlabel='chromosome',
               ylabel='-log10(p-value)',
               lines= [],
               top1 = 50,
               top2 = 50,
               colors = colors)
plt.savefig('ChicagoPlot_BaselineRegOut+-PEER10.png', dpi=2000)





Li00=[Baseline_Log2_NoReg,Baseline_Log2_RegOut,Baseline_PEER10_NoReg,Baseline_PEER10_RegOut,Ischemia_Log2_NoReg,Ischemia_Log2_RegOut,Ischemia_PEER10_NoReg,Ischemia_PEER10_RegOut];
for x in Li00:
    mpl.rcParams['figure.dpi']=300
    mpl.rcParams['savefig.dpi']=300
    mpl.rcParams['figure.figsize']=5.375, 5.375
    qqplot([x['f3']], 
           [str(x)], 
           color=['b'], 
           fill_dens=[0.2], 
           error_type='theoretical', 
           distribution='beta',
           title='QQ plot of "x" eQTLs')
    plt.savefig('qq_theoretical_error.png_str(x)', dpi=300)


get_lambda(State_Baseline['f3'], definition = 'median') #(Lambda=1.5253357472223925)
get_lambda(State_Ischemia['f3'], definition = 'median') #(Lambda=1.505432069571327)

mpl.rcParams['figure.dpi']=300
mpl.rcParams['savefig.dpi']=300
mpl.rcParams['figure.figsize']=5.375, 5.375
qqplot([State_Baseline['f3']], 
       ['State_Baseline'], 
       color=['b'], 
       fill_dens=[0.2], 
       error_type='theoretical', 
       distribution='beta',
       title='QQ plot of State_Baseline eQTLs')
plt.savefig('qq_StateBaseline_theoretical_error.png', dpi=300)

qqplot([State_Ischemia['f3']], 
       ['State_Ischemia'], 
       color=['b'], 
       fill_dens=[0.2], 
       error_type='theoretical', 
       distribution='beta',
       title='QQ plot of State_Ischemia eQTLs')
plt.savefig('qq_StateIschemia_theoretical_error.png', dpi=300)

mpl.rcParams['figure.dpi']=300
mpl.rcParams['savefig.dpi']=300
mpl.rcParams['figure.figsize']=5.375, 5.375
qqplot([State_Baseline['f3'], State_Ischemia['f3']], 
       ['Baseline', 'Ischemia'], 
       color=['b','r'], 
       fill_dens=[0.2,0.2], 
       error_type='experimental', 
       distribution='beta',
       title='QQ Plot of State-specific eQTLs')
plt.savefig('qq_two_StateQTLs.png', dpi=300)



Baseline_NoReg_PEER10_FDR1_nonRedundantEgenes=np.genfromtxt('Baseline_NoReg_PEER10_FDR1.nonRedundantEgenes.txt', dtype=None);
print(get_lambda(Baseline_NoReg_PEER10_FDR1_nonRedundantEgenes['f3'], definition = 'median'));
mpl.rcParams['figure.dpi']=300
mpl.rcParams['savefig.dpi']=300
mpl.rcParams['figure.figsize']=5.375, 5.375
qqplot([Baseline_NoReg_Log2_FDR1_nonRedundantEgenes['f3'], Baseline_NoReg_PEER10_FDR1_nonRedundantEgenes['f3']], 
       ['Log2(1.22)', 'PEER10(1.36)'], 
       color=['b','r'], 
       fill_dens=[0.2,0.2], 
       error_type='experimental', 
       distribution='beta',
       title='QQ Plot of Baseline_NoReg eQTLs')
plt.savefig('QQ2_Baseline_NoReg.nonRedundantEgenes.png', dpi=300)


Baseline_RegOut_Log2_FDR1_nonRedundantEgenes=np.genfromtxt('Baseline_RegOut_Log2_FDR1.nonRedundantEgenes.txt', dtype=None);
print(get_lambda(Baseline_RegOut_Log2_FDR1_nonRedundantEgenes['f3'], definition = 'median'));

Baseline_RegOut_PEER10_FDR1_nonRedundantEgenes=np.genfromtxt('Baseline_RegOut_PEER10_FDR1.nonRedundantEgenes.txt', dtype=None);
print(get_lambda(Baseline_RegOut_PEER10_FDR1_nonRedundantEgenes['f3'], definition = 'median'));
mpl.rcParams['figure.dpi']=300
mpl.rcParams['savefig.dpi']=300
mpl.rcParams['figure.figsize']=5.375, 5.375
qqplot([Baseline_RegOut_Log2_FDR1_nonRedundantEgenes['f3'], Baseline_RegOut_PEER10_FDR1_nonRedundantEgenes['f3']], 
       ['Log2(1.26)', 'PEER10(1.36)'], 
       color=['b','r'], 
       fill_dens=[0.2,0.2], 
       error_type='experimental', 
       distribution='beta',
       title='QQ Plot of Baseline_RegOut eQTLs')
plt.savefig('QQ2_Baseline_RegOut.nonRedundantEgenes.png', dpi=300)

Ischemia_NoReg_PEER10_FDR1_nonRedundantEgenes=np.genfromtxt('Ischemia_NoReg_PEER10_FDR1.nonRedundantEgenes.txt', dtype=None);
print(get_lambda(Ischemia_NoReg_PEER10_FDR1_nonRedundantEgenes['f3'], definition = 'median'));
mpl.rcParams['figure.dpi']=300
mpl.rcParams['savefig.dpi']=300
mpl.rcParams['figure.figsize']=5.375, 5.375
qqplot([Ischemia_NoReg_Log2_FDR1_nonRedundantEgenes['f3'], Ischemia_NoReg_PEER10_FDR1_nonRedundantEgenes['f3']], 
       ['Log2(1.22)', 'PEER10(1.36)'], 
       color=['b','r'], 
       fill_dens=[0.2,0.2], 
       error_type='experimental', 
       distribution='beta',
       title='QQ Plot of Ischemia_NoReg eQTLs')
plt.savefig('QQ2_Ischemia_NoReg.nonRedundantEgenes.png', dpi=300)



Ischemia_RegOut_Log2_FDR1_nonRedundantEgenes=np.genfromtxt('Ischemia_RegOut_Log2_FDR1.nonRedundantEgenes.txt', dtype=None);
print(get_lambda(Ischemia_RegOut_Log2_FDR1_nonRedundantEgenes['f3'], definition = 'median'));

Ischemia_RegOut_PEER10_FDR1_nonRedundantEgenes=np.genfromtxt('Ischemia_RegOut_PEER10_FDR1.nonRedundantEgenes.txt', dtype=None);
print(get_lambda(Ischemia_RegOut_PEER10_FDR1_nonRedundantEgenes['f3'], definition = 'median'));
mpl.rcParams['figure.dpi']=300
mpl.rcParams['savefig.dpi']=300
mpl.rcParams['figure.figsize']=5.375, 5.375
qqplot([Ischemia_RegOut_Log2_FDR1_nonRedundantEgenes['f3'], Ischemia_RegOut_PEER10_FDR1_nonRedundantEgenes['f3']], 
       ['Log2(1.26)', 'PEER10(1.36)'], 
       color=['b','r'], 
       fill_dens=[0.2,0.2], 
       error_type='experimental', 
       distribution='beta',
       title='QQ Plot of Ischemia_RegOut eQTLs')
plt.savefig('QQ2_Ischemia_RegOut.nonRedundantEgenes.png', dpi=300)

PatrikTestChr1 = np.genfromtxt('MatrixEQTL.cis.result.PatrikTestDataChr1.txt', dtype=None)
print(get_lambda(PatrikTestChr1['f3'], definition = 'median'))
1.06192544166
mpl.rcParams['figure.dpi']=300
mpl.rcParams['savefig.dpi']=300
mpl.rcParams['figure.figsize']=5.375, 5.375
qqplot([PatrikTestChr1['f3']], 
       ['PatrikTestChr1(lambda=1.062)'], 
       color=['b'], 
       fill_dens=[0.2], 
       error_type='theoretical', 
       distribution='beta',
       title='QQ plot of PatrikTestChr1 eQTLs')
plt.savefig('qq_PatrikTestChr1_theoretical_error.png', dpi=300)



if __name__=="__main__":
    main()
