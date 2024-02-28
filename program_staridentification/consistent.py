import numpy
import pandas as pd


A=pd.read_csv('./Result_cg.csv')
B=pd.read_csv('./Result_wcg.csv')
etrsh=1
RA_err=numpy.zeros(len(A))
dec_err=numpy.zeros(len(A))
roll_err=numpy.zeros(len(A))

for i in range(len(A)):
    diff1=abs(A.RA[i]-B.RA[i])
    if diff1>etrsh:
        RA_err[i]=1
    else:
        RA_err[i]=0

    diff2=abs(A.dec[i]-B.dec[i])
    if diff2>etrsh:
        dec_err[i]=1
    else:
        dec_err[i]=0

    diff3=abs(A.roll[i]-B.roll[i])
    if diff3>etrsh:
        roll_err[i]=1
    else:
        roll_err[i]=0

df=pd.DataFrame({'RA_error':RA_err, 'dec_error':dec_err, 'roll_error':roll_err})
#print(df)


cg_res=pd.DataFrame({'No':[],'RA':[],'dec':[],'roll':[]})
wcg_res=pd.DataFrame({'No':[],'RA':[],'dec':[],'roll':[]})

for j in range(len(A)):
    if RA_err[j]==1:
        cg_res.loc[len(cg_res)]={'No':j+1,'RA':A.RA[j],'dec':A.dec[j],'roll':A.roll[j]}
        wcg_res.loc[len(wcg_res)]={'No':j+1,'RA':B.RA[j],'dec':B.dec[j],'roll':B.roll[j]}
        
# Assuming cg_res and wcg_res are your DataFrames
cg_res.to_csv('cg_result.csv', index=False)
wcg_res.to_csv('wcg_result.csv', index=False)

