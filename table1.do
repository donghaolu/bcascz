**
** TABLE 1 EPIDEMIOLOGICAL ASSOCIATIONS
**

* created by Donghao Lu, donghao.lu@ki.se


**
** INPATIENT DX ONLY
**


*****************************
********SCZ TO BCA **********
*****************************

use SczBCaCohort20180319, clear

*mean age
gen age=(ref_date-birthdate)/365.25
summarize age if case==1, d

*code schizophrenia
gen scz=0
replace scz=1 if scz_inpdate<ref_date

*code parity: 0, 1-2, >=3
gen para=0
replace para=1 if inlist(parity,1,2)
replace para=2 if parity>=3

*code previous psy dx
gen prepsy=0
replace prepsy=1 if psy1_date<ref_date & psy1_date<scz_inpdate

*N of individuals/events
tab scz case, matcell(freq) 

*OR - basic model
clogit case scz ib3.edu ib2.region, group(setno) or 

*OR - second model
clogit case scz ib3.edu ib2.region ib1.para prepsy obesity, group(setno) or 

  
  


*****************************
********BCA TO SCZ **********
*****************************

use SczBCaCohort20180319, clear

*code schizophrenia
gen new_exit=min(scz_inpdate, exit, date("12/31/2010","MDY"))
gen scz=0
replace scz=1 if new_exit==scz_inpdate
replace new_exit=new_exit+1 if new_exit==ref_date

*age at dx
gen age=(scz_inpdate-birthdate)/365.25 if scz==1
summarize age if ref_date<=new_exit, d

*code parity: 0, 1-2, >=3
gen para=0
replace para=1 if inlist(parity,1,2)
replace para=2 if parity>=3

*code previous psy dx
gen prepsy=0
replace prepsy=1 if psy1_date<ref_date

*set time
stset new_exit, failure(scz==1) enter(ref_date) exit(time new_exit) origin(ref_date) scale(365.25) 

*N of individuals/events
tab _d case if _st==1, matcell(freq) 

  
*HR - basic model
stcox case ib3.edu ib2.region, strata(setno)
  
*HR - second model
stcox case ib3.edu ib2.region ib1.para prepsy obesity, strata(setno)




**
** INPATIENT & OUTPATIENT DX
**

*****************************
********SCZ TO BCA **********
*****************************

use SczBCaCohort20180319, clear

*mean age
gen age=(ref_date-birthdate)/365.25
summarize age if case==1, d

*code schizophrenia
gen scz=0
replace scz=1 if scz_date<ref_date


*median follow-up from scz to bca
gen t=(ref_date-scz_date)/365.25 if scz==1
summarize t, d
drop t

*code parity: 0, 1-2, >=3
gen para=0
replace para=1 if inlist(parity,1,2)
replace para=2 if parity>=3

*code previous psy dx
gen prepsy=0
replace prepsy=1 if psy1_date<ref_date & psy1_date<scz_date

*N of individuals/events
tab scz case, matcell(freq) 

*OR - basic model
clogit case scz ib3.edu ib2.region, group(setno) or 

*OR - second model
clogit case scz ib3.edu ib2.region ib1.para prepsy obesity, group(setno) or 

  
  
*****************************
********BCA TO SCZ **********
*****************************

use SczBCaCohort20180319, clear

*code schizophrenia
gen new_exit=min(scz_date, exit, date("12/31/2010","MDY"))
gen scz=0
replace scz=1 if new_exit==scz_date
replace new_exit=new_exit+1 if new_exit==ref_date

*code parity: 0, 1-2, >=3
gen para=0
replace para=1 if inlist(parity,1,2)
replace para=2 if parity>=3

*code previous psy dx
gen prepsy=0
replace prepsy=1 if psy1_date<ref_date

*set time
stset new_exit, failure(scz==1) enter(ref_date) exit(time new_exit) origin(ref_date) scale(365.25) 

*N of individuals/events
tab _d case if _st==1, matcell(freq) 

*HR - basic model
stcox case ib3.edu ib2.region, strata(setno)

*HR - second model
stcox case ib3.edu ib2.region ib1.para prepsy obesity, strata(setno)

  

**
** END OF FILE
**






