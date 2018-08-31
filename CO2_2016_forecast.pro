pro CO2_2016_forecast

; Script to forecast 2015/2016 CO2 rise.
; CDJ. 18 Nov. 2015
;
; checked and reviewed for archive as requried by the journal.
; data files and this IDL script to be stored together. Output is the timeseries
; of 12 monthly values for 2016 forecast as publsihed by Betts et al., 2016
;
; CDJ. 14 June 2018
;
; Internally reviewed by Andy Wiltshire, 30 Aug 2018
;
;
; 1. Read in data, and choose which station/index to use
;

; 1.1. CO2 data
;
; Mauna Loa data downloaded from here: http://scrippsco2.ucsd.edu/data/atmospheric_co2
; on 1 Dec 2015
; 
; this has "raw" and "filled" data. Use the filled. (column 9)
;

; Mauna Loa
a1=0.
a2=0.
a3=0.
yco2_mon=0.
a5=0.
a6=0.
a7=0.
a8=0.
co2_mon=0.
x=dc_read_free("monthly_mlo.csv",nskip=57,delim=',',a1,a2,a3,yco2_mon,a5,a6,a7,a8,co2_mon,/col,$
                re=[1,2,3,4,5,6,7,8,9])

; create annual means from 1959-2014 (inclusive)
yco2=1959
co2=avg(co2_mon(where(fix(yco2_mon) eq 1959)))
for i=1960,2014 do begin
 yco2=[yco2,i]
 co2=[co2,avg(co2_mon(where(fix(yco2_mon) eq i)))]
endfor

y0=min(yco2)
ny=2015-y0
print,y0,ny

;
;------------------------------------------------------------------------
;

; 1.2. emissions data
; downloaded emissions from GCP data here: http://cdiac.ornl.gov/GCP/
; on 10 Nov 2015
;
; I stripped out just the fossil and land-use emissions to this text file:
x=dc_read_free("GCP_Emiss.txt",yemiss,eff,elu,/col,re=[1,2,3])
emiss=eff+elu

; just get years according to which CO2 data we use
i0=where(yemiss eq y0)
yemiss = yemiss(i0:*)
eff    = eff(i0:*)
elu    = elu(i0:*)
emiss  = emiss(i0:*)

;
;------------------------------------------------------------------------
;

; 1.3. Nino3/3.4 data
; 

; 1.3.1 Historical from HadISST1.1
;
file='HadISST1.1_sst_1870on_1dg_anm6190_ninos_save'
restore,file
date_nino=findgen(n_elements(nino3))/12.+1870+1./24.

; create annuals from y0 onwards
ynino = findgen(ny) + y0 + 10./12.
nino3_ann  = fltarr(ny)
nino34_ann = fltarr(ny)

for i=0,ny-1 do begin 
  i0=min(where(fix(date_nino) eq y0+i)) 
;   print,i,i0,y0+i,date_nino(i0+3),date_nino(i0+14) 
  nino3_ann(i)  = avg(nino3(i0+3:i0+14)) 
  nino34_ann(i) = avg(nino34(i0+3:i0+14))
endfor

; 1.3.2 projected, 40 realisations
;
file='nino34_combined_realisations.dat'

n34_proj=fltarr(16,40)
yn34_proj=findgen(16)/12.+2015.+1./24.
openr,1,file
readf,1,n34_proj
close,1


; 1.3.3 read in HadISST3 ensemble of historical realisations and create
; timeseries of median annual values.

file_n34='Nino_34_1870.txt'
status = dc_read_free(file_n34, /col, yy)
nm = n_elements(yy)
mm = fltarr(nm)
nino = fltarr(100, nm)
nino0 = fltarr(100)
;
openr, unit, file_n34, /get_lun
for i=0,nm-1 do begin 
  readf, unit, yy0, mm0, nino0, format='(i4,x,i2,x,100(f7.3,x))' 
  yy(i) = yy0 
  mm(i) = mm0 
  nino(*,i) = nino0 
endfor
close, unit
free_lun, unit

index=where(yy ge min(yco2))
yy_n34 = yy(index)
mm_n34 = mm(index)
nino34_ens = nino(*,index) 
date_nino_ens=yy_n34+(mm_n34-0.5)/12.0

; create annuals from y0 onwards
ynino_ens = findgen(ny) + y0 + 10./12.
nino34_ann_ens  = fltarr(100,ny)

for j=0,99 do begin 
  for i=0,ny-1 do begin 
    i0=min(where(fix(date_nino_ens) eq y0+i)) 
    nino34_ann_ens(j,i)  = avg(nino34_ens(j,i0+3:i0+14)) 
  endfor 
endfor

; calculate median annual values:
nino34_ann_median=fltarr(ny)
for i=0,ny-1 do begin 
  xx=nino34_ann_ens(*,i) 
  zz=sort(xx) 
  yy=xx(zz) 
  nino34_ann_median(i)=avg(yy(49:50))
endfor

;
;------------------------------------------------------------------------
;
; 2. Calculate regression coeffs and work out reconstructed delta-CO2s
;

; 2.1 observed delta-CO2, centred on 1 Jan of each year
;
; label this from the later of the 2 years - i.e. the delta CO2 for 1 Jan 2015
; means the increment from 2014 annual mean to 2015 annual mean.
ydco2=yco2(1:*)
dco2=co2(1:*)-co2

;------------------------------------------------------------------------
; 2.2 do a regression onto nino and emissions
; use Nino3.4, and without 1992/3 Pinatubo years
;
index=where((ydco2 ne 1992) and (ydco2 ne 1993))
X34_novolc = [TRANSPOSE(nino34_ann(index)), TRANSPOSE(emiss(index))]
regr34_novolc = REGRESS(X34_novolc, dco2(index), SIGMA=sigma34_novolc, CONST=const34_novolc)

; PLUS an extra one using the median HadSST3 data
;
X34_novolc_median = [TRANSPOSE(nino34_ann_median(index)), TRANSPOSE(emiss(index))]
regr34_novolc_median = REGRESS(X34_novolc_median, dco2(index), SIGMA=sigma34_novolc_median, CONST=const34_novolc_median)

dco2_recon34_novolc = const34_novolc + regr34_novolc(0)*nino34_ann + $
                   regr34_novolc(1)*emiss
dco2_recon34_novolc_median = const34_novolc_median + regr34_novolc_median(0)*nino34_ann_median + $
                   regr34_novolc_median(1)*emiss

;------------------------------------------------------------------------
; 2.3 now predict increments and annual CO2 for 2015/2016 using (a) obs Nino3.4
; for 2015 and projected Nino3.4 for 2016. The latter will have 40 realisations

; 2015
;
; according to GCP budget, Eff=9.8, Elu=1.1
dco2_recon_2015 = const34_novolc + regr34_novolc(0)*nino34_ann(ny-1) + $
                   regr34_novolc(1)*10.9
dco2_recon_2015_median = const34_novolc_median + regr34_novolc_median(0)*nino34_ann_median(ny-1) + $
                   regr34_novolc_median(1)*10.9

; 2016
; according to GCP budget, Eff=9.2, Elu=1.1
n34_proj_ann=fltarr(40)
for i=0,39 do n34_proj_ann(i)=avg(n34_proj(3:14,i))

dco2_recon_2016 = const34_novolc + regr34_novolc(0)*n34_proj_ann + $
                   regr34_novolc(1)*10.3
dco2_recon_2016_median = const34_novolc_median + regr34_novolc_median(0)*n34_proj_ann + $
                   regr34_novolc_median(1)*10.3

; use new values for SST forecast. These are mean and +- 2 sigma
; on the SST forecast:
dco2_recon_2016_Jeff = const34_novolc_median + $
                       regr34_novolc_median(0)*[1.79,2.02,2.25] + $
                       regr34_novolc_median(1)*10.3


;------------------------------------------------------------------------
; 2.4 now predict monthly CO2 itself using a mean S.C. from 2010-14
;
ind=where(yco2_mon ge 2010 and yco2_mon le 2014.99)
ysc=yco2_mon(ind)
csc=co2_mon(ind)

; calculate a linear trend from Jan 2015 - Jan 2010 divided by 60 months
tr=(co2_mon(ind(59)+1)-co2_mon(ind(0)))/60.
csc=csc-findgen(60)*tr

ysc=ysc(0:11)-2010
for i=0,11 do csc(i)=avg(csc(indgen(5)*12+i))
csc=csc(0:11)-avg(csc(0:11))

; now add this on to the 2015/16 annual numbers, and re-introduce the linear trend
y_mon_recon   = [2015+ysc,2016+ysc]

co2_2015=co2(ny-1)+dco2_recon_2015
co2_2016=co2_2015+avg(dco2_recon_2016)
co2_mon_recon = [co2_2015+csc+(indgen(12)-5.5)*tr, $
                 co2_2016+csc+(indgen(12)-5.5)*tr]

co2_2015_median=co2(ny-1)+dco2_recon_2015_median
co2_2016_median=co2_2015_median+dco2_recon_2016_jeff(1)
co2_mon_recon_median = [co2_2015_median+csc+(indgen(12)-5.5)*tr, $
                 co2_2016_median+csc+(indgen(12)-5.5)*tr]

print,'co2_2015_median = ', co2_2015_median
print,'co2_2016_median = ', co2_2016_median

print,'2016 Monthly Mauna Loa predicted CO2'
print,'(co2_mon_recon_median)'
for i=0,11 do begin
  print,i+1, co2_mon_recon_median(i+12)
endfor

;------------------------------------------------------------------------
; 2.5 Error/uncertainty estimates
;

print,'--'
print,'uncertainty terms'
print

n34_proj_ann_ens=avg(n34_proj,0)
n34_proj_ann_sd=stdev(n34_proj_ann_ens)

print,'Uncertainty in N34 projection: (1-sigma) ', n34_proj_ann_sd
print,'Uncertainty in D-CO2 due to N34 projection: (1-sigma) ', regr34_novolc_median(0)*n34_proj_ann_sd

d_dco2=dco2(index)-dco2_recon34_novolc_median(index)
dco2_error=stdev(d_dco2)
print
print,'Error (1-sigma) in past D-CO2 reconstruction ', dco2_error

print
print,'combined uncertainty (sum-of-squares) = ', sqrt (dco2_error^2 + $
        (regr34_novolc_median(0)*n34_proj_ann_sd)^2 )

end
