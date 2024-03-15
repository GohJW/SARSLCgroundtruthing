
######################################################################
# Author :
# Date written : 19-10-2024
# Date last modified : 19-10-2024
# Purpose : to downgrade resolution/resample input complex SAR image  
# Assume range is row dimension and azimuth is col dimension
# If no squint, applying on slant or geometrically corrected (gc) images is equivalent
# If squint, apply on slc images using this code and do gc afterwards 
# input and output defined as weighted and both on slant plane or ground plane
######################################################################


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal.windows import taylor,hamming,hann
from math import log10
import os 
import sys
sys.path.append(os.path.join(sys.path[0],".."))

from SLCtogeotiff.adjustSLC import compressImg

#from SARvisualize import dispimg 

def reswtfactor(win_type='taylor'):  
    if win_type.upper()=='HAMMING':
        win_factor=1.3
    elif win_type.upper()=='TAYLOR':
        win_factor=1.2
    elif win_type.upper()=='RECT':
        win_factor=0.89
    elif win_type.upper()=='HANNING':
        win_factor=1.4466   
    else: 
        print('Invalid window type.')
        win_factor=[]

    return win_factor         

def dr_coherent_multiple_inwtoutwt_preserveintensity(img_in,iprowres,iprowpixsize,ipcolres,ipcolpixsize,oprowres,oprowpixsize,opcolres,opcolpixsize,win_typein,win_typeout,nofftshift,norm_fac,numsubaperture,numsubband):

    win_factorin=reswtfactor(win_typein)
    win_factorout=reswtfactor(win_typeout)


    if (oprowres/win_factorout)<(iprowres/win_factorin):     #Eg in some cases of iceeye  
       print(f'oprowres {oprowres}, win_factorout {win_factorout}, iprowres {iprowres}, win_factorin {win_factorin}')
       print('Reqd row res better than available res')
       oprowres=iprowres/win_factorin*win_factorout

    if (opcolres/win_factorout)<(ipcolres/win_factorin):     #Eg in some cases of iceeye 
       print(f'oprowres {oprowres}, win_factorout {win_factorout}, iprowres {iprowres}, win_factorin {win_factorin}')
       print('Reqd col res better than available res')
       opcolres=ipcolres/win_factorin*win_factorout


    row,col=np.shape(img_in) 
    if (ipcolres/win_factorin)<ipcolpixsize:    #Eg in some cases of iceye
       ipcolres=ipcolpixsize*win_factorin
  
    validrow=int(np.rint(iprowpixsize/(iprowres/win_factorin)*row))
    startvalidrow=int(np.rint(row/2-validrow/2))
    if startvalidrow<0:
        startvalidrow=0
    endvalidrow=startvalidrow+validrow-1
    validcol=int(np.rint(ipcolpixsize/(ipcolres/win_factorin)*col))
    startvalidcol=int(np.rint(col/2-validcol/2))
    if startvalidcol<0:
        startvalidcol=0
    endvalidcol=startvalidcol+validcol-1

    #print(validrow,validcol)
    #print(startvalidrow,endvalidrow,row,startvalidcol,endvalidcol,col)

    img_f=np.fft.ifft2(img_in); 
    #check phhistory is at centre; if not at centre, set fftshiftchoice to True 
    #plt.imshow(np.abs(img_f))
    #plt.title('Orig Phase History')
    #plt.show()

    if nofftshift==False:
        img_f=np.fft.fftshift(img_f)

    if startvalidrow<0:
        startvalidrow=0
    if endvalidrow>row:
        endvalidrow=row
    if startvalidcol<1:
        startvalidcol=0
    if endvalidcol>col:
        endvalidcol=col

    
    img_fvalid=img_f[startvalidrow:endvalidrow+1,startvalidcol:endvalidcol+1]
    print('startvalidrow',startvalidrow,'endvalidrow',endvalidrow,'startvalidcol',startvalidcol,'endvalidcol',endvalidcol)

    
    winin=np.zeros((row,col)) #default rectangular apodization ie no apodization
    if win_typein.upper()=='TAYLOR':
       wtg1drow=taylor(validrow, nbar=5, sll=35).reshape(1,-1)
       wtg1dcol=taylor(validcol, nbar=5, sll=35).reshape(-1,1)
    elif win_typein.upper()=='HANNING':
       wtg1drow=hann(validrow).reshape(1,-1)
       wtg1dcol=hann(validcol).reshape(-1,1)    
    elif win_typein.upper()=='HAMMING':
       wtg1drow=hamming(validrow).reshape(1,-1)
       wtg1dcol=hamming(validcol).reshape(-1,1)
    elif win_typein.upper()=='RECT':
       wtg1drow=np.ones(validrow).reshape(1,-1)
       wtg1dcol=np.ones(validcol).reshape(-1,1)         

    winin[startvalidrow:endvalidrow+1,startvalidcol:endvalidcol+1]=np.outer(wtg1drow,wtg1dcol)
    img_fvalid=np.divide(img_fvalid,winin[startvalidrow:endvalidrow+1,startvalidcol:endvalidcol+1])
    
    oldrow,oldcol=np.shape(img_fvalid)

    newrow=int(np.rint((iprowres/win_factorin)/(oprowres/win_factorout)*oldrow))

    #print(ipcolres,win_factorin,opcolres,win_factorout)
    newcol=int(np.rint((ipcolres/win_factorin)/(opcolres/win_factorout)*oldcol))

    print('newrow',newrow,'newcol',newcol)
    
    if numsubaperture>1:
      deltarow=int(np.floor((oldrow-newrow)/(numsubaperture-1)))
    if numsubband>1:  
      deltacol=int(np.floor((oldcol-newcol)/(numsubband-1)))

    zeropaddedrows=int(np.rint(iprowres/win_factorin*oldrow/oprowpixsize))
    zeropaddedcols=int(np.rint(ipcolres/win_factorin*oldcol/opcolpixsize))

    #remove this line after more checking, zeropaddedrows=int(np.rint(iprowpixsize*row/oprowpixsize))
    #remove this line after more checking, zeropaddedcols=int(np.rint(ipcolpixsize*col/opcolpixsize))  

    #print('zeropaddedrows',zeropaddedrows,'zeropaddedcols',zeropaddedcols)

    img_outlist=[[np.zeros((zeropaddedrows,zeropaddedcols),dtype=np.complex64) for _ in range(numsubband)] for _ in range(numsubaperture)]

    for countaperture in range(numsubaperture):
        if numsubaperture>1:
           centrerow=int(np.rint(newrow/2+countaperture*deltarow))
        else:
           centrerow=int(np.rint(oldrow/2))
    
        for countband in range(numsubband):
            if numsubband>1:
                centrecol=int(np.rint(newcol/2+countband*deltacol))                
            else:
                centrecol=int(np.rint(oldcol/2))

            if (newrow<oldrow):
                startnewrow=int(np.rint(centrerow-newrow/2))
                if startnewrow<0:
                    startnewrow=0
                endnewrow=startnewrow+newrow-1
            elif (newrow==oldrow):
                startnewrow=0
                endnewrow=oldrow-1
            elif (newrow>oldrow):
                print('invalid')
                break

            if (newcol<oldcol):
                startnewcol=int(np.rint(centrecol-newcol/2))
                if startnewcol<0:
                    startnewcol=0
                endnewcol=startnewcol+newcol-1
            elif (newcol==oldcol):
                startnewcol=0
                endnewcol=oldcol-1
            elif (newcol>oldcol):
                print('invalid')
                break

            
            
            
            img_f=img_fvalid[startnewrow:endnewrow+1,startnewcol:endnewcol+1]
            rowf,colf=np.shape(img_f) 
            winout=np.zeros((rowf,colf)) #default rectangular apodization ie no apodization
            if win_typeout.upper()=='TAYLOR':
              wtg1drow=taylor(rowf, nbar=5, sll=35).reshape(1,-1)
              wtg1dcol=taylor(colf, nbar=5, sll=35).reshape(-1,1)
            elif win_typeout.upper()=='HANNING':
              wtg1drow=hann(rowf).reshape(1,-1)
              wtg1dcol=hann(colf).reshape(-1,1)   
            elif win_typeout.upper()=='HAMMING':
              wtg1drow=hamming(rowf).reshape(1,-1)
              wtg1dcol=hamming(colf).reshape(-1,1)   
            elif win_typeout.upper()=='RECT':
              wtg1drow=np.ones(rowf).reshape(1,-1)
              wtg1dcol=np.ones(colf).reshape(-1,1)       
            winout=np.outer(wtg1drow,wtg1dcol)
            img_f=np.multiply(img_f,winout)

            

            img_fzeropadded=np.zeros((zeropaddedrows,zeropaddedcols),dtype=np.complex64)
            startzeropaddedrow=int(np.rint(zeropaddedrows/2-newrow/2))
            if startzeropaddedrow<0:
                startzeropaddedrow=0
            endzeropaddedrow=startzeropaddedrow+rowf-1
            startzeropaddedcol=int(np.rint(zeropaddedcols/2-newcol/2))
            if startzeropaddedcol<0:
                startzeropaddedcol=0
            endzeropaddedcol=startzeropaddedcol+colf-1
            #print(startzeropaddedrow,endzeropaddedrow,startzeropaddedcol,endzeropaddedcol)
            #print(np.shape(img_f))
            #print(np.shape(img_fzeropadded))
            #print(np.shape(img_fzeropadded[startzeropaddedrow:endzeropaddedrow+1,startzeropaddedcol:endzeropaddedcol+1]))
            #print(startzeropaddedrow,endzeropaddedrow+1,startzeropaddedcol,endzeropaddedcol+1)

            img_fzeropadded[startzeropaddedrow:endzeropaddedrow+1,startzeropaddedcol:endzeropaddedcol+1]=img_f

            
            if nofftshift==False:
                img_fzeropadded=np.fft.ifftshift(img_fzeropadded)

            img_out=np.fft.fft2(img_fzeropadded)
            
            #imgpmg,mindpmg,maxdpmg = dispimg.dispimg_modified(img_out,0.05)

            #plt.subplot(1,2,1)
            #plt.imshow(imgpmg,cmap='gray')
            #plt.title('degraded res slc new cmd')
            #plt.show()
            #del imgpmg, mindpmg, maxdpmg
            norm_fac_value = 0
            if norm_fac==True:
                

               tempphhist=np.zeros((row,col))
               tempphhist[startvalidrow:endvalidrow+1,startvalidcol:endvalidcol+1]=winin[startvalidrow:endvalidrow+1,startvalidcol:endvalidcol+1]
               tempimgin=np.fft.fft2(tempphhist)
            
               tempimg_f=np.fft.ifft2(tempimgin)


            
               tempimg_fvalid=tempimg_f[startvalidrow:endvalidrow+1,startvalidcol:endvalidcol+1]
               tempimg_fvalid=np.divide(tempimg_fvalid,winin[startvalidrow:endvalidrow+1,startvalidcol:endvalidcol+1])

               tempimg_f=tempimg_fvalid[startnewrow:endnewrow+1,startnewcol:endnewcol+1]
               tempimg_f=np.multiply(tempimg_f,winout)

               tempimg_fzeropadded=np.zeros((zeropaddedrows,zeropaddedcols),dtype=np.complex64)
            
               tempimg_fzeropadded[startzeropaddedrow:endzeropaddedrow+1,startzeropaddedcol:endzeropaddedcol+1]=tempimg_f

               tempimg_out=np.fft.fft2(tempimg_fzeropadded)

               
            
               norm_fac_value=np.max(np.abs(tempimgin))/np.max(np.abs(tempimg_out))
               img_out=np.multiply(img_out,norm_fac_value)     
            
        
            #imgpmg,mindpmg,maxdpmg = dispimg.dispimg_modified(img_out,0.01)

            #plt.subplot(1,2,1)
            #plt.imshow(imgpmg,cmap='gray')
            #plt.title('degraded res slc')
            #plt.show()
            #del imgpmg, mindpmg, maxdpmg

            img_outlist[countaperture][countband]=img_out

                


    oprowres=row/newrow*iprowpixsize*win_factorout
    opcolres=col/newcol*ipcolpixsize*win_factorout
    oprowpixsize=row/zeropaddedrows*iprowpixsize
    opcolpixsize=col/zeropaddedcols*ipcolpixsize

    return   img_outlist,oprowres,oprowpixsize,opcolres,opcolpixsize,norm_fac_value  