import sys,os
#sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
#import tools
from texttable import txttableclass
from astropy.io import ascii
from astropy.table import Table, Column
from astropy.coordinates import Angle
import astropy.units as un
import re
import pickle
from scipy.interpolate import UnivariateSpline as spline

class use_texttableclass:

   def __init__(self):
      self.filename=''
      self.surv_name=''
      self.path2TRANSFfile=os.environ.get('PIPE_PS1TRANSFDIR')+'/'


   def read_band(self,filename,surv_name):
       b_list=[]
       txt=use_texttableclass()
       rt=txttableclass()
       rt.loadfile(filename)

       keys1 = rt.CUT_inrange('Surv',surv_name,surv_name)

       for i in keys1:

          b_list.append(rt.getentry(i,'Surv_f'))

       return b_list

   def read_slope_off(self,filename,surv_name,filt):
       slope=0.0
       off=0.0
       txt=use_texttableclass()
       rt=txttableclass()

       rt.loadfile(filename)
       keys1 = rt.CUT_inrange('Surv',surv_name,surv_name)


       for i in keys1:
           band= rt.getentry(i,'Surv_f')

           if band == filt:

               print('>>>>>> cutting %s. Taking only %s and filter %s.' % (filename,surv_name, filt))
               rt.printtxttable(cols=['Surv','PS1_f','Surv_f','Slope_data','Off_data'],keys=i)
               slope= rt.getentry(i,'Slope_data')
               off= rt.getentry(i,'Off_data')


               return slope,off,band

   def read_slope_off_filts(self,filename,surv_name,filt):
       slope=0.0
       off=0.0
       txt=use_texttableclass()
       rt=txttableclass()

       rt.loadfile(filename)
       keys1 = rt.CUT_inrange('Surv',surv_name,surv_name)


       for i in keys1:
           band= rt.getentry(i,'Surv_f')

           if band == filt:

               print('>>>>>> cutting %s. Taking only %s and filter %s.' % (filename,surv_name, filt))
               rt.printtxttable(cols=['Surv','PS1_f','Surv_f','Slope_data','Off_data','PS1_f1','PS1_f2'],keys=i)
               slope= rt.getentry(i,'Slope_data')
               off= rt.getentry(i,'Off_data')
               filt1= rt.getentry(i,'PS1_f1')
               filt2= rt.getentry(i,'PS1_f2')

               return slope,off,band,filt1,filt2


   def read_cat(self,outfilename,sexagesimal,filterlist=['g','r','i','z','y']):
      mag=dict()
      err=dict()
      filt_num=dict()


      data = ascii.read(outfilename)

      filt_num={'g':2,'r':4,'i':6,'z':8,'y':10}

      for line in data:
         mag.setdefault('ra',[]).append(line[0])
         mag.setdefault('dec',[]).append(line[1])




         for filt in filterlist:
             # no idea what this was..
             #mag.setdefault(filt,[]).append("{:.4f}".format(float(line[filt_num[filt]])))
             #err.setdefault(filt,[]).append("{:.4f}".format(float(line[filt_num[filt]+1])))
             mag.setdefault(filt,[]).append("{:.4f}".format(float(line[filt])))
             err.setdefault(filt,[]).append("{:.4f}".format(float(line['d'+filt])))



      return mag,err,filterlist#,colormag,colorerr

   def save_cat(self,mag,tr_mag,err,filterlist,outfilename,path2name,b_list,surv_name,selected_band):
      var=[]
      header=['#ra','dec']
      var.append(mag['ra'])
      var.append(mag['dec'])


      for filt in selected_band:
        if filt in b_list:
            err_filt=['d'+filt]

            print('>>>>>> Saving transformed %s filter in %s.' % (filt,outfilename))
            var.append(tr_mag[filt])
            if filt == 'B':var.append(err['g'])
            elif filt == 'V':var.append(err['g'])
            elif filt == 'u':var.append(err['g'])
            else: var.append(err[filt])
            header.extend(filt)
            header.extend(err_filt)
        else:
           print('>>>>>> WARNING: Not saving %s filter in %s because it is not in %s bands list: %s.' % (filt,outfilename,path2name,b_list))
           #print header
        data=Table(var,names=header)
        # Hack to get rid of bad data in catalog
        masknames = header[2:]
        for name in masknames:
          data = data[[float(d[name]) < 30.0 for d in data]]
          data = data[[float(d[name]) > 0.0 for d in data]]
      ascii.write(data,outfilename,overwrite=True)


   def load_spline(self,filt,surv_name):
      if surv_name == 'DES0':
         if filt == 'u': name_spline=self.path2TRANSFfile+'PS12DES_u_s1.393_db0.15_spline.pkl'
         if filt == 'g': name_spline=self.path2TRANSFfile+'PS12DES_g_s0.0003_db0.15_spline.pkl'
         if filt == 'r': name_spline=self.path2TRANSFfile+'PS12DES_r_s0.000325_db0.15_spline.pkl'
         if filt == 'i': name_spline=self.path2TRANSFfile+'PS12DES_i_s0.0008_db0.15_spline.pkl'
         if filt == 'z': name_spline=self.path2TRANSFfile+'PS12DES_z_s0.00635_db0.15_spline.pkl'
      if surv_name == 'SNLS0':
         if filt == 'u': name_spline=self.path2TRANSFfile+'PS12SNLS_u_s0.88_db0.15_spline.pkl'
         if filt == 'g': name_spline=self.path2TRANSFfile+'PS12SNLS_g_s0.00018_db0.15_spline.pkl'
         if filt == 'r': name_spline=self.path2TRANSFfile+'PS12SNLS_r_s0.0001_db0.15_spline.pkl'
         if filt == 'i': name_spline=self.path2TRANSFfile+'PS12SNLS_i_s0.00015_db0.15_spline.pkl'
         if filt == 'z': name_spline=self.path2TRANSFfile+'PS12SNLS_z_s0.0005_db0.15_spline.pkl'
      if surv_name == 'CSP1':
         if filt == 'u': name_spline=self.path2TRANSFfile+'PS12CSP_u_s2.1_db0.1_spline.pkl'
         if filt == 'B': name_spline=self.path2TRANSFfile+'PS12CSP_B_s0.09_db0.15_spline.pkl'
         if filt == 'V': name_spline=self.path2TRANSFfile+'PS12CSP_V_s0.0101_db0.15_spline.pkl'
         if filt == 'g': name_spline=self.path2TRANSFfile+'PS12CSP_g_s0.0018_db0.15_spline.pkl'
         if filt == 'r': name_spline=self.path2TRANSFfile+'PS12CSP_r_s0.0001_db0.15_spline.pkl'
         if filt == 'i': name_spline=self.path2TRANSFfile+'PS12CSP_i_s0.0001_db0.15_spline.pkl'



      print('> Loading spline structure from: %s' % name_spline)
      file2 = open(name_spline, 'rb')
      # Python 2/3 compatibility
      try:
        ld_spline = pickle.load(file2)
      except:
        u = pickle._Unpickler(file2)
        u.encoding = 'latin1'
        ld_spline = u.load()
      file2.close()

      return(ld_spline)

class PS1toINSTclass:


    def convert(self,path2name,surv_name,outfilename,selected_band,sexagesimal,spline,filterlist=['g','r','i','z','y']):

       texttable=use_texttableclass()


       b_list= texttable.read_band(path2name,surv_name)
       if not spline:
          print('>>>>>> reading survay: ',surv_name)
          print('>>>>>> reading file: ',path2name)
          print('>>>>>> outfile name: ',outfilename)
       # AR: you need to check if survey exists, and if not throw an error!!!
          if len(b_list)==0:
             raise RuntimeError("survey %s does not exist!" % (surv_name))

          mag,err,filterlist=texttable.read_cat(outfilename,sexagesimal,filterlist=filterlist)#,coloroutfilename,colorfilterlist,band)
          tr_mag=dict()


       if not spline:
          for filt in selected_band:
             if filt =='u' and filt in b_list:
                print('>>>>>> working on filter u ')
                print('>>>>>> yes! %s filter is in %s filters list: %s' % ('u',path2name,b_list))

                slope,off,band=texttable.read_slope_off(path2name,surv_name,'u')

                slope=float(slope)

                off=float(off)
                print('>>>>>> converting %s PS1 fiter to %s %s filter. slope = %s, off = %s derived from: PS1g - %su = slope*(PS1g- PS1i) + off' % ('g','u',surv_name,slope,off,surv_name))

                print('>>>>>> PS1g-slope*(PS1g-PS1i)-off=%su(PS1g)' %(surv_name))
                for t in range(len(mag['g'])):
                   tr_mag.setdefault('u',[]).append("{:.4f}".format(float(mag['g'][t])-(slope)*(float(mag['g'][t])-float(mag['i'][t]))-(off)) )
                #print '%f-(%f)*(%f-%f)-(%f)=%f' %(float(mag['g'][t]),slope,float(mag['g'][t]),float(mag['i'][t]),off,float(tr_mag['u'][t]))

             elif filt =='B' and filt in b_list:
                print('>>>>>> working on filter B ')
                print('>>>>>> yes! %s filter is in %s filters list: %s' % ('B',path2name,b_list))

                slope,off,band=texttable.read_slope_off(path2name,surv_name,'B')

                slope=float(slope)

                off=float(off)
                print('>>>>>> converting %s PS1 fiter to %s %s filter. slope = %s, off = %s derived from: PS1g - %sB = slope*(PS1g- PS1i) + off' % ('g','B',surv_name,slope,off,surv_name))

                print('>>>>>> PS1g-slope*(PS1g-PS1i)-off=%sB(PS1g)' %(surv_name))
                for t in range(len(mag['g'])):
                   tr_mag.setdefault('B',[]).append("{:.4f}".format(float(mag['g'][t])-(slope)*(float(mag['g'][t])-float(mag['i'][t]))-(off)) )
                #print '%f-(%f)*(%f-%f)-(%f)=%f' %(float(mag['g'][t]),slope,float(mag['g'][t]),float(mag['i'][t]),off,float(tr_mag['u'][t]))

             elif filt =='V' and filt in b_list:
                print('>>>>>> working on filter V ')
                print('>>>>>> yes! %s filter is in %s filters list: %s' % ('V',path2name,b_list))

                slope,off,band=texttable.read_slope_off(path2name,surv_name,'V')

                slope=float(slope)

                off=float(off)
                print('>>>>>> converting %s PS1 fiter to %s %s filter. slope = %s, off = %s derived from: PS1r - %sV = slope*(PS1g- PS1i) + off' % ('r','V',surv_name,slope,off,surv_name))

                print('>>>>>> PS1r-slope*(PS1g-PS1i)-off=%sV(PS1g)' %(surv_name))
                for t in range(len(mag['r'])):
                   tr_mag.setdefault('V',[]).append("{:.4f}".format(float(mag['r'][t])-(slope)*(float(mag['g'][t])-float(mag['i'][t]))-(off)) )
                #print '%f-(%f)*(%f-%f)-(%f)=%f' %(float(mag['g'][t]),slope,float(mag['g'][t]),float(mag['i'][t]),off,float(tr_mag['u'][t]))


             elif filt != 'u' and filt != 'B' and filt != 'V' and filt in b_list:


                print('>>>>>> working on filter %s ' % (filt))
                print('>>>>>> yes! %s filter is in %s filters list: %s' % (filt,path2name,b_list))
                #slope,off,band=texttable.read_slope_off(path2name,surv_name,filt)
                slope,off,band,filt1,filt2=texttable.read_slope_off_filts(path2name,surv_name,filt)

                print(f'>>>>>> converting {band} PS1 fiter to {filt} {surv_name} filter. slope = {slope}, off = {off} derived from: PS1{filt} -{surv_name}{filt} = slope*(PS1{filt1}- PS1{filt2})+off')
                slope=float(slope)
                off=float(off)
                print('>>>>>> PS1{filt}-slope*(PS1{filt1}-PS1{filt2})-off={surv_name}{filt}(PS1{filt})')
                for t in range(len(mag[filt])):
                   tr_mag.setdefault(filt,[]).append("{:.4f}".format(float(mag[filt][t])-(slope)*(float(mag[filt1][t])-float(mag[filt2][t]))-(off)) )
                #print '%f-(%f)*(%f-%f)-(%f)=%f' %(float(mag[filt][t]),slope,float(mag['g'][t]),float(mag['i'][t]),off,float(tr_mag[filt][t]))
             else:

                print('>>>>>> WARNING: Skipping %s filter because it is not in %s bands list: %s.' % (filt,path2name,b_list))



       elif spline:
          print('Spline!')
          mag,err,filterlist=texttable.read_cat(outfilename,sexagesimal)#,coloroutfilename,colorfilterlist,band)
          tr_mag=dict()
          for filt in selected_band:
             if filt =='B' and filt in b_list:
                print('>>>>>> working on filter B ')
                ld_spline=texttable.load_spline(filt,surv_name)
                for t in range(len(mag['g'])):
                   tr_mag.setdefault('B',[]).append("{:.4f}".format(float(mag['g'][t])- ld_spline.__call__(float(mag['g'][t])-float(mag['i'][t]),ext=0)))
             elif filt =='V' and filt in b_list:
                print('>>>>>> working on filter V ')
                ld_spline=texttable.load_spline(filt,surv_name)
                for t in range(len(mag['r'])):
                   tr_mag.setdefault('V',[]).append("{:.4f}".format(float(mag['r'][t])- ld_spline.__call__(float(mag['g'][t])-float(mag['i'][t]),ext=0)))
             elif filt =='u' and filt in b_list:
                print('>>>>>> working on filter u ')
                ld_spline=texttable.load_spline(filt,surv_name)
                for t in range(len(mag['g'])):
                   tr_mag.setdefault('u',[]).append("{:.4f}".format(float(mag['g'][t])- ld_spline.__call__(float(mag['g'][t])-float(mag['i'][t]),ext=0)))
             elif filt != 'u'  and filt != 'B' and filt != 'V' and filt in b_list:
                print('>>>>>> working on filter %s ' % (filt))
                ld_spline=texttable.load_spline(filt,surv_name)
                for t in range(len(mag[filt])):
                   tr_mag.setdefault(filt,[]).append("{:.4f}".format(float(mag[filt][t])- ld_spline.__call__(float(mag['g'][t])-float(mag['i'][t]),ext=0)))
             else:

                print('>>>>>> WARNING: Skipping %s filter because it is not in %s bands list: %s.' % (filt,path2name,b_list))

       texttable.save_cat(mag,tr_mag,err,filterlist,outfilename,path2name,b_list,surv_name,selected_band)
