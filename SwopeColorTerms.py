#!/usr/bin/env python
# D. Jones - 3/3/24
# Get Swope color offsets compared to Pan-STARRS

import numpy as np
import astropy.table as at
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import least_squares
from scipy.stats import binned_statistic
import glob
import os
import pylab as plt
plt.ion()

_max_dm = 0.02
_min_sep_arcsec = 2
_sigmaclip=3

#### Color Terms ####
# u: g-r
# B: g-r
# V: g-r
# g: g-r
# r: g-r
# i: r-i

swope2ps = {
    #'u':['g','g','r'],
    'B':['g','g','r'],
    'V':['g','g','r'],
    'g':['g','g','r'],
    'r':['r','g','r'],
    'i':['i','r','i']
}

colorcuts = {
    'B':[0.25,0.85],
    'V':[0.3,1.4],
    'g':[0.2,1.3],
    'r':[0.2,1.3],
    'i':[0.1,1.5]
}

### Steps ###
# 1. Read in swope fluxes: calspec/fluxfilesswope/
# 2. Get PS1 catalogs for each
# 3a. Solve for avg ZPT (no color terms) &
# 3b. Plot swope mags versus PS1 colors
# 4. Get best-fit slopes
# 5a. Play around with cuts: airmass, seeing, maglim(?), date &
# 5b. Figure out if those params make any significant difference in slopes
# 5c. Decide on color cuts for ZPT stars
# 6. Figure out how to put into photpipe

def errfnc(x):
    return(np.std(x)/np.sqrt(len(x)))

def getradecbox(ra,dec,size):

    minDec = dec-0.5*size
    if minDec<=-90.0:minDec=-89.9
    maxDec = dec+0.5*size
    if maxDec>=90.0:maxDec=89.9

    invcosdec = max(1.0/np.cos(dec*np.pi/180.0),
                    1.0/np.cos(minDec  *np.pi/180.0),
                    1.0/np.cos(maxDec  *np.pi/180.0))

    # get the (conservative) boxlimits
    ramin = ra-0.5*size*invcosdec
    ramax = ra+0.5*size*invcosdec
    decmin = dec-0.5*size
    decmax = dec+0.5*size

    if (ra-ramin)<-180:
        ramin-=360.0
        ramax-=360.0
    elif (ra-ramin)>180:
        ramin+=360.0
        ramax+=360.0

    return ramin,ramax,decmin,decmax

class SwopeColors:

    def getPScats(self,field_info):

        for f in field_info:
            outfile = f'PS1cats/all/{f[2]}.PS1_3piDR2.photcat'
            if not os.path.exists(outfile):
                cmd = f'python getPS1DR2cat.py {f[0]} {f[1]} --dMmax {_max_dm} --size 0.75 --outfile {outfile} -gri --nocutsinfile'
                print(f'running: {cmd}')
                os.system(cmd)

        # should maybe cut out stars near SN location but punt for now

    def MatchCats(self,swope_files,swope_filt):

        matched_table = at.Table(
            [np.array([]),np.array([]),np.array([]),np.array([]),
             np.array([]),np.array([]),np.array([]),np.array([]),
             np.array([]),np.array([]),np.array([])],
            names=('cat_id','ra','dec',f'{swope_filt}_Swope_flux',
                   f'{swope_filt}_Swope_fluxerr','g_PS','gerr_PS',
                   'r_PS','rerr_PS','i_PS','ierr_PS')
        )
        
        for j,sf in enumerate(swope_files):
            objname = sf.split('/')[-1].split('.')[0]
            pscat = f'PS1cats/all/{objname}.PS1_3piDR2.photcat'

            ps = at.Table.read(pscat,format='ascii')
            ps = ps[(ps['g'] > -998) & (ps['r'] > -998) & (ps['i'] > -998)]
            sw = at.Table.read(sf,format='ascii',names=('ra','dec','flux','fluxerr'))

            try:
                sw = sw[~sw['fluxerr'].mask]
            except AttributeError:
                pass
            
            with open(sf) as fin:
                for line in fin:
                    if line.startswith('#airmass'):
                        airmass = line.split()[1].replace('\n','')
                    elif line.startswith('#mjd'):
                        mjd = line.split()[1].replace('\n','')
                    elif line.startswith('#seeing'):
                        seeing = line.split()[1].replace('\n','')
                    elif line.startswith('#m5sigma'):
                        m5sigma = line.split()[1].replace('\n','')
                    elif line.startswith('#exptime'):
                        exptime = line.split()[1].replace('\n','')
                    elif line.startswith('#ra'):
                        break

            for i,s in enumerate(sw):
                if s['flux'] < 0 or s['flux'] != s['flux']:
                    continue
                if 1.0857*s['fluxerr']/s['flux'] > _max_dm:
                    continue
                scs = SkyCoord(s['ra'],s['dec'],unit=u.deg)
                scp = SkyCoord(ps['ra'],ps['dec'],unit=u.deg)
                sep = scs.separation(scp).arcsec
                if min(sep) < _min_sep_arcsec:
                    iMatched = np.where(sep == np.min(sep))[0]
                    if len(iMatched) == 1:
                        ps_color = \
                            ps[swope2ps[swope_filt][1]][iMatched][0] - ps[swope2ps[swope_filt][2]][iMatched][0]
                        if len(colorcuts[swope_filt]) == 0:
                            matched_table.add_row(
                                (
                                    j,s['ra'],s['dec'],s['flux'],s['fluxerr'],
                                    ps['g'][iMatched],ps['dg'][iMatched],
                                    ps['r'][iMatched],ps['dr'][iMatched],
                                    ps['i'][iMatched],ps['di'][iMatched],
                                )
                            )
                        elif (ps_color > colorcuts[swope_filt][0]) & \
                           (ps_color < colorcuts[swope_filt][1]):
                            matched_table.add_row(
                                (
                                    j,s['ra'],s['dec'],s['flux'],s['fluxerr'],
                                    ps['g'][iMatched],ps['dg'][iMatched],
                                    ps['r'][iMatched],ps['dr'][iMatched],
                                    ps['i'][iMatched],ps['di'][iMatched],
                                )
                            )

        matched_table.write(f'MatchedCats/{swope_filt}_matched.txt',overwrite=True,format='ascii')

        return matched_table

    def fit_zpts(self,matched_table,swope_filt,npars=2):

        inst_mags = -2.5*np.log10(matched_table[f"{swope_filt}_Swope_flux"])
        inst_magerrs = 1.0857*matched_table[f"{swope_filt}_Swope_fluxerr"]/matched_table[f"{swope_filt}_Swope_flux"]
        unq_ids = np.unique(matched_table["cat_id"])

        x_full = matched_table[f"{swope2ps[swope_filt][1]}_PS"] - matched_table[f"{swope2ps[swope_filt][2]}_PS"]
        y_full = inst_mags[:] - matched_table[f"{swope2ps[swope_filt][0]}_PS"]
        yerr_full = inst_magerrs[:]

        def line_fit_model(p,x):
            model = 0
            for i,f in enumerate(p[::-1]):
                model += f*x**i

            return model
                
        def line_fit_noclip(p,x,y):

            model = 0
            for i,f in enumerate(p[::-1]):
                model += f*x**i

            return y - model


        iNoClip = np.ones(len(y_full),dtype=bool)
        for i in range(5):
            ### ZPTs for each file individually
            for c in unq_ids:
                # first we need zeropoints for everything
                
                iNoClip_single = np.where((matched_table["cat_id"] == c) & (iNoClip))[0]
                x = x_full[iNoClip_single]
                y = y_full[iNoClip_single]
                yerr = yerr_full[iNoClip_single]

                
                def line_fit(p):
                    
                    model = 0
                    for i,f in enumerate(p[::-1]):
                        model += f*x**i
                
                    return (y - model) #/yerr

                guess = [1]*npars #+ [25]*len(unq_ids)
                result = least_squares(line_fit,guess)
                # difference between Swope and PS at median color should be zero
                inst_mags[matched_table["cat_id"] == c] -= line_fit_model(result.x,np.median(x))

                #plt.plot(x,y,'.')
                #magbins = np.linspace(-2,3,15)
                #plt.plot(magbins,line_fit_model(result.x[0:len(unq_ids)],magbins),color='k')
                #binned_mags = binned_statistic(x,y,bins=magbins,statistic='median').statistic
                #plt.plot((magbins[1:]+magbins[:-1])/2.,binned_mags,'o-')

            y_full = inst_mags[:] - matched_table[f"{swope2ps[swope_filt][0]}_PS"]

            # crazy outliers are causing problems
            if i == 0 and swope_filt == 'B':
                iNoClip = (y_full > -0.5) & (y_full < 0.5)
            def line_fit_all(p):
                 
                model = 0
                for i,f in enumerate(p[::-1]):
                    model += f*x_full[iNoClip]**i
                
                return (y_full[iNoClip] - model) #/yerr_full[iNoClip]
            guess = [1]*npars
            result_full = least_squares(line_fit_all,guess)


            resids = line_fit_noclip(result_full.x,x_full,y_full)
            iNoClip = np.abs(resids) < _sigmaclip*np.std(resids) + np.mean(resids)
            iClip = np.abs(resids) > _sigmaclip*np.std(resids) + np.mean(resids)

        plt.plot(x_full[iNoClip],y_full[iNoClip],'.',label='not clipped',ms=1)
        plt.plot(x_full[iClip],y_full[iClip],'.',label='clipped',ms=1)
        #plt.scatter(x_full[iNoClip],y_full[iNoClip],c=yerr_full[iNoClip],zorder=100,s=1)
        magbins = np.linspace(colorcuts[swope_filt][0],colorcuts[swope_filt][1],10)
        try:
            binned_mags = binned_statistic(x_full[iNoClip],y_full[iNoClip],bins=magbins,statistic='median').statistic
            binned_magerrs = binned_statistic(x_full[iNoClip],y_full[iNoClip],bins=magbins,statistic=errfnc).statistic
            plt.errorbar((magbins[1:]+magbins[:-1])/2.,binned_mags,yerr=binned_magerrs,fmt='o-',label='median bins')
        except ValueError:
            print('everything clipped???')
        plt.ylabel(f'${swope_filt}_{{Swope}}$ - ${swope2ps[swope_filt][0]}_{{PS1}}$',fontsize=15)
        plt.xlabel(f'${swope2ps[swope_filt][1]}_{{PS1}} - {swope2ps[swope_filt][2]}_{{PS1}}$',fontsize=15)

        if npars == 2:
            plt.plot(magbins,line_fit_model(result_full.x,magbins),color='k',
                     label=f'y = {result_full.x[0]:.3f}*x + {result_full.x[1]:.3f}',zorder=101)
            print(f'y = {result_full.x[0]:.3f}*x + {result_full.x[1]:.3f}')
        elif npars == 3:
            plt.plot(magbins,line_fit_model(result_full.x,magbins),color='k',
                     label=f'y = {result_full.x[0]:.3f}*x^2 + {result_full.x[1]:.3f}*x + {result_full.x[2]:.3f}',zorder=101)
            print(f'y = {result_full.x[0]:.3f}*x^2 + {result_full.x[1]:.3f}*x + {result_full.x[2]:.3f}')
        else:
            plt.plot(magbins,line_fit_model(result_full.x,magbins),color='k',
                     label=f'y = {result_full.x[0]:.3f}*x^3 + {result_full.x[1]:.3f}*x^2 + {result_full.x[2]:.3f}*x + {result_full.x[3]:.3f}',zorder=101)
        plt.legend()

        
        magbins_pred = line_fit_model(result_full.x,(magbins[1:]+magbins[:-1])/2.)
        
        try:
            plt.title(f"$\sigma$ = {np.std(resids[iNoClip]):.3f}, max binned dev = {np.max(np.abs(magbins_pred[binned_mags == binned_mags]-binned_mags[binned_mags == binned_mags])):.3f}")
        except:
            pass
            
        import pdb; pdb.set_trace()
        
    def main(self,do_matching=True):
        for sfilt,pfilts in zip(
                swope2ps.keys(),swope2ps.values()
                ):

            if do_matching:
                # swope files
                swope_files = glob.glob(f'SwopeCats/{sfilt}/*dat')
                # median ra/dec for each
                field_info = []
                for s in swope_files:
                    data = at.Table.read(
                        s,format='ascii',names=('ra','dec','flux','fluxerr')
                    )
                    med_ra,med_dec = np.median(data['ra']),np.median(data['dec'])
                    field_info += [[med_ra,med_dec,s.split('/')[-1].split('.')[0]],]

                # get PS cats
                self.getPScats(field_info)

                # do the swope/PS catalog matching
                # write out the matched files to save time later(?)
                matched_table = self.MatchCats(swope_files,sfilt)
            else:
                matched_table = at.Table.read(f'MatchedCats/{sfilt}_matched.txt',format='ascii')
                
            # let's fit for the zeropoints and slopes
            self.fit_zpts(matched_table,sfilt)
            
            # now do some plotting
            
            
if __name__ == "__main__":
    sc = SwopeColors()
    sc.main()
