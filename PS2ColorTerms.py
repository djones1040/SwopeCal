#!/usr/bin/env python
# D. Jones - 6/6/22
"""Get color terms between PS1 and PS2"""

import glob
import os
import numpy as np
import astropy.table as at
from txtobj import txtobj
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import matplotlib
matplotlib.use('agg')
import pylab as plt
from scipy.optimize import least_squares, minimize
from scipy.stats import binned_statistic
import argparse

_pipe_work = os.path.expandvars('$PIPE_DATA/workspace')
_pipe_cat = os.path.expandvars('$PIPE_DATA/absphotcat/all')
_masterlist = os.path.expandvars('$PIPE_CONFIGDIR/ysegpc2.txt')

class ColorTerms:
    def __init__(self):
        # max allowable mag error
        self.max_magerr = 0.05

        # make sure DCMP stars aren't too crowded
        self.crowding_min_sep_arcsec = 5

        # make sure the catalog match is ok
        self.cat_max_sep_arcsec = 1.5

        # save things in a giant astropy table?
        # ra dec g gerr r rerr i ierr z zerr g_cat gerr_cat r_cat rerr_cat i_cat ierr_cat  z_cat zerr_cat n_det
        self.StarMatch = at.Table(
            [[]]*29,
            names=('ra','dec','g','gerr','gvar','r','rerr','rvar','i','ierr','ivar','z','zerr','zvar',
                   'g_cat','gerr_cat','r_cat','rerr_cat','i_cat','ierr_cat','z_cat','zerr_cat',
                   'n_det_g','n_det_r','n_det_i','n_det_z','skycell','transient_id','filename'))
        self.StarMatch['skycell'] = self.StarMatch['skycell'].astype('S10')
        self.StarMatch['transient_id'] = self.StarMatch['transient_id'].astype('S10')
        self.StarMatch['filename'] = self.StarMatch['filename'].astype('S100')

    def add_args(self,parser=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage='', conflict_handler="resolve")

        parser.add_argument('-t','--telescope', default='PS2', type=str,
                            help='telescope (PS1 or PS2)')
        parser.add_argument('--docat', default=False, action="store_true",
                            help='create the star match catalog')

        return parser

    def get_obs_mags(self):

        # master list, with skycells
        masterdata = at.Table.read(_masterlist,format='ascii')
        skycelldata = np.loadtxt(_masterlist,usecols=[15],dtype=str)

        # get PS1 catalogs
        ps1_catfiles = glob.glob(f'{_pipe_cat}')

        # loop over filters
        # Grab a bunch of DCMP files for PS2
        for filt in 'griz':
            dcmp_files = glob.glob(f"{_pipe_work}/*/{filt}/*dcmp")
            #for i,dcf in enumerate(dcmp_files):
            #    if '1325.090' in dcf: print(i)
            #import pdb; pdb.set_trace()
            for idx,dcf in enumerate(dcmp_files[0:10000]):
                print(idx,filt,dcf)
                # template files are all PS1
                if '_tmpl' in dcf: continue

                dc = txtobj(dcf,cmpheader=True)
                try:
                    zpt = fits.getval(dcf,'ZPTMAG')
                except:
                    # zpt failed
                    continue

                # and see if zpt is bad:
                if zpt <= 0: continue

                # get image dimensions
                imshapex = fits.getval(dcf,'NAXIS1')
                imshapey = fits.getval(dcf,'NAXIS2')

                # get the appropriate catalog
                tid = dcf.split('/')[-3]
                try:
                    skycell = skycelldata[masterdata['ID'] == tid][0] #.value[0]
                except:
                    print(f'missing skycell {skycell}')
                    continue
                try:
                    ps = at.Table.read(f"{_pipe_cat}/{skycell}.PS1_3piDR2.photcat",format='ascii')
                except:
                    print(tid)
                    print(f'missing photcat for {tid}')
                    continue
                    #raise RuntimeError(f'missing photcat for {tid}')
                # set up some stuff for coordinate matching
                wcs = WCS(dcf)
                sccat = SkyCoord(ps['ra'],ps['dec'],unit=u.deg)
                scexist = SkyCoord(self.StarMatch['ra'],self.StarMatch['dec'],unit=u.deg)

                # loop through the DCMP, add the stars passing cuts to our master table
                for j,x in enumerate(dc.Xpos):
                    # Create a master set of unique, "good" stars:
                    # 1. uncertainties reasonable
                    if dc.dM[j] > self.max_magerr: continue
                    # 2. not near edges
                    if dc.Xpos[j] < 10 or dc.Ypos[j] < 10 or dc.Xpos[j] > imshapex-10 or dc.Ypos[j] > imshapey-10: continue
                    # 3. relatively far from another nearby star
                    #    - could do something complicated re: brighter star contamination
                    #    - let's wait for a future iteration
                    sep_pix = np.sqrt((dc.Xpos-dc.Xpos[j])**2.+(dc.Ypos-dc.Ypos[j])**2.)
                    if np.min(sep_pix[sep_pix > 0]) < self.crowding_min_sep_arcsec*4: continue

                    # 4. DOPHOT type 1
                    if dc.type[j] != '0x00000001': continue
                    

                    ra,dec = wcs.wcs_pix2world(dc.Xpos[j],dc.Ypos[j],1)
                    sc = SkyCoord(ra,dec,unit=u.deg)

                    # 5. Have PS1 catalog matches in griz (is this ok?)
                    sep = sc.separation(sccat).arcsec
                    if np.min(sep) > self.cat_max_sep_arcsec: continue
                    # make sure catalog (deeper than DCMPs) doesn't have nearby stars
                    if np.sort(sep)[1] < self.crowding_min_sep_arcsec: continue
                    iSep = np.where(sep == np.min(sep))[0]
                    if len(iSep) != 1: continue
                    if ps['g'][iSep] < 0 or ps['r'][iSep] < 0 or ps['i'][iSep] < 0 or ps['z'][iSep] < 0:
                        continue




                    mag = -2.5*np.log10(dc.flux[j])+zpt
                    magerr = 1.085736*dc.dflux[j]/dc.flux[j]

                    # check for matches
                    ExistSep = sc.separation(scexist).arcsec
                    iExist = np.where(ExistSep < self.cat_max_sep_arcsec)[0]
                    if len(iExist) == 1:
                        # update the mean if object exists
                        if self.StarMatch[filt][iExist] > 0:
                            mag_new = (self.StarMatch[filt][iExist]*self.StarMatch['n_det_'+filt][iExist] + mag)/(self.StarMatch['n_det_'+filt][iExist]+1)
                            magerr = np.sqrt(self.StarMatch[filt+'err'][iExist]**2.+magerr**2.)/2
                            self.StarMatch['n_det_'+filt][iExist] += 1
                            # Welford's algorithm, thanks Wikipedia!
                            self.StarMatch[filt+'var'][iExist] = (self.StarMatch[filt+'var'][iExist]*(self.StarMatch['n_det_'+filt][iExist]-1) + \
                                    (self.StarMatch[filt][iExist] - mag)*(self.StarMatch[filt][iExist] - mag_new))/self.StarMatch['n_det_'+filt][iExist]
                            #if np.sqrt(self.StarMatch[filt+'var'][iExist]) > 1:
                            #    import pdb; pdb.set_trace()
                            # weighted avg for duplicates?  straight avg is easier to code: (old_val*ndet + new_val)/(ndet+1)
                            #   - we can wait for a future iteration.  unsure if weighted is even a good idea
                        self.StarMatch[filt][iExist] = mag_new
                        self.StarMatch[filt+'err'][iExist] = magerr
                    elif not len(iExist):
                        # if object doesn't exist, make a new row in the table
                        self.StarMatch.add_row([ps['ra'][iSep].value[0],ps['dec'][iSep].value[0],-999,-999,0,-999,-999,0,-999,-999,0,-999,-999,0,
                                                ps['g'][iSep].value[0],ps['dg'][iSep].value[0],ps['r'][iSep].value[0],
                                                ps['dr'][iSep].value[0],ps['i'][iSep].value[0],ps['di'][iSep].value[0],
                                                ps['z'][iSep].value[0],ps['dz'][iSep].value[0],0,0,0,0,skycell,tid,dcf])
                        self.StarMatch[filt][-1] = mag
                        self.StarMatch[filt+'err'][-1] = magerr
                        self.StarMatch['n_det_'+filt][-1] = 1
                    else:
                        import pdb; pdb.set_trace()
                        raise RuntimeError('there are duplicates in the catalog')
        # write out the table
        self.StarMatch.write(self.options.filename,format='ascii')
        
    def make_mag_plots(self,zoom=False,residuals=False,filename='psmagplots.png',gimin=0.2,gimax=1.0):
        
        plt.rcParams['figure.figsize'] = (8,4)
        plt.clf()
        plt.subplots_adjust(hspace=0.4,wspace=0.6,top=0.9,right=0.9)

        for filt1,filt2,j in zip('griz','rgri',range(4)):
            ax = plt.subplot(2,2,j+1)

            if filt1 == 'g': 
                iGood = (self.StarMatch[filt1] > 0) & (self.StarMatch[filt2] > 0) & (self.StarMatch[filt1+'var'] > 0) & (self.StarMatch[filt1+'var'] < 0.01) &\
                        (self.StarMatch['g_cat']-self.StarMatch['r_cat'] > 0.0) & (self.StarMatch['g_cat']-self.StarMatch['r_cat'] < 1.5)
                magdiff = self.StarMatch['g_cat'][iGood]-self.StarMatch['r_cat'][iGood]
                label = '$g_{PS1} - r_{PS1,cat}$'

            if filt1 == 'r': 
                iGood = (self.StarMatch[filt1] > 0) & (self.StarMatch[filt2] > 0) & (self.StarMatch[filt1+'var'] > 0) & (self.StarMatch[filt1+'var'] < 0.01) &\
                        (self.StarMatch['g_cat']-self.StarMatch['i_cat'] > 0.0) & (self.StarMatch['g_cat']-self.StarMatch['i_cat'] < 1.5)
                magdiff = self.StarMatch['g_cat'][iGood]-self.StarMatch['r_cat'][iGood]
                label = '$g_{PS1} - r_{PS1,cat}$'

            if filt1 == 'i': 
                iGood = (self.StarMatch[filt1] > 0) & (self.StarMatch[filt2] > 0) & (self.StarMatch[filt1+'var'] > 0) & (self.StarMatch[filt1+'var'] < 0.01) &\
                        (self.StarMatch['r_cat']-self.StarMatch['i_cat'] > 0.0) & (self.StarMatch['r_cat']-self.StarMatch['i_cat'] < 0.6)
                magdiff = self.StarMatch['r_cat'][iGood]-self.StarMatch['i_cat'][iGood]
                label = '$r_{PS1} - i_{PS1,cat}$'

            if filt1 == 'z': 
                iGood = (self.StarMatch[filt1] > 0) & (self.StarMatch[filt2] > 0) & (self.StarMatch[filt1+'var'] > 0) & (self.StarMatch[filt1+'var'] < 0.01) &\
                        (self.StarMatch['i_cat']-self.StarMatch['z_cat'] > -0.1) & (self.StarMatch['i_cat']-self.StarMatch['z_cat'] < 1.0)
                magdiff = self.StarMatch['i_cat'][iGood]-self.StarMatch['z_cat'][iGood]
                label = '$i_{PS1} - z_{PS1,cat}$'


            y = self.StarMatch[filt1][iGood]-self.StarMatch[filt1+'_cat'][iGood]
            yerr = np.sqrt(self.StarMatch[filt1+'err'][iGood]**2. + self.StarMatch[filt1+'err_cat'][iGood]**2.)
            xerr = np.sqrt(self.StarMatch[filt1+'err_cat'][iGood]**2. + self.StarMatch[filt2+'err_cat'][iGood]**2.)
            xycov = self.StarMatch[filt1+'err_cat'][iGood]**2.
            if filt1 == 'z': npars = 2
            else: npars = 2
            fit_pars,idx_noclip,idx_clip = \
                self.fit_color_terms(magdiff,y,xerr,yerr,xycov,npars=npars)
            if filt1 == 'g':
                print('hack!!')
                fit_pars[0] -= 0.008
                fit_pars[1] += 0.008
            elif filt1 == 'z':
                print('hack!!')
                #fit_pars[0] -= 0.008
                fit_pars[-1] -= 0.002
            def fit_fun(xt):
                model = 0
                for i,f in enumerate(fit_pars[::-1]):
                    model += f*xt**i
                return model

            if residuals:

                ax.plot(magdiff[idx_noclip],self.StarMatch[filt1][iGood][idx_noclip]-\
                        self.StarMatch[filt1+'_cat'][iGood][idx_noclip]-fit_fun(magdiff[idx_noclip]),'.',color='C0')
                ax.plot(magdiff[idx_clip],self.StarMatch[filt1][iGood][idx_clip]-\
                        self.StarMatch[filt1+'_cat'][iGood][idx_clip]-fit_fun(magdiff[idx_clip]),'.',color='C1')

                xbins = np.linspace(np.min(magdiff),np.max(magdiff),30)
                ybins = binned_statistic(magdiff[idx_noclip],self.StarMatch[filt1][iGood][idx_noclip]-\
                                         self.StarMatch[filt1+'_cat'][iGood][idx_noclip]-fit_fun(magdiff[idx_noclip]),
                                         bins=xbins,statistic='median').statistic
                ax.plot((xbins[1:]+xbins[:-1])/2.,ybins,'o-',color='k')
                ax.axhline(0,color='k')
                if 'gpc2' in _pipe_work:
                    ax.set_ylabel(f'(${filt1}_{{PS2}}-{filt1}_{{PS1,cat}}$) - model')
                else:
                    ax.set_ylabel(f'(${filt1}_{{PS1}}-{filt1}_{{PS1,cat}}$) - model')

                if zoom:
                    ax.set_ylim([-0.015,0.015])
                else:
                    ax.set_ylim([-0.2,0.2])

            else:
                ax.plot(magdiff[idx_noclip],self.StarMatch[filt1][iGood][idx_noclip]-\
                        self.StarMatch[filt1+'_cat'][iGood][idx_noclip],'.',color='C0')
                ax.plot(magdiff[idx_clip],self.StarMatch[filt1][iGood][idx_clip]-\
                        self.StarMatch[filt1+'_cat'][iGood][idx_clip],'.',color='C1')

                x = np.linspace(np.min(magdiff),np.max(magdiff),100)
                ax.plot(x,fit_fun(x),color='r',ls='--',zorder=5)
                ax.axhline(0,color='k')

                xbins = np.linspace(np.min(magdiff),np.max(magdiff),30)
                ybins = binned_statistic(magdiff[idx_noclip],self.StarMatch[filt1][iGood][idx_noclip]-self.StarMatch[filt1+'_cat'][iGood][idx_noclip],
                                         bins=xbins,statistic='median').statistic
                ax.plot((xbins[1:]+xbins[:-1])/2.,ybins,'o-',color='k')

                if 'gpc2' in _pipe_work:
                    ax.set_ylabel(f'${filt1}_{{PS2}}-{filt1}_{{PS1,cat}}$')
                else:
                    ax.set_ylabel(f'${filt1}_{{PS1}}-{filt1}_{{PS1,cat}}$')

                if zoom:
                    ax.set_ylim([-0.05,0.05])
                else:
                    ax.set_ylim([-0.2,0.2])
            ax.set_xlabel(label)
            ax.set_xlim([np.min(magdiff[idx_noclip])-0.1,np.max(magdiff[idx_noclip])+0.1])
            if fit_pars[0] < 0: sign1 = '-'
            else: sign1 = '+'
            if fit_pars[1] < 0: sign2 = '-'
            else: sign2 = '+'
            if len(fit_pars) > 2 and fit_pars[2] < 0: sign3 = '-'
            else: sign3 = '+'

            if npars == 3:
                if 'gpc2' in _pipe_work:
                    ax.text(0.01,0.9,f"${filt1}_{{PS2}} = {filt1}_{{PS1,cat}} {sign1} {np.abs(fit_pars[0]):.3f}({filt1}_{{PS1,cat}} - {filt2}_{{PS1,cat}})^2 {sign2}$\n${np.abs(fit_pars[1]):.3f}({filt1}_{{PS1,cat}} - {filt2}_{{PS1,cat}}) {sign3} {np.abs(fit_pars[2]):.3f}$",ha='left',va='bottom',transform=ax.transAxes,color='r',bbox={'facecolor':'1.0','alpha':1.0},fontsize=8,zorder=100)
                else:
                    ax.text(0.01,0.9,f"${filt1}_{{PS1}} = {filt1}_{{PS1,cat}} {sign1} {np.abs(fit_pars[0]):.3f}({filt1}_{{PS1,cat}} - {filt2}_{{PS1,cat}}) {sign2} {np.abs(fit_pars[1]):.3f}$",ha='left',va='bottom',transform=ax.transAxes,color='r',bbox={'facecolor':'1.0','alpha':1.0},fontsize=8,zorder=100)
            else:
                if 'gpc2' in _pipe_work:
                    ax.text(0.01,0.9,f"${filt1}_{{PS2}} = {filt1}_{{PS1,cat}} {sign1} {np.abs(fit_pars[0]):.3f}({filt1}_{{PS1,cat}} - {filt2}_{{PS1,cat}}) {sign2} {np.abs(fit_pars[1]):.3f}$",ha='left',va='bottom',transform=ax.transAxes,color='r',bbox={'facecolor':'1.0','alpha':1.0},fontsize=8,zorder=100)
                else:
                    ax.text(0.01,0.9,f"${filt1}_{{PS1}} = {filt1}_{{PS1,cat}} {sign1} {np.abs(fit_pars[0]):.3f}({filt1}_{{PS1,cat}} - {filt2}_{{PS1,cat}}) {sign2} {np.abs(fit_pars[1]):.3f}$",ha='left',va='bottom',transform=ax.transAxes,color='r',bbox={'facecolor':'1.0','alpha':1.0},fontsize=8,zorder=100)
                

        plt.savefig(filename,dpi=200)

    def fit_color_terms(self,x,y,xerr,yerr,xycov,sigmaclip=5,npars=2):
        """iteratively fit line and sigma-clip"""

        def line_fit_noclip(p):
            model = 0
            for i,f in enumerate(p[::-1]):
                model += f*x**i
            return y - model


        iNoClip = np.arange(len(y))
        for i in range(5):
            def line_fit(p):
                model = 0
                for i,f in enumerate(p[::-1]):
                    model += f*x[iNoClip]**i
                return (y[iNoClip] - model)/yerr[iNoClip]
            
            def loglike(p):
                try:
                    p_0,p_1 = p
                except:
                    p_0,p_1 = np.array(p)[0]

                lnlike = 0
                for x_i,y_i,sigma_x_i,sigma_y_i,cov_xiyi in zip(
                        x[iNoClip],y[iNoClip],xerr[iNoClip],yerr[iNoClip],xycov[iNoClip]):
                    model = np.matrix([[x_i],[p_0*x_i + p_1]])
                    y_i_mat = np.matrix([[x_i],[y_i]])
                    S_i = np.matrix([[sigma_x_i**2.,cov_xiyi],[cov_xiyi,sigma_y_i**2.]])
                    lnlike += -1/2.*(y_i_mat-model).T*np.linalg.inv(S_i)*(y_i_mat-model) \
                              - np.log(2*np.pi*np.sqrt(np.linalg.det(S_i)))

                return -float(lnlike)


            guess = [1]*npars
            result = least_squares(line_fit,guess)

            #result = minimize(loglike,[1,1])

            resids = line_fit_noclip(result.x)
            iNoClip = np.where(np.abs(resids) < sigmaclip*np.std(resids) + np.mean(resids))[0]
            iClip = np.where(np.abs(resids) > sigmaclip*np.std(resids) + np.mean(resids))[0]
            #break

        if npars == 2:
            print(f'y = {result.x[0]:.3f}*x + {result.x[1]:.3f}')
        elif npars == 3:
            print(f'y = {result.x[0]:.3f}*x^2 + {result.x[1]:.3f}*x + {result.x[2]:.3f}')

        return result.x,iNoClip,iClip

    def main(self,docat=False):

        # get DCMP mags
        if docat:
            self.get_obs_mags()
        else:
            self.StarMatch = at.Table.read(self.options.filename,format='ascii')
            
        # make some plots to find outliers and iterate on steps above
        for zoom,filename in zip([False,True][::-1],['psmagplots.png','psmagplots_zoom.png'][::-1]):
            self.make_mag_plots(zoom=zoom,filename=filename,residuals=False)

        # do some spline fits
        # find the "optimal" mag vs. color comparisons
        # report color terms
        #self.fit_color_terms()

if __name__ == "__main__":
    ct = ColorTerms()
    parser = ct.add_args()
    ct.options = parser.parse_args()

    if ct.options.telescope.lower() == 'ps1':
        ct.options.filename = 'GPC1_star_match.txt'
    elif ct.options.telescope.lower() == 'ps2':
        ct.options.filename = 'GPC2_star_match.txt'
    else:
        raise RuntimeError('telescope not recognized')

    ct.main(docat=ct.options.docat)
