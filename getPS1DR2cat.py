#!/usr/bin/env python
from __future__ import print_function
from astropy.io import fits
import sys, os,re,types,string,math,random
import optparse
import numpy as np
import requests
import time
import casjobs
import mastcasjobs
from PS1toINST import PS1toINSTclass
import astropy.table as at
from astropy.coordinates import SkyCoord
import astropy.units as u


def makepath(path,raiseError=1):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)

    return(0)




def makepath4file(filename,raiseError=1):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path,raiseError=raiseError))
    else:
        return(0)

def deg2sex(degrees, ra=False, outputformatRA='%02d:%02d:%06.3f',outputformatDEC='%1s%02d:%02d:%05.2f'):
    if type(degrees) is  bytes:
        degrees=float(degrees)
    if ra:
        # a.k.a. indeg and outhours
        if degrees < 0.0:
            while degrees<0:degrees+=360.
        if degrees > 360.0:
            while degrees>360.0:degrees-=360.
        degrees /= 15.

    if degrees < 0:
        sign = '-'
    else:
        sign = '+'

    degrees = abs(degrees)

    d1  = (degrees - (degrees % 1))
    rem = (degrees % 1) * 60
    d2  = rem - (rem % 1)
    srem = (rem % 1) * 60
    d3 = srem

    if ra:
      return outputformatRA % (d1, d2, d3)
    else:
      return outputformatDEC % (sign, d1, d2, d3)

### Converts sexigesimal notation to decimal degrees or decimal hours (if option 'ra=True' invoked)
def sex2deg(sexigecimal, ra=False):
    ### 2005/12/02 - AR: make sure it is in sexagesimal format!
    # is it a string? if not check if it is None
    if not (type(sexigecimal) is bytes):
        if type(sexigecimal) == None:
            raise RuntimeError("ERROR: sex2deg cannot handle 'None'!")
        return sexigecimal
    # Does it have ':' or ' '? If not, it must be a float in string format, just convert it and return
    if re.search('\:|\s',sexigecimal) == None:
        return(string.atof(sexigecimal))

    s1, s2, s3 = list(map(string.atof, re.split('[DHMSdhms:\s]', sexigecimal.strip())[0:3]))

    # Get the sign
    if re.search('-', sexigecimal):
        sign = -1
    else:
        sign = 1

    deg = abs(s1) + s2 / 60. + s3 / 3600.

    deg *= sign

    if ra:
        deg *= 15.

    return deg

# Returns the passed in RA in decimal degrees
# input RA can be in 'HH:MM:SS.ss', 'HH MM SS.ss' or in decimal degrees 
def RaInDeg(Ra):
    import types
    if type(Ra)==bytes:
        if re.search('[DHMShms: ]',Ra.strip()):
            return(sex2deg(Ra,ra=True))
    return(float(Ra))

# Returns the passed in Dec in decimal degrees
# input Dec can be in 'DD:MM:SS.ss', 'DD MM SS.ss' or in decimal degrees 
def DecInDeg(Dec):
    import types
    if type(Dec)==bytes:
        if re.search('[DHMShms: ]',Dec.strip()):
            return(sex2deg(Dec,ra=False))
    return(float(Dec))

def rmfile(filename,raiseError=1,gzip=False):
    " if file exists, remove it "
    if os.path.lexists(filename):
        os.remove(filename)
        if os.path.isfile(filename):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename)
            else:
                return(1)
    if gzip and os.path.lexists(filename+'.gz'):
        os.remove(filename+'.gz')
        if os.path.isfile(filename+'.gz'):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename+'.gz')
            else:
                return(2)
    return(0)

def rmfiles(filenames,raiseError=1,gzip=False):
    if not (type(filenames) is list):
        raise RuntimeError("List type expected as input to rmfiles")
    errorflag = 0
    for filename in filenames:
        errorflag |= rmfile(filename,raiseError=raiseError,gzip=gzip)
    return(errorflag)

class getPS1DR2catclass:
    def __init__(self):
        self.ra = None
        self.dec = None

        self.alldata={}

        self.filename = None

        self.filterlist=[]
        self.colorfilterlist=[]
        self.errorcolsFlag= True
        self.Nonestring = 'nan'

        # Eddies columns are named 'median(0)' and 'err(0)' for g, etc
        self.filt2index = {'g':0,'r':1,'i':2,'z':3,'y':4}
        self.supercal={'g':0.020,'r':0.033,'i':0.024,'z':0.028,'y':0.011}

        self.inst=''


    def autooutfilename(self):
        if self.ra==None or self.dec==None:
            raise RuntimeError('No RA/Dec defined to be used for auto filename!')

        filename = '%011.7f_%011.7f.PS1.' % (self.ra,self.dec)
        filename += ''.join(self.filterlist)
        filename += '.cat'
        
        return(filename)

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('-v', '--verbose', action="count", dest="verbose",default=0)
        parser.add_option('-d','--debug', default=False, action="store_true",
                          help='Debugging output')
        parser.add_option('--cathtml'  , default='http://faun.rc.fas.harvard.edu/eschlafly/wise/fitscat_qy', type="string",
                          help='root html for Eddies cats (default=%default)')
        parser.add_option('--tmpdir'  , default='/tmp', type="string",
                          help='Directory for temporary files (default=%default)')
        parser.add_option('--avoid_sn_region'  , default=False, action="store_true",
                          help='avoid sn region if set (default=%default)')
        parser.add_option('--eventfile'  , default=None, type="string",
                          help='event file, not required (default=%default)')
        parser.add_option('--eventid'  , default=None, type="string",
                          help='event ID, not required (default=%default)')
        parser.add_option('--event_sep_arcmin'  , default=1, type="float",
                          help='don\'t get stars w/i this many arcmin of our event, used if eventfile is not None (default=%default)')
        parser.add_option('-B',  default=False, action="store_true",
                          help='B band')
        parser.add_option('-V',  default=False, action="store_true",
                          help='V band')
        parser.add_option('-u',  default=False, action="store_true",
                          help='u band')
        parser.add_option('-g',  default=False, action="store_true",
                          help='g band')
        parser.add_option('-r',  default=False, action="store_true",
                          help='r band')
        parser.add_option('-i',  default=False, action="store_true",
                          help='i band')
        parser.add_option('-z',  default=False, action="store_true",
                          help='z band')
        parser.add_option('-y',  default=False, action="store_true",
                          help='y band')
        parser.add_option('--sexagesimal',  default=False, action="store_true",
                          help='output RA/Dec are in sexagesimal format')
        parser.add_option('-s','--show',  default=False, action="store_true",
                          help='print catalog to screen')
        parser.add_option('--keepnans',  default=False, action="store_true",
                          help='keep entries with NaN values')
        parser.add_option('--requiregriz',  default=False, action="store_true",
                          help='only keep objects that have griz measurements')
        parser.add_option('--requiregri',  default=False, action="store_true",
                          help='only keep objects that have griz measurements')
        parser.add_option('-o','--outfile'  , default=None , type="string",
                          help='file name of output file. If not specified, then automatic filename depending on RA,Dec and filters (default=%default)')
        parser.add_option('-c','--clobber',  default=False, action="store_true",
                          help='overwrite file if it already exists')
        parser.add_option('-s','--size'  , default=None , type="string",
                          help='sidelength of box in degree. If format AxB, then rectangle with sidelength A,B along Ra, Dec, respectively(default=%default)')
        parser.add_option('--skipsave',  default=False, action="store_true",
                          help='Don\'t save the catalog')
        parser.add_option('--Mmax'  , default=None , type="float",
                          help='cut all objects with M<Mmin for specified filters (default=%default)')
        parser.add_option('--Mmin'  , default=None , type="float",
                          help='cut all objects with M>Mmax for specified filters (default=%default)')
        parser.add_option('--dMmax'  , default=None , type="float",
                          help='cut all objects with dM<dMmin for specified filters (default=%default)')
        parser.add_option('--dMmin'  , default=None , type="float",
                          help='cut all objects with dM>dMmax for specified filters (default=%default)')
        parser.add_option('--decam'  , default=False , action="store_true",
                          help='convert the PS1 mag to DECAM(DES)')
        parser.add_option('--megacam'  , default=False , action="store_true",
                          help='convert the PS1 mag to MEGACAM(SNLS)')
        parser.add_option('--swope'  , default=False , action="store_true",
                          help='convert the PS1 mag to SWOPE(CSP)')
        parser.add_option('--ps2'  , default=False , action="store_true",
                          help='convert the PS1 mag to PS2')
        parser.add_option('--transfdir'  , default=False , type="string",
                          help='path to the PS1transformations dir (default=%default)')
        parser.add_option('--skip_sc'  , default=False ,action="store_true",
                          help='Skip the supercal correction to the PS1 mag (default=%default)')
        parser.add_option('--spline'  , default=False ,action="store_true",
                          help='Use spline fit (default=linear)')
        parser.add_option('--maxtries'  , default=2 ,type="int",
                          help='number of times to try again if download fails (default=linear)')
        parser.add_option('--nocutsinfile'  , default=False , action="store_true",
                          help='for multiple filters, can be helpful to ignore the dMmax cuts in the ASCII file trimming stage')

        return(parser)

    def check4TRANSFdir(self):

        if self.options.transfdir: self.path2name=self.options.transfdir
        else:
            if 'PIPE_PS1TRANSFDIR' in os.environ: 
                if os.path.isdir(os.getenv('PIPE_PS1TRANSFDIR')): self.path2name=os.getenv('PIPE_PS1TRANSFDIR')+'/PS1trans_linear.txt'
                else:
                    print('WARNING: %s does not exist, looking into %s' %(os.getenv('PIPE_PS1TRANSFDIR'),os.getcwd()+'/PS1transformations'))
                    if os.path.isdir('./PS1transformations'): 
                        self.path2name='./PS1transformations/PS1trans_linear.txt' 
                    else: sys.exit('ERROR: Unable to find PS1transformations dir using the current path: %s. Try again with --transfdir option specifing the new path for the dir' % os.getcwd())
            else:
                print('WARNING: setenv PIPE_PS1TRANSFDIR does not exist, looking into %s' %os.getcwd()+'/PS1transformations')
                if os.path.isdir('./PS1transformations'): self.path2name='./PS1transformations/PS1trans_linear.txt' 
                else: sys.exit('ERROR:  Unable to find PS1transformations dir using the current path: %s. Try again  with --transfdir option specifing the new path for the dir' % os.getcwd())


    def getfilterlist_from_options(self):

        if self.options.decam:  
            self.inst='DES0'
        if self.options.megacam:  
            self.inst='SNLS0'
            self.options.spline=True
            self.options.skip_sc=True
        if self.options.swope:
            self.inst='CSP1'
        if self.options.ps2:
            self.inst='PS2'

        if self.options.u and self.inst == None: 
            print('> PS1 has no u band. Use --decam or --decam option.')
            sys.exit(0)
       
        if not self.inst or self.options.ps2:
            if self.options.g: self.filterlist.append('g')
            if self.options.r: self.filterlist.append('r')
            if self.options.i: self.filterlist.append('i')
            if self.options.z: self.filterlist.append('z')
            if self.options.y: self.filterlist.append('y')
        self.filterlist_output = self.filterlist[:]

        if self.options.requiregriz and 'y' in self.filterlist:
            self.filterlist = ['g','r','i','z','y']
        elif self.options.requiregriz:
            self.filterlist += ['g','r','i','z']
            self.filterlist = np.unique(self.filterlist)
        elif self.options.requiregri:
            self.filterlist += ['g','r','i']
            self.filterlist = np.unique(self.filterlist)
        elif self.filterlist == []:
            self.filterlist = ['g','r','i','z','y']
            self.filterlist_output = self.filterlist[:]
        if self.options.decam or self.options.megacam or self.options.swope or self.options.ps2:
            self.filterlist_output = self.filterlist[:]

    def add2list(self,alldata,col1,newdata,col2):
        Nrowsnew = len(newdata[col2])
        if col1 in alldata:
            Nrows1 = len(alldata[col1])
            if Nrowsnew>0:
                tmp=np.zeros((Nrows1+Nrowsnew))
                tmp[:Nrows1]=alldata[col1]
                tmp[Nrows1:]=newdata[col2]
                alldata[col1]=tmp
        else:
            if Nrowsnew>0:
                alldata[col1]=newdata[col2]
        return(0)




    def fluxcolname(self,filt):
        return('median(%d)' % self.filt2index[filt])

    def errcolname(self,filt):
        return('err(%d)' % self.filt2index[filt])    
 
    def getdataforRADEC(self,ra,dec,alldata=None,urldonelist=None,ramin=None,ramax=None,decmin=None,decmax=None,tmpfiledir='.',dMmax=0.05,filt=None):
        tmpfiledir = re.sub('\/$','',os.path.abspath(os.path.expanduser(tmpfiledir)))

        if alldata==None:
            alldata=self.alldata

        if ra<0.0:ra+=360.0
        if ra>=360.0:ra-=360.0

        if ramin!=None:
            if (ra-ramin)<-180:
                ramin-=360.0
                ramax-=360.0
            elif (ra-ramin)>180:
                ramin+=360.0
                ramax+=360.0

        if self.options.avoid_sn_region:
            import pipeclasses
            params = pipeclasses.paramfileclass()
            params.loadfile(os.environ['PIPE_PARAMS'])
            params.loadfile(os.environ['EXTRAPARAMFILE'],addflag=True)

            eventdata = at.Table.read(params.get('ELC_EVENTLIST'),format='ascii')
            # I'm gonna pay for this line later...
            eventid = self.options.outfile.split('/')[-1].split('.')[0]
            iEvt = eventdata['ID'] == eventid
            eventra,eventdec = eventdata['ra'][iEvt][0],eventdata['dec'][iEvt][0]
            scsn = SkyCoord(eventra,eventdec,unit=u.deg)

        # root http
        http_addr=re.sub('\/$','',self.options.cathtml)

        # just in case that ra,dec are integers, add some wiggle room...
        epsil = 0.000001
        # name of the catalog
        web_name='cat.l='+repr(int(math.floor(ra)))+'_'+repr(int(math.ceil(ra+epsil)))+'.b='+repr(int(math.floor(dec)))+'_'+repr(int(math.ceil(dec+epsil)))+'.fits'

        if self.options.verbose>1:
            if ramin!=None:
                print('RA=%.6f Dec=%.6f (ramin=%.6f ramax=%.6f): %s' % (ra,dec,ramin,ramax,web_name))
            else:
                print('RA=%.6f Dec=%.6f: %s' % (ra,dec,web_name))
                
        # put it all together into a url 
        url = '%s/%s.gz' % (http_addr,web_name)

        # Did we already add this url to alldata?
        if (type(urldonelist) is list) and (url in urldonelist):
            if self.options.verbose:
                print('%s already downloaded and added to cat, skipping...' % (web_name+'.gz'))
            return(0,url)

        
        tmpdir = re.sub('\/$','',self.options.tmpdir)
        if self.options.debug:
            tmpfile = '%s/%s.delme.txt' % (tmpfiledir,web_name) 
        else:
            tmpfile = '%s/%s.%d.delme.txt' % (tmpfiledir,web_name,random.randint(0,1000000)) 
        
        if self.options.verbose>1:
            print('txt file temporarely saved into %s' % tmpfile)

        makepath4file(tmpfile)

        if self.options.debug and os.path.isfile(tmpfile):
            print('DEBUG: skipping downloading %s into %s, it already exists!' % (url,tmpfile))
        else:
            rmfile(tmpfile)
            makepath4file(tmpfile)

            # casjobs query here

            query = """select o.objID, o.nDetections, o.raMean, o.decMean,
m.gMeanPSFMag, m.gMeanPSFMagErr, m.rMeanPSFMag, m.rMeanPSFMagErr,
m.iMeanPSFMag, m.iMeanPSFMagErr, m.iMeanKronMag, m.iMeanKronMagErr, m.zMeanPSFMag, m.zMeanPSFMagErr,
m.yMeanPSFMag, m.yMeanPSFMagErr

from ObjectThin o
inner join MeanObject m on o.objID=m.objID
where
    o.raMean between %.7f and %.7f
    and o.decMean between %.7f and %.7f
    and o.nDetections>1 and m.iMeanPSFMag - m.iMeanKronMag < 0.05 and m.iMeanPSFMagErr < %.2f
"""%(ramin,ramax,decmin,decmax,dMmax)

            if self.options.verbose:
                print('running quick query:')
                print(query)
            jobs = casjobs.CasJobs(userid='892987546',password='BossTent1',base_url="http://mastweb.stsci.edu/ps1casjobs/services/jobs.asmx")
            #jobs = mastcasjobs.MastCasJobs(username='djones1040',password='BossTent1',context='PanSTARRS_DR2')

            quicksuccess = False
            try:
                job_output = jobs.quick(query,context='PanSTARRS_DR2',task_name='PS1cat_ra%.7f_dec%.7f'%(ra,dec))
                quicksuccess = True
            except:
                print('quick job failed!!  submitting job to queue.  this might be slow...')

                dropquery = """DROP TABLE [PS1cat_ra%i_dec%i]
go"""%(ra,np.abs(dec))

                try:
                    droptable = jobs.quick(dropquery,context="MYDB",task_name="drop_table")
                except: pass

                mktablequery = """CREATE TABLE PS1cat_ra%i_dec%i(
objID bigint, nDetections float, raMean float, decMean float,
gMeanPSFMag float, gMeanPSFMagErr float, rMeanPSFMag float, rMeanPSFMagErr float,
iMeanPSFMag float, iMeanPSFMagErr float, iMeanKronMag float, iMeanKronMagErr float, zMeanPSFMag float, zMeanPSFMagErr float,
yMeanPSFMag float, yMeanPSFMagErr float)"""%(ra,np.abs(dec))
                
                query = """INSERT INTO MYDB.PS1cat_ra%i_dec%i select o.objID, o.nDetections, o.raMean, o.decMean,
m.gMeanPSFMag, m.gMeanPSFMagErr, m.rMeanPSFMag, m.rMeanPSFMagErr,
m.iMeanPSFMag, m.iMeanPSFMagErr, m.iMeanKronMag, m.iMeanKronMagErr, m.zMeanPSFMag, m.zMeanPSFMagErr,
m.yMeanPSFMag, m.yMeanPSFMagErr

from ObjectThin o
inner join MeanObject m on o.objID=m.objID
where
    o.raMean between %.7f and %.7f
    and o.decMean between %.7f and %.7f
    and o.nDetections>1 and m.iMeanPSFMag - m.iMeanKronMag < 0.05 and m.iMeanPSFMagErr < %.2f
"""%(ra,np.abs(dec),ramin,ramax,decmin,decmax,dMmax)
                
                quicktable = jobs.quick(mktablequery,context="MYDB",task_name="create_table")
                job_id = jobs.submit(query,context='PanSTARRS_DR2',task_name='PS1cat_ra%.7f_dec%.7f'%(ra,dec))
                status = jobs.monitor(job_id)
                tstart = time.time()
                while status[1] != 'finished' and status[1] != 'failed':
                    print('time elapsed: %i seconds'%(time.time()-tstart))
                    print('Job is %s...'%status[1])
                    time.sleep(30)
                if status[1] == 'failed':
                    raise RuntimeError('casjobs query failed!!')
                job_id = jobs.request_output('PS1cat_ra%i_dec%i'%(ra,np.abs(dec)),'CSV')
                status = jobs.monitor(job_id)
                tstart = time.time()
                while status[1] != 'finished' and status[1] != 'failed':
                    print('time elapsed: %i seconds'%(time.time()-tstart))
                    print('Job is %s...'%status[1])
                    time.sleep(30)
                if status[1] == 'failed':
                    raise RuntimeError('casjobs query failed!!')

                tmpcasroot,ext = os.path.splitext(tmpfile)
                tmpcasfile = '%s_casjobs.csv'%tmpcasroot
                jobs.get_output(job_id,tmpcasfile)
                print('saving job to %s'%tmpcasfile)

                droptable = jobs.quick(dropquery,context="MyDB",task_name="drop_table")

            if quicksuccess: job_outputlines = job_output.split('\n')
            else:
                job_outputlines = open(tmpcasfile).readlines()
                os.system('rm %s'%tmpcasfile)

            fout = open(tmpfile,'w')
            headerline = '#       ra         dec'
            for f in self.filterlist: headerline += '        %s      d%s'%(f,f)
            print(headerline,file=fout)
            for line in job_outputlines:
                line = line.replace('\n','')
                if quicksuccess and line.startswith('['):
                    cols = np.array(line.split(','))
                elif not quicksuccess and line.startswith('objID'):
                    cols = np.array(line.split(','))
                else:
                    if not line: continue
                    data = np.array(line.split(',')).astype(float)
                    if quicksuccess:
                        outline = '%.7f  %.7f'%(data[cols == '[raMean]:Float'][0],data[cols == '[decMean]:Float'][0])
                        for f in self.filterlist:
                            outline += '  %.4f  %.4f'%(data[cols == '[%sMeanPSFMag]:Float'%f][0],data[cols == '[%sMeanPSFMagErr]:Float'%f][0])
                    else:
                        outline = '%.7f  %.7f'%(data[cols == 'raMean'][0],data[cols == 'decMean'][0])
                        for f in self.filterlist:
                            outline += '  %.4f  %.4f'%(data[cols == '%sMeanPSFMag'%f][0],data[cols == '%sMeanPSFMagErr'%f][0])

                    if self.options.avoid_sn_region:
                        # we need to make sure ZPT stars are not within some distance
                        # from our SN.  This is a clumsy way to make sure that we're 
                        # avoiding complex backgrounds
                        sc = SkyCoord(float(outline.split()[0]),float(outline.split()[1]),unit=u.deg)
                        if scsn.separation(sc).arcmin > self.options.event_sep_arcmin:
                            print(outline,file=fout)
                    else:
                        print(outline,file=fout)
            fout.close()


        return(0,tmpfile)

    def mkcuts(self,tmpfile,dMmax=None):
        fin = open(tmpfile,'r')
        fout = open('%s.2'%tmpfile,'w')
        for line in fin:
            if line.startswith('#'):
                cols = np.array(line.split()[1:])
                print(line.replace('\n',''),file=fout)
            else:
                printline = True
                for f in self.filterlist:
                    data = np.array(line.split()).astype(float)
                    mag,magerr = data[cols == f][0],data[cols == 'd%s'%f][0]
                    if (mag < 0 or magerr > dMmax or magerr < 0) and not self.options.nocutsinfile: printline=False
                if printline:
                    print(line.replace('\n',''),file=fout)
        fin.close()
        fout.close()
        os.remove(tmpfile)
        os.rename('%s.2'%tmpfile,tmpfile)

    def checkCoords(self,tmpfile):
        inImCount = 0

        fin = open(tmpfile,'r')
        for line in fin:
            if line.startswith('#'):
                cols = np.array(line.replace('#','').split())
            else:
                data = np.array(line.split()).astype(float)
                ra,dec = data[cols == 'ra'][0],data[cols == 'dec'][0]
                if (ra > self.ramin and ra < self.ramax and dec > self.decmin and dec < self.decmax) or \
                   (ra > self.ramin+360 and ra < self.ramax+360 and dec > self.decmin and dec < self.decmax):
                    inImCount += 1
        fin.close()

        if inImCount > 1:
            return(1)
        else:
            return(0)

    def getradecbox(self, ra, dec):
        # keep track of RA,Dec
        self.ra  = RaInDeg(ra)
        self.dec = DecInDeg(dec)
        if self.options.verbose>1:
            print('RA: %.7f, Dec:%.7f' % (self.ra,self.dec))
        RAboxsize = DECboxsize = None
        if self.options.size!=None:
            # get the boxsize, and ra/dec limits
            if re.search('\w+x\w+',self.options.size):
                m=re.search('(^\S+)x(\S+$)',self.options.size)
                if m!=None:
                    RAboxsize = float(m.groups()[0])
                    DECboxsize = float(m.groups()[1])
                    print('box sizes: ',RAboxsize,DECboxsize)
                else: 
                    raise RuntimeError('Could not parse size %s for' % self.options.size)
            else:
                RAboxsize = DECboxsize = float(self.options.size)

            # get the maximum 1.0/cos(DEC) term: used for RA cut
            minDec = self.dec-0.5*DECboxsize
            if minDec<=-90.0:minDec=-89.9
            maxDec = self.dec+0.5*DECboxsize
            if maxDec>=90.0:maxDec=89.9

            invcosdec = max(1.0/math.cos(self.dec*math.pi/180.0),
                            1.0/math.cos(minDec  *math.pi/180.0),
                            1.0/math.cos(maxDec  *math.pi/180.0))

            # get the (conservative) boxlimits
            ramin = self.ra-0.5*RAboxsize*invcosdec
            ramax = self.ra+0.5*RAboxsize*invcosdec
            decmin = self.dec-0.5*DECboxsize
            decmax = self.dec+0.5*DECboxsize
            # check for url for center and all 4 corners...
            radeclist =[(self.ra,self.dec),
                        (ramin,decmin), 
                        (ramin,decmax), 
                        (ramax,decmin), 
                        (ramax,decmax)] 
            
            
            if ramax-ramin > 1 and decmax-decmin > 1:
                radeclist=[]
                for n in range(int(ramax-ramin)+1):
                    app_ra=ramin+n
                    for l in range(int(decmax-decmin)+1):
                        app_dec=decmin+l
                        radeclist.extend([(app_ra,app_dec)])
                    
            elif ramax-ramin > 1 and decmax-decmin < 1:
                radeclist=[]
                for n in range(int(ramax-ramin)+1):
                    app_ra=ramin+n
                    for dec in [decmin,decmax]:
                        app_dec=dec
                        radeclist.extend([(app_ra,app_dec)])

            elif ramax-ramin < 1 and decmax-decmin > 1:
                radeclist=[]
                for ra in [ramin,ramax]:
                    app_ra=ra
                    for l in range(int(decmax-decmin)+1):
                        app_dec=decmin+l
                        radeclist.extend([(app_ra,app_dec)])


        else:
            ramin = ramax = decmin = decmax = None
            radeclist = [(self.ra,self.dec)]

        self.ramin,self.ramax,self.decmin,self.decmax = ramin,ramax,decmin,decmax

    def getcatalog(self, ra, dec,switch, tmpfiledir = '.'):
        # keep track of RA,Dec
        self.ra  = RaInDeg(ra)
        self.dec = DecInDeg(dec)
        if self.options.verbose>1:
            print('RA: %.7f, Dec:%.7f' % (self.ra,self.dec))

        # Which filters?
        self.getfilterlist_from_options()
        if self.options.verbose>1:
            print('Filters:',self.filterlist)


        # list of urls for which the data is already parsed. so they don't have to be done again...
        urldonelist=[]
        #print(radeclist)

        print('### getting catalog at %f %f'% (self.ra,self.dec))
        (errorflag,tmpfile) = self.getdataforRADEC(self.ra,self.dec,urldonelist=urldonelist,ramin=self.ramin,ramax=self.ramax,decmin=self.decmin,decmax=self.decmax,tmpfiledir=tmpfiledir,dMmax=self.options.dMmax)

        return(0,tmpfile)

    def get_output_format_info(self,list_of_filters):
        # get cols ...
        cols=['ra','dec']
        colsformat = ['%11.7f','%11.7f']
        header = '#%11s %11s' % ('ra','dec')
        for filt in list_of_filters:
            colsformat.append('%7.4f')
            cols.append(filt)
            header += ' %7s' % (filt)
            if self.errorcolsFlag:
                colsformat.append('%7.4f')
                cols.append('d'+filt)
                header += ' %7s' % ('d'+filt)

        return(cols,colsformat,header)

    def getformattedcatalog(self,cols,colsformat, alldata=None, indices = None, Nonestring=None, addreturnFlag=False):
        if alldata==None:
            alldata=self.alldata

        if indices == None:
            indices = range(len(alldata['ra']))

        if Nonestring==None:
            Nonestring = self.Nonestring

        cs = range(len(cols))
        
        lines=[]
        for i in indices:
            line = ''
            nanflag = False
            for c in cs:                    
                #print len(alldata[cols[c]]),i
                if math.isnan(alldata[cols[c]][i]) or math.isinf(alldata[cols[c]][i]):
                    nanflag=True
                    line += ' %7s' %  Nonestring
                else:
                    if cols[c] in ['ra','dec'] and self.options.sexagesimal:
                        line += ' '+deg2sex(alldata[cols[c]][i],ra=(cols[c]=='ra'))
                    else:
                        line += ' '+colsformat[c] % alldata[cols[c]][i] 
                    
            # skip bad lines if not --keepnans
            if nanflag and (not self.options.keepnans):
                continue

            if addreturnFlag:
                line += '\n'

            lines.append(line)

        return(lines)
            

    def savecatalog(self, outfilename, keys = None, cols = None, errorcolsFlag=True, Nonestring='NaN',clobber=False):
        makepath4file(outfilename)

        # file already exists?
        if os.path.isfile(outfilename):
            if clobber:
                if self.options.verbose>0:
                    print('file %s exists, clobbering it...' % outfilename)
                rmfile(outfilename)
            else:
                print('file %s already exists, stopping (if you want to overwrite file, use --clobber)' % outfilename)
                return(1)

        if self.options.verbose>1:
            print('Getting formatted catalog...')

        (cols,colsformat,header) = self.get_output_format_info(self.filterlist)

        (lines) = self.getformattedcatalog(cols,colsformat,Nonestring=Nonestring,addreturnFlag=True)

        if self.options.verbose:
            print('Saving %d entries into %s' % (len(lines),outfilename))
        f=open(outfilename,'w')
        f.writelines([header+'\n'])
        f.writelines(lines)
        f.close()

        # keep track of filename
        self.filename = outfilename

        if not os.path.isfile(outfilename):
            raise RuntimeError('Could not write to %s' % outfilename)

        return(0)

    def showcatalog(self, keys = None, cols = None, errorcolsFlag=True, Nonestring='NaN'):
        (cols,colsformat,header) = self.get_output_format_info()
        (lines) = self.getformattedcatalog(cols,colsformat,Nonestring=Nonestring)
        print(header)
        for line in lines: print(line)
        return(0)
    
    def converter(self,outfilename):
        #added by gio 01/04/2016 convertion from PS1 to DECAM
        #surv_name=self.inst
        selected_band=[]
        #print self.filterlist
        if self.options.B: 
            selected_band.append('B')
            #print selected_band
        elif self.options.V: 
            selected_band.append('V')
            #print selected_band
        elif self.options.u: 
            selected_band.append('u')
            #print selected_band
        elif self.options.g: 
            selected_band.append('g')
            #print selected_band
        elif self.options.r: 
            selected_band.append('r')
            #print selected_band
        elif self.options.i: 
            selected_band.append('i')
            #print selected_band
        elif self.options.z: 
            selected_band.append('z')
            #print selected_band
        elif self.options.y: 
            selected_band.append('y')
            #print selected_band
        else:
            selected_band.append('u')
            selected_band.extend(self.filterlist)


        # hack: no conversion
        #PS1toINST=PS1toINSTclass()
        #PS1toINST.convert(self.path2name,self.inst,outfilename,selected_band,self.options.sexagesimal,self.options.spline,filterlist=self.filterlist)

    def make_output_catalog(self,tmpfile,outfilename):

        header_to_write = ['ra','dec']+np.concatenate([[f,'d'+f] for f in self.filterlist_output]).tolist()
        
        with open(tmpfile) as fin, open(outfilename,'w') as fout:
            print('#'+' '.join(header_to_write),file=fout)

            for line in fin:
                line = line.replace('\n','')
                if line.startswith('#'):
                    vals = np.array(line.replace('#','').split())
                else:
                    outline = ''
                    lineparts = np.array(line.split())
                    for hk in header_to_write:
                        outline += lineparts[vals == hk][0] + ' '
                    print(outline[:-1],file=fout)

if __name__=='__main__':

    getPS1cat=getPS1DR2catclass()
    parser = getPS1cat.add_options(usage='getPS1cat.py RA Dec')

    
    if len(sys.argv)<3:
        options,  args = parser.parse_args(args=['-h'])        
        sys.exit(0)

    
    # this is a hack to allow negative numbers as arguments (The one bug in optparse). 
    # this means that all options must be after the arguments.
    #getPS1cat.options,  args = parser.parse_args()
    getPS1cat.options,  dummy = parser.parse_args(args=sys.argv[3:])
    args = sys.argv[1:3]
    if len(dummy)!=0:
        print('ERROR: extra arguments!',dummy)
        sys.exit(0)

    (ra,dec)=args[:2]

    # hack -- no filter transformations please
    #getPS1cat.check4TRANSFdir()

    # outfilename
    if getPS1cat.options.outfile!=None:
        outfilename = getPS1cat.options.outfile
    else:
        outfilename = getPS1cat.autooutfilename()

    getPS1cat.getradecbox(ra,dec)
    if not getPS1cat.options.skipsave and not getPS1cat.options.clobber and os.path.exists(outfilename):
        success = getPS1cat.checkCoords(outfilename)
        if success:
            print('photcat %s exists and matches the image!!!  Not clobbering'%outfilename)
            print('SUCCESS %s' % os.path.basename(sys.argv[0]))
            sys.exit(0)
        else:
            print('photcat %s exists but coords don\'t match the image!!!  clobbering'%outfilename)
    # directory for temporary files
    tmpfiledir = os.path.dirname(os.path.abspath(os.path.expanduser(outfilename)))

    # get the catalog from the web
    err,tmpfile = getPS1cat.getcatalog(ra,dec,switch=0,tmpfiledir=tmpfiledir)

    getPS1cat.mkcuts(tmpfile,getPS1cat.options.dMmax)

    # make sure things are correct
    success = getPS1cat.checkCoords(tmpfile)
    count = 1
    while not success and count < getPS1cat.options.maxtries:
        # if it didn't work, try again
        err,tmpfile = getPS1cat.getcatalog(ra,dec,switch=0,tmpfiledir=tmpfiledir)
        getPS1cat.mkcuts(tmpfile,getPS1cat.options.dMmax)
        success = getPS1cat.checkCoords(tmpfile)
        count += 1
    if not success:
        os.system('rm %s'%tmpfile)
        raise RuntimeError('Error : getPS1DR2cat.py downloaded the wrong catalog??')

    # show the catalog?
    #if getPS1cat.options.show:
    #    getPS1cat.showcatalog()

    # save the catalog?
    if not getPS1cat.options.skipsave:
        if len(getPS1cat.filterlist) == len(getPS1cat.filterlist_output):
            os.system('mv %s %s'%(tmpfile,outfilename))
        else:
            getPS1cat.make_output_catalog(tmpfile,outfilename)
        #getPS1cat.savecatalog(outfilename,clobber=getPS1cat.options.clobber)    

    #move PS1cat to Surv?   
    if getPS1cat.options.decam or getPS1cat.options.megacam or getPS1cat.options.swope or getPS1cat.options.ps2:
        getPS1cat.converter(outfilename)
    
    print('SUCCESS %s' % os.path.basename(sys.argv[0]))
