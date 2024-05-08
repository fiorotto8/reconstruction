#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
import uproot
import numpy as np
import math
import debug_code.tools_lib as tl

class cameraGeometry:
    def __init__(self,params):
        self.pixelwidth = params['pixelwidth']
        self.npixx = params['npixx']
        self.name = params['name']
        self.vignette = params['vignette']
        self.xmin = params['xmin'] if 'xmin' in params else 0
        self.xmax = params['xmax'] if 'xmax' in params else self.npixx
        self.ymin = params['ymin'] if 'ymin' in params else 0
        self.ymax = params['ymax'] if 'ymax' in params else self.npixx
        
        
class cameraTools:
    def __init__(self,geometry):
        self.geometry = geometry
        # attach to a dict to make it persistent
        # the matrix is the max size possible, still ok if rebinned (because it is redone from the TH2D when it is readout)
        self.vignetteMap = { self.geometry.name : np.zeros((int(self.geometry.npixx),int(self.geometry.npixx))) }

    def pedsub(self,img,pedarr):
        return img - pedarr

    def satur_corr(self,img):
        e = 1.60217662e-7 # electric charge in C
        d2 = 0.015625 # mm^2
        omega = 0.00018 #solid  angle  covered  by  thephotocamera
        alpha = 0.08 # coefficient for varied gas mixture
        sigma0 = 2.5 # small variable of density of charge, we set it as limit for different fit parameters, unit of pC/mm^2
        a0 = 0.1855 # 0.2368 
        #fit function is y=ax^2+bx+c 
        a = a0*e/(d2*alpha*omega) # from fit of charge densities functions, converted to photon units
        b = (1.- 2*a0*sigma0)
        c= a0*sigma0*sigma0*(d2*alpha*omega)/e
        img0 = sigma0*(d2*alpha*omega)/e  # minimal desity of photons, like sigma0 but in photon unit 
        img_sat_corr = np.where(img > img0, (a*img*img + b*img + c), img)
        return img_sat_corr
    
    def zsfullres(self,img_sub,noisearr,nsigma=1):
        img_zs = np.where(img_sub > nsigma * noisearr, img_sub, 0)
        return img_zs
        
    def arrrebin(self,img,rebin):
        newshape = int(self.geometry.npixx/rebin)
        img_rebin = tl.rebin(img,(newshape,newshape))
        return img_rebin

    def acceptance(self,img,xmin,xmax,ymin,ymax):
        img[:xmin,:]=0
        img[xmax:,:]=0
        img[:,:ymin]=0
        img[:,ymax:]=0
        return img
        
    # this returns x,y,z before the zero suppression
    def getImage(self,th2):
        return np.array(th2)

    def noisearray(self,th2):
        noisearr = np.zeros( (th2.GetNbinsX(),th2.GetNbinsY()) )
        for ix in range(th2.GetNbinsX()):
            for iy in range(th2.GetNbinsY()):
                noisearr[ix][iy] = th2.GetBinError(ix+1,iy+1)
        return noisearr

    def getRestrictedImage(self,th2,xmin,xmax,ymin,ymax):
        nx = th2.GetNbinsX(); ny = th2.GetNbinsY();
        nxp = xmax-xmin
        nyp = ymax-ymin
        th2_rs = ROOT.TH2D(th2.GetName()+'_rs',th2.GetName()+'_rs',nxp,xmin,xmax,nyp,ymin,ymax)
        th2_rs.SetDirectory(None)
        for ix,x in enumerate(range(xmin,xmax)):
            for iy,y in enumerate(range(ymin,ymax)):
                orig_ixb = th2.GetXaxis().FindBin(x)
                orig_iyb = th2.GetYaxis().FindBin(y)
                z = th2.GetBinContent(orig_ixb,orig_iyb)
                th2_rs.SetBinContent(ix+1,iy+1,z)
        return th2_rs

    def loadVignettingMap(self):
        print ("Loading vignette map from: {vf}...".format(vf=self.geometry.vignette))
        det = self.geometry.name
        if det == 'lemon': # not implemented (we were taking the efficienct region within the FC)
            return self.vignetteMap[det]
        elif det == 'lime' or 'Mango_full':
            if not self.vignetteMap[det].any():
                tf = uproot.open(self.geometry.vignette)
                namehmap = 'normmap_'+self.geometry.name
                if det == 'Mango_full':
                       namehmap = 'normmap_lime'		#in vignette_runs... there is no mango_full, so lime is used
                vignetteMapRebinned = tf[namehmap].values()
                print("py ndim = ",vignetteMapRebinned.ndim)
                rebinx = int(self.geometry.npixx/vignetteMapRebinned.shape[0])
                rebiny = int(self.geometry.npixx/vignetteMapRebinned.shape[1])
                macroPixel = np.zeros((rebinx,rebiny))
                print ("Macro-pixel of the vignetting map has size = ",macroPixel.shape)
                for ibx in range(vignetteMapRebinned.shape[0]):
                    for iby in range(vignetteMapRebinned.shape[1]):
                        macroPixel[:,:] = 1./vignetteMapRebinned[iby,ibx]
                        (self.vignetteMap[det])[int(ibx*rebinx):int((ibx+1)*rebinx),int(iby*rebiny):int((iby+1)*rebiny)] = macroPixel
                return self.vignetteMap[det]
            else:
                return self.vignetteMap[det]
        else:
            print ('WARNING! Geometry ',det,' not foreseen. Return correction 1')
            return self.vignetteMap[det]

    def vignette_corr(self,img,vignette):
        return np.multiply(img,vignette)


