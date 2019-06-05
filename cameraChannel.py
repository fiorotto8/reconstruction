#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
import numpy as np
import math

class cameraGeometry:
    def __init__(self):
        # 240 mm = major axis of the ellipse = 2048 pixels
        self.pixelwidth = 117E-3 # mm
        
class cameraTools:
    def __init__(self):
        pass
        
    def isGoodChannelSafe(self,pedmap,ix,iy):
        pedvals = []; pedrmss = []
        for x in range(ix-1,ix+2):
            for y in range(iy-1,iy+2):
                pedvals.append(pedmap.GetBinContent(x,y))
                pedrmss.append(pedmap.GetBinError(x,y))
        if any([p>110 for p in pedvals]): return False
        if any([r<0.2 for r in pedvals]): return False
        #if pedrms > 5: return False
        return True

    def isGoodChannelFast(self,pedval,pedrms):
        if pedval > 110: return False
        if pedrms < 0.2: return False
        if pedrms > 5: return False
        return True
    
    def zs(self,th2,pedmap,nsigma=3,plot=False):
        #print "zero suppressing..."
        nx = th2.GetNbinsX(); ny = th2.GetNbinsY();
        xmin,xmax=(th2.GetXaxis().GetXmin(),th2.GetXaxis().GetXmax())
        ymin,ymax=(th2.GetYaxis().GetXmin(),th2.GetYaxis().GetXmax())
        th2_zs = ROOT.TH2D(th2.GetName()+'_zs',th2.GetName()+'_zs',nx,xmin,xmax,ny,ymin,ymax)
        th2_unzs = th2_zs.Clone(th2.GetName()+'_unzs')
        th2_zs.SetDirectory(None); th2_unzs.SetDirectory(None);
        for ix in range(1,nx+1):
            for iy in range(1,ny+1):
                ped   = pedmap.GetBinContent(ix,iy)
                noise = pedmap.GetBinError(ix,iy)
                if not self.isGoodChannelFast(ped,noise): continue
                z = th2.GetBinContent(ix,iy)-ped
                th2_unzs.SetBinContent(ix,iy,z)
                #print("x,y,z=",ix," ",iy," ",z,"   noise = ",noise , "ped-ib = ", ped_ixb," - ", ped_ixb)
                if z>nsigma*noise:
                    th2_zs.SetBinContent(ix,iy,z)
                    #print "x,y,z=",ix," ",iy," ",z,"   noise = ",noise
        #th2_zs.GetZaxis().SetRangeUser(0,1)
        if plot:
            canv = ROOT.TCanvas('zs','',600,600)
            th2_zs.Draw('colz')
            for ext in ['pdf','root']:
                canv.SaveAs('{name}.{ext}'.format(name=th2.GetName()+'_zs',ext=ext))
        return (th2_zs,th2_unzs)

    def zshits(self,hits,pedmap,nsigma=3,plot=False):
        print("zero suppressing a subset of ",len(hits)," hits...")
        ret=[]
        for h in hits:
            x,y,z=h[0],h[1],h[2]
            # warning: this works only for FR coordinates, otherwise need the findbin etc...
            ped = pedmap.GetBinContent(x+1,y+1)
            noise = pedmap.GetBinError(x+1,y+1)
            z_ps = z-ped
            if self.isGoodChannelFast(ped,noise) and z_ps>nsigma*noise and z_ps<30:
                ret.append(x,y,z_ps)
            else:
                ret.append(x,y,0)
        return ret
    
    def getData(self,th2):
        if not th2.InheritsFrom("TH2"):
            print("ERROR! The input object should be a TH2")
            return []
        x_bins = th2.GetNbinsX()
        y_bins = th2.GetNbinsY()
        bins = np.zeros((x_bins,y_bins))
        for y_bin in range(y_bins): 
            for x_bin in range(x_bins): 
                z = th2.GetBinContent(x_bin + 1,y_bin + 1)
                bins[x_bin,y_bin] = z
        return bins

    def getActiveCoords(self,th2):
        ret = []
        if not th2.InheritsFrom("TH2"):
            print("ERROR! The input object should be a TH2")
            return ret
        x_bins = th2.GetNbinsX()
        y_bins = th2.GetNbinsY()
        for y_bin in range(y_bins): 
            for x_bin in range(x_bins): 
                x = th2.GetXaxis().GetBinCenter(x_bin+1)
                y = th2.GetYaxis().GetBinCenter(y_bin+1)
                z = th2.GetBinContent(x_bin + 1,y_bin + 1)
                if z>0:
                    ret.append((x,y,math.log(z)))
        return np.array(ret)

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
        # canv = ROOT.TCanvas('zs','',600,600)
        # th2_rs.Draw('colz')
        # for ext in ['png','pdf','root']:
        #     canv.SaveAs('{name}_rs.{ext}'.format(name=th2.GetName(),ext=ext))
        return th2_rs
