import ROOT
import uproot
import numpy as np
import pandas as pd
import glob
import re
from tqdm import tqdm
from multiprocessing import Pool
import argparse
import dask
import dask.array as da
from dask import delayed
import concurrent.futures

def get_numbers_from_filename(filename):
    return int(re.search(r'\d+', filename).group(0))
def grapherr(x,y,ex,ey,x_string, y_string,name=None, color=9, markerstyle=33, markersize=2,write=True):
        plot = ROOT.TGraphErrors(len(x),  np.array(x  ,dtype="d")  ,   np.array(y  ,dtype="d") , np.array(   ex   ,dtype="d"),np.array( ey   ,dtype="d"))
        if name is None: plot.SetNameTitle(y_string+" vs "+x_string,y_string+" vs "+x_string)
        else: plot.SetNameTitle(name, name)
        plot.GetXaxis().SetTitle(x_string)
        plot.GetYaxis().SetTitle(y_string)
        plot.SetMarkerColor(color)#blue
        plot.SetMarkerStyle(markerstyle)
        plot.SetMarkerSize(markersize)
        if write==True: plot.Write()
        return plot
def graph(x,y,x_string, y_string,name=None, color=9, markerstyle=33, markersize=2,write=True):
        plot = ROOT.TGraph(len(x),  np.array(x  ,dtype="d")  ,   np.array(y  ,dtype="d"))
        if name is None: plot.SetNameTitle(y_string+" vs "+x_string,y_string+" vs "+x_string)
        else: plot.SetNameTitle(name, name)
        plot.GetXaxis().SetTitle(x_string)
        plot.GetYaxis().SetTitle(y_string)
        plot.SetMarkerColor(color)#blue
        plot.SetMarkerStyle(markerstyle)
        plot.SetMarkerSize(markersize)
        if write==True: plot.Write()
        return plot
def fill_h(histo_name, array):
    for x in range (len(array)):
        histo_name.Fill((np.array(array[x] ,dtype="d")))
def hist(list, x_name, channels=1000, linecolor=4, linewidth=4,write=True,startZero=False,BinSpacing=None):
    array=np.array(list ,dtype="d")
    if BinSpacing is None: ch=channels
    else: ch = int((np.max(array) - np.min(array)) / BinSpacing)
    if startZero==False: hist=ROOT.TH1D(x_name,x_name,ch,0.99*np.min(array),1.01*np.max(array))
    else: hist=ROOT.TH1D(x_name,x_name,ch,0,1.01*np.max(array))
    fill_h(hist,array)
    hist.SetLineColor(linecolor)
    hist.SetLineWidth(linewidth)
    hist.GetXaxis().SetTitle(x_name)
    hist.GetYaxis().SetTitle("Entries")
    if write==True: hist.Write()
    #hist.SetStats(False)
    hist.GetYaxis().SetMaxDigits(3);
    hist.GetXaxis().SetMaxDigits(3);
    return hist
def nparr(list):
    return np.array(list, dtype="d")

#Parser
parser = argparse.ArgumentParser(description='Analyze pedestal', epilog='Version: 1.0')
parser.add_argument('-j','--cores',help='select # cores', action='store', type=int, default=1)
parser.add_argument('-t','--test',help='select test mode', action='store', type=bool, default=False)
args = parser.parse_args()
# Number of concurrent processes
num_processes = args.cores
#check all pedestals
if args.test==True: files=glob.glob("ped_test/*.root")
else: files=glob.glob("pedestals/*.root")
files.reverse()
#create Pool

"""
averages,averages_transposed,run_num=[],[],[]
mean_left,mean_right,middle_mean=[],[],[]
for file in tqdm(files[:12]):
    run_num.append(get_numbers_from_filename(file))
    root_file = uproot.open(file)
    mean_h2 = root_file["pedmap"]
    proj_means = mean_h2.to_numpy()[0]
    transposed_proj_means = np.transpose(proj_means)

    avg_proj_means = np.mean(proj_means, axis=1)
    avg_transposed_proj_means = np.mean(transposed_proj_means, axis=1)

    averages.append(avg_proj_means)
    averages_transposed.append(avg_transposed_proj_means)
    mean_left.append(np.mean(avg_proj_means[0:10]))
    mean_right.append(np.mean(avg_proj_means[-11:-1]))
    middle_mean.append(np.mean(avg_proj_means[1000:1304]))
"""
def process_file(file):
    run_num = get_numbers_from_filename(file)
    root_file = uproot.open(file)
    mean_h2 = root_file["pedmap"]
    proj_means = mean_h2.to_numpy()[0]
    transposed_proj_means = np.transpose(proj_means)

    avg_proj_means = np.mean(proj_means, axis=1)
    avg_transposed_proj_means = np.mean(transposed_proj_means, axis=1)

    return {
        'run_num': run_num,
        'avg_proj_means': avg_proj_means,
        'avg_transposed_proj_means': avg_transposed_proj_means,
        'mean_left': np.mean(avg_proj_means[:10]),
        'mean_right': np.mean(avg_proj_means[-11:]),
        'middle_mean': np.mean(avg_proj_means[1000:1304]),
    }

num_threads = args.cores  # Adjust the number of threads as needed

averages, averages_transposed, run_num = [], [], []
mean_left, mean_right, middle_mean = [], [], []

with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
    results = list(tqdm(executor.map(process_file, files), total=len(files)))

# Extracting results
for result in results:
    run_num.append(result['run_num'])
    averages.append(result['avg_proj_means'])
    averages_transposed.append(result['avg_transposed_proj_means'])
    mean_left.append(result['mean_left'])
    mean_right.append(result['mean_right'])
    middle_mean.append(result['middle_mean'])



main=ROOT.TFile("anal_ped_XY.root","RECREATE")
main.mkdir("Raw")
main.cd("Raw")
expo_left=ROOT.TF1("expo_left", "[0]*exp([1]*x)+[2]",0,400)
expo_right=ROOT.TF1("expo_right", "[0]*exp([1]*x)+[2]",1900,2304)
max_left, slope_left,max_right, slope_right=np.empty(len(averages)),np.empty(len(averages)),np.empty(len(averages)),np.empty(len(averages))
for i in range(len(averages)):
    temp=graph(np.arange(len(averages[i])), averages_transposed[i], "X position", "average noise", "left X projection "+str(run_num[i]),write=False)
    temp1=graph(np.arange(len(averages[i])), averages_transposed[i], "X position", "average noise", "right X projection "+str(run_num[i]),write=False)
    expo_left.SetParameters(averages[i][0],0,averages[i][1000])
    temp.Fit("expo_left","RQ","r")
    temp.Write()
    expo_right.SetParameters(1E-9,-1*expo_left.GetParameter(1),averages[i][1000])
    temp1.Fit("expo_right","RQ","r")
    temp1.Write()
    max_left[i],slope_left[i]=(expo_left.GetParameter(0)+expo_left.GetParameter(2)),expo_left.GetParameter(1)
    max_right[i],slope_right[i]=((expo_right.GetParameter(0)*np.exp(expo_right.GetParameter(1)*2304))+expo_right.GetParameter(2)),expo_right.GetParameter(1)





main.cd()
#CHECK VARIABLES
def canvas_history(x,y,x_string,y_string,name):
    graph=grapherr(x,y,1E-20*np.ones(len(run_num)),1E-20*np.ones(len(run_num)),x_string,y_string,name)
    can1=ROOT.TCanvas(name, name, 1200   ,1200)
    can1.SetFillColor(0)
    can1.SetBorderMode(0)
    can1.SetBorderSize(2)
    can1.SetLeftMargin(0.18)
    can1.SetRightMargin(0.02)
    can1.SetTopMargin(0.1)
    can1.SetBottomMargin(0.1)
    can1.SetFrameBorderMode(0)
    can1.SetFrameBorderMode(0)
    can1.SetFixedAspectRatio()
    #can1.cd()
    graph.Draw("AP")
    #graph.GetXaxis().SetDecimals(1)
    #graph.GetXaxis().SetMaxDigits(2)

    can1.Update()
    ymax=ROOT.gPad.GetUymax()
    ymin=ROOT.gPad.GetUymin()
    xmax=ROOT.gPad.GetUxmax()
    xmin=ROOT.gPad.GetUxmin()

    #lines 23802,24328,AmBe Campaign with also GEM Off
    line1=ROOT.TLine(23.802,ymin,23.802,ymax)
    line1.SetLineColor(40)
    line1.SetLineWidth(2)
    line1.Draw("SAME")
    line2=ROOT.TLine(25.427,ymin,25.427,ymax)
    line2.SetLineColor(40)
    line2.SetLineWidth(2)
    line2.Draw("SAME")
    pave=ROOT.TPaveText(23.802,ymax-0.02*(ymax-ymin),25.427,ymax-0.02*(ymax-ymin))
    pave.SetBorderSize(2)
    pave.SetFillColor(0)
    pave.SetTextSize(0.02)
    pave.SetTextColor(40)
    pave.AddText("AmBe Campaign")
    pave.Draw("SAME")

    #lines 27858,30831,Ba Campaign with also collimator off
    line3=ROOT.TLine(27.858,ymin,27.858,ymax)
    line3.SetLineColor(41)
    line3.SetLineWidth(2)
    line3.Draw("SAME")
    line4=ROOT.TLine(30.831,ymin,30.831,ymax)
    line4.SetLineColor(41)
    line4.SetLineWidth(2)
    line4.Draw("SAME")
    pave1=ROOT.TPaveText(27.858,ymax-0.02*(ymax-ymin),30.831,ymax-0.02*(ymax-ymin))
    pave1.SetBorderSize(2)
    pave1.SetFillColor(0)
    pave1.SetTextSize(0.02)
    pave1.SetTextColor(41)
    pave1.AddText("Ba Campaign")
    pave1.Draw("SAME")

    #lines 30845,36663,Ba Campaign with also collimator off
    line5=ROOT.TLine(30.845,ymin,30.845,ymax)
    line5.SetLineColor(42)
    line5.SetLineWidth(2)
    line5.Draw("SAME")
    line6=ROOT.TLine(36.663,ymin,36.663,ymax)
    line6.SetLineColor(42)
    line6.SetLineWidth(2)
    line6.Draw("SAME")
    pave2=ROOT.TPaveText(30.845,ymax-0.02*(ymax-ymin),36.663,ymax-0.02*(ymax-ymin))
    pave2.SetBorderSize(2)
    pave2.SetFillColor(0)
    pave2.SetTextSize(0.02)
    pave2.SetTextColor(42)
    pave2.AddText("Eu Campaign")
    pave2.Draw("SAME")

    #lines 37955,38292,Am Campaign
    line7=ROOT.TLine(37.955,ymin,37.955,ymax)
    line7.SetLineColor(43)
    line7.SetLineWidth(2)
    line7.Draw("SAME")
    line8=ROOT.TLine(38.292,ymin,38.292,ymax)
    line8.SetLineColor(43)
    line8.SetLineWidth(2)
    line8.Draw("SAME")
    pave3=ROOT.TPaveText(37.955,ymax-0.02*(ymax-ymin),38.292,ymax-0.02*(ymax-ymin))
    pave3.SetBorderSize(2)
    pave3.SetFillColor(0)
    pave3.SetTextSize(0.02)
    pave3.SetTextColor(43)
    pave3.AddText("Am Campaign")
    pave3.Draw("SAME")

    can1.Write()
    can1.SaveAs("output_plots/"+name+".png")
run_num=nparr(run_num)
canvas_history(run_num/1000,max_left,"Pedestal Run Number (#times 10^{3})","max from left FIT","Max noise left part FIT")
canvas_history(run_num/1000,max_right,"Pedestal Run Number (#times 10^{3})","max from left FIT","Max noise right part  FIT")
canvas_history(run_num/1000,slope_left,"Pedestal Run Number (#times 10^{3})","slope from left FIT","Exp slope left part FIT")
canvas_history(run_num/1000,slope_right,"Pedestal Run Number (#times 10^{3})","slope from left FIT","Exp slope right part FIT")
canvas_history(run_num/1000,mean_right,"Pedestal Run Number (#times 10^{3})","max from left 10 points","Max noise right part")
canvas_history(run_num/1000,mean_left,"Pedestal Run Number (#times 10^{3})","slope from left 10 points","Max noise left part")
canvas_history(run_num/1000,middle_mean,"Pedestal Run Number (#times 10^{3})","average in middle 300 points","noise in middle")
