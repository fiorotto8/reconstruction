import pandas as pd
import ROOT
import numpy as np

# Function to create a TGraph object
def graph(x, y, x_string, y_string, name=None, color=46, markerstyle=33, markersize=2, write=True):
    # Create a TGraph with x and y data
    plot = ROOT.TGraph(len(x), np.array(x, dtype="d"), np.array(y, dtype="d"))

    # Set name and title of the TGraph
    if name is None:
        plot.SetNameTitle(y_string + " vs " + x_string, y_string + " vs " + x_string)
    else:
        plot.SetNameTitle(name, name)

    # Set X and Y axis labels
    plot.GetXaxis().SetTitle(x_string)
    plot.GetYaxis().SetTitle(y_string)

    # Set marker color, style, and size
    plot.SetMarkerColor(color)  # blue
    plot.SetMarkerStyle(markerstyle)
    plot.SetMarkerSize(markersize)

    # Write the TGraph to file if specified
    if write:
        plot.Write()

    return plot


# Read CSV files into pandas DataFrames
hum_df = pd.read_csv("Humidity-data-2023-12-01 11_01_56.csv")
run_df = pd.read_csv("run number-data-2023-12-01 11_02_32.csv")

# Extract the second columns and apply a filter
hum_vals = hum_df.iloc[:, 1]
run_vals = run_df.iloc[:, 1]
filter_condition = hum_vals > -10
hum_vals = hum_vals[filter_condition]
run_vals = run_vals[filter_condition]

# Create a TGraph for normalized humidity
hum_plot = graph(run_vals / 1000, hum_vals / np.max(hum_vals), "Run Num", "Normalized Humidity", write=False)

# Read a TFile containing a TGraphErrors for normalized mean
root_file = ROOT.TFile("anal_ped.root", "READ")
mean_plot = root_file.Get("Normalized Mean from TH2")
# Retrieve the TGraphErrors from the TCanvas
mean_plot = mean_plot.GetPrimitive("Normalized Mean from TH2")  # Replace "GraphErrors_Name" with the actual name of your TGraphErrors
mean_plot_clone=mean_plot.Clone()
root_file.Close()
# Create a TFile for storing output
main = ROOT.TFile("hum_ped.root", "RECREATE")



# Create a canvas
canvas = ROOT.TCanvas("canvas", "Mean from TH2 and Humidity", 1000, 1000)
canvas.SetLeftMargin(0.15)

# Create the first pad (top)
pad1 = ROOT.TPad("pad1", "pad1", 0, 0.5, 1, 1)
pad1.Draw()
pad1.cd()
hum_plot.Draw("AP")
# Customize pad1 or graph1 as needed

# Create the second pad (bottom)
canvas.cd()
pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.5)
pad2.Draw()
pad2.cd()
mean_plot.Draw("AP")
# Customize pad2 or graph2 as needed










"""
# Draw the TGraph for normalized humidity on the left Y-axis
hum_plot.Draw("AP")

# Draw the TGraphErrors for normalized mean on the right Y-axis
mean_plot.Draw("P SAME")
"""

# Set legend if needed
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(hum_plot, "Normalized Humidity", "p")
legend.AddEntry(mean_plot, "Normalized Mean from TH2", "p")
legend.Draw()

# Show the canvas
canvas.Update()
canvas.Write()
canvas.SaveAs("test.png")
#canvas.WaitPrimitive()


###########################################################################################


#io ho hum_vals and run_vals
# Get the number of points in the TGraph
num_points = mean_plot_clone.GetN()
# Initialize empty lists for x and y values
x_mean,y_mean=[],[]
# Loop over each point in the TGraph and append x, y values to the lists
for i in range(num_points):
    x_mean.append(mean_plot.GetX()[i])
    y_mean.append(mean_plot.GetY()[i])

means,hums=[],[]
for j in range(len(run_vals)):
    if j!=826 and j!=1677:
        if (run_vals[j]/1000) in x_mean:
            #print(run_vals[j])
            if run_vals[j]<27880 or run_vals[j]>34790:
                if run_vals[j]<37870 or run_vals[j]>37930:
                    #print(run_vals[j], x_mean[(x_mean.index(run_vals[j]/1000))])
                    means.append(y_mean[(x_mean.index(run_vals[j]/1000))])
                    hums.append(hum_vals[j])
print(len(hums), len(means))
graph(hums, means, "Humidity(%)", "Mean from TH2",color=9)


