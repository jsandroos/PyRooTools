### define utility functions ###

def histmaker(name, title, nbins, xmin, xmax, xlab, ylab, color, style, width):
    histogram = TH1D(name, title, nbins, xmin, xmax)
    histogram.GetXaxis().SetTitle(xlab)
    histogram.GetYaxis().SetTitle(ylab)
    histogram.SetLineColor(color)
    histogram.SetLineStyle(style)
    histogram.SetLineWidth(width)
    histogram.Sumw2()
    histogram.SetStats(0)
    return histogram

def histmaker2D(name, title, xbins, xmin, xmax, ybins, ymin, ymax, xlab, ylab, zlab):
    hist = TH2D(name, title, xbins, xmin, xmax, ybins, ymin, ymax)
    hist.SetXTitle(xlab)
    hist.SetYTitle(ylab)
    hist.SetZTitle(zlab)
    hist.Sumw2()
    hist.SetStats(0)
    return hist

# divide two histograms
def divideHistograms(h1,h2,name="h_division",debug=False):
    h1clone=h1.Clone()
    h2clone=h2.Clone()
    h1clone.Divide(h1clone,h2clone)
    h1clone.SetName(name)
    h1clone.SetTitle(name)
    return h1clone

def get_bin_centers(edges):
    """Get the bin centers for a given set of bin edges.
       This works even if bins don't have equal width."""
    edges = np.array(edges, dtype=np.float)
    if is_logarithmic(edges):
        return np.sqrt(edges[:-1]*edges[1:])
    else:
        return (edges[:-1] + edges[1:])/2.

# write a 2D Histogram Filler
def D2HistFiller(input): #maked root 2d histogram out of 2d numpy array
    hist = TH2D("h1", "2d generic histogram", len(input), -1, 1, len(input[0]), true_bin_edges[0][0], true_bin_edges[0][len(true_bin_edges[0])-1])
    # fill by cos, energy
    for i in range(0, len(input)):
        for k in range(0, len(input[0])):
            val = input[i][k]
            
            #print bin, val
        #print " bin: ", i, val
        #print "step 5"
            hist.SetBinContent(i+1, k+1, val)
    return hist

def flatten_m_matrix(m_matrix):
    from itertools import product
    shape = m_matrix.shape
    dims = len(shape)
    assert dims % 2 == 0

    tr_dimensions = shape[:dims/2]
    re_dimensions = shape[dims/2:]



    flat_m_matrix = []
    r_indx = 0

    for bin_pdf in product(*map(range, reversed(tr_dimensions))):
        flat_m_matrix.append([])
        bin_pdf = tuple(reversed(bin_pdf))
        
        for bin_v in product(*map(range, reversed(re_dimensions))):
            bin_v = tuple(reversed(bin_v))
            flat_m_matrix[r_indx].append(m_matrix[bin_pdf][bin_v])

        r_indx += 1

    flat_m_matrix = np.array(flat_m_matrix)

    return flat_m_matrix

def SetDivisionError(h1, h2, rh): #the two histograms, and the rate histogram
    for i in range(0, h1.GetSize()):
        b1 = h1.GetBinContent(i)
        b2 = h2.GetBinContent(i)
        e1 = h1.GetBinError(i)
        e2 = h2.GetBinError(i)
        norm = (np.power(b1*b2,2))
        if norm < 1.0:
            rh.SetBinError(i, 0)
        else:
            err = np.sqrt((np.power(e1*b2,2)+np.power(e2*b1,2)) / norm)
            rh.SetBinError(i, err)
    return rh

def PlotPRratio(hist1, hist2, canv, label1, label2, labelratio, NonZero=0):
    c4 = ROOT.TCanvas("c4","c4 ExpUnfold results", 0, 0, 1900, 1000)
    #c4.Divide(0,2)
    pad1 = ROOT.TPad("pad1","",0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)# Upper and lower plot are joined
    pad1.SetGridx()#        // Vertical grid
    pad1.SetGridy()
    pad1.Draw()#          // Draw the upper pad: pad1
    pad1.cd()#               // pad1 becomes the current pad

    xtitle = hist1.GetXaxis().GetTitle()
    ytitle = hist1.GetYaxis().GetTitle()
    Lcolor = hist2.GetLineColor()
    Mcolor = hist2.GetMarkerColor()
    TLcolor = hist1.GetLineColor()
    TMcolor = hist1.GetMarkerColor()
    
    print "colors: ", Lcolor, Mcolor, TLcolor, TMcolor
    
    if NonZero==1:
        d1 = []
        d2 = []
        for i in range(0, hist1.GetSize()):
            bin1 = hist1.GetBinContent(i)
            bin2 = hist2.GetBinContent(i)
            if bin1 > 0.000001:
                d1.append([bin1, hist1.GetBinError(i)])
                d2.append([bin2, hist2.GetBinError(i)])
        print "lengths: ", len(d1), len(d2)
        h1 = TH1D("H1", hist1.GetTitle(), len(d1), 0, len(d1))
        h2 = TH1D("H2", hist2.GetTitle(), len(d1), 0, len(d1))
        for k in range(0, len(d1)):
            h1.Fill(k, d1[k][0])
            h2.Fill(k, d2[k][0])
            h1.SetBinError(k, d1[k][1])
            h2.SetBinError(k, d2[k][1])
    else:
        h1 = hist1.Clone()
        h2 = hist2.Clone()

    ROOT.gStyle.SetOptStat(0)
    h1.GetYaxis().CenterTitle()
    h1.GetYaxis().SetTitleSize(0.037)
    h1.GetYaxis().SetTitleOffset(0.87)
    h1.SetLineWidth(2)
    h1.SetLineColor(TLcolor)
    h1.SetMarkerColor(TMcolor)
    h1.GetYaxis().SetTitle(ytitle)
    h2.SetMarkerSize(1.0)
    h2.SetMarkerStyle(8)
    h2.SetFillColor(Lcolor)
    h2.SetFillStyle(3008)
    h2.SetLineColor(Lcolor)
    h2.SetMarkerColor(Mcolor)
    h1.Draw("hist E1")
    h2.Draw("same E1")
    pad1.Modified()
    pad1.Update()
    c4.Modified()
    c4.Update()


    c4.cd()#          // Go back to the main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2",0, 0.02, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.2075319)
    pad2.SetGridx()
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()
    pad2.Range(-125, -0.5211422, 1125.718,1.99)
    ratio = h2.Clone()
    ratio.Sumw2()
    ratio.Divide(h1)
    ratio.SetLineColor(38)
    ratio.SetLineWidth(3)
    ratio.GetYaxis().SetRangeUser(0.4, 1.59)
    ratio.SetXTitle(xtitle)
    ratio.GetXaxis().SetTitleSize(0.1)
    #ratio.GetXaxis().SetTitleOffset(1.2)
    ratio.GetXaxis().SetLabelSize(0.1)
    ratio.SetTitle("")
    ratio.GetYaxis().SetLabelSize(0.08)
    ratio.GetYaxis().SetTitleSize(0.11)
    ratio.GetYaxis().SetTitleOffset(0.3)
    ratio.SetYTitle("Ratio")
    ratio.GetYaxis().CenterTitle()
    ratio = SetDivisionError(h2, h1, ratio)
    ratio.SetFillColor(38)
    ratio.SetFillStyle(3008)
    ratio.SetMarkerColor(ROOT.kBlue-2)
    ratio.SetMarkerStyle(7)
    ratio.Draw("E0")

    pad2.Modified()
    pad2.Update()
    c4.cd()
    c4.Modified()
    c4.cd()
    c4.Update()
    
    c4.cd()
    pad1.cd()
    legUnf = ROOT.TLegend(0.7, 0.75, 0.95, 0.92)
    legUnf.AddEntry(h1, label1[0], label1[1])
    legUnf.AddEntry(h2, label2[0], label2[1])
    legUnf.AddEntry(ratio, labelratio[0], labelratio[1])
    legUnf.Draw()
    pad1.Modified()
    #pad1.Update()
    c4.Modified()
    c4.cd()
    c4.Update()
    c4.SetSelected(c4)
    
    sleep(30)
