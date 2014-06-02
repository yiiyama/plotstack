import ROOT

def drawPave(pave, onRight = True):
    if onRight:
        pave.SetX1NDC(0.72)
        pave.SetX2NDC(0.9)
    else:
        pave.SetX1NDC(0.1)
        pave.SetX2NDC(0.28)

    pave.Draw()
    pave.SetDrawOption(' ')


class Group(object):
    OBSERVED = 0
    BACKGROUND = 1
    SIGNAL = 2

    def __init__(self, name, title, color, category):
        if category not in [Group.OBSERVED, Group.BACKGROUND, Group.SIGNAL]:
            raise RuntimeError('Invalid category')
            
        self.name = name
        self.title = title
        self.color = color
        self.category = category


class Stack(object):

    USESIGNIFICANCE = False
    
    def __init__(self, hdef):
        self.name = hdef.name
        self.hdef = hdef
        self.groups = []
        self.histograms = {}

        self.obsHistogram = self.hdef.generate('obs')
        self.obsHistogram.SetLineColor(ROOT.kBlack)
        self.obsHistogram.SetMarkerStyle(8)
        self.obsHistogram.SetMarkerSize(0.4)
        self.obsHistogram.SetMarkerColor(ROOT.kBlack)
        self.obsHistogram.SetFillStyle(0)

        self.bkgHistogram = self.hdef.generate('bkg')
        self.bkgHistogram.SetLineWidth(0)
        self.bkgHistogram.SetMarkerSize(0)
        self.bkgHistogram.SetMarkerStyle(0)
        self.bkgHistogram.SetFillStyle(3003)
        self.bkgHistogram.SetFillColor(ROOT.kBlack)

    def addGroup(self, group):
        self.groups.append(group)
        histogram = self.hdef.generate(group.name)
        histogram.SetLineColor(group.color)
        if group.category == Group.OBSERVED:
            histogram.SetMarkerStyle(8)
            histogram.SetMarkerSize(0.4)
            histogram.SetMarkerColor(group.color)
            histogram.SetFillStyle(0)
        elif group.category == Group.BACKGROUND:
            histogram.SetFillColor(group.color)
            histogram.SetFillStyle(1001)
        if group.category == Group.SIGNAL:
            histogram.SetLineWidth(2)
            histogram.SetFillStyle(0)

        self.histograms[group] = histogram
    
    def add(self, group, h):
        self.histograms[group].Add(h)
        if group.category == Group.OBSERVED:
            self.obsHistogram.Add(h)
        elif group.category == Group.BACKGROUND:
            self.bkgHistogram.Add(h)

    def draw(self, plotsDir, texts = [], arbitraryUnit = False, maskObserved = False):
        if self.hdef.dimension == 1:
            self.draw1D(plotsDir, texts, arbitraryUnit, maskObserved)
        else:
            self.draw2D(plotsDir, texts, arbitraryUnit, maskObserved)

    def draw1D(self, plotsDir, texts, arbitraryUnit, maskObserved):
        bkgGroups = filter(lambda x: x.category == Group.BACKGROUND, self.groups)
        sigGroups = filter(lambda x: x.category == Group.SIGNAL, self.groups)

        stack = ROOT.THStack(self.hdef.name, self.hdef.title)

        for group in bkgGroups:
            stack.Add(self.histograms[group])

        canvas = ROOT.TCanvas('stack', 'stack')
        canvas.Divide(1, 2, 0., 0.)

        distPad = canvas.GetPad(1)
        ratioPad = canvas.GetPad(2)

        distPad.SetPad(0., 0.25, 1., 1.)
        distPad.SetLeftMargin(0.07)
        distPad.SetRightMargin(0.07)
        distPad.SetTopMargin(0.1)
        distPad.SetBottomMargin(0.02)
        
        ratioPad.SetPad(0., 0., 1., 0.25)
        ratioPad.SetGridy(2)
        ratioPad.SetLeftMargin(0.07)
        ratioPad.SetRightMargin(0.07)
        ratioPad.SetTopMargin(0.02)
        ratioPad.SetBottomMargin(0.25)

        distPad.SetLogy(self.hdef.logscale)

        legend = ROOT.TLegend(0., 0., 0., 0.9)
        legend.SetFillStyle(4000)
        legend.SetBorderSize(0)
        legend.SetTextFont(62)
        legend.SetTextSize(0.03)
        legend.SetTextAlign(12)
        legend.ConvertNDCtoPad()

        if self.obsHistogram.GetSumOfWeights() > 0.:
            legend.AddEntry(self.obsHistogram, 'Observed', 'LP')

        for group in reversed(bkgGroups):
            legend.AddEntry(self.histograms[group], group.title, 'F')

        for group in sigGroups:
            legend.AddEntry(self.histograms[group], group.title, 'L')

        y = 0.9 - legend.GetNRows() * 0.04
        legend.SetY1NDC(y)
        for text in texts:
            height = text.GetY2NDC() - text.GetY1NDC()
            text.SetY2NDC(y)
            y -= height
            text.SetY1NDC(y)

        if self.hdef.logscale:
            minContent = self.bkgHistogram.GetMaximum()
            for iX in range(1, self.bkgHistogram.GetNbinsX() + 1):
                cont = self.bkgHistogram.GetBinContent(iX)
                if cont != 0. and cont < minContent:
                    minContent = cont
            stack.SetMinimum(minContent * 0.2)
        else:
            stack.SetMinimum(0.)
        stack.SetMaximum(self.bkgHistogram.GetMaximum())

        distPad.cd()

        stack.Draw('HIST')
        self.bkgHistogram.Draw('E2 SAME')
        for group in sigGroups:
            self.histograms[group].Draw('HIST SAME')

        self.obsHistogram.Draw('EP SAME')

        distPad.Update()

        stack.GetXaxis().SetTitle(self.hdef.xtitle)
        stack.GetYaxis().SetTitle(self.hdef.ytitle)

        stack.GetXaxis().SetLabelSize(0.)
        stack.GetXaxis().SetTitleSize(0.)
        stack.GetYaxis().SetTitleOffset(0.8)
        if arbitraryUnit: stack.GetYaxis().SetLabelSize(0)

        if stack.GetMaximum() < self.obsHistogram.GetMaximum():
            if distPad.GetLogy():
                stack.SetMaximum(self.obsHistogram.GetMaximum() * 5.)
            else:
                stack.SetMaximum(self.obsHistogram.GetMaximum() * 1.2)

        if maskObserved and self.hdef.maskedRegion:
            if self.hdef.maskedRegion[0] == '-inf':
                x1 = self.obsHistogram.GetXaxis().GetXmin()
            else:
                x1 = self.hdef.maskedRegion[0]
            if self.hdef.maskedRegion[1] == 'inf':
                x2 = self.obsHistogram.GetXaxis().GetXmax()
            else:
                x2 = self.hdef.maskedRegion[1]

            box = ROOT.TPave(x1, 0., x2, 0., 0, '')
            box.ConvertNDCtoPad()
            box.SetY1NDC(0.02)
            box.SetY2NDC(0.9)
            box.SetFillColor(ROOT.kGray)
            box.SetFillStyle(3944)
            box.Draw()

        onRight = self.bkgHistogram.GetMean() < self.hdef.xedges[-1] / 2.
        drawPave(legend, onRight)
        for text in texts:
            drawPave(text, onRight)

        ratioPad.cd()

        if Stack.USESIGNIFICANCE:
            hDiff = self.obsHistogram.Clone(self.name + '_diff')

            for iX in range(1, self.hdef.nx + 1):
                bkg = self.bkgHistogram.GetBinContent(iX)
                bkgErr = self.bkgHistogram.GetBinError(iX)
                if bkgErr == 0.:
                    bkgErr = 1. / self.bkgHistogram.GetXaxis().GetBinWidth(iX)

                hDiff.SetBinError(iX, self.obsHistogram.GetBinError(iX) / bkgErr)
                hDiff.SetBinContent(iX, (self.obsHistogram.GetBinContent(iX) - bkg) / bkgErr)

            hDiff.Draw('EP')
            hDiff.GetYaxis().SetRangeUser(-3., 3.)
            hDiff.GetYaxis().SetTitle('(obs - bkg) / #delta_{bkg}')
            hFrame = hDiff

            line = ROOT.TLine(hDiff.GetXaxis().GetXmin(), 0., hDiff.GetXaxis().GetXmax(), 0.)
        else:
            hUncert = self.bkgHistogram.Clone(self.name + '_uncert')
            for iX in range(1, self.hdef.nx + 1):
                bkg = self.bkgHistogram.GetBinContent(iX)
                if bkg == 0.: continue
                hUncert.SetBinError(iX, self.bkgHistogram.GetBinError(iX) / bkg)
                hUncert.SetBinContent(iX, 1.)

            hUncert.Draw('E2')
            hUncert.GetYaxis().SetRangeUser(0.5, 1.5)
            hUncert.GetYaxis().SetTitle('obs / bkg')
            hFrame = hUncert
            
            hRatio = self.obsHistogram.Clone(self.name + '_ratio')

            for iX in range(1, self.hdef.nx + 1):
                bkg = self.bkgHistogram.GetBinContent(iX)
                bkgErr = self.bkgHistogram.GetBinError(iX)
                if bkg == 0.: continue

                hRatio.SetBinError(iX, self.obsHistogram.GetBinError(iX) / bkg)
                hRatio.SetBinContent(iX, self.obsHistogram.GetBinContent(iX) / bkg)

            hRatio.Draw('EP SAME')

            line = ROOT.TLine(hUncert.GetXaxis().GetXmin(), 1., hUncert.GetXaxis().GetXmax(), 1.)

        line.Draw()

        hFrame.SetTitle("")
        hFrame.GetXaxis().SetLabelSize(0.1)
        hFrame.GetXaxis().SetTitleSize(0.1)
        hFrame.GetYaxis().SetTitleOffset(0.2)
        hFrame.GetYaxis().SetLabelSize(0.1)
        hFrame.GetYaxis().SetTitleSize(0.1)
        hFrame.GetYaxis().SetNdivisions(110)

        canvas.SaveAs(plotsDir + "/" + self.name + ".pdf")

    def draw2D(self, plotsDir, texts, arbitraryUnit, maskObserved):
        drawOption = self.hdef.drawOption.upper()
        if not drawOption:
            drawOption = 'COLZ'

        canvas = ROOT.TCanvas('stack', 'stack')
        canvas.SetLeftMargin(0.07)
        canvas.SetRightMargin(0.07)
        canvas.SetTopMargin(0.1)
        canvas.SetBottomMargin(0.1)

        canvas.SetLogz(self.hdef.logscale)

        y = 0.9
        for text in texts:
            height = text.GetY2NDC() - text.GetY1NDC()
            text.SetY2NDC(y)
            y -= height
            text.SetY1NDC(y)

        sigGroups = filter(lambda x: x.category == Group.SIGNAL, self.groups)

        if 'BOX' in drawOption and 'TEXT' in drawOption:
            drawOption = drawOption.replace('TEXT', '')

            self.bkgHistogram.Draw(drawOption)
            for text in texts:
                drawPave(text, True)

            legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.3)
            legend.SetFillStyle(4000)
            legend.SetBorderSize(0)
            legend.SetTextFont(62)
            legend.SetTextSize(0.03)
            legend.SetTextAlign(12)
            legend.ConvertNDCtoPad()            
    
            self.bkgHistogram.SetFillStyle(0)
            self.bkgHistogram.SetLineColor(ROOT.kBlue)
            legend.AddEntry(self.bkgHistogram, 'Estimated', 'L')

            self.obsHistogram.Draw(drawOption + ' SAME')
            legend.AddEntry(self.obsHistogram, 'Observed', 'L')
    
            legend.Draw()

            text = ROOT.TLatex()
            text.SetTextSize(0.02 * self.bkgHistogram.GetMarkerSize())
            text.SetTextAlign(22)
            for iX in range(1, self.hdef.nx + 1):
                x = self.bkgHistogram.GetXaxis().GetBinCenter(iX)
                for iY in range(1, len(self.hdef.yedges)):
                    bin = self.bkgHistogram.GetBin(iX, iY)
                    y = self.bkgHistogram.GetYaxis().GetBinLowEdge(iY)
                    offset = self.bkgHistogram.GetYaxis().GetBinWidth(iY) / 3.

                    text.DrawLatex(x, y + offset, '%.1f#pm%.1f' % (self.bkgHistogram.GetBinContent(bin), self.bkgHistogram.GetBinError(bin))).SetTextColor(ROOT.kBlue)
                    text.DrawLatex(x, y + offset * 2., '%.1f#pm%.1f' % (self.obsHistogram.GetBinContent(bin), self.obsHistogram.GetBinError(bin))).SetTextColor(ROOT.kBlack)
                    
            canvas.SaveAs(plotsDir + "/" + self.name + ".pdf")
    
        else:
            self.bkgHistogram.Draw(drawOption)

            if arbitraryUnit: self.bkgHistogram.GetZaxis().SetLabelSize(0)
            
            for text in texts:
                drawPave(text, True)

            canvas.SaveAs(plotsDir + '/' + self.bkgHistogram.GetName() + '.pdf')

            self.obsHistogram.Draw(drawOption)

            if arbitraryUnit: self.obsHistogram.GetZaxis().SetLabelSize(0)
            
            for text in texts:
                drawPave(text, True)

            canvas.SaveAs(plotsDir + '/' + self.obsHistogram.GetName() + '.pdf')

        for group in sigGroups:
            histo = self.histograms[group]
            if histo.GetEntries() == 0.: continue
            histo.Draw(drawOption)

            if arbitraryUnit: histo.GetZaxis().SetLabelSize(0)
            
            for text in texts:
                drawPave(text, True)

            canvas.SaveAs(plotsDir + "/" + histo.GetName() + ".pdf")


class StackConfig(object):
    def __init__(self, plotMaker):
        self.plotMaker = plotMaker
        self.sampleSets = {}
        self.floatFixer = None
        self.specialPlotMakers = {}
