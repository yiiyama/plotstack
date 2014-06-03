import math
import array
import ROOT

class Histogram:
    def __init__(self, hdef, directory = None, suffix = ''):
        self.hdef = hdef

        if suffix:
            self.name = self.hdef.name + '_' + suffix
        else:
            self.name = self.hdef.name

        if not directory:
            self.hWeighted = hdef.generate(suffix)
            self.hScaleUp = hdef.generate(suffix + '_ScaleUp')
            self.hScaleDown = hdef.generate(suffix + '_ScaleDown')
            self.hRaw = hdef.generate(suffix + '_Raw')
        else:
            self.hWeighted = directory.Get(self.name)
            self.hScaleUp = directory.Get(self.name + '_ScaleUp')
            self.hScaleDown = directory.Get(self.name + '_ScaleDown')
            self.hRaw = directory.Get(self.name + '_Raw')

        self._core = ROOT.Histogram(
            self.hWeighted,
            self.hScaleUp,
            self.hScaleDown,
            self.hRaw
        )
        self._core.overflowable = self.hdef.overflowable

    def add(self, other, factor = 1.):
        self.hWeighted.Add(other.hWeighted, factor)
        self.hScaleUp.Add(other.hScaleUp, factor)
        self.hScaleDown.Add(other.hScaleDown, factor)
        self.hRaw.Add(other.hRaw, factor)

    def scale(self, scale):
        self.hWeighted.Scale(scale)
        self.hScaleUp.Scale(scale)
        self.hScaleDown.Scale(scale)

    def setError(self):
        for iX in range(1, self.hdef.nx + 1):
            for iY in range(1, self.hdef.ny + 1):
                bin = self.hWeighted.GetBin(iX, iY)
                cont = self.hWeighted.GetBinContent(bin)
                stat = self.hWeighted.GetBinError(bin)
                scal = max(self.hScaleUp.GetBinContent(bin) - cont, cont - self.hScaleDown.GetBinContent(bin))
                self.hWeighted.SetBinError(bin, math.sqrt(stat * stat + scal * scal))

    def postFill(self, applyMask = False):
        if self.hdef.maskedRegion and applyMask:
            maskedBins = []
            if self.hdef.maskedRegion[0] == '-inf':
                maskLowX = 0
            else:
                maskLowX = self.hRaw.GetXaxis().FindFixBin(self.hdef.maskedRegion[0])
            if self.hdef.maskedRegion[1] == 'inf':
                maskHighX = self.hdef.nx + 1
            else:
                maskHighX = self.hRaw.GetXaxis().FindFixBin(self.hdef.maskedRegion[1])

            if self.hdef.dimension == 1:
                maskedBins = [bin for bin in range(maskLowX, maskHighX + 1)]
            else:
                if self.hdef.maskedRegion[2] == '-inf':
                    maskLowY = 0
                else:
                    maskLowY = self.hRaw.GetYaxis().FindFixBin(self.hdef.maskedRegion[2])
                if self.hdef.maskedRegion[3] == 'inf':
                    maskHighY = self.hdef.ny + 1
                else:
                    maskHighY = self.hRaw.GetYaxis().FindFixBin(self.hdef.maskedRegion[2])

                for iY in range(maskLowY, maskHighY + 1):
                    for iX in range(maskLowX, maskHighX + 1):
                        maskedBins.append(self.hRaw.GetBin(iX, iY))

            for bin in set(maskedBins):
                entries = self.hRaw.GetBinContent(bin)
                for h in [self.hWeighted, self.hScaleUp, self.hScaleDown, self.hRaw]:
                    h.SetBinContent(bin, 0.)
                    h.SetBinError(bin, 0.)
                    h.SetEntries(h.GetEntries() - entries)


class HDef(object):
    def __init__(self, name, title, binning, xtitle = '', ytitle = '', overflowable = False, mask = None, logscale = False, drawOption = ''):
        self.name = name
        self.title = title
        
        if len(binning) == 1 or type(binning) == list:
            if type(binning) == tuple:
                self.xedges = list(binning[0])
            else:
                self.xedges = binning
            self.dimension = 1
            
        elif len(binning) == 3:
            nx, xmin, xmax = binning
            dx = (xmax - xmin) / nx
            self.xedges = [xmin + ix * dx for ix in range(nx + 1)]
            self.dimension = 1
            
        elif len(binning) == 2:
            self.xedges = binning[0]
            self.yedges = binning[1]
            self.dimension = 2
            
        elif len(binning) == 6:
            nx, xmin, xmax, ny, ymin, ymax = binning
            dx = (xmax - xmin) / nx
            dy = (ymax - ymin) / ny
            self.xedges = [xmin + ix * dx for ix in range(nx + 1)]
            self.yedges = [ymin + iy * dy for iy in range(ny + 1)]
            self.dimension = 2

        self.overflowable = overflowable

        if type(xtitle) == tuple:
            self.xtitle = xtitle[0] + '(' + xtitle[1] + ')'
            self.xunit = xtitle[1]
        else:
            self.xtitle = xtitle
            self.xunit = ''

        if self.dimension == 1:
            self.ytitle = 'events'
            if self.xunit: self.ytitle += ' / ' + self.xunit
            else: self.ytitle += ' / unit'
            if self.overflowable: self.ytitle += ' (last bin: overflow events)'
        else:
            if type(ytitle) == tuple:
                self.ytitle = ytitle[0] + '(' + ytitle[1] + ')'
                self.yunit = ytitle[1]
            else:
                self.ytitle = ytitle
                self.yunit = ''

        if self.overflowable:
            self.nx = len(self.xedges)
            if self.dimension == 2:
                self.ny = len(self.yedges)
            else:
                self.ny = 1
        else:
            self.nx = len(self.xedges) - 1
            if self.dimension == 2:
                self.ny = len(self.yedges) - 1
            else:
                self.ny = 1

        self.maskedRegion = mask

        self.logscale = logscale

        self.drawOption = drawOption

        self.xlabels = []
        self.ylabels = []

    def generate(self, suffix = ''):
        if suffix:
            name = self.name + '_' + suffix
        else:
            name = self.name

        if self.overflowable:
            xedges = array.array('d', self.xedges + [self.xedges[-1] + (self.xedges[-1] - self.xedges[0]) / 40.])
        else:
            xedges = array.array('d', self.xedges)

        if self.dimension == 1:
            h = ROOT.TH1F(name, self.title, len(xedges) - 1, xedges)

        else:
            if self.overflowable:
                yedges = array.array('d', self.yedges + [self.yedges[-1] + (self.yedges[-1] - self.yedges[0]) / 40.])
            else:
                yedges = array.array('d', self.yedges)

            h = ROOT.TH2F(name, self.title, len(xedges) - 1, xedges, len(yedges) - 1, yedges)
        
        h.GetXaxis().SetTitle(self.xtitle)
        h.GetYaxis().SetTitle(self.ytitle)

        if len(self.xlabels) == len(self.xedges) - 1:
            iX = 1
            for label in self.xlabels:
                h.GetXaxis().SetBinLabel(iX, label)
                iX += 1

        if self.dimension == 2:
            if len(self.ylabels) == len(self.yedges) - 1:
                iY = 1
                for label in self.ylabels:
                    h.GetYaxis().SetBinLabel(iY, label)
                    iY += 1
            
            title = 'events'
            if self.xunit and self.yunit:
                if self.xunit == self.yunit: unit = self.xunit + '^2'
                else: unit = self.xunit + '#cdot' + self.yunit
            elif self.xunit:
                unit = self.xunit
            elif self.yunit:
                unit = self.yunit
            else:
                unit = 'unit area'
            title += ' / ' + unit
    
            h.GetZaxis().SetTitle(title)

        return h


class HistogramContainer(object):
    def __init__(self, name):
        self.name = name
        self.histograms = []

    def bookHistograms(self, hdefs, outputFile):
        self.histograms = []
        for hdef in hdefs:
            errIgn = ROOT.gErrorIgnoreLevel
            ROOT.gErrorIgnoreLevel = 4000
            if not outputFile.cd(hdef.name):
                outputFile.mkdir(hdef.name).cd()
            ROOT.gErrorIgnoreLevel = errIgn
                
            self.histograms.append(Histogram(hdef, suffix = self.name))

    def loadHistograms(self, hdefs, sourceFile):
        self.histograms = []
        for hdef in hdefs:
            directory = sourceFile.GetDirectory(hdef.name)
            self.histograms.append(Histogram(hdef, suffix = self.name, directory = directory))

    def getHistogram(self, name):
        return next(h for h in self.histograms if h.name == name + '_' + self.name)

    def setHistogramErrors(self):
        for h in self.histograms:
            h.setError()

    def scaleHistograms(self, scale):
        for h in self.histograms:
            h.scale(scale)

    def postFill(self, applyMask):
        for h in self.histograms:
            h.postFill(applyMask)
