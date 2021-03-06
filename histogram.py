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

    def add(self, other, scale = 1.):
        self.hWeighted.Add(other.hWeighted, scale)
        self.hScaleUp.Add(other.hScaleUp, scale)
        self.hScaleDown.Add(other.hScaleDown, scale)
        self.hRaw.Add(other.hRaw, scale)

    def scale(self, scale):
        self.hWeighted.Scale(scale)
        self.hScaleUp.Scale(scale)
        self.hScaleDown.Scale(scale)

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
    def __init__(self, name, binning, xtitle = '', xlabels = [], ytitle = '', ylabels = [], cond = [], overflowable = False, mask = None, logscale = False, xlog = False, vrange = None, drawOption = ''):
        self.name = name

        if len(binning) == 1 or type(binning) == list:
            if type(binning) == tuple:
                self.xedges = list(binning[0])
            else:
                self.xedges = list(binning)
            self.dimension = 1
            
        elif len(binning) == 3:
            nx, xmin, xmax = binning
            dx = (xmax - xmin) / nx
            self.xedges = [xmin + ix * dx for ix in range(nx + 1)]
            self.dimension = 1
            
        elif len(binning) == 2:
            self.xedges = list(binning[0])
            self.yedges = list(binning[1])
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
            self.xunit = xtitle[1]
            if self.xunit != 'NoUnit':
                self.xtitle = xtitle[0] + ' (' + xtitle[1] + ')'
            else:
                self.xtitle = xtitle[0]
        else:
            self.xunit = ''
            self.xtitle = xtitle

        if self.dimension == 1:
            self.ytitle = 'events'
            if self.xunit == 'NoUnit': pass
            elif self.xunit: self.ytitle += ' / ' + self.xunit
            else: self.ytitle += ' / unit'

            if self.overflowable: self.ytitle += ' (last bin: overflow events)'
            self.yunit = ''
        else:
            if type(ytitle) == tuple:
                self.ytitle = ytitle[0] + ' (' + ytitle[1] + ')'
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

        if type(cond) == list or type(cond) == tuple:
            self.conditions = list(cond)
        elif type(cond) == str:
            self.conditions = [cond]

        self.maskedRegion = mask

        self.logscale = logscale
        self.xlog = xlog

        self.vrange = vrange

        self.drawOption = drawOption

        self.xlabels = xlabels
        self.ylabels = ylabels

    def clone(self):
        if self.dimension == 1:
            clone = HDef(self.name, self.xedges)
            clone.dimension = 1
        else:
            clone = HDef(self.name, (self.xedges, self.yedges))
            clone.dimension = 2

        clone.conditions = list(self.conditions)
        clone.overflowable = self.overflowable
        clone.xtitle = self.xtitle
        clone.xunit = self.xunit
        clone.ytitle = self.ytitle
        clone.yunit = self.yunit
        clone.nx = self.nx
        clone.ny = self.ny
        clone.maskedRegion = self.maskedRegion
        clone.logscale = self.logscale
        clone.xlog = self.xlog
        clone.vrange = self.vrange
        clone.drawOption = self.drawOption
        clone.xlabels = self.xlabels
        clone.ylabels = self.ylabels

        return clone

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
            h = ROOT.TH1D(name, '', len(xedges) - 1, xedges)

        else:
            if self.overflowable:
                yedges = array.array('d', self.yedges + [self.yedges[-1] + (self.yedges[-1] - self.yedges[0]) / 40.])
            else:
                yedges = array.array('d', self.yedges)

            h = ROOT.TH2D(name, '', len(xedges) - 1, xedges, len(yedges) - 1, yedges)
        
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
            xunit = self.xunit
            if xunit == 'NoUnit': xunit = ''

            if xunit and self.yunit:
                if xunit == self.yunit: unit = xunit + '^2'
                else: unit = xunit + '#cdot' + self.yunit
            elif xunit:
                unit = xunit
            elif self.yunit:
                unit = self.yunit
            else:
                unit = 'unit area'
            title += ' / ' + unit
    
            h.GetZaxis().SetTitle(title)

        return h

