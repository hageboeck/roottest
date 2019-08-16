# RooFit example of D0 mass fit with toy data.

from ROOT import RooWorkspace, RooArgSet, RooDataSet, RooFit, RooAbsData, TCanvas, TColor, TTree, TPad, gStyle, TLine, gSystem
import ROOT
import time
import sys

runPartial = False

def main():
    ws = RooWorkspace()
    
    import_model(ws)
    model = ws.pdf('model')
    data = generate_data(ws, model)
    do_mass_fit(ws, data)

def import_model(workspace):
    print('Make D0 model.')

    # Fit parameters
    workspace.factory('prod::sigma(frac[0.5,0,1],s[10,8,12])')
    workspace.factory('f[1.7, 0, 3]')
    workspace.factory('prod::sigma1(sigma,f)')
    workspace.factory('a1[0, -0.2, 0.2]')
    workspace.factory('a2[0, -0.2, 0.2]')
    workspace.factory('alpha[3]')

    workspace.factory('Gaussian::sig2(D0_LoKi_DTF_M[1805,1925],mean[1865, 1850, 1870],sigma1)')
    workspace.factory('CBShape::sig1(D0_LoKi_DTF_M,mean, sigma, alpha, n[3])')

    workspace.factory('Chebychev::bkg(D0_LoKi_DTF_M,{a1,a2})')

    workspace.factory('SUM::sig(frac*sig1,sig2)')
    if runPartial:
        workspace.factory('SUM::model(nsig[1e5, 0.5e5, 1.5e5]*sig,nbkg[0.5e5, 0.2e5, 1.5e5]*bkg)')
    else:
        workspace.factory('SUM::model(nsig[1e6, 0.5e6, 1.5e6]*sig,nbkg[0.5e6, 0.2e6, 1.5e6]*bkg)')
    # ~ workspace.Print()
    print('Making model finished.')

def generate_data(workspace, pdf):
    # Get what we need from workspace
    D0_LoKi_DTF_M = workspace.var('D0_LoKi_DTF_M')
    myset = RooArgSet(D0_LoKi_DTF_M)

    # Generate toy events
    dataset = pdf.generate(myset, 1500000 if not runPartial else 150000)
    return dataset
    
def kick_parameters(workspace, parameters):
    for param in parameters:
        param.setVal(param.getVal() * 0.9 + 0.1)
    workspace.saveSnapshot("params_kicked", parameters)

def do_mass_fit(workspace, dataset):
    print('Do fit.')

    # Get what we need from workspace
    D0_LoKi_DTF_M = workspace.var('D0_LoKi_DTF_M')

    model = workspace.pdf('model')
    bkg = workspace.pdf('bkg')
    sig = workspace.pdf('sig')
    sig1 = workspace.pdf('sig1')
    sig2 = workspace.pdf('sig2')
    
    parameters = model.getParameters(dataset)
    workspace.saveSnapshot('params_initial', parameters)
    kick_parameters(workspace, parameters)

    # Run the scalar reference:
    workspace.loadSnapshot('params_kicked')
    print('\n\nStart scalar fit:')
    beginScalarFit = time.clock()
    fitResultScalar = model.fitTo(dataset,
        RooFit.Extended(),
        RooFit.Save(),
        RooFit.PrintLevel(-1))
    endScalarFit = time.clock()
    
    fitResultScalar.Print()
    if fitResultScalar.status() != 0:
        print('Scalar fit failed')
        sys.exit(2)
    
    
    # Reset parameters and run batch
    workspace.loadSnapshot('params_kicked')
    print('\n\nStart batch fit:')
    beginBatchFit = time.clock()
    fitResultBatch = model.fitTo(dataset,
        RooFit.Extended(),
        RooFit.Save(),
        RooFit.PrintLevel(-1),
        RooFit.BatchMode(True),
        RooFit.Optimize(1))
    endBatchFit = time.clock()
    
    fitResultBatch.Print()
    if fitResultBatch.status() != 0:
        print('Batch fit failed')
        sys.exit(2)

    
    print('Timer batch  fit:', endBatchFit-beginBatchFit, 's')
    print('Timer scalar fit:', endScalarFit-beginScalarFit, 's')
    
    # This will print to stdout if something is wrong. Pick it up using
    # FAILREGEX
    fitResultBatch.isIdentical(fitResultScalar, 2.E-4)


     # Some plottin and chi2 tests
    frame = D0_LoKi_DTF_M.frame()
    dataset.plotOn(frame, RooFit.DataError(RooAbsData.SumW2))
    model.plotOn(frame, RooFit.Name("model"))

    npars = len(fitResultBatch.floatParsFinal())
    chi2_ndof = frame.chiSquare(npars)
    ndof = D0_LoKi_DTF_M.numBins() - npars
    chi2 = chi2_ndof * ndof
    print('Fit has chi2/ndof = ', chi2, '/', ndof, '=', chi2_ndof, '\tP=',
        ROOT.TMath.Prob(chi2, ndof))
    if not runPartial and not (1.1 <= chi2_ndof and chi2_ndof <= 1.13):
        print('Chi2/ndof is off')
        sys.exit(3)
    
    # Calculate pulls
    frame_p = D0_LoKi_DTF_M.frame(RooFit.Title(''))
    frame_p.addPlotable(frame.pullHist(), 'P')

    model.plotOn(frame, RooFit.Components('bkg'), RooFit.LineColor(419), RooFit.LineStyle(5),
        RooFit.Name('bkg'))
    model.plotOn(frame, RooFit.Components('sig1'), RooFit.LineColor(874), RooFit.LineStyle(2),
        RooFit.Name('sig1'))
    model.plotOn(frame, RooFit.Components('sig2'), RooFit.LineColor(886), RooFit.LineStyle(4),
        RooFit.Name('sig2'))
    model.paramOn(frame, RooFit.Format('NEU', RooFit.AutoPrecision(1)), RooFit.Layout(0.7,0.9,0.9))
    frame.getAttText().SetTextSize(0.03)
    frame.getAttText().SetTextFont(132)
    frame.GetXaxis().SetLabelOffset(999) # Remove X axis labels
    frame.GetYaxis().SetLabelFont(132)
    frame.GetYaxis().SetTitleFont(132)
    
    legend = ROOT.TLegend(0.15, 0.9, 0.4, 0.6)
    legend.AddEntry(frame.findObject("model"), 'Full model', 'l')
    legend.AddEntry(frame.findObject("bkg"), 'bkg', 'l')
    legend.AddEntry(frame.findObject("sig1"), 'sig1', 'l')
    legend.AddEntry(frame.findObject("sig2"), 'sig2', 'l')
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Pull frame
    frame_p.SetMinimum(-5)
    frame_p.SetMaximum(5)
    frame_p.GetXaxis().SetTitle('#font[12]{m}(#font[12]{D^{0}} ) [MeV]')
    frame_p.GetXaxis().SetTitleOffset(1.)
    frame_p.GetXaxis().SetLabelFont(132)
    frame_p.GetXaxis().SetTitleFont(132)
    frame_p.SetLabelSize(0.06, 'XY')
    frame_p.GetYaxis().SetTitle('(data - model)/#font[12]{#sigma}_{data}')
    frame_p.GetYaxis().SetTitleOffset(0.5)
    frame_p.GetYaxis().SetNdivisions(209)
    frame_p.GetYaxis().SetLabelFont(132)
    frame_p.GetYaxis().SetTitleFont(132)
    frame_p.SetTitleSize(0.06, 'XY')

    # Create Canvas and Pads and Draw plot
    c1 = TCanvas('', '', 600, 600)
    upperPad = TPad('upperPad', 'upperPad', .005, .3475, .995, .995)
    lowerPad = TPad('lowerPad', 'lowerPad', .005, .005, .995, .3475)
    upperPad.SetBottomMargin(0.03)
    lowerPad.SetTopMargin(0.02)
    lowerPad.SetBottomMargin(0.3)

    upperPad.Draw()
    lowerPad.Draw()
    upperPad.cd()
    frame.Draw()
    legend.Draw()
    legend.Print()
    lowerPad.cd()
    gStyle.SetOptTitle(0)
    line = TLine(1805,0,1925,0)
    line.SetLineColor(1)
    frame_p.Draw()
    line.Draw()
    c1.SaveAs('RooFit_LHCb001_example.pdf')


if __name__ == '__main__':
    main()
