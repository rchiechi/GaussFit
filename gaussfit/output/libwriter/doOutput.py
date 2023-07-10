def doOutput(writer):
    ''' Run through all the methods for writing output files.'''
    writer.logger.info("Writing files...")
    writer.WriteParseInfo()
    writer.WriteSummary()
    # TODO Vtrans looses its mind if you only parse positive or negative biases
    try:
        writer.WriteVtrans()
    except ValueError:
        print("Vtrans data are too broken to output ¯\_(ツ)_/¯")  # noqa #pylint: disable=W1401
    writer.WriteGNUplot('Vtransplot')
    writer.WriteFN()
    writer.WriteVT()
    writer.WriteSLM()
    writer.WriteGNUplot('VTplot')
    writer.WriteGauss()
    writer.WriteSegmentedGauss()
    writer.WriteSegmentedGauss('nofirst')
    writer.WriteFilteredGauss()
    writer.WriteGNUplot('JVplot')
    writer.WriteGNUplot('NDCplot')
    writer.WriteData()
    writer.WriteDJDV()
    writer.WriteNDC()
    writer.WriteFiltered()
    writer.WriteData(True)
    writer.WriteLag()
    writer.WriteRData()
    writer.WriteGNUplot('Rplot')
    try:
        writer.WriteHistograms()
        writer.WriteFilteredHistograms()
        writer.WriteGNUplot('JVhistplot')
    except IndexError as msg:
        print("Error outputting histrograms %s" % str(msg))
    try:
        writer.WriteGHistogram()
    except IndexError:
        print("Error outputting Ghistrograms")
    try:
        writer.WriteGMatrix('GMatrix')
        writer.WriteGNUplot('GMatrix', ['parula.pal'])
        writer.WriteGMatrix('NDCMatrix')
        writer.WriteGNUplot('NDCMatrix', ['parula.pal'])
        writer.WriteGMatrix('LogJMatrix')
        writer.WriteGNUplot('LogJMatrix', ['parula.pal'])
    except IndexError:
        print("Error outputting GMatrix")
    writer.logger.info("Done!")
