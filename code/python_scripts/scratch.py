for sim in groupedDirs:
    for layer in groupedDirs[sim]:
        modelNames = []
        modelDirs = []
        for PNGfile in allModelFigs:
            if sim[4:] in PNGfile and layer.replace('-', '_') in PNGfile:
                modeldir, filename = os.path.split(PNGfile)
                modelNames.append(filename)
                modelDirs.append(modeldir)
            else:
                print('model not added')






##
parametersLateX, parametersLateXPath = create_parameter_table(dataPath)


##

parametersIm = convert_from_bytes(open(pathPDF, 'rb').read())[0]
pathPNG = os.path.join(dataPath, 'sim_shared_parameters.png')
# Crop image to its minimum size
parametersIm = parametersIm.crop((100, 500, 1000, 1400))
parametersIm.show()
parametersIm.save(pathPNG, "PNG", quality=95)



##
import pylatex as pl

doc = pl.Document()
doc.packages.append(pl.Package('booktabs'))

# Difference to the other answer:
with doc.create(pl.Table(position='htbp')) as table:
    table.add_caption('Simulation parameters')
    table.append(pl.Command('centering'))
    table.append(pl.NoEscape(parametersLateX))
#
doc.generate_pdf(os.path.join(dataPath, 'sim_shared_parameters'), clean_tex=False)


