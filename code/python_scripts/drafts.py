##
datapath = '/Volumes/Transcend/derivatives_test_sim_control_noiselevels'
simNames, layerNames, groupedDirs = get_sim_and_layer_names(datapath)

# Load parameters

parametersPNG_dir = os.path.join(datapath, 'sim_shared_parameters.png')
parametersImage = Image.open(parametersPNG_dir)

# Order figures
namesToSort = ['Fig-nanmean via Fisher z (CS separate)',
               'Fig-CI95 via Fisher z (CS separate)',
               'Fig-nanmean via Fisher z (CS averaged)',
               'Fig-CI95 via Fisher z (CS averaged)',
               'Fig-Mean vs CI95 (CS separate)',
               'Fig-Mean vs CI95 (CS averaged)',
               'Fig-Time-resolved correlation (C-S averaged) - all time series',
               'Fig-Time-resolved correlation (C-S averaged) - averaged time series_',
               'Fig-Time-resolved correlation (C-S averaged) - averaged time series -',
               'Fig-Time-resolved correlation (C-S averaged) - averaged time series (diagonal)']

# Iterate across layers and simulations
allIms = {}
for sim in groupedDirs:
    allIms[sim] = {}
    for layer in groupedDirs[sim]:

        currDirs = groupedDirs[sim][layer]

        indices = []
        for name in namesToSort:
            # print(name)
            idx = [i for i, s in enumerate(currDirs) if name in s]
            # print(idx)
            indices.append(idx[:])
        # print(indices)

        # Delete empty elements
        tidyIndices = list(filter(None, indices))
        tidyIndices = [idx[0] for idx in tidyIndices]
        # print(tidyIndices)

        # Reorder image array
        orderedPaths = []
        for idx in tidyIndices:
            orderedPaths.append(currDirs[idx])
        # print(orderedPaths)

        # Open images and get dimensions

        images = [Image.open(x) for x in orderedPaths]

        allIms[sim][layer] = images


##
from PIL import Image, ImageFont, ImageDraw

plot = allIms[sim][layer]

#
font = ImageFont.truetype("Times New Roman.ttf", 40, encoding="unic")


#

textIm = Image.new(mode = "RGB", size = (2500,100), color = "white")
drawText = ImageDraw.Draw(textIm)
drawText.text((100,20), sim+'_'+layer, font=font, fill='black')
#
concatPlot = get_concat_tile_resize([[textIm],
                                     [parametersImage],
                                     [images[0], images[1], images[2], images[3] ],
                                     [images[4], images[5], images[6], images[7] , images[8]]])

concatPlot.show()

plotNoTitle = get_concat_tile_resize([[parametersImage],
                                     [images[0], images[1], images[2], images[3]],
                                     [images[4], images[5], images[6], images[7], images[8]]])
plotNoTitle.show()
##
from fpdf import FPDF
pdf = FPDF(orientation = 'L', unit = 'mm', format='A4')

for sim in allPlots:
    for layer in allPlots[sim]:
        currIm = allPlots[sim][layer]
        pdf.add_page()
        pdf.image(currIm)

##
path, filename = os.path.split(datapath)
pdfName = os.path.join(datapath, filename+'_report.pdf')

frontPage = Image.new(mode="RGB", size=(2480, 3508), color="white")
drawText = ImageDraw.Draw(frontPage)
font = ImageFont.truetype("Times New Roman.ttf", 100, encoding="unic")
drawText.text((100, 20), 'Results: '+filename, font=font, fill='black')
frontPage.save(pdfName, "PDF", quality = 95)


##
for sim in allPlots:
    allPlots[sim] = collections.OrderedDict(sorted(allPlots[sim].items()))

