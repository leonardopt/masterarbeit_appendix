## Imports
import collections
import glob
from PIL import Image, ImageDraw, ImageFont
import os
import pandas as pd
import re
from pdf2image import convert_from_bytes
import pylatex as pl

# Define some functions
# From https://github.com/nkmk/python-snippets/blob/4e232ef06628025ef6d3c4ed7775f5f4e25ebe19/notebook/pillow_concat.py#L66-L99

def create_parameter_table(dataPath):
    # Get and clean up parameters
    parameters = pd.read_csv(os.path.join(dataPath, 'sim_shared_parameters.csv'), index_col=False)
    parameters = parameters.drop(columns = ['layer_descr', 'model_names_1', 'model_names_2', 'model_names_3'])

    # Clean up

    parameters = parameters.rename(columns={"taskweights_1":"taskweights_aeh",
                               "taskweights_2":"taskweights_kon",
                               "taskweights_3":"taskweights_sum",
                               "phaseweights_1": "phaseweights_cue",
                               "phaseweights_2": "phaseweights_trial",
                               "phaseweights_3":"phaseweights_rest",
                               "response_1": "response_cue",
                               "response_2": "response_trial",
                               "response_3":"response_rest"
                               })


    # Edit parameter table

    # Save parameters as LateX table
    dataPath = '/Users/leonardo.pettini/Desktop/Masterarbeit_BCCN/ROI_task_positive_network/group'
    # Format floats
    pd.options.display.float_format = lambda x: '{:.0f}'.format(x) if int(x) == x else '{:,.3f}'.format(x)
    parameters = parameters.set_index('layer_number').T
    parametersLateX = parameters.to_latex()
    tableFileName = 'sim_shared_parameters.tex'
    parametersLateXPath =  os.path.join(dataPath, tableFileName)

    with open(parametersLateXPath,'w') as tf:
        tf.write(parametersLateX)

    print(parametersLateX)

    return parametersLateX, parametersLateXPath

def get_concat_h(im1, im2):
    newIm = Image.new('RGB', (im1.width + im2.width, im1.height))
    newIm.paste(im1, (0, 0))
    newIm.paste(im2, (im1.width, 0))
    return newIm

def get_concat_v_blank(im1, im2, color=(0, 0, 0)):
    dst = Image.new('RGB', (max(im1.width, im2.width), im1.height + im2.height), color)
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst

def get_concat_h_resize(im1, im2, resample=Image.BICUBIC, resize_big_image=True):
    if im1.height == im2.height:
        _im1 = im1
        _im2 = im2
    elif (((im1.height > im2.height) and resize_big_image) or
          ((im1.height < im2.height) and not resize_big_image)):
        _im1 = im1.resize((int(im1.width * im2.height / im1.height), im2.height), resample=resample)
        _im2 = im2
    else:
        _im1 = im1
        _im2 = im2.resize((int(im2.width * im1.height / im2.height), im1.height), resample=resample)
    dst = Image.new('RGB', (_im1.width + _im2.width, _im1.height))
    dst.paste(_im1, (0, 0))
    dst.paste(_im2, (_im1.width, 0))
    return dst

def get_concat_v_resize(im1, im2, resample=Image.BICUBIC, resize_big_image=True):
    if im1.width == im2.width:
        _im1 = im1
        _im2 = im2
    elif (((im1.width > im2.width) and resize_big_image) or
          ((im1.width < im2.width) and not resize_big_image)):
        _im1 = im1.resize((im2.width, int(im1.height * im2.width / im1.width)), resample=resample)
        _im2 = im2
    else:
        _im1 = im1
        _im2 = im2.resize((im1.width, int(im2.height * im1.width / im2.width)), resample=resample)
    dst = Image.new('RGB', (_im1.width, _im1.height + _im2.height))
    dst.paste(_im1, (0, 0))
    dst.paste(_im2, (0, _im1.height))
    return dst

# def get_concat_h_multi_blank(im_list):
#     _im = im_list.pop(0)
#     for im in im_list:
#         _im = get_concat_h_blank(_im, im)
#     return _im

def get_concat_h_multi_resize(im_list, resample=Image.BICUBIC):
    min_height = min(im.height for im in im_list)
    im_list_resize = [im.resize((int(im.width * min_height / im.height), min_height),resample=resample)
                      for im in im_list]
    total_width = sum(im.width for im in im_list_resize)
    dst = Image.new('RGB', (total_width, min_height))
    pos_x = 0
    for im in im_list_resize:
        dst.paste(im, (pos_x, 0))
        pos_x += im.width
    return dst

def get_concat_v_multi_resize(im_list, resample=Image.BICUBIC):
    min_width = min(im.width for im in im_list)
    im_list_resize = [im.resize((min_width, int(im.height * min_width / im.width)),resample=resample)
                      for im in im_list]
    total_height = sum(im.height for im in im_list_resize)
    dst = Image.new('RGB', (min_width, total_height))
    pos_y = 0
    for im in im_list_resize:
        dst.paste(im, (0, pos_y))
        pos_y += im.height
    return dst

def get_concat_tile_resize(im_list_2d, resample=Image.BICUBIC):
    im_list_v = [get_concat_h_multi_resize(im_list_h, resample=resample) for im_list_h in im_list_2d]
    return get_concat_v_multi_resize(im_list_v, resample=resample)

def get_sim_and_layer_names(dataPath):
    # Get data
    allFilePaths = []
    for file in os.listdir(dataPath):
        if file.endswith(".png") and not file.startswith('.'):
            allFilePaths.append(os.path.join(dataPath, file))

    # Get directory names
    allFileNames = []
    for dir in allFilePaths:
        drive, path = os.path.splitdrive(dir)
        path, filename = os.path.split(path)
        allFileNames.append(filename)

    # Get simulation names
    simNames = []
    for name in allFileNames:
        currName = re.search(r".+?(?=layer)", name)
        if currName is not None:
            currName = currName.group(0)
            if currName[0] != '.':
                simNames.append(currName[:-1])

    # Names of the simulations
    simNames = list(set(simNames))

    # Do same thing for layers
    layerNames = []
    for name in allFileNames:
        currName = re.search(r"(?=layer).[a-z]*-[0-9]", name)
        if currName is not None:
            currName = currName.group(0)
            layerNames.append(currName)

    # Names of the layers
    layerNames = list(set(layerNames))


    # Group paths by names
    groupedDirs = {}
    for sim in simNames:
        groupedDirs[sim] = {}
        for layer in layerNames:
            groupedDirs[sim][layer] = []
            for dir in allFilePaths:
                if sim + '_layer' in dir and layer in dir:
                    groupedDirs[sim][layer].append(dir)

    return simNames, layerNames, groupedDirs

def create_report(dataPath, groupedDirs, appendixN=None):

    # # Load data
    # # Get data
    # allFilePaths = []
    # for file in os.listdir(dataPath):
    #     if file.endswith(".png"):
    #         allFilePaths.append(os.path.join(dataPath, file))

   # # Get paths for comp model plots
    modelPath = glob.glob(os.path.join(dataPath, 'all_model_figures'))[0]

    # Load parameters
    parametersPNG_dir = os.path.join(dataPath, 'sim_shared_parameters.png')
    parametersImage = Image.open(parametersPNG_dir)

    # Order figures
    namesToSort = [
    'fig-similarity_matrix_mean_rho_complex_simple_separate',
    'fig-similarity_matrix_CI95_complex_simple_separate',
    'fig-similarity_matrix_mean_rho_complex_simple_averaged',
    'fig-similarity_matrix_CI95_complex_simple_averaged',
    'fig-mean_vs_CI95_complex_simple_separate',
    'fig-mean_vs_CI95_complex_simple_averaged',
    'fig-time_resolved_correlation_all_time_series_complex_simple_averaged',
    'fig-time_resolved_correlation_horizontal_averaged_time_series_complex_simple_averaged',
    'fig-time_resolved_correlation_diagonal_averaged_time_series_complex_simple_averaged']

    # Iterate across layers and simulations
    allPlots = {}
    for nSim, sim in enumerate(groupedDirs):
        allPlots[sim] = {}

        for layer in groupedDirs[sim]:

            currDirs = groupedDirs[sim][layer]

            indices = []
            for name in namesToSort:
                #print(name)
                idx = [i for i, s in enumerate(currDirs) if name in s]
                #print(idx)
                indices.append(idx[:])
            #print(indices)

            # Delete empty elements
            tidyIndices = list(filter(None, indices))
            tidyIndices = [idx[0] for idx in tidyIndices]
            #print(tidyIndices)

            # Reorder image array
            orderedPaths = []
            for idx in tidyIndices:
                orderedPaths.append(currDirs[idx])
            #print(orderedPaths)


            images = [Image.open(x) for x in orderedPaths]
            if images is None:
                print(f'Paths: {orderedPaths}')
                print(f'Simulation: {sim}, {layer}')

            widths, heights = zip(*(i.size for i in images))

            # Open images and get dimensions
            allModelFigs = []
            for file in os.listdir(modelPath):
                if file.endswith(".png") and file[0]!='.':
                    allModelFigs.append(os.path.join(modelPath, file))

            # Open model figures
            modelNames = []
            modelDirs = []
            for PNGfile in allModelFigs:
                if sim[4:] in PNGfile and layer.replace('-', '_') in PNGfile:
                    modeldir, filename = os.path.split(PNGfile)
                    modelNames.append(filename)
                    modelDirs.append(modeldir)

            # Give sorted names
            sortedModelNames = ['models_before_convolution.png',
                                'models_before_convolution_with_voxel_noise.png',
                                'models_after_convolution_with_signal_noise.png']
            modelIndices = []
            for name in sortedModelNames:
                idx = [i for i, s in enumerate(modelNames) if name in s]
                # print(idx)
                modelIndices.append(idx[:])

            # Reorder image array according to indices
            orderedModelFigs = []
            orderedModelFigDirs = []
            for idx in modelIndices:
                orderedModelFigs.append(modelNames[idx[0]])
                orderedModelFigDirs.append(os.path.join(modelDirs[idx[0]], modelNames[idx[0]]))

            # Load images and create text files with title
            font = ImageFont.truetype("Times New Roman.ttf", 60, encoding="unic")
            modelFigs = []
            titleIms = []
            for n, dir in enumerate(orderedModelFigDirs):
                currIm = Image.open(dir)
                modelFigs.append(currIm.crop((300, 500, 4600, 1850)))
                currTitleIm = Image.new(mode="RGB", size=(currIm.width, 100), color="white")
                drawText = ImageDraw.Draw(currTitleIm)
                drawText.text((500, 20), sortedModelNames[n][:-4], font=font, fill='black')
                titleIms.append(currTitleIm)

            # Initialise A4 page
            width, height = int(8.27 * 300), int(11.7 * 300)  # A4 at 300dpi
            page = Image.new('RGB', (width, height), 'white')

            # Add title
            font = ImageFont.truetype("Times New Roman.ttf", 40, encoding="unic")
            textIm = Image.new(mode="RGB", size=(2500, 100), color="white")
            drawText = ImageDraw.Draw(textIm)
            drawText.text((10, 10), f'Appendix {appendixN}.{nSim+1}.{layer[-1]}: {sim}_{layer}', font=font, fill='black')
            page.paste(textIm, box=(450, 150))

            # Models
            # Merge model figures
            modelsIm = get_concat_tile_resize([ [titleIms[0]],
                                                [modelFigs[0]],
                                                [titleIms[1]],
                                                [modelFigs[1]],
                                                [titleIms[2]],
                                                [modelFigs[2]]])
            parsz = tuple([round(z * 0.8) for z in parametersImage.size])
            modsz = tuple([round(z * 0.13) for z in modelsIm.size])
            parIm = parametersImage.resize(parsz)
            modIm = modelsIm.resize(modsz)
            page.paste(parIm, box=(350, 260))
            page.paste(modIm, box=(1200, 280))

            # Plots
            # New sizes
            simMatsSz = (1000, 970)
            meanciSz = (450, 800)
            timeserSz = (860, 810)
            meanSepCS = images[0].resize(simMatsSz)
            CISepCS = images[1].resize(simMatsSz)
            meanAvCS = images[2].resize(simMatsSz)
            CIAvCS = images[3].resize(simMatsSz)
            meanciSepCS = images[4].resize(meanciSz)
            meanciAvCS = images[5].resize(meanciSz)
            allCurves = images[6].resize(timeserSz)
            horizCurve = images[7].resize(timeserSz)
            diagCurve = images[8].resize(timeserSz)
            simMatsSz = tuple([round(z * 0.4) for z in images[0].size])
            meanciSz = tuple([round(z * 0.42) for z in images[4].size])
            timeserSz = tuple([round(z * 0.4) for z in images[6].size])
            meanSepCS = images[0].resize(simMatsSz)
            CISepCS = images[1].resize(simMatsSz)
            meanAvCS = images[2].resize(simMatsSz)
            CIAvCS = images[3].resize(simMatsSz)
            meanciSepCS = images[4].resize(meanciSz)
            meanciAvCS = images[5].resize(meanciSz)
            allCurves = images[6].resize(timeserSz)
            horizCurve = images[7].resize(timeserSz)
            diagCurve = images[8].resize(timeserSz)
            # Determine offsets and paste
            # Horiz off for all
            yall = 300
            # Sim matrices
            xmat1 = 300
            xmat2 = xmat1 + 800
            ymat1 = 600 + yall
            ymat2 = ymat1 + 700
            page.paste(meanSepCS, box=(xmat1, ymat1))
            page.paste(CISepCS, box=(xmat2, ymat1))
            page.paste(meanAvCS, box=(xmat1, ymat2))
            page.paste(CIAvCS, box=(xmat2, ymat2))
            # mean vs. CI
            xmn1 = xmat2 + 850
            ymn = 20
            page.paste(meanciSepCS, box=(xmn1, ymat1 + ymn))
            page.paste(meanciAvCS, box=(xmn1, ymat2 + ymn))
            # Time resolved corr
            xts1 = 315
            xoff = 700
            xts2 = xts1 + xoff
            xts3 = xts2 + xoff
            yts = 2100 + yall
            page.paste(allCurves, box=(xts1, yts))
            page.paste(horizCurve, box=(xts2, yts))
            page.paste(diagCurve, box=(xts3, yts))

            # Add directory
            rootDir = '/analysis/share/corinna_kai/afx/'
            simulDir = modelPath.split('/')[3]
            simDir = os.path.join(simulDir, 'similarity-analysis_from-FIR1_'+ sim, 'group', layer)
            sourceDir = os.path.join(rootDir, simDir)
            # Add to page
            font = ImageFont.truetype("Times New Roman.ttf", 25, encoding="unic")
            dirIm = Image.new(mode="RGB", size=(2500, 100), color="white")
            drawDirIm = ImageDraw.Draw(dirIm)
            drawDirIm.text((10, 10), sourceDir, font=font, fill='grey')
            page.paste(dirIm, box=(450, 3200))


            # Store in dictionary (with model figures)
            allPlots[sim][layer] = page

            print(sim+"_"+layer+" - done")

    return allPlots

def generate_pdf(allPlots, dataPath, fileName=None, PDFtitle=None):

    if not fileName:
        path, filename = os.path.split(dataPath)
        pdfName = os.path.join(dataPath, filename+'_report.pdf')
    else:
        filename = fileName
        pdfName = os.path.join(dataPath, fileName)

    # Initialise A4 page
    width, height = int(8.27 * 300), int(11.7 * 300)  # A4 at 300dpi
    frontPage = Image.new(mode="RGB", size=(width, height), color="white")
    drawText = ImageDraw.Draw(frontPage)
    font = ImageFont.truetype("Times New Roman.ttf", 72, encoding="unic")
    drawText.text((400, 350), PDFtitle, font=font, fill='black')
    frontPage.save(pdfName, "PDF", quality=95)

    allPlots = collections.OrderedDict(sorted(allPlots.items()))

    for sim, items in allPlots.items():
        allPlots[sim] = collections.OrderedDict(sorted(allPlots[sim].items()))
        # simPage = Image.new(mode="RGB", size=(2500, 1510), color="white")
        # simText = ImageDraw.Draw(simPage)
        # simText.text((100, 20), sim, font=font, fill='black')
        # simPage.save(pdfName, "PDF", quality=95, append=True)

        for layer, items in allPlots[sim].items():
            # currPlotTitle = allPlots[sim][layer][0]
            # currPlotAnalysis = allPlots[sim][layer][1]
            # currPlotTitle.resize((currPlotAnalysis.height, currPlotAnalysis.width)).save(pdfName, "PDF", quality = 95, append=True)
            allPlots[sim][layer].save(pdfName, "PDF", quality = 95, append=True)

    return pdfName





## Create report
# Determine base directory
basedir = '/Volumes/Transcend/masterarbeit_figures'
simulationsToDo = [
                    'derivatives_test_sim_control_noiselevels',
                    'derivatives_test_sim_control_weightlevels',
                    'derivatives_test_sim_control_mixedeffectvoxels5',
                    'derivatives_test_sim_control_mixedeffectvoxels10'
                    ]

allDatapaths = [os.path.join(basedir, sim) for sim in simulationsToDo]

##
for n, dataPath in enumerate(allDatapaths):

    ## Create parameter table
    parametersLateX, parametersLateXPath = create_parameter_table(dataPath)
    print("LateX code for parameters table created")

    ##
    # Create PDF with LaTeX table
    doc = pl.Document()
    doc.packages.append(pl.Package('booktabs'))

    # Difference to the other answer:
    with doc.create(pl.Table(position='htbp')) as table:
        table.add_caption('Simulation parameters')
        table.append(pl.Command('centering'))
        table.append(pl.NoEscape(parametersLateX))
    #
    doc.generate_pdf(os.path.join(dataPath, 'sim_shared_parameters'), clean_tex=False)

    ## Convert PDF into PNG
    pathPDF = os.path.join(dataPath, 'sim_shared_parameters.pdf')
    parametersIm = convert_from_bytes(open(pathPDF, 'rb').read())[0]
    pathPNG = os.path.join(dataPath, 'sim_shared_parameters.png')
    # Crop image to its minimum size
    parametersIm = parametersIm.crop((400, 350, 1300, 1100))
    parametersIm.save(pathPNG, "PNG", quality=95)
    print("Parameters file converted to PNG")

    ## Get grouped directories
    simNames, layerNames, groupedDirs = get_sim_and_layer_names(dataPath)
    # Order dictionary
    groupedDirs = collections.OrderedDict(sorted(groupedDirs.items()))

    ## Create pages for the report
    allPlots = create_report(dataPath, groupedDirs, appendixN=f'C.{n+1}')
    print("All figures collected")

    ## Generate PDF
    pdf = generate_pdf(allPlots, basedir, fileName=simulationsToDo[n] + '.pdf', PDFtitle=f'Appendix C.{n+1}: simulated data - {simulationsToDo[n][29:]}')
    print("PDF generated")



