## Imports
import collections
import glob
from PIL import Image, ImageDraw, ImageFont
import os
import pandas as pd
import re
from pdf2image import convert_from_bytes
import pylatex as pl
import re
import ocrmypdf


def load_images(currDirs):

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

    images = [Image.open(x) for x in orderedPaths]
    widths, heights = zip(*(i.size for i in images))

    return images

def add_figures_to_page(images, roiname=False, appendixN=None):

    # Initialise A4 page
    width, height = int(8.27 * 300), int(11.7 * 300)  # A4 at 300dpi
    page = Image.new('RGB', (width, height), 'white')
    # Add title
    font = ImageFont.truetype("Times New Roman.ttf", 40, encoding="unic")
    titleIm = Image.new(mode="RGB", size=(height, 100), color="white")
    drawText = ImageDraw.Draw(titleIm)

    if appendixN:
        if roiname:
            drawText.text((10, 10), f'Appendix {appendixN}: {roiname}', font=font, fill='black')
        else:
            drawText.text((10, 10), f'Appendix {appendixN}: whole brain', font=font, fill='black')
    else:
        drawText.text((10, 10), '', font=font, fill='black')

    page.paste(titleIm, box=(450, 150))

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
    rootDir = '/analysis/share/corinna_kai/afx/derivatives/similarity-analysis_from-FIR1'
    if roiname:
        rootDir = rootDir + '_ROI'
        sourceDir = os.path.join(rootDir, roiname, 'group')
    else:
        sourceDir = os.path.join(rootDir, 'group')

    # Add to page
    font = ImageFont.truetype("Times New Roman.ttf", 25, encoding="unic")
    dirIm = Image.new(mode="RGB", size=(2500, 100), color="white")
    drawDirIm = ImageDraw.Draw(dirIm)
    drawDirIm.text((10, 10), sourceDir, font=font, fill='grey')
    page.paste(dirIm, box=(450, 3200))

    return page

def initialise_PDF_report(PDFtitle, PDFoutputDir):

    # Initialise A4 page
    width, height = int(8.27 * 300), int(11.7 * 300)  # A4 at 300dpi
    frontPage = Image.new(mode="RGB", size=(width, height), color="white")
    drawText = ImageDraw.Draw(frontPage)
    font = ImageFont.truetype("Times New Roman.ttf", 72, encoding="unic")
    drawText.text((400, 350), PDFtitle, font=font, fill='black')
    frontPage.save(PDFoutputDir, "PDF", quality=95)

def get_concat_h(im1, im2, im3):
    newIm = Image.new('RGB', (im1.width + im2.width + im3.width, im1.height))
    newIm.paste(im1, (0, 0))
    newIm.paste(im2, (im1.width, 0))
    newIm.paste(im3, (im1.width+im2.width, 0))
    return newIm


## Base directory
basedir = '/Volumes/Transcend/masterarbeit_figures'

## Make report - whole brain

# Get directory
wholebrainDir = os.path.join(basedir, 'empirical_whole_brain')
wholebrainFileDirs = glob.glob(os.path.join(wholebrainDir, '*.png'))
# Output file name
pdfWholebraindir = os.path.join(basedir, 'empirical_whole_brain_report.pdf')
initialise_PDF_report(f'Appendix A: empirical data - whole brain', pdfWholebraindir)

# Load images
images = load_images(wholebrainFileDirs)
page = add_figures_to_page(images, appendixN=f'A.{1}')
# Open brain image
currBrainImages = [Image.open(x).crop((300, 150, 2000, 1200)) for x in glob.glob(os.path.join(wholebrainDir, 'brain_fig', '*.png'))]
# Concatenate all brain pics
brainPic = get_concat_h(currBrainImages[0], currBrainImages[1], currBrainImages[2])
sizeMask = tuple([round(z * 0.4) for z in brainPic.size])
brainPic = brainPic.resize(sizeMask)
# Paste brain pic to the page
page.paste(brainPic, box=(315, 300))
# Save in dir
page.save(pdfWholebraindir, "PDF", quality=95, append=True)
# page.show()

## Make report - ROIs

# Get directory
ROIDir = os.path.join(basedir, 'empirical_ROIs')
# Output file name
pdfROIdir = os.path.join(basedir, 'empirical_ROIs_report.pdf')
initialise_PDF_report(f'Appendix B: empirical data - Regions of Interest (ROIs)', pdfROIdir)
# Load mask renderings
maskDir = os.path.join(ROIDir, 'masks')
# Get all mask png files
maskFilesAll = glob.glob(os.path.join(maskDir, '*.png'))
# Open images
roiFileDirs = glob.glob(os.path.join(ROIDir, '*.png'))
# Get ROI names
ROInames = list(set([re.search("(?<=ROI_riw)(.*)(?=_fig-)", x).group(0) for x in roiFileDirs if re.search("ROI_riw(.*?)_fig-", x) is not None]))


##
# Load images
for n, roi in enumerate(ROInames):
    # Get figures for current ROI
    currDir = [dir for dir in roiFileDirs if roi in dir]
    # Load images
    images = load_images(currDir)
    # Add page
    page = add_figures_to_page(images, roiname=roi, appendixN=f'B.{n+1}')
    # Open mask files
    currMaskFiles = [x for x in maskFilesAll if roi in x]
    currMaskFiles.sort()
    currMaskImages = [Image.open(x).crop((300, 150, 2000, 1200)) for x in currMaskFiles]
    # Concatenate all mask pics
    maskPic = get_concat_h(currMaskImages[0], currMaskImages[1], currMaskImages[2])
    sizeMask = tuple([round(z * 0.4) for z in maskPic.size])
    maskPic = maskPic.resize(sizeMask)
    # Paste mask pic to the page
    page.paste(maskPic, box=(315, 300))
    page.save(pdfROIdir, "PDF", quality=95, append=True)



##
page.show()

