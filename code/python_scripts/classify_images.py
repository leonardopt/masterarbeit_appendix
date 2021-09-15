import re
import os
from PIL import Image

def get_sim_and_layer_names(datapath):

    ## Get data
    allFilePaths = []
    for file in os.listdir(datapath):
        if file.endswith(".png"):
            allFilePaths.append(os.path.join(datapath, file))

    ## Get directory names
    allFileNames = []
    for dir in allFilePaths:
        drive, path = os.path.splitdrive(dir)
        path, filename = os.path.split(path)
        allFileNames.append(filename)

    ## Get simulation names

    simNames = []
    for name in allFileNames:
        currName = re.search(r".+?(?=layer)", name)
        if currName is not None:
            currName = currName.group(0)
            if currName[0]!= '.':
                simNames.append(currName[:-1])

    # Names of the simulations
    simNames = set(simNames)


    ## Do same thing for layers
    layerNames = []
    for name in allFileNames:
        currName = re.search(r"(?=layer).[a-z]*-[0-9]", name)
        if currName is not None:
            currName = currName.group(0)
            layerNames.append(currName)

    # Names of the layers
    layerNames = set(layerNames)

    return simNames, layerNames