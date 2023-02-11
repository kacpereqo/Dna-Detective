from indigo import *
from indigo_renderer import *

indigo = Indigo()
renderer = IndigoRenderer(indigo)

mol = indigo.loadMolecule("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
mol.layout()  # if not called, will be done automatically by the renderer
indigo.setOption("render-output-format", "png")
indigo.setOption("render-comment", "Caffeine")
indigo.setOption("render-comment-position", "top")
indigo.setOption("render-image-size", 200, 250)
indigo.setOption("render-background-color", 1.0, 1.0, 1.0)
renderer.renderToFile(mol, "caffeine.png")
