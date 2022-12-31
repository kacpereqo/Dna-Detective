from app.rna_translation import routes as rna_translation
from app.visualization_2d import routes as visualization
from app.mol_properties import routes as mol_properties

#|------------------------------------------------------------------------------|#

def include_routes(app):
    app.include_router(rna_translation.router)
    app.include_router(visualization.router)
    app.include_router(mol_properties.router)

#|------------------------------------------------------------------------------|#
