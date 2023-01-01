from app.rna_translation import routes as rna_translation
from app.visualization_2d import routes as visualization
from app.propeties import routes as propeties
from app.hydrophobicity import routes as hydrophobicity
from app.charge import routes as charge
from app.flexibility import routes as flexibility
from app.aromacity import routes as aromacity 

#|------------------------------------------------------------------------------|#

def include_routes(app):
    app.include_router(rna_translation.router)
    app.include_router(visualization.router)
    app.include_router(propeties.router)
    app.include_router(hydrophobicity.router)
    app.include_router(charge.router)
    app.include_router(flexibility.router)
    app.include_router(aromacity.router)

#|------------------------------------------------------------------------------|#
