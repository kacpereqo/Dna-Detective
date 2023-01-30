from backend.services.rna_translation import routes as rna_translation
from backend.services.visualization_2d import routes as visualization
from backend.services.propeties import routes as propeties
from backend.services.hydrophobicity import routes as hydrophobicity
from backend.services.charge import routes as charge
from backend.services.flexibility import routes as flexibility
from backend.services.aromacity import routes as aromacity
from backend.services.instability import routes as instability
from backend.services.statistics import routes as statics
from backend.services.caching import routers as caching
from backend.services.facts import routers as facts

# |------------------------------------------------------------------------------|#


def include_routes(app):
    app.include_router(caching.router)
    app.include_router(rna_translation.router)
    app.include_router(visualization.router)
    app.include_router(propeties.router)
    app.include_router(hydrophobicity.router)
    app.include_router(charge.router)
    app.include_router(flexibility.router)
    app.include_router(aromacity.router)
    app.include_router(instability.router)
    app.include_router(statics.router)
    app.include_router(facts.router)

# |------------------------------------------------------------------------------|#
