from services.rna_translation import routes as rna_translation
from services.visualization_2d import routes as visualization
from services.propeties import routes as propeties
from services.hydrophobicity import routes as hydrophobicity
from services.charge import routes as charge
from services.flexibility import routes as flexibility
from services.aromacity import routes as aromacity
from services.instability import routes as instability
from services.statistics import routes as statics
from services.caching import routers as caching
from services.facts import routers as facts
from services.polarity import routes as polarity

# |------------------------------------------------------------------------------|#


def include_routes(app):
    app.include_router(polarity.router)
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
