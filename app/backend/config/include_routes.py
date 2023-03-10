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
from services.auth import routes as auth
from services.profile import routes as profile
from services.bulkiness import routes as bulkiness
from services.recognition import routes as recognition
from services.structure import routes as structure

# |------------------------------------------------------------------------------|#


def include_routes(app):
    app.include_router(structure.router)
    app.include_router(bulkiness.router)
    app.include_router(recognition.router)
    app.include_router(profile.router)
    app.include_router(auth.router)
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
