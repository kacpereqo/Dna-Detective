from app.rna_translation import routes as rna_translation
from app.visualization_2d import routes as visualization

def include_routes(app):
    app.include_router(rna_translation.router)
    app.include_router(visualization.router)