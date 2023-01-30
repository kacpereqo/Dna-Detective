from fastapi import APIRouter
from random import choice

router = APIRouter()


@router.get("/fact")
async def get_fact():
    facts = [
        "Słonia afrykańskiego można rozpoznać po unikalnym uzębieniu - każde zęby jest inne.",
        "Kręgi na ogonie węża rosną razem z nim i nie odrastają po utracie.",
        "Serce narwala może ważyć nawet nieco ponad 20 kg i mieć do 1,5 m średnicy.",
        "Wiele gatunków ryb potrafi zmieniać kolor, aby się ukryć lub przestraszyć drapieżnika.",
        "Żyrafy mają najdłuższe szyje ze wszystkich ssaków proporcjonalnie do ich ciała.",
        "Drzewa komunikują się ze sobą przez wymianę chemikaliów i dźwięków.",
        "Ślimaki posiadają dwa pary ramion, jedną do poruszania się, a drugą do jedzenia.",
        "Kraby kameleona potrafią zmieniać kolor, aby pasować do otoczenia lub wywołać strach u drapieżnika.",
    ]
    return {"fact": choice(facts)}
