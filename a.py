import random
from superwires import games, color
```python

games.init(screen_width=640, screen_height=480, fps=50)


class Pan(games.Sprite):
    """Patelnia"""
    image = games.load_image("patelnia.bmp")

    def __init__(self):
        super().__init__(image=Pan.image, x=games.mouse.x, bottom=games.screen.height)
        self.score = games.Text(value=0, size=25, color=color.black,
                                top=5, right=games.screen.width - 20)
        games.screen.add(self.score)

    def update(self):
        self.x = games.mouse.x
        if self.left < 0:
            self.left = 0
        if self.right > games.screen.width:
            self.width = games.screen.width
        self.check_catch()

    def check_catch(self):
        for pizza in self.overlapping_sprites:
            self.score.value += 10
            self.score.right = games.screen.width - 20
            pizza.handle_caught()


class Pizza(games.Sprite):
    image = games.load_image("pizza.bmp")
    speed = 1

    def __init__(self, x, y=90, speed=1):
        super().__init__(image=Pizza.image, x=x, y=y, dy=speed)

    def update(self):
        if self.bottom > games.screen.height:
            self.end_game()
            self.destroy()

    def handle_caught(self):
        self.destroy()

    def end_game(self):
        end_message = games.Message(value="KONIEC GRY!", size=90, color=color.red, x=games.screen.width / 2,
                                    y=games.screen.height / 2, lifetime=5 * games.screen.fps, after_death=games.screen.quit)
        games.screen.add(end_message)


class Chef(games.Sprite):
    """Szef kuchni"""
    image = games.load_image("kucharz.bmp")

    def __init__(self, y=60, speed=2, odds_change=200):
        super().__init__(image=Chef.image, x=games.screen.width / 2, y=y, dx=speed)
        self.odds_change = odds_change
        self.time_til_drop = 0
        self.dropped_pizzas = 0

    def update(self):
        if self.left < 0 or self.right > games.screen.width:
            self.dx = -self.dx
        elif random.randrange(self.odds_change) == 0:
            self.dx = -self.dx
        self.check_drop()

    def check_drop(self):
        if self.time_til_drop > 0:
            self.time_til_drop -= 1
        else:
            new_pizza = Pizza(x=self.x, speed=2 ** self.dropped_pizzas)
            games.screen.add(new_pizza)
            self.time_til_drop = int(new_pizza.height * 1.5 / Pizza.speed) + 1


def main():
    wall_image = games.load_image("sciana.jpg", transparent=False)
    games.screen.background = wall_image

    the_chef = Chef()
    games.screen.add(the_chef)

    the_pan = Pan()
    games.screen.add(the_pan)

    games.mouse.is_visible = False
    games.screen.event_grab = True

    games.screen.mainloop()


main()
