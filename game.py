# -*- coding: utf-8 -*-


from astrobox.space_field import SpaceField
from dmitrieva import DmitrievaDrone


if __name__ == '__main__':
    asteroids_count = 15
    scene = SpaceField(
        speed=3,
        asteroids_count=asteroids_count,
    )
    drones = [DmitrievaDrone(asteroids_count=asteroids_count) for _ in range(5)]
    scene.go()


# Второй этап: зачёт!
