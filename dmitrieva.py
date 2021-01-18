# -*- coding: utf-8 -*-
from astrobox.core import Drone
from robogame_engine.geometry import Point, Vector
from robogame_engine.theme import theme
from astrobox.themes.default import MOTHERSHIP_HEALING_DISTANCE


class DmitrievaDrone(Drone):
    """
    'middle_harvester' - сборщик элериума с центра поля
    'harvester' - сборщик элериума со всех объектов
    'base_destroyer' - дрон, разрушающий близжайшую базу врага
    'defender' - защитник базы
    'terminator' - активно атакующий дрон
    """
    list_of_drones = []
    all_complete = set()
    index = 0
    role = None
    not_completely_full_range_of_flight, empty_range_of_flight, full_range_of_flight = 0, 0, 0

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.list_of_drones.append(self)

    def on_born(self):
        self.index = self.list_of_drones.index(self)
        self.get_role(first_step=True)

    def on_stop_at_asteroid(self, asteroid):
        if self.second_step():
            self.get_role()
        elif (self.role == 'harvester' or self.role == 'middle_harvester') and self.target in self.asteroids:
            if self.payload + self.target.payload < 100 and self._get_nearest_asteroids(index=1):
                next_target = self._get_nearest_asteroids(index=1)
            else:
                next_target = self.mothership
            self.load_from(asteroid)
            self.turn_to(next_target)

    def on_load_complete(self):
        if not self.is_full:
            if self._get_nearest_asteroids():
                self.target = self._get_nearest_asteroids()
            elif self.near_enemy_to_load_from():
                self.target = self.get_enemy(for_steal=True)
        else:
            self.target = self.my_mothership
        self.move_at(self.target)

    def on_stop_at_mothership(self, mothership):
        if mothership == self.my_mothership:
            if self.is_empty and self._all_objects_empty():
                self.all_complete.add(self)
                self.check_game_over()
            else:
                self.unload_to(mothership)
        else:
            self.load_from(mothership)

    def on_unload_complete(self):
        if self.role == 'middle_harvester':
            self.move_at(self.get_position_circle())
        if self._all_objects_empty():
            self.all_complete.add(self)
            self.check_game_over()

    def on_wake_up(self):
        if self.second_step():
            self.get_role()
        elif self.role == 'middle_harvester':
            self.middle_harvester_act()
        elif self.role == 'defender':
            self.defender_act()
        elif self.role == 'base_destroyer':
            self.base_destroyer_act()
        elif self.role == 'terminator':
            self.terminator_act()
        else:
            self.harvester_act()

    def on_heartbeat(self):
        self.check_drone_health()
        if self.role == 'harvester' and self.near_enemy_to_load_from():
            self.load_from(self.get_enemy(for_steal=True))
        elif self._all_objects_empty():
            self.choose_target_and_move_at(self.my_mothership)

    def get_role(self, first_step=False, last_step=False, fourth_step=False):
        if first_step:
            self.role = 'middle_harvester'
        elif fourth_step:
            self.role = 'terminator'
        elif last_step:
            self.role = 'harvester'
        else:  # third step
            if self.can_shoot(target=self.get_enemies_mothership()):
                self.role = 'base_destroyer'
            else:
                self.role = 'defender'

    def middle_harvester_act(self):
        """

        На первом этапе сражения все дроны получают эту роль. Далее идет более конкретная разбивка
        """
        self.choose_target_and_move_at(self._get_nearest_asteroids(index=self.index, middle=True))

    def base_destroyer_act(self):
        """

        После распределения ролей, дрон по возможности становится base_destroyer, так как враг без базы не представляет
        опасности.
        """
        if not self.near(self.get_position_circle()):
            self.move_at(self.get_position_circle())
        else:
            if self.get_enemies_mothership() and self.can_shoot(self.get_enemies_mothership()):
                self.target = self.get_enemies_mothership()
                self.turn_to_and_shot_target()
            else:
                self.get_role()

    def defender_act(self):
        """

        Дроны распределяются вокруг базы с размахом sweep в радиусе действия защитного поля и бьют дронов врага.
        """
        if not self.near(self.get_position_circle()):
            self.move_at(self.get_position_circle())
        else:
            self.target = self.get_enemy()
            self.shoot_enemy()
            if self.can_shoot(self.get_enemies_mothership()):
                self.get_role()
            elif self.is_safety_to_attack():
                self.get_role(fourth_step=True)
            elif not self.get_enemy():
                self.get_role(last_step=True)

    def terminator_act(self):
        """

        Дрон выходит за пределы защитного поля своей базы и делает всё, чтобы противники умерли.
        """
        self.target = self.get_enemy(common_target=True)
        victim = self.target.my_mothership.coord if self.target and self.target.my_mothership.is_alive else None
        go_aggressive_attack = self.is_safety_to_attack(one_enemy=True) and victim
        self.target = victim if go_aggressive_attack else self.get_enemy(common_target=True)
        distance = self.gun.shot_distance if go_aggressive_attack else MOTHERSHIP_HEALING_DISTANCE * 1.5
        start_point = victim if go_aggressive_attack else None
        position = self.get_position_circle(distance=distance, start_point=start_point)

        if not self.target:
            self.get_role(last_step=True)
        elif not self.is_safety_to_attack():
            self.get_role()
        elif self.near(position):
            self.turn_to_and_shot_target()
        else:
            self.move_at(position)

    def harvester_act(self):
        """

        На последнем этапе сражения все дроны получают эту роль. Задача - сбор элериума со всех объектов
        """
        if self.get_enemies_mothership():
            self.target = self.get_enemies_mothership()
            position = self.get_position_circle(start_point=self.target.coord, distance=self.gun.shot_distance)
            if self.near(position):
                self.turn_to_and_shot_target()
            else:
                self.move_at(position)
        elif self.get_enemy(for_steal=True):  # хоть 1 умерший враг
            self.choose_target_and_move_at(self.get_enemy(for_steal=True))
        elif self._get_nearest_asteroids(index=len(self.harvesters_distribution())):  # каждому по астероиду
            self.choose_target_and_move_at(self._get_nearest_asteroids(index=self.get_personal_asteroid_index()))
        elif self._get_nearest_asteroids():  # хоть 1 астероид
            self.choose_target_and_move_at(self._get_nearest_asteroids())
        elif self.get_enemies_mothership(for_steal=True):  # хоть 1 база врага
            self.choose_target_and_move_at(self.get_enemies_mothership(for_steal=True))
        else:
            self.choose_target_and_move_at(self.my_mothership)

    def choose_target_and_move_at(self, target):
        self.target = target
        self.move_at(self.target)

    def get_position_circle(
            self, distance=MOTHERSHIP_HEALING_DISTANCE * 0.9, start_point=None, end_point=None, sweep=45):
        position = Circle(self.my_mothership, self.list_of_drones, self.index)
        return position.get_position_circle(distance, start_point, end_point, sweep)

    def get_enemies_mothership(self, for_destroy=True, for_steal=False):
        all_motherships = []
        for enemies_mothership in self.scene.motherships:
            if for_destroy and enemies_mothership != self.my_mothership and enemies_mothership.is_alive:
                all_motherships.append([enemies_mothership, self.distance_to(enemies_mothership)])

            elif for_steal and not enemies_mothership.is_alive and not enemies_mothership.is_empty:
                all_motherships.append([enemies_mothership, self.distance_to(enemies_mothership)])
        if all_motherships:
            all_motherships.sort(key=lambda x: x[1])
            return all_motherships[0][0]
        return None

    def get_enemy(self, number=0, common_target=False, for_kill=True, for_steal=False):
        enemies = []
        starting_point = self.my_mothership if common_target else self
        for enemy in self.scene.drones:
            if for_kill and enemy and enemy.team != self.team and enemy.is_alive:
                enemies.append((enemy, starting_point.distance_to(enemy)))
            elif for_steal and enemy and enemy.team != self.team and not enemy.is_alive and not enemy.is_empty:
                enemies.append((enemy, starting_point.distance_to(enemy)))
        if number >= len(enemies):
            return None
        if enemies:
            enemies.sort(key=lambda x: x[1])
            return enemies[number][0]
        return None

    def is_safety_to_attack(self, one_enemy=False):
        enemies = []
        for enemy in self.scene.drones:
            if enemy.team != self.team and enemy.is_alive:
                enemies.append((enemy, self.my_mothership.distance_to(enemy)))
        if one_enemy:
            return len(enemies) == 1 and len(self.list_of_drones) >= 3
        return len(enemies) + 1 < len(self.list_of_drones)

    def shoot_enemy(self):
        number = 1
        while not self.can_shoot(target=self.target):
            if self.get_enemy(number):
                self.target = self.get_enemy(number)
                number += 1
            else:
                self.target = None
                break
        if self.can_shoot(self.target):
            self.turn_to_and_shot_target()

    def turn_to_and_shot_target(self):
        self.turn_to(self.target)
        self.gun.shot(self.target)

    def _has_decision_linear(self, x0, y0, x, y, xm, ym, r):
        b = (x * y0 - x0 * y) / (x - x0)
        k = (y - y0) / (x - x0)
        disc = 4 * (b * k - k * ym - xm) ** 2 - 4 * (1 + k ** 2) * (b ** 2 - 2 * b * ym + xm ** 2 + ym ** 2 - r ** 2)
        return disc > 0

    def _has_decision_straight(self, x, xm, r):
        return r ** 2 >= (x - xm) ** 2

    def _sufficient_distance(self, distance_points):
        return distance_points <= self.gun.shot_distance

    def _no_teammates_at_gunpoint(self, x0, y0, x, y, r):
        drones, permission = [drone for drone in self.list_of_drones if drone != self], []
        for drone in drones:
            xm, ym = drone.coord.x, drone.coord.y
            if x == x0:
                permission.append(not self._has_decision_straight(x, xm, r))
            elif y == y0:
                permission.append(not self._has_decision_straight(y, ym, r))
            else:
                permission.append(not self._has_decision_linear(x0, y0, x, y, xm, ym, r))
        return all(permission)

    def _sufficient_distance(self, distance_points):
        return distance_points <= self.gun.shot_distance

    def can_shoot(self, target, save_distance=50):
        if not target:
            return False
        coordinates_from = self.coord
        x0, y0, x, y = coordinates_from.x, coordinates_from.y, target.x, target.y
        r = save_distance
        return self._sufficient_distance(self.distance_to(target)) and self._no_teammates_at_gunpoint(x0, y0, x, y, r)

    def second_step(self):
        return self.role == 'middle_harvester' and self.near(self.get_position_circle())

    def harvesters_distribution(self):
        amount, drone_distribution = 0, {}
        for drone in self.list_of_drones:
            if drone.role == 'harvester':
                drone_distribution[drone.index] = amount
                amount += 1
        return drone_distribution

    def get_personal_asteroid_index(self):
        drone_distribution = self.harvesters_distribution()
        return drone_distribution[self.index]

    def check_drone_health(self, safe_share_of_health=0.6):
        if self.health == 0 and self in self.list_of_drones:
            self.list_of_drones.remove(self)
            for drone in self.list_of_drones:
                drone.index = self.list_of_drones.index(drone)
        elif 0 < self.health < theme.DRONE_MAX_SHIELD * safe_share_of_health:
            self.choose_target_and_move_at(self.my_mothership)

    def near_enemy_to_load_from(self):
        return self.get_enemy(for_steal=True) and self.distance_to(self.get_enemy(for_steal=True)) < 1

    def check_game_over(self):
        if len(self.all_complete) == len(self.list_of_drones):
            self._show_statistics()

    def _all_objects_empty(self):
        objects_empty = [
            self._get_nearest_asteroids() is None,
            self.get_enemies_mothership(for_steal=True) is None,
            self.get_enemy(for_steal=True) is None
        ]
        return all(objects_empty)

    def _get_nearest_asteroids(self, index=0, middle=False, position=False):
        distances_to_asteroids = {}
        for num, asteroid in enumerate(self.asteroids):
            distances_to_asteroids[asteroid] = self.distance_to(asteroid)

        distances_to_asteroids_sorted = sorted(distances_to_asteroids.items(), key=lambda item: item[1])
        if position:
            return distances_to_asteroids_sorted[index][0]
        fulls = [asteroid[0] for asteroid in distances_to_asteroids_sorted if asteroid[0].payload >= 100]
        not_empty = [asteroid[0] for asteroid in distances_to_asteroids_sorted if asteroid[0].payload > 0]

        if middle:
            middle_index = len(self.asteroids) // 2
            middle_index -= len(self.list_of_drones) // 2
            fulls = fulls[middle_index:middle_index + len(self.list_of_drones)]
        if len(fulls) >= index + 1:
            return fulls[index]
        elif len(not_empty) >= index + 1:
            return not_empty[index]
        return None

    def _show_statistics(self):
        self.stop()
        not_completely_full_range_of_flight, empty_range_of_flight, full_range_of_flight = 0, 0, 0
        for drone in self.list_of_drones:
            not_completely_full_range_of_flight += drone.not_completely_full_range_of_flight
            empty_range_of_flight += drone.empty_range_of_flight
            full_range_of_flight += drone.full_range_of_flight
        common = empty_range_of_flight + not_completely_full_range_of_flight + full_range_of_flight
        pr1 = round(empty_range_of_flight * 100 / common, 1)
        pr2 = round(not_completely_full_range_of_flight * 100 / common, 1)
        pr3 = round(full_range_of_flight * 100 / common, 1)
        print(f'дальность полета не загруженными {pr1} %'
              f'\nдальность  полета загруженными не полностью {pr2} %'
              f'\nдальность  полета загруженными полностью {pr3} %')

    def move_at(self, target, speed=None):
        super().move_at(target=target, speed=speed)
        if 0 < self.payload < 100:
            self.not_completely_full_range_of_flight += self.distance_to(target)
        elif self.payload == 0:
            self.empty_range_of_flight += self.distance_to(target)
        else:
            self.full_range_of_flight += self.distance_to(target)


class Circle:
    positions_circle = {}
    edges = []

    def __init__(self, my_mothership, list_of_drones, index):
        self.my_mothership = my_mothership
        self.list_of_drones = list_of_drones
        self.index = index

    def _append_edges(self, amount_of_drones, sweep):
        amount_of_parts = amount_of_drones // 2
        edge_step = sweep // amount_of_parts
        index, from_ = 1, 0 if amount_of_drones % 2 != 0 else edge_step
        for _ in range(2):
            for i in range(from_, sweep + 1, edge_step):
                self.edges.append(i * index)
            index *= -1
            from_ = edge_step

    def _find_position(self, vec_position, amount_of_drones, drone, sweep):
        if not self.edges:
            self._append_edges(amount_of_drones, sweep)
        final_vector = vec_position
        final_vector.rotate(self.edges[drone.index])
        return final_vector

    def _find_new_positions_circle(self, start_point, end_point, distance, sweep):
        start_dot = start_point if start_point else self.my_mothership.coord
        end_dot = end_point if end_point else Point(theme.FIELD_WIDTH // 2, theme.FIELD_HEIGHT // 2)
        amount_of_drones = len(self.list_of_drones)
        vec = Vector.from_points(start_dot, end_dot)
        norm_vec = Vector(vec.x / vec.module, vec.y / vec.module)
        vec_position = Vector(norm_vec.x * distance, norm_vec.y * distance)
        for drone in self.list_of_drones:
            final_vector = self._find_position(vec_position, amount_of_drones, drone, sweep)
            final_point = (start_dot.x + final_vector.x, start_dot.y + final_vector.y)
            point = Point(final_point[0], final_point[1])
            self.positions_circle[distance].append(point)
            vec_position.rotate(-self.edges[drone.index])

    def get_position_circle(self, distance, start_point, end_point, sweep):
        amount_of_drones = len(self.list_of_drones)
        if distance not in self.positions_circle.keys() or len(self.positions_circle[distance]) >= amount_of_drones:
            self.positions_circle[distance], self.edges = [], []
            self._find_new_positions_circle(start_point, end_point, distance, sweep)
        return self.positions_circle[distance][self.index]