from messages import Messages

class Square:
    # This will make squares with the sides of length 2.
    # It is built by specifying the origin.

    def __init__(self, x, y):
        self.x = x
        self.y = y


    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def translate(self, xval, yval):
        newx = self.x + xval
        newy = self.y + yval
        return Square(newx, newy)

    def rot_90(self):
         return Square(-self.y, self.x)

    def flip(self):
        return Square(-self.x, self.y)
    
    def __str__(self) -> str:
        return f'Vector2({self.x//2},{-1*self.y//2})'

class FakePolyomino:
    def __init__(self, squares, shapecode = {"translation": (0, 0),"flipped": False,  "rotation": 0}):
        self.squares = squares
        self.shapecode = shapecode

    def rot_90(self):
        new_squares = []
        for square in self.squares:
            new_squares.append(square.rot_90())

        new_translation = (-self.shapecode["translation"][1],
        self.shapecode["translation"][0])

        new_flipped = self.shapecode["flipped"]
        new_rotation = (self.shapecode["rotation"]+90) % 360

        new_shapecode = {
        "translation": new_translation,
        "flipped": new_flipped,
        "rotation": new_rotation
        }

        return FakePolyomino(new_squares, shapecode = new_shapecode)

class Polyomino:
    # These will be polyominoes

    def __init__(self, squares, priority = [], shapecode = {"translation": (0, 0),"flipped": False,  "rotation": 0}, collision_data = []):

        self.squares = squares
        self.priority = priority
        self.shapecode = shapecode
        self.collision_data = self.collisions() if collision_data == [] else collision_data

    def __eq__(self, other):
        # we assume here that there is the tile is not mirror symmetric or rotationtionally symmetric
        return self.shapecode == other.shapecode

    def orientations(self, fake = False):
        orientation_list =[
        self,
        self.rot_90(fake = fake),
        self.rot_90(fake = fake).rot_90(),
        self.rot_90(fake = fake).rot_90().rot_90(),
        self.flip(fake = fake),
        self.flip(fake = fake).rot_90(),
        self.flip(fake = fake).rot_90().rot_90(),
        self.flip(fake = fake).rot_90().rot_90().rot_90(),
        ]

        new_orientations_list = []
        for orientation in orientation_list:
            if orientation not in new_orientations_list:
                new_orientations_list.append(orientation)
        return new_orientations_list

    def collisions(self):
        collision_dict = {
        0: set(),
        1: set(),
        2: set(),
        3: set(),
        4: set(),
        5: set(),
        6: set(),
        7: set(),
        }

        key_dict = {
        "0F": 0,
        "90F": 1,
        "180F": 2,
        "270F": 3,
        "0T": 4,
        "90T": 5,
        "180T": 6,
        "270T": 7,
        }

        for orientation in self.orientations(fake = True):
            pre_key = f"""{orientation.shapecode["rotation"]}{"T" if orientation.shapecode["flipped"] else "F"}"""
            key = key_dict[pre_key]
            for or_square in orientation.squares:
                for base_square in self.squares:
                    coord = (
                    base_square.x - or_square.x,
                    base_square.y - or_square.y,
                    )
                    collision_dict[key].add(coord)
        return collision_dict

    def translate(self, xval, yval):
        new_squares = []
        for square in self.squares:
            new_squares.append(square.translate(xval, yval))

        new_translation = (self.shapecode["translation"][0] + xval,
        self.shapecode["translation"][1] + yval)

        new_flipped = self.shapecode["flipped"]
        new_rotation = self.shapecode["rotation"]

        new_shapecode = {
        "translation": new_translation,
        "flipped": new_flipped,
        "rotation": new_rotation
        }

        new_collision_data = {}
        for bob in range(8):
            new_collision_data[bob] = {(coord[0] + xval, coord[1] + yval) for coord in self.collision_data[bob]}

        return Polyomino(new_squares, shapecode = new_shapecode, collision_data = new_collision_data)

    def translate_rel(self, square1, square2):
        return self.translate(square1.x - square2.x, square1.y - square2.y)

    def rot_90(self, fake = False):
        new_squares = []
        for square in self.squares:
            new_squares.append(square.rot_90())

        new_translation = (-self.shapecode["translation"][1],
        self.shapecode["translation"][0])

        new_flipped = self.shapecode["flipped"]
        new_rotation = (self.shapecode["rotation"]+90) % 360

        new_shapecode = {
        "translation": new_translation,
        "flipped": new_flipped,
        "rotation": new_rotation
        }

        if fake:
            return FakePolyomino(new_squares, shapecode = new_shapecode)
        else:
            new_collision_data = {}
            for i in range(8):
                if i < 4:
                    new_collision_data[i] = {(-coord[1], coord[0]) for coord in self.collision_data[ (i+3) % 4]}
                elif i <8:
                    new_collision_data[i] = {(-coord[1], coord[0]) for coord in self.collision_data[ ((i+3) % 4)+4]}

            return Polyomino(new_squares, shapecode = new_shapecode, collision_data = new_collision_data)

    def flip(self, fake = False):
        new_squares = []
        for square in self.squares:
            new_squares.append(square.flip())

        new_translation = (
        -self.shapecode["translation"][0],
        self.shapecode["translation"][1])

        new_flipped = not self.shapecode["flipped"]
        new_rotation = (360 - self.shapecode["rotation"]) % 360

        new_shapecode = {
        "translation": new_translation,
        "flipped": new_flipped,
        "rotation": new_rotation
        }

        if fake:
            return FakePolyomino(new_squares, shapecode = new_shapecode)
        else:
            new_collision_data = {
            0: {(-coord[0], coord[1]) for coord in self.collision_data[ 4 ] },
            1: {(-coord[0], coord[1]) for coord in self.collision_data[ 7 ] },
            2: {(-coord[0], coord[1]) for coord in self.collision_data[ 6 ] },
            3: {(-coord[0], coord[1]) for coord in self.collision_data[ 5 ] },
            4: {(-coord[0], coord[1]) for coord in self.collision_data[ 0 ] },
            5: {(-coord[0], coord[1]) for coord in self.collision_data[ 3 ] },
            6: {(-coord[0], coord[1]) for coord in self.collision_data[ 2 ] },
            7: {(-coord[0], coord[1]) for coord in self.collision_data[ 1 ] },
            }


            return Polyomino(new_squares, shapecode = new_shapecode, collision_data = new_collision_data)


    def config_collision(self, coord, config, key):
        if coord in self.collision_data[key]:
            return False
        for shape in config:
            if coord in shape.collision_data[key]:
                return False
        return True


    def corona_maker(self, base_orientations, heesch=False, printing=True):

        def not_occupied_in(elem, config, extra = False):
            config_squares = self.squares.copy() if extra else []

            for shape in config:
                config_squares.extend(shape.squares)
            for square in config_squares:
                if elem == square:
                    return False

            return True


        def check_combinations(base_orientations, config, outs_square):
            new_possible_config = []

            key_dict = {
            "0F": 0,
            "90F": 1,
            "180F": 2,
            "270F": 3,
            "0T": 4,
            "90T": 5,
            "180T": 6,
            "270T": 7,
            }

            for index in range(len(base_orientations)):
                orientation = base_orientations[index]
                pre_key = f"""{orientation.shapecode["rotation"]}{"T" if orientation.shapecode["flipped"] else "F"}"""
                key = key_dict[pre_key]
                for ns_square in orientation.squares:
                    coord = (outs_square.x - ns_square.x, outs_square.y - ns_square.y)
                    if self.config_collision(coord, config, key):
                        new_config = config.copy()
                        new_config.append(orientation.translate(*coord))
                        new_possible_config.append(new_config)
                    else:
                        continue
            return new_possible_config


        possible_configs = []
        outside_list = self.outside()
        for i in range(len(outside_list)):
            if heesch:
                message = f" \r we are now at {int((i+1)/len(outside_list)*100)}% "
                print(message, end="")
            outs_square = outside_list[i]
            if len(possible_configs) == 0:
                possible_configs.extend(check_combinations(base_orientations, [], outs_square))
                if printing:
                    print(len(possible_configs))

            else:
                new_possible_configs = []
                for config in possible_configs:
                    # if the current boundary tile is already occupied in the current configuration, 
                    # then we add it to the next configurations, otherwise we do need to find a new validTile to occupy it
                    if not_occupied_in(outs_square, config):
                        new_possible_configs.extend(check_combinations(base_orientations, config, outs_square))
                    else:
                        new_possible_configs.append(config)
                if new_possible_configs == []:
                    return []
                possible_configs = new_possible_configs.copy()
                if printing:
                    print(len(possible_configs))
        return [[config] for config in possible_configs] if heesch else possible_configs

    def outside(self):

        bigsquarelist = []
        for vert in self.squares:
            bigsquarelist.extend(bigsquare_maker(vert.x, vert.y))
        res = self.priority.copy()
        for square in bigsquarelist:
            if square not in res:
                res.append(square)
        res = [square for square in res if square not in self.squares]
        return res

    def heesch_corona(self, possible_configs, coronalist, c_index):
        def sec_corona_maker(conf_corona_config):

            # We now translate, rotate and flip the whole config so that the start_tile
            # is the same as the base tile.

            def transform(st_tile, configuration):

                def flip(flipped, tile):

                    # A function that flips a tile if flipping is true,
                    # otherwise it just returns the tile

                    if flipped:
                        return tile.flip()
                    else:
                        return tile

                def rotate(rotation, tile):

                    # A function that can rotate a tile
                    # by 90 degrees multiple times

                    if rotation == 0:
                        return tile
                    else:
                        return rotate(rotation - 90, tile.rot_90())

                def translate(translation, tile):

                    # A function that translates a tile
                    # just to keep things more organised

                    return tile.translate(*translation)

                sc = st_tile.shapecode
                translation = ( -sc["translation"][0], -sc["translation"][1])
                flipped = sc["flipped"]
                rotation = 360 - sc["rotation"]

                transformed_config = [flip(flipped, rotate(rotation, translate(translation, tile))) for tile in configuration]

                return transformed_config

            
            # Now we take a for loop over all the coronas in possible_configs
            # When this corona fits we append it to the list of the new possible configurations
            # A corona fits if all of the tiles that do collide with fixed configuration
            # are actually part of that configuration.

            def collides_with(tile, configuration):

                # In here we check if a tile collides with one of the tiles in a configuration

                key_dict = {
                "0F": 0,
                "90F": 1,
                "180F": 2,
                "270F": 3,
                "0T": 4,
                "90T": 5,
                "180T": 6,
                "270T": 7,
                }

                coord = tile.shapecode["translation"]
                pre_key = f"""{tile.shapecode["rotation"]}{"T" if tile.shapecode["flipped"] else "F"}"""
                key = key_dict[pre_key]

                for conf_tile in configuration:
                    if coord in conf_tile.collision_data[key]:
                        return True
                return False

            def retransform(st_tile, corona):

                def flip(flipped, tile):

                    # A function that flips a tile if flipping is true,
                    # otherwise it just returns the tile

                    if flipped:
                        return tile.flip()
                    else:
                        return tile

                def rotate(rotation, tile):

                    # A function that can rotate a tile
                    # by 90 degrees multiple times

                    if rotation == 0:
                        return tile
                    else:
                        return rotate(rotation - 90, tile.rot_90())

                def translate(translation, tile):

                    # A function that translates a tile
                    # just to keep things more organised

                    return tile.translate(*translation)

                sc = st_tile.shapecode
                translation = sc["translation"]
                flipped = sc["flipped"]
                rotation = sc["rotation"]

                retransformed_corona = [translate(translation, rotate(rotation, flip(flipped, tile))) for tile in corona]

                return retransformed_corona

            
            # new_c_configs now contains the all the possible coronas
            # around the start_tile that fit with the fixed configuration.
            

            
            # In the same way as we did with the corona_maker function
            # we now want to combine these with the old config to get a list of possible_c_configs
            # for every tile in the corona of the config we will now change
            # the possible_c_configs to be all of the currently possible_c_configs

            conf_corona = conf_corona_config[c_index]
            pre_possible_c_configs = conf_corona_config.copy()
            pre_possible_c_configs.append([])
            possible_c_configs = [ pre_possible_c_configs ]

            # So here we loop over the tiles in the corona of the config
            # where we start at tile 1, since we already have done tile 0

            for start_num in range(len(conf_corona)):
                start_tile = conf_corona[start_num]

                # We now go in a loop over all the c_configs in possible_c_configs """
                new_possible_c_configs = []
                for c_config in possible_c_configs:
                    
                    # Since we now want to work with c_config with corona structure
                    # we need to define something new that we want to translate
                    # which we call the absolute c config

                    abs_c_config = [self]
                    for abs_index in range(c_index+2):
                        abs_c_config.extend(c_config[abs_index])
                    transformed_abs_c_config = transform(start_tile, abs_c_config)

                    # For this recentered absoltue c_config we will append all of the
                    # retransformed corona that fit with the transformed c_config

                    new_c_configs = []
                    for corona in coronalist:
                        if all(c_tile in transformed_abs_c_config for c_tile in corona if collides_with(c_tile, transformed_abs_c_config)):
                            retransformed_corona = [ tile for tile in retransform( start_tile, corona ) if tile not in abs_c_config ]
                            
                            # We first copy the current c_config in corona structure
                            # and add the new retransformed_corona to it
                            # then we append this to new_c_configs, this then contains all the new_possible_c_configs at the end
                            
                            c_config_with_added_corona = c_config.copy()
                            c_config_with_added_corona[c_index+1] = c_config_with_added_corona[c_index+1] + retransformed_corona
                            new_c_configs.append( c_config_with_added_corona )

                    new_possible_c_configs.extend( new_c_configs )

                # If we no corona fits anymore for every possible c_config,
                # then the original config is dead

                if new_possible_c_configs == []:
                    possible_c_configs = []
                    break

                possible_c_configs = new_possible_c_configs.copy()

            return possible_c_configs

        second_configs = []
        for conf_corona_index in range(len(possible_configs)):
            Messages.coronaCount(conf_corona_index, len(possible_configs))
            # the conf_corona is the last corona in the possible_config
            conf_corona_config = possible_configs[conf_corona_index]


            second_configs.extend(sec_corona_maker(conf_corona_config))
        
        Messages.totalNewCoronaAmount(len(second_configs))
        return second_configs

    def heesch_computer(self):
        def has_holes(config, output = False ):
            total_squares = set()
            for shape in config+[self]:
                total_squares = total_squares.union({(square.x, square.y) for square in shape.squares})
            extended_squares = set()
            for square in total_squares:
                extended_squares = extended_squares.union({
                square,
                (square[0], square[1]-2),
                (square[0]+2, square[1]-2),
                (square[0]+2, square[1]),
                (square[0]+2, square[1]+2),
                (square[0], square[1]+2),
                (square[0]-2, square[1]+2),
                (square[0]-2, square[1]),
                (square[0]-2, square[1]-2),
                })
            without_inside = extended_squares.difference(total_squares)

            # getting a starting position
            xmin = min({coord[0] for coord in without_inside})
            ymin = min({coord[1] for coord in without_inside if coord[0] == xmin})
            mincoord = (xmin, ymin)

            # computing the connected component and removing them
            coords = { mincoord }
            while coords != set():
                without_inside = without_inside.difference(coords)
                new_coords = set()
                for coord in coords:
                    nextcoords = {
                    (coord[0], coord[1]-2),
                    (coord[0]+2, coord[1]),
                    (coord[0], coord[1]+2),
                    (coord[0]-2, coord[1]),
                    }
                    nextcoords = {n_coord for n_coord in nextcoords if n_coord in without_inside}
                    new_coords = new_coords.union(nextcoords)
                coords = new_coords

            return without_inside if output else (without_inside != set())
        
        # We first build the first coronas
        Messages.firstCorona()
        coronalist = self.corona_maker(self.orientations(), printing=True)

        # The case when the tile has Heesch number 0
        if coronalist == []:
            return []
        

        no_hole_coronalist = [corona for corona in coronalist if not has_holes(corona)]
        possible_configs = [ [corona] for corona in no_hole_coronalist]

        # Now we build the higher order coronas
        i = 0
        while True:
            Messages.laterCorona(i)
            new_possible_configs = self.heesch_corona(possible_configs, no_hole_coronalist, i)
            if new_possible_configs == []:
                Messages.heeschNumber(i)
                return possible_configs
            else:
                possible_configs = new_possible_configs.copy()
                i += 1

def bigsquare_maker(x, y):
    squares = [
    Square(x, y),
    Square(x, y-2),
    Square(x+2, y-2),
    Square(x+2, y),
    Square(x+2, y+2),
    Square(x, y+2),
    Square(x-2, y+2),
    Square(x-2, y),
    Square(x-2, y-2),
    ]
    return squares

def test_check_combinations(tile, base_orientations, config, outs_square):
    new_possible_config = []

    key_dict = {
    "0F": 0,
    "90F": 1,
    "180F": 2,
    "270F": 3,
    "0T": 4,
    "90T": 5,
    "180T": 6,
    "270T": 7,
    }

    for index in range(len(base_orientations)):
        orientation = base_orientations[index]
        pre_key = f"""{orientation.shapecode["rotation"]}{"T" if orientation.shapecode["flipped"] else "F"}"""
        key = key_dict[pre_key]
        for ns_square in orientation.squares:
            coord = (outs_square.x - ns_square.x, outs_square.y - ns_square.y)
            if tile.config_collision(coord, config, key):
                new_possible_config.append(orientation.translate(*coord))
            else:
                continue
    return new_possible_config