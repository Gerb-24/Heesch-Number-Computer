class Vector2:
    def __init__(self, x, y) -> None:
        self.x = x
        self.y = y
    
    def equals(self, other):
        return (self.x == other.x and self.y == other.y)
    
    def rotate(self):
        return Vector2(-self.y, self.x)
    
    def flip(self):
        return Vector2(-self.x, self.y)
    
    def add(self, other):
        return Vector2(self.x+other.x, self.y+other.y)
    
    def subtract(self, other):
        return Vector2(self.x-other.x, self.y-other.y)
    
    def __eq__(self, other):
        return isinstance(other, Vector2) and self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))
    
    def __str__(self) -> str:
        return f'Vector2({self.x},{self.y})'

class Tile:
    def __init__(self, coords) -> None:
        self.coords = coords

    def translate(self, translation):
        return Tile(tuple(coord.add(translation) for coord in self.coords))
    
    # computes if two tiles intersect with eachother
    def doIntersect(self, other):
        for coord in self.coords:
            for _coord in other.coords:
                if coord.equals(_coord):
                    return True
        return False
    
    def getBoundary(self, unitVectors):
        coordSet = set(self.coords)
        boundarySet = { coord.add(unitVector) for coord in self.coords for unitVector in unitVectors}
        return tuple(boundarySet.difference(coordSet))

class CollisionData:
    def __init__(self, collisions=None) -> None:
        self.collisions = collisions if collisions else createNewCollisions()

    def getCollisionsFromBaseOrientations(self, baseTile: "Tile", baseOrientations ):
        for tileType, baseOrientation in baseOrientations.items():
            for orientationCoord in baseOrientation.coords:
                for baseTileCoord in baseTile.coords:
                    translation = baseTileCoord.subtract(orientationCoord)
                    self.collisions[tileType].add( translation )
    
    # we still need to know the collisions for a different shapecode then base shapecode.
    def getFlippedCollisions(self):
        newCollisions = createNewCollisions()
        # when we flip the basetile, then the flipped collisions are for the flipped orientation
        for orientationIndex in self.collisions:
            newCollisions[flipIndex(orientationIndex)] = set(flipCoords(tuple(self.collisions[orientationIndex])))
        return CollisionData(newCollisions)
    
    def getRotatedCollisions(self):
        newCollisions = createNewCollisions()
        # when we rotate the basetile, then the rotated collisions are for the rotated orientation
        for orientationIndex in self.collisions:
            newCollisions[rotateIndex(orientationIndex)] = set(rotateCoords(tuple(self.collisions[orientationIndex])))
        return CollisionData(newCollisions)

    def getTranslatedCollisions(self, translation):
        newCollisions = createNewCollisions()
        for orientationIndex in self.collisions:
            newCollisions[orientationIndex] = set(translateCoords(tuple(self.collisions[orientationIndex]), translation))
        return CollisionData(newCollisions)

    def getCollisionsFromShapeCode(self, orientationIndex, translation):
        # This works but is not very DRY
        if orientationIndex == 0:
            return self.getTranslatedCollisions(translation)
        elif orientationIndex == 1:
            return self.getRotatedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 2:
            return self.getRotatedCollisions().getRotatedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 3:
            return self.getRotatedCollisions().getRotatedCollisions().getRotatedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 4:
            return self.getFlippedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 5:
            return self.getRotatedCollisions().getFlippedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 6:
            return self.getRotatedCollisions().getRotatedCollisions().getFlippedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 7:
            return self.getRotatedCollisions().getRotatedCollisions().getRotatedCollisions().getFlippedCollisions().getTranslatedCollisions(translation)

class ShapeCode:
    def __init__(self, orientationIndex, translation) -> None:
        self.orientation = orientationIndex
        self.translation = translation

    def __str__(self) -> str:
        return f'Shapecode({self.orientation}, {self.translation})'
    
    def printMotionCanvas(self) -> str:
        return f'[[{self.translation.x},{self.translation.y}],{enumNames[self.orientation]}],'

class Configuration:
    def __init__(self, configuration, boundary, index=0) -> None:
        self.configuration = configuration
        self.boundary = boundary

    def getNextConfigurations(self, baseOrientations):
        startCoord = self.boundary[0]
        validNewTiles =  getValidNewTiles(baseOrientations, self.configuration, startCoord)
        nextConfigurations = set()
        for shapeCode in validNewTiles:
            # If we build up the whole configuration structure as a tree, then the previous configuration is redudant information.
            newConfiguration = Configuration( [*self.configuration, shapeCode ], getNewBoundaryForValidNewTile(self.boundary, baseOrientations, shapeCode) )
            nextConfigurations.add(newConfiguration)
        return nextConfigurations

    def dfs(self, baseOrientations):
        if not self.boundary:
            print(self.printMotionCanvas())
            return self
        # print(f'still searching {len(self.boundary)}')
        nextConfigurations = self.getNextConfigurations(baseOrientations)
        return [ nextConfiguration.dfs(baseOrientations) for nextConfiguration in nextConfigurations ]
    
    def printMotionCanvas(self):
        baseTileText = self.configuration[0].printMotionCanvas()
        coronaText = '\n'.join([conf.printMotionCanvas() for conf in self.configuration[1:]])
        outText = f'''
  yield* coronas[0].addTiles([
    {baseTileText}
  ], 1)

  const output: [PossibleVector2, TileType][] =  [
    {coronaText}
]
'''
        return outText
        

def rotateCoords( coordsAsVectors ):
    return tuple(vec.rotate() for vec in coordsAsVectors)

def flipCoords( coordsAsVectors ):
    return tuple(vec.flip() for vec in coordsAsVectors)

def translateCoords( coordsAsVectors, translate ):
    return tuple(vec.add(translate) for vec in coordsAsVectors)

def rotateIndex( orientationIndex ):
    if orientationIndex <= 3:
        return (orientationIndex + 1) % 4
    # Apperently we cycle through the rotations the other way when flipped.
    # This should get fixed later?
    return (((orientationIndex - 4) - 1) % 4 ) + 4

def flipIndex( orientationIndex ):
    return (orientationIndex + 4)%8

def createNewCollisions():
    return { index: set() for index in range(amountOfOrientations) }


amountOfOrientations = 8
allOrientations = ((0,0), (1,0), (2,0), (3,0),
                   (0,1), (1,1), (2,1), (3,1))
enumNames = ('TileType.F0', 'TileType.F90', 'TileType.F180', 'TileType.F270',
             'TileType.T0', 'TileType.T90', 'TileType.T180', 'TileType.T270')

unitVectors = (Vector2(1,0), Vector2(-1,0), Vector2(0,1), Vector2(0, -1), Vector2(1,1), Vector2(-1,-1), Vector2(1,-1), Vector2(-1,1))
coords =     [[0, 0],
		[1, 0], [1,1], [1,2], [1,3],
         [2,2] ]


coordsAsVectors = tuple(Vector2(*coord) for coord in coords)


# Here we construct the tiles with the different rotations and amount of flips.
baseOrientations = {
    0:      Tile(coordsAsVectors),
    1:      Tile(rotateCoords(coordsAsVectors)),
    2:      Tile(rotateCoords(rotateCoords(coordsAsVectors))),
    3:      Tile(rotateCoords(rotateCoords(rotateCoords(coordsAsVectors)))),
    4:      Tile(flipCoords(coordsAsVectors)),
    5:      Tile(flipCoords(rotateCoords(coordsAsVectors))),
    6:      Tile(flipCoords(rotateCoords(rotateCoords(coordsAsVectors)))),
    7:      Tile(flipCoords(rotateCoords(rotateCoords(rotateCoords(coordsAsVectors)))))
}

baseTile = baseOrientations[0]
baseBoundary = baseTile.getBoundary(unitVectors)

# using these we compute the different translations in which these intersect with out base tile. 
# We can then forget about the tiles coordinates, as this is the only data that matters.

collisionData = CollisionData()
collisionData.getCollisionsFromBaseOrientations( baseTile, baseOrientations )


# For the first corona we take loop over the coordinates in the boundary and want to fill them in.
# So we want to loop over all the coordinates of all the baseOrientations, 
# and see if the translation to the corresponding square in the boundary 
# is not in the set of collisions for all of the tiles in the current configuration set. 
# If it is not, then we make a deepcopy of our configuration set, where we add this new tile.

# The baseConfiguration is a list of ShapeCodes,
# we then want to find all the new ShapeCodes so that they cover the startCoord but not the  tiles in the baseConfiguration

def getValidNewTiles(baseOrientations, baseConfiguration, startCoord: Vector2):
    validNewTiles = set()
    for orientationIndex, baseOrientation in baseOrientations.items():
        for coord in baseOrientation.coords:
            translation = startCoord.subtract(coord)
            doesNotCollide = True
            for shapeCode in baseConfiguration:
                # recomputing the collisionData for every new tile makes no sense
                if translation in collisionData.getCollisionsFromShapeCode(shapeCode.orientation, shapeCode.translation).collisions[orientationIndex]:
                    doesNotCollide = False
                    break
            if doesNotCollide:
                validNewTiles.add( ShapeCode(orientationIndex, translation) )
    return validNewTiles

# We now have the valid new tiles, so we add these to get a new configuration.
# To go on with adding tiles, we need to find the new startingCoord.
# Going to the next tile in the baseBoundary does not work, our added tile might already cover it.
# So we need to update our boundary and continue to use this

def getTilesFromShapeCode(baseOrientations, shapeCode: ShapeCode):
    baseOrientation = baseOrientations[ shapeCode.orientation ]
    tiles = translateCoords(baseOrientation.coords, shapeCode.translation)
    return tiles

def getNewBoundaryForValidNewTile(currentBoundary, baseOrientations, shapeCode: ShapeCode):
    newTiles = getTilesFromShapeCode(baseOrientations, shapeCode)
    newBoundary =  tuple(set(currentBoundary).difference(set(newTiles)))
    return newBoundary

baseConfiguration = Configuration( [ShapeCode(0, Vector2(0,0))], baseBoundary)
baseConfiguration.dfs(baseOrientations)